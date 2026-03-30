//
// LE created by chenghau
//
#include "util.h"
#include "LEInputData.h"
#include "LENodeData.h"
#include "LEBCSetup.h"
#include "SSLEEquation.h"
#include "OctToPhysical.h"
#include <algorithm>
#include <array>
#include <sfcTreeLoop_matvec_io.h>
#include <petscvec.h>
#include "CalcStress.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
using namespace PETSc;

namespace
{
void writeMaterialPropertyField(DA *octDA,
                                const std::vector<TREENODE> &treePart,
                                const DomainExtents &domainExtents,
                                const SSLEEquation &equation)
{
  if (!octDA->isActive())
  {
    return;
  }

  const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
  const DENDRITE_UINT sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();

  std::vector<double> nodeCoords(static_cast<std::size_t>(nPe) * DIM, 0.0);
  std::vector<double> elemData(static_cast<std::size_t>(octDA->getLocalElementSz()) * 2, 0.0);
  OctToPhysical octToPhysical(domainExtents);

  DENDRITE_UINT elemIdx = 0;
  ot::MatvecBaseCoords<DIM> loop(sz, octDA->getElementOrder(), false, 0, tnCoords,
                                 &(*treePart.cbegin()), treePart.size(),
                                 *partFront, *partBack);

  while (!loop.isFinished())
  {
    if (loop.isPre() && loop.subtreeInfo().isLeaf())
    {
      if (!nodeCoords.empty())
      {
        const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
        std::copy(nodeCoordsFlat, nodeCoordsFlat + nodeCoords.size(), nodeCoords.begin());
        octToPhysical.convertCoordsToPhys(nodeCoords.data(), nPe);

        std::array<double, DIM> centroid{};
        for (DENDRITE_UINT node = 0; node < nPe; ++node)
        {
          for (int d = 0; d < DIM; ++d)
          {
            centroid[d] += nodeCoords[node * DIM + d];
          }
        }
        const double inv = 1.0 / static_cast<double>(nPe);
        for (auto &value : centroid)
        {
          value *= inv;
        }

        const auto matProps = equation.queryMaterialAt(centroid[0], (DIM >= 2) ? centroid[1] : 0.0);
        elemData[elemIdx * 2 + 0] = matProps.E;
        elemData[elemIdx * 2 + 1] = matProps.nu;
      }
      ++elemIdx;
      loop.next();
    }
    else
    {
      loop.step();
    }
  }

  static const char *varNames[]{"young_modulus", "poisson_ratio"};
  IO::writeVecTopVtu(octDA, treePart, elemData.data(), "MaterialProp", "material_props",
                     varNames, domainExtents, true, false, 2);
}
}

int main(int argc, char *argv[])
{

  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  LEInputData inputData;

  if (!(inputData.ReadFromFile()))
  {
    if (!rank)
    {
      throw TALYFEMLIB::TALYException() << "Can't read the config file \n";
    }
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  inputData.PrintInputData();
  ///------------------------------ Command line option to restart from a checkpoint -------------------------------////
  bool resume_from_checkpoint = false;
  {
    PetscBool resume = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-resume_from_checkpoint", &resume, nullptr);
    resume_from_checkpoint = (resume == PETSC_TRUE);
  }
  Checkpointer checkpointer(inputData.CheckpointNumbackup, "CheckPoint");
  bool restart_sum = false;
  {
    PetscBool restart_sum_ = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-restart_sum", &restart_sum_, nullptr);
    restart_sum = (restart_sum_ == PETSC_TRUE);
  }
  /// --------------------------------------------------------------------------------------------------------------////

  const DENDRITE_UINT eleOrder = inputData.elemOrder;
  const DENDRITE_UINT levelBase = inputData.mesh_def.refine_lvl_base;
  const bool mfree = inputData.ifMatrixFree;

// Linear elasticity
#if (DIM == 2)
  static const char *varname[]{"UX", "UY"};
#endif
#if (DIM == 3)
  static const char *varname[]{"UX", "UY", "UZ"};
#endif

  ///------------------------------------------Creation/loading of mesh----------------------------------------------///
  DomainInfo cubeDomain, physDomain;
  cubeDomain.min.fill(0);
  cubeDomain.max.fill(1);
  physDomain.min.fill(0);
  physDomain.max.fill(1);

  DomainExtents domainExtents(cubeDomain, physDomain);

  DA *octDA = nullptr;
  DistTREE dTree;
  SubDomain subDomain(domainExtents, resume_from_checkpoint);

  std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
    return (subDomain.functionToRetain(physCoords, physSize));
  };


//// time the mesh construction
    double start_time = MPI_Wtime();
    octDA = createSubDA(dTree, functionToRetain, levelBase, eleOrder);
    // subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);

    double end_time = MPI_Wtime();
    PrintStatus("Time to create Mesh = ",end_time-start_time);

    TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());

    /// --------------------------------------------------------------------------------------------------------------////
    const auto &treePartition = dTree.getTreePartFiltered();
    IO::writeBoundaryElements(octDA, dTree.getTreePartFiltered(), "boundary", "subDA", domainExtents);

    SubDomainBoundary boundary(&subDomain, octDA, domainExtents);

    LEBCSetup LEBC(&boundary, &inputData);

    /// Gridfield setup
  TalyMesh<LENodeData> talyMesh(octDA->getElementOrder());

  // domain boundary
  SubDomainBoundary *subDomainBoundary = nullptr;
  subDomainBoundary = &boundary;

  // ndof
  static const DENDRITE_UINT ndof = LENodeData::LE_DOF;

  // time info
  //// dummy variable for output span
  std::vector<int> OutputSpan;
  TimeInfo ti(0.0, inputData.dt, inputData.totalT,OutputSpan);

  // for SC

  auto leEq = new TalyEquation<SSLEEquation, LENodeData>(&talyMesh,octDA, dTree.getTreePartFiltered(), domainExtents, ndof,
                                                         &ti, false, nullptr, &inputData);


  LinearSolver *leSolver = setLinearSolver(leEq, octDA,dTree,ndof, mfree,false,true);


  inputData.solverOptionsLE.apply_to_petsc_options("-le_");
  {
    KSP SSLEEquation_ksp = leSolver->ksp();
    KSPSetOptionsPrefix(SSLEEquation_ksp, "le_");
    KSPSetFromOptions(SSLEEquation_ksp);
  }

  writeMaterialPropertyField(octDA, treePartition, domainExtents, *leEq->equation());

  //// now what we want to do is loop around all the three BCTypes
  std::vector<BCTYPE> bctypes ={e100,e001,e010};
  double Ceff[3][3] = {{0}};   // columns correspond to load cases
  double eps0 = 1.0;           // or 1e-3
  double gam0 = 1.0;           // or 1e-3

  auto fillColumn = [&](int col, const double sig[3], double denom) {
    Ceff[0][col] = sig[0] / denom;
    Ceff[1][col] = sig[1] / denom;
    Ceff[2][col] = sig[2] / denom;
  };

  int col = 0;
  for (BCTYPE bctype : bctypes)
  {
    PrintStatus("Setting BC for Type",bctype);
    leSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary
                                            {
                                              Boundary b;
                                              LEBC.setBoundaryConditions(b, pos,bctype);
                                              return b; });

    leSolver->solve();
    std::string prefix = "displacement"+std::to_string(bctype);
    petscVectopvtu(octDA, dTree.getTreePartFiltered(), leSolver->getCurrentSolution(), "Displacements",
                   prefix.c_str(),varname, domainExtents, false, false, ndof);



    // CalcStress calcStress(octDA, dTree.getTreePartFiltered(), {VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)}, domainExtents,
    //                       &subDomain, &inputData, ti);
    // CalcHomogenizedResponse(DA *octDA,
    //                         const std::vector<TREENODE> &treePart,
    //                         const VecInfo &v,
    //                         const DomainExtents &domain,
    //                         const SubDomain *subDomain,
    //                         LEInputData *idata,
    //                         const TimeInfo ti)

    Vec solution = leSolver->getCurrentSolution();
    VecInfo solutionInfo(solution, LENodeData::LE_DOF, 0);

    CalcHomogenizedResponse calc_homogenized_response(octDA, dTree.getTreePartFiltered(), solutionInfo,
      domainExtents,&subDomain,&inputData, ti, leEq->equation());

    auto avg = calc_homogenized_response.getAverages();


    double sig[3] = {avg.sxx, avg.syy, avg.txy};

    if (bctype == e100) fillColumn(0, sig, eps0);
    if (bctype == e010) fillColumn(1, sig, eps0);
    if (bctype == e001) fillColumn(2, sig, gam0);


  }
  double C11 = Ceff[0][0];
  double C12 = Ceff[0][1];
  double C66 = Ceff[2][2];

  double nu_eff = 0.0;
  if (inputData.caseType == CaseType::PLANESTRESS) {
    nu_eff = C12 / C11;
  } else { // PLANESTRAIN
    nu_eff = C12 / (C11 + C12);
  }

  double E_eff = 2.0 * C66 * (1.0 + nu_eff);

  if (rank == 0) {
    auto rel = [](double a, double b) {
      return std::abs(a - b) / (std::max({1.0, std::abs(a), std::abs(b)}));
    };

    const double tol = 1e-2;

    const double rC12 = rel(Ceff[0][1], Ceff[1][0]);
    const double rC11 = rel(Ceff[0][0], Ceff[1][1]);

    const double C13 = Ceff[0][2];
    const double C23 = Ceff[1][2];

    const bool sym_ok   = (rC12 < tol);
    const bool diag_ok  = (rC11 < tol);
    const bool shear_ok = (std::abs(C13) < tol * std::abs(Ceff[0][0])) &&
                          (std::abs(C23) < tol * std::abs(Ceff[0][0]));

    const bool iso_ok = sym_ok && diag_ok && shear_ok;

    std::cout << "\n================ Homogenization Summary ================\n";
    std::cout << "Effective moduli:\n";
    std::cout << "  E_eff  = " << E_eff  << "\n";
    std::cout << "  nu_eff = " << nu_eff << "\n\n";

    std::cout << "Material symmetry:\n";
    std::cout << "  " << (iso_ok ? "Approximately isotropic" : "Anisotropic") << "\n\n";

    std::cout << "Diagnostics (tol = " << tol << "):\n";
    std::cout << "  rel(C12 - C21) = " << rC12 << "\n";
    std::cout << "  rel(C11 - C22) = " << rC11 << "\n";
    std::cout << "  C13            = " << C13  << "\n";
    std::cout << "  C23            = " << C23  << "\n";
    std::cout << "========================================================\n\n";
  }




  delete leEq;
  delete leSolver;

  dendrite_finalize(octDA);
}
