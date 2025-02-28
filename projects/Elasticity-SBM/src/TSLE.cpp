//
// LE created by chenghau
//
#include "util.h"
#include "LEInputData.h"
#include "LENodeData.h"
#include "LEBCSetup.h"
#include "TSLEEquation.h"
#include <petscvec.h>
#include <IMGA/IMGASolverUtils.h>

using namespace PETSc;

// [fix] const
void GetVel(Vec prev1LEvel, Vec prev1LEacc, const Vec LEacc, const double dtime, const double beta1)
{
  // [bug fix] (y,a,x) VecAYPX: y=a*y+x / VecAXPY:y=a*x+y
  VecAXPY(prev1LEvel, dtime * (1 - beta1), prev1LEacc); // prev1LEvel=prev1LEvel+dtime*(1-beta1)*prev1LEacc
  VecAXPY(prev1LEvel, dtime * beta1, LEacc);            // prev1LEvel=prev1LEvel+dtime*beta1*LEacc
}

void GetDisp(Vec prev1LEdisp, Vec prev1LEvel, Vec prev1LEacc, const Vec LEacc, const double dtime, const double beta2)
{
  // [bug fix] (y,a,x) VecAYPX: y=a*y+x / VecAXPY:y=a*x+y
  VecAXPY(prev1LEdisp, dtime, prev1LEvel);                             // prev1LEdisp=prev1LEdisp+dtime*prev1LEvel
  VecAXPY(prev1LEdisp, 0.5 * dtime * dtime * (1 - beta2), prev1LEacc); // prev1LEdisp=prev1LEdisp+0.5*dtime*dtime*(1-beta2)*prev1LEacc
  VecAXPY(prev1LEdisp, 0.5 * dtime * dtime * beta2, LEacc);            // prev1LEdisp=prev1LEdisp+0.5*dtime*dtime*beta2*LEacc
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  dendrite_init(argc, argv);
  int rank = TALYFEMLIB::GetMPIRank();
  int nProc = TALYFEMLIB::GetMPISize();
  MPI_Status status;

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

  const DENDRITE_UINT eleOrder = inputData.elemOrder;
  const DENDRITE_UINT levelBase = inputData.mesh_def.refine_lvl_base;
  const bool mfree = inputData.ifMatrixFree;

  // Linear elasticity
  static const char *varname[]{"UX", "UY"};

  /// for dynamic LE
  TimeInfo ti(0.0, inputData.dt, inputData.totalT);
  Vec LEacc, prev1LEacc;
  Vec LEvel, prev1LEvel;
  Vec LEdisp, prev1LEdisp;

  ///------------------------------------------Creation/loading of mesh----------------------------------------------///
  DomainExtents domainExtents(inputData.mesh_def.fullDADomain, inputData.mesh_def.physDomain);
  DA *octDA = nullptr;
  DistTREE dTree;
  SubDomain subDomain(domainExtents);

  /// immersed boundary
  std::vector<GEOMETRY::STL *> stls; // 3D
  std::vector<GEOMETRY::MSH *> mshs; // 2D

  /// IBM setup
  IMGA imga(domainExtents);
  std::vector<GEOMETRY::Geometry *> ibm_geoms;
  std::vector<GeomRefinement> ibm_refinements;
  std::vector<std::vector<double>> node_geom;

  for (const auto &geom_def : inputData.ibm_geom_def)
  {
#if (DIM == 2)
    if (geom_def.type == IBMGeomDef::Type::MESHOBJECT_2D)
    {
      std::vector<double> node;
      util_funcs::read_msh(geom_def.mesh_path, node);
      node_geom.emplace_back(node);

      mshs.push_back(new GEOMETRY::MSH(geom_def.mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D));
      std::array<DENDRITE_REAL, DIM> shift;

      shift[0] = geom_def.InitialDisplacement[0];
      shift[1] = geom_def.InitialDisplacement[1];

      auto geom_retain_side = RetainSide::OUT;
      if (geom_def.outer_boundary)
      {
        geom_retain_side = RetainSide::IN;
      }

      ibm_geoms.emplace_back(new GEOMETRY::Geometry(mshs.back(), Point<DIM>(shift), geom_retain_side));
      ibm_refinements.emplace_back(geom_def.geomRefine);
      GeomRefinement geomRefine;
    }
#endif

#if (DIM == 3)
    if (geom_def.type == IBMGeomDef::Type::MESHOBJECT)
    {
      stls.push_back(new GEOMETRY::STL(geom_def.mesh_path, GEOMETRY::InOutTest::RAY_TRACING));
      std::array<DENDRITE_REAL, DIM> shift;
      shift[0] = geom_def.InitialDisplacement[0];
      shift[1] = geom_def.InitialDisplacement[1];
      shift[2] = geom_def.InitialDisplacement[2];
      ibm_geoms.emplace_back(new GEOMETRY::Geometry(stls.back(), Point<DIM>(shift)));
      ibm_refinements.emplace_back(geom_def.geomRefine);
    }
#endif
  }

  /// carving subda
  std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords, double scale)
  {
    return (subDomain.functionToRetain(octCoords, scale));
  };

  for (int i = 0; i < ibm_geoms.size(); i++)
  {
    imga.addGeometry(ibm_geoms.at(i), ibm_refinements.at(i));
  }

  octDA = createSubDA(dTree, functionToRetain, levelBase, eleOrder);
  subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);

  util_funcs::performRefinementSubDAIBM(octDA, dTree.getTreePartFiltered(), domainExtents, dTree, inputData, &subDomain, ibm_geoms);
  const auto &treePartition = dTree.getTreePartFiltered();
  IO::writeBoundaryElements(octDA, treePartition, "boundary", "subDA", domainExtents);
  SubDomainBoundary boundary(&subDomain, octDA, domainExtents);

  TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());

  /// --------------------------------------------------------------------------------------------------------------////

  LEBCSetup LEBC(&boundary, &inputData);

  /// Gridfield setup
  TalyMesh<LENodeData> talyMesh(octDA->getElementOrder());

  // domain boundary
  SubDomainBoundary *subDomainBoundary = nullptr;
  subDomainBoundary = &boundary;

  //ndof
  static const DENDRITE_UINT ndof = LENodeData::valueno();

  imga.initIMGAComputation(octDA, treePartition);
  Marker elementMarker(octDA, treePartition, domainExtents, &imga, MarkerType::GAUSS_POINT);
  //printGaussPointsToFile("GaussPoints.txt", imga.getSurfaceGaussPoints());

  std::vector<PetscInt> dirichletNodes;
  getIBMDirichletNodes(octDA, treePartition, dirichletNodes, &elementMarker);

  auto leEq =
      new TalyEquation<TSLEEquation, LENodeData>(&talyMesh, octDA, treePartition,
                                               subDomain.domainExtents(), ndof,&ti,true, subDomainBoundary,
                                               &inputData);

  leEq->assignIBMConstructs(&imga, elementMarker.getMarkers().data()); // change to "le"

  LinearSolver *leSolver = setLinearSolver(leEq, octDA, ndof, mfree);
  leSolver->setIBMDirichletNodes(dirichletNodes);
  leEq->setTime(&ti);
  inputData.solverOptionsLE.apply_to_petsc_options("-le_");
  {
    KSP TSLEEquation_ksp = leSolver->ksp();
    KSPSetOptionsPrefix(TSLEEquation_ksp, "le_");
    KSPSetFromOptions(TSLEEquation_ksp);
  }

  PrintStatus("Setting BC for LE");
  leSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary
                                          {
                                            Boundary b;
                                            LEBC.setBoundaryConditions(b, pos);
                                            return b;
                                          });

  ///------------------------------------------------- dynamic LE -----------------------------------------------------------///

  octDA->petscCreateVector(LEacc, false, false, LENodeData::LE_DOF);
  octDA->petscCreateVector(prev1LEacc, false, false, LENodeData::LE_DOF);
  octDA->petscCreateVector(LEvel, false, false, LENodeData::LE_DOF);
  octDA->petscCreateVector(prev1LEvel, false, false, LENodeData::LE_DOF);
  octDA->petscCreateVector(LEdisp, false, false, LENodeData::LE_DOF);
  octDA->petscCreateVector(prev1LEdisp, false, false, LENodeData::LE_DOF);

  /// Set initial conditions of all vectors
  VecSet(LEacc, 0.0);
  VecSet(prev1LEacc, 0.0);
  VecSet(LEvel, 0.0);
  VecSet(prev1LEvel, 0.0);
  VecSet(LEdisp, 0.0);
  VecSet(prev1LEdisp, 0.0);

  /// bridge btw main and nodedata (dendrite <-> taly)
  leEq->setVectors({VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::AX),
                    VecInfo(LEvel, LENodeData::LE_DOF, LENodeData::VX),
                    VecInfo(LEdisp, LENodeData::LE_DOF, LENodeData::UX),
                    VecInfo(prev1LEacc, LENodeData::LE_DOF, LENodeData::AX_PRE1),
                    VecInfo(prev1LEvel, LENodeData::LE_DOF, LENodeData::VX_PRE1),
                    VecInfo(prev1LEdisp, LENodeData::LE_DOF, LENodeData::UX_PRE1)},
                   SYNC_TYPE::VECTOR_ONLY);

  util_funcs::save_timestep(octDA, treePartition, LEdisp, ndof, ti, subDomain, "le", varname);

  while (ti.getCurrentTime() < ti.getEndTime())
  {
    ti.print();

    VecCopy(LEacc, prev1LEacc);
    VecCopy(LEvel, prev1LEvel);
    VecCopy(LEdisp, prev1LEdisp);

    leSolver->solve();

    double dtime = inputData.dt[0];

    VecCopy(leSolver->getCurrentSolution(), LEacc);
    double beta1 = 0.5, beta2 = 0.5;

    // [bug fix] using veccopy !
    GetVel(LEvel, prev1LEacc, LEacc, dtime, beta1);
    GetDisp(LEdisp, prev1LEvel, prev1LEacc, LEacc, dtime, beta2);

    ti.increment();
    double scaleFactor = inputData.scaleFactor;

    /// ----------------------------------------------file writing---------------------------------------------------///
    int count_geo = 0;
    if ((ti.getCurrentTime() >= inputData.OutputStartTime) && (ti.getTimeStepNumber() % inputData.OutputInterval == 0))
    {
      util_funcs::save_timestep(octDA, treePartition, LEdisp, ndof, ti, subDomain, "le", varname); 
    }
  }
  ///------------------------------------------------- dynamic LE:end -----------------------------------------------------------///

  delete leEq;
  delete leSolver;
  VecDestroy(&LEacc);
  VecDestroy(&prev1LEacc);
  VecDestroy(&LEvel);
  VecDestroy(&prev1LEvel);
  VecDestroy(&LEdisp);
  VecDestroy(&prev1LEdisp);
  dendrite_finalize(octDA);
}
