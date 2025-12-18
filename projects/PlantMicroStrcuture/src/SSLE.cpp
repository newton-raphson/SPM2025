//
// LE created by chenghau
//
#include "util.h"
#include "LEInputData.h"
#include "LENodeData.h"
#include "LEBCSetup.h"
#include "SSLEEquation.h"
#include <petscvec.h>
#include "CalcStress.h"
#ifdef DEEPTRACE
#include <onnxruntime_cxx_api.h>
#endif


using namespace PETSc;

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

#ifndef PROFILING
  /// this is for quick testing through different refine lvl
  if (argc == 2 and inputData.BaselvlFromArgument)
  {
    inputData.mesh_def.refine_lvl_base = std::atoi(argv[1]);
    TALYFEMLIB::PrintStatus("---------------------------------------");
    TALYFEMLIB::PrintStatus("Refine level base from argument = ", inputData.mesh_def.refine_lvl_base);
    TALYFEMLIB::PrintStatus("---------------------------------------");
  }
  /// this is for quick testing through different refine lvl and lambda
  if (argc == 3 and inputData.BaselvlFromArgument)
  {
    inputData.mesh_def.refine_lvl_base = std::atoi(argv[1]);
    inputData.RatioGPSBM = (double)std::atoi(argv[2]) / 100;
    TALYFEMLIB::PrintStatus("---------------------------------------");
    TALYFEMLIB::PrintStatus("Refine level base from argument = ", inputData.mesh_def.refine_lvl_base);
    TALYFEMLIB::PrintStatus("Lambda from argument = ", inputData.RatioGPSBM);
    TALYFEMLIB::PrintStatus("---------------------------------------");
  }
#else
  TALYFEMLIB::PrintStatus("---------------------------------------");
  TALYFEMLIB::PrintStatus("please put -log_view after the run to print out the time");
  TALYFEMLIB::PrintStatus("---------------------------------------");
#endif

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
  DomainExtents domainExtents(inputData.mesh_def.fullDADomain, inputData.mesh_def.physDomain);
  DA *octDA = nullptr;
  DistTREE dTree;
  SubDomain subDomain(domainExtents, resume_from_checkpoint);


//// time the mesh construction
    double start_time = MPI_Wtime();
    octDA = createSubDA(dTree, functionToRetain, levelBase, eleOrder);
    subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);

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


  LinearSolver *leSolver = setLinearSolver(leEq, octDA,dTree,ndof, mfree);


  inputData.solverOptionsLE.apply_to_petsc_options("-le_");
  {
    KSP SSLEEquation_ksp = leSolver->ksp();
    KSPSetOptionsPrefix(SSLEEquation_ksp, "le_");
    KSPSetFromOptions(SSLEEquation_ksp);
  }

  PrintStatus("Setting BC for LE");
  leSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary
                                          {
                                            Boundary b;

                                            LEBC.setBoundaryConditions(b, pos);
                                            return b; });


  /// Calculate Cmatrix here
  inputData.Cmatrix.resize(3 * (DIM - 1));
  util_funcs::CalcCmatrix(&inputData, inputData.Cmatrix);

  /// Solve
  TimerGroup<MPITimer> timers;
  std::vector<std::string> timer_labels = {"Solving"};
  std::map<std::string, int> timer_tags;
  for (int i = 0; i < timer_labels.size(); i++)
  {
    timer_tags.insert(std::pair<std::string, int>(timer_labels[i], i));
    timers.AddTimer(timer_labels[i]);
  }
  PrintStatus("before solving!");
  timers.Start(timer_tags["Solving"]);
  leSolver->solve();
  timers.Stop(timer_tags["Solving"]);
  PrintStatus("solved!");
  timers.PrintTotalTimeSeconds();

    petscVectopvtu(octDA, dTree.getTreePartFiltered(), leSolver->getCurrentSolution(), "results",
                   "le",varname, domainExtents, false, false, ndof);


  CalcStress calcStress(octDA, dTree.getTreePartFiltered(), {VecInfo(leSolver->getCurrentSolution(), LENodeData::LE_DOF, LENodeData::UX)}, domainExtents,
                        &subDomain, &inputData, ti);
#if (DIM ==3)
  static const char *stress_varname[]{"strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_xz", "strain_yz",
                                      "stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_xz", "stress_yz", "vonMises"};
#endif

#if (DIM ==2)
    static const char *stress_varname[]{"strain_xx", "strain_yy", "strain_xy",
                                      "stress_xx", "stress_yy", "stress_xy", "vonMises"};
#endif

  std::vector<double> StressVectorPerElement;
  calcStress.getElementalstress(StressVectorPerElement);
  IO::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), StressVectorPerElement.data(), "results",
                     "Stress", stress_varname,
                     domainExtents, true, false, 6 * (DIM - 1) + 1);


  Vec U_le = leSolver->getCurrentSolution();



  static const char *varname2[]{"Displacment"};



#if (DIM == 3)

  IS x_is;
  IS y_is;
  IS z_is;
  Vec U_mag = util_funcs::GetMag3D(U_le, x_is, y_is, z_is);

//  util_funcs::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), U_mag, "results",
//                             "U_mag", varname2,
//                             domainExtents, false, false, ndof);

    petscVectopvtu(octDA, dTree.getTreePartFiltered(),U_mag, "results",
                   "U_mag",varname2, domainExtents, false, false, ndof);

#endif
#if (DIM == 2)
  IS x_is, y_is;
  Vec U_mag = util_funcs::GetMag(U_le, x_is, y_is);
  //
  // util_funcs::save_timestep(octDA, treePartition, U_mag, 1, ti, subDomain, "leMag", varname2);
  VecInfo v(U_mag, 1, 0);
  // Analytic LEAnalytic(octDA, treePartition, v, analytic_sol, subDomain.domainExtents());
  //  LEAnalytic.getL2error();
#endif
#if (DIM == 2)
  Vec U_x, U_y;
  util_funcs::GetVec(U_le, x_is, y_is, U_x, U_y);
  VecInfo vx(U_x, 1, 0);
  VecInfo vy(U_y, 1, 0);
#endif
  delete leEq;
  delete leSolver;

  dendrite_finalize(octDA);
}
