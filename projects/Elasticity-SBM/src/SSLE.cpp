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


#pragma mark OptSug
//#include "CheckSurface.h"
#include "DACoarse.h"
#include "DARefine.h"



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


    ///------------------------------------------Specific thing to do with some cases ----------------------------------------------///
    if (inputData.bccaseType == CSV_FORCE) {

        // create k-d tree for traction points reading from CSV
         util_funcs::readTXT(inputData.CsvForce.filename,inputData.traction_position_,inputData.traction_vector_);
         util_funcs::removeDuplicates(inputData.traction_position_,inputData.traction_vector_);
        inputData.kdTree_ = new my_kd_tree_t(3 /*dim*/, inputData.traction_position_, 10/* max leaf */);

    }

  ///------------------------------------------Creation/loading of mesh----------------------------------------------///
  DomainExtents domainExtents(inputData.mesh_def.fullDADomain, inputData.mesh_def.physDomain);
  DA *octDA = nullptr;
  DistTREE dTree;
  SubDomain subDomain(domainExtents, resume_from_checkpoint);




    if(inputData.ibm_geom_def[0].type != IBMGeomDef::Type::DeepTrace){
        throw TALYFEMLIB::TALYException() << "The model is not a deeptrace model";
    }

    const auto &deepTrace = inputData.ibm_geom_def[0];

    //    USING DEEP TRACE FOR GETTING THE MESH ###
    // Call the function you want to measure

    // Initialize the ONNX Runtime environment
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "test");

    // Set up session options
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);

    // Create a session with the loaded model
    Ort::Session session(env, inputData.model_path, session_options);

    // Define names of the input and output tensors
    const char *input_name = "input"; // Replace with the actual input name obtained from the model
    const char *output_name = "output"; // Replace with the actual output name obtained from the model



    const auto geomDef = [&](const double *coords) {


        float x = coords[0]*deepTrace.scale;
        float y = coords[1]*deepTrace.scale;
        float z = 0.0; // Initialize z with 0.0 by default
#if (DIM == 3)
        z = coords[2]*deepTrace.scale;
#endif

        // Prepare input tensor
        std::vector<float> input_tensor_values = {x, y, z};
        std::vector<int64_t> input_tensor_shape = {1, 3}; // Assuming the model expects a [1,3] shape tensor

        // Create memory info to describe the allocation info
        Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

        // Create input tensor
        Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
                memory_info, input_tensor_values.data(), input_tensor_values.size(), input_tensor_shape.data(),
                input_tensor_shape.size());

        // Perform the inference
        auto output_tensors = session.Run(Ort::RunOptions{nullptr}, &input_name, &input_tensor, 1, &output_name, 1);

        float *floatarr = output_tensors.front().GetTensorMutableData<float>();

        float dist = 0;

        dist = floatarr[0]/deepTrace.scale;


        if (dist > 0) {
            return ibm::Partition::IN;
        }
        return ibm::Partition::OUT;
    };


    subDomain.addObject(geomDef);


    /// Pre Distance Calculation
    PointCloud<double> CenterPts;

    std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords, double scale) {
        return subDomain.functionToRetain(octCoords, scale);
    };

    //// this is just here to maintain the consistency of the code

    IMGA *imga = new IMGA(domainExtents,  IBM_METHOD::SBM);



    my_kd_tree_t kd_tree(3 /*dim*/, CenterPts, {10 /* max leaf */});


    octDA = createSubDA(dTree, functionToRetain, levelBase, eleOrder);
    subDomain.finalize(octDA, dTree.getTreePartFiltered(), domainExtents);

    int no_refine = util_funcs::performRefinementSubDA(octDA, domainExtents, dTree, inputData, &subDomain);

    PrintStatus("Number of refinement = ", no_refine);

    TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());

    /// --------------------------------------------------------------------------------------------------------------////
    const auto &treePartition = dTree.getTreePartFiltered();
    IO::writeBoundaryElements(octDA, dTree.getTreePartFiltered(), "boundary", "subDA", domainExtents);

    imga->initIMGAComputation(octDA, treePartition);
    Marker elementMarker(octDA, treePartition, domainExtents, imga, MarkerType::GAUSS_POINT);

    Marker elementMarker_nodeBase(octDA, treePartition, domainExtents, imga, MarkerType::ELEMENT_NODES);
    elementMarker_nodeBase.printMarker();


    SubDomainBoundary boundary(&subDomain, octDA, domainExtents);

    LEBCSetup LEBC(&boundary, &inputData);


    Vec nodalFalseElement;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> elementMarkers_withFalseIntercepted;
    util_funcs::generateNewMarkers(octDA, dTree, domainExtents, subDomain, &boundary, imga, elementMarkers_withFalseIntercepted, nodalFalseElement, &inputData);


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
                                                         &ti, true, subDomainBoundary, &inputData, IBM_METHOD::SBM, imga);

  leEq->equation()->setSbmCalc(&kd_tree);

  std::vector<PetscInt> dirichletNodes;
    LinearSolver *leSolver = setLinearSolver(leEq, octDA,dTree,ndof, mfree); // for IBM

  leSolver->setIBMDirichletNodes(dirichletNodes); // it does not set anything (dirichletNodes do not have anything), but it is here to "prevent segmentation fault" !


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

  leEq->setVectors({VecInfo(nodalFalseElement, 1, LENodeData::NODE_ID)});
  leEq->assignIBMConstructs(imga, elementMarkers_withFalseIntercepted.data(), LENodeData::NODE_ID);

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
//  util_funcs::writeVecTopVtu(octDA, dTree.getTreePartFiltered(), leSolver->getCurrentSolution(), , varname,
//                             subDomain.domainExtents(), false, false, );

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

  delete leEq;
  delete leSolver;

  dendrite_finalize(octDA);
}
