//
// Created by maksbh on 5/21/20.
//


#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <NSNodeData.h>
#include <NSEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
#include <NSUtils.h>
#include <NSBoundaryConditions.h>
#include <Checkpoint/Checkpointer.h>

#ifdef IBM
#include <IMGA/Marker.h>
#endif

#include "DARefine.h"
#include "DACoarsen.h"


using namespace PETSc;


int main(int argc, char *argv[]) {

    dendrite_init(argc, argv);
    int rank = TALYFEMLIB::GetMPIRank();

    NSInputData inputData;
    if (!(inputData.ReadFromFile())) {
        if (!rank) {
            throw TALYFEMLIB::TALYException() << "Can't read the config file \n";
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (!inputData.CheckInputData()) {
        throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
    TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
    TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

    if (argc == 4) {
        inputData.basisFunction = static_cast<TALYFEMLIB::kBasisFunction>(std::atoi(argv[1]));
        inputData.meshDef.baseLevel = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
        inputData.mfree = static_cast<bool>(std::atoi(argv[3]));
    }

    TALYFEMLIB::PrintStatus("eleOrder ", inputData.basisFunction);
    TALYFEMLIB::PrintStatus("Level ", inputData.meshDef.baseLevel);
    TALYFEMLIB::PrintStatus("Mfree ", inputData.mfree);
    TALYFEMLIB::PrintStatus("DIM =  ", DIM);


    ///------------------------------------------- Things for periodic -------------------------------------------///
    std::array<DENDRITE_UINT, DIM> PERIODS;
    for (int dim = 0; dim < DIM; dim++) {
        PERIODS[dim] = ((1u << m_uiMaxDepth) * inputData.meshDef.periodicScale[dim]);
    }
    periodic::PCoord<DENDRITE_UINT, DIM>::periods(PERIODS);


    ///---------------------------------------------- Time construct ------------------------------------------------////
    std::vector<DENDRITE_REAL> dt = inputData.dt;
    std::vector<DENDRITE_REAL> totalT = inputData.totalT;
    std::vector<int> OutputSpan_V = inputData.OutputSpan_V;
    TimeInfo ti(0.0, dt, totalT, OutputSpan_V);

    ////--------------------------------------------------- Vectors---------------------------------------------------////
    Vec prev1Solution, prev2Solution;
    Vec averaged_soln; // time averaged velocity




    ///------------------------------ Command line option to restart from a checkpoint -------------------------------////
    bool resume_from_checkpoint = false;
    {
        PetscBool resume = PETSC_FALSE;
        PetscOptionsGetBool(nullptr, nullptr, "-resume_from_checkpoint", &resume, nullptr);
        resume_from_checkpoint = (resume == PETSC_TRUE);
    }
    Checkpointer checkpointer(inputData.numberOfBackups, "CheckPoint");
    /// --------------------------------------------------------------------------------------------------------------////

    DomainInfo cubeDomain, physDomain;
    cubeDomain.min = inputData.meshDef.min;
    cubeDomain.max = inputData.meshDef.max;
    physDomain.min = inputData.physDomain.min;
    physDomain.max = inputData.physDomain.max;

    DomainExtents domain(cubeDomain, physDomain);
    SubDomain subDomain(domain);

    PrintStatus("fullDADomain: ", inputData.meshDef.max[0], ", ", inputData.meshDef.max[1], ", ",
                inputData.meshDef.max[2]);
    PrintStatus("physicalDomain: ", physDomain.max[0], ", ", physDomain.max[1], ", ",
                physDomain.max[2]);

    DistTREE distTree;
    DA *octDA = nullptr;

//    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
//        return (subDomain.functionToRetain(physCoords, physSize));
//    };
    std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords,
                                                                                 double scale) {
        return subDomain.functionToRetain(octCoords, scale);
    };

#ifdef IBM
    ///------------------------------------------Loading geometries ----------------------------------------------///
    GeomInfo geomInfo;
    IMGA *imga = new IMGA(domain, (inputData.ImmersedMethodType == NSInputData::ShiftedBoundary) ? IBM_METHOD::SBM
                                                                                                 : IBM_METHOD::NITSCHE);

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
        float x = coords[0];
        float y = coords[1];
        float z = 0.0; // Initialize z with 0.0 by default
#if (DIM == 3)
        z = coords[2];
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

        if (inputData.stlRetainSide == RetainSide::IN) {
            dist = floatarr[0];
        } else {
            dist = -floatarr[0];
        }

        if (dist > 0) {
            return ibm::Partition::IN;
        }
        return ibm::Partition::OUT;
    };


    subDomain.addObject(geomDef);
#endif

    PrintStatus("we are here");
    if (!resume_from_checkpoint) {
        int baseLevel;
        if (inputData.lvl.empty() or inputData.InletBCType != NSInputData::LDC) {
            baseLevel = inputData.meshDef.baseLevel;
        } else {
            baseLevel = inputData.lvl[0];
            for (int idx = inputData.totalT.size() - 1; idx > 0; idx--) {
                if (ti.getCurrentTime() >= inputData.totalT.at(idx - 1)) {
                    baseLevel = inputData.lvl[idx];
                    break;
                }
            }
        }
        baseLevel=7;
        PrintStatus(baseLevel);
        /// Create DA
        octDA = createSubDA(distTree, functionToRetain, baseLevel, inputData.basisFunction);
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domain);
//        performNoChangeRefinement(octDA, distTree, domain, 5, subDomain, inputData);
        PrintStatus("we are here");
//        performRefinement(octDA, distTree, domain, baseLevel, subDomain, inputData, ti);
        IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "mesh", "mesh", domain);
        PrintStatus("we are here");
        exit(0);

        //////////////// Preventing losing elements for uniform mesh case /////////////////
        if (inputData.region_refine.empty()) {
            bool doRefine = inputData.carved_out_geoms_def.empty() ? true : (inputData.meshDef.baseLevel ==
                                                                             inputData.carved_out_geoms_def[0].refine_lvl);
            if (doRefine) {
                while (true) {
                    DARefine refine(octDA, distTree.getTreePartFiltered(), domain, baseLevel, false);

                    DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,
                                                      ot::RemeshPartition::SurrogateOutByIn);
                    if (newDA == nullptr) {
                        break;
                    }
                    std::swap(newDA, octDA);
                    delete newDA;

                    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domain);
                }

                while (true) {
                    DARefine refine(octDA, distTree.getTreePartFiltered(), domain, baseLevel + 1, true);

                    DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,
                                                      ot::RemeshPartition::SurrogateOutByIn);
                    if (newDA == nullptr) {
                        break;
                    }
                    std::swap(newDA, octDA);
                    delete newDA;

                    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domain);
                }

                {
                    while (true) {
                        DACoarse refine(octDA, distTree.getTreePartFiltered(), domain, baseLevel, false);
                        DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,
                                                          ot::RemeshPartition::SurrogateInByOut);
                        if (newDA == nullptr) {
                            break;
                        }
                        std::swap(newDA, octDA);
                        delete newDA;

                        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domain);
                    }
                }
            }
        }


    } else {
        /// after loading vectors are at t+dt;
        TALYFEMLIB::PrintStatus("Loading checkpoint");
        std::vector<VecInfo> vecs;
        checkpointer.loadFromCheckpoint(octDA, distTree, functionToRetain, vecs, domain, &ti, false);
//        perfromNoChangeRefinement(octDA, distTree, domain, &inputData, &subDomain, vecs, ti);

        if (octDA->isActive()) {
            prev1Solution = vecs[0].v;
            prev2Solution = vecs[1].v;
            averaged_soln = vecs[2].v; // time averaged velocity
        }
        TALYFEMLIB::PrintStatus("Checkpoint Loaded");

    }


    std::vector<PetscInt> dirichletNodes;
#ifdef IBM

    imga->initIMGAComputation(octDA, distTree.getTreePartFiltered());

    SubDomainBoundary subDomainBoundary(&subDomain, octDA, domain);

    Vec nodalFalseElement;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> elementMarkers_withFalseIntercepted;
    generateNewMarkers(octDA, distTree, domain, subDomain, &subDomainBoundary, imga,
                       elementMarkers_withFalseIntercepted, nodalFalseElement, &inputData);

    Marker *elementMarker_withoutFalseIntercept = new Marker(octDA, distTree.getTreePartFiltered(), domain, imga,
                                                             MarkerType::GAUSS_POINT);
    /// dirichlet nodes mean the nodal point that do not belong to any "active" elements (we do not infill volume GPs on that elements)
    dirichletNodes = GetIBMDirichletNodes(octDA, distTree, inputData, elementMarker_withoutFalseIntercept);
#endif

    const auto treePart = distTree.getTreePartFiltered();

    SubDomainBoundary boundary(&subDomain, octDA, domain);

    if (resume_from_checkpoint) {
        PrintStatus("Refining and Intergridding");

        int baseLevel;
        if (inputData.lvl.empty() or inputData.InletBCType != NSInputData::LDC) {
            baseLevel = inputData.meshDef.baseLevel;
        } else {
            baseLevel = inputData.lvl[0];
            for (int idx = inputData.totalT.size() - 1; idx > 0; idx--) {
                if (ti.getCurrentTime() >= inputData.totalT.at(idx - 1)) {
                    baseLevel = inputData.lvl[idx];
                    break;
                }
            }
        }

        PrintStatus("Refine lvl = ", baseLevel);
        int no_refine_return_resume = 0;
        refineAndIntergrid(octDA, distTree, domain, inputData,
                           prev1Solution, prev2Solution, averaged_soln, subDomain, false,
                           baseLevel, ti, no_refine_return_resume);
    }

    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "mesh", "mesh", domain);


    TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());

    ///------------------------------------------Pre Distance Calculation (k-d tree)-----------------------------------------///
//    std::vector<PointCloud<double>> CenterPts = SetupPointCloud(geomInfo);
//
//    std::vector<my_kd_tree_ptr> kd_trees(CenterPts.size());
//
//    for (int i = 0; i < CenterPts.size();i++) {
//        auto kd_tree = std::make_unique<my_kd_tree_t>(3 /*dim*/, CenterPts[i], 10 /* max leaf */);
//        kd_trees[i] = std::move(kd_tree);
//    }

    /// creating vectors to compue time averaged velocity
    PetscScalar *averaged_soln_array, *soln_array;

    ////--------------------------------------------------- Equation---------------------------------------------------////
    TalyMesh<NSNodeData> talyMesh(octDA->getElementOrder());
    BoundaryDistanceData boundaryDistanceData(octDA, DIM);

    NSBoundaryConditions nsBC(&inputData, &ti, &boundary);

    std::unordered_map<std::string, std::vector<double>> queryToIndexMap_vector;
#ifdef IBM
    auto nsEq = InitializeEquations(inputData, talyMesh, octDA, &subDomainBoundary, distTree,
                                    subDomain, ti, imga, boundaryDistanceData, queryToIndexMap_vector);
#else
    auto nsEq = new TalyEquation<NSEquation, NSNodeData>(octDA, treePart, domain, NSNodeData::NS_DOF, &ti, false,
                                                         nullptr, &inputData);
#endif


    NonlinearSolver *nsSolver = setNonLinearSolver(nsEq, octDA, distTree, NSNodeData::NS_DOF, inputData.mfree);
    inputData.solverOptionsNS.apply_to_petsc_options("-");
    SNES m_snes = nsSolver->snes();
    SNESSetFromOptions(m_snes);
    if (inputData.ifMMS) {
        nsBC.setAnalyticalFunction(AnalyticalSolution);
    }
#ifdef IBM
    nsSolver->setIBMDirichletNodes(dirichletNodes);
#endif

    nsSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
        Boundary b;
        nsBC.getMomentumBoundaryCondition(b, pos);
        return b;
    });

#if(DIM == 3)
    static const char *varname[]{"u", "v", "w", "p"};
#endif
#if(DIM == 2)
    static const char *varname[]{"u", "v", "p"};
#endif

    if (!resume_from_checkpoint) {
        if (octDA->isActive()) {
            octDA->petscCreateVector(prev1Solution, false, false, NSNodeData::NS_DOF);
            octDA->petscCreateVector(prev2Solution, false, false, NSNodeData::NS_DOF);
            octDA->petscCreateVector(averaged_soln, false, false, NSNodeData::NS_DOF);

            setInitialCondition(octDA, inputData, prev1Solution);
        }
    }

    const auto resetDAVecs = [&] {
        /// ----------------------------------------Sync vectors: For NS-------------------------------------------///
#ifdef IBM
        nsEq->setVectors({VecInfo(nsSolver->getCurrentSolution(), NSNodeData::NS_DOF, NSNodeData::VEL_X),
                          VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1),
                          VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2),
                          VecInfo(nodalFalseElement, 1, NSNodeData::NODE_ID)
                         });

#else
        nsEq->setVectors({VecInfo(nsSolver->getCurrentSolution(), NSNodeData::NS_DOF, NSNodeData::VEL_X),
                          VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1),
                          VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2)
                         });
#endif

        if (octDA->isActive()) {
            VecCopy(prev1Solution, nsSolver->getCurrentSolution());

            SNES m_snes = nsSolver->snes();
            SNESSetFromOptions(m_snes);

        }
    };

    ///Vectors are synced here
    resetDAVecs();


    VecCopy(prev1Solution, nsSolver->getCurrentSolution());
    petscVectopvtu(octDA, treePart, nsSolver->getCurrentSolution(), "initial", "ns_init", varname, domain, false, false,
                   NSNodeData::NS_DOF);

    const int order = inputData.RelativeOrderIntercepted;
    nsEq->setInterceptedRelativeOrder(&order);

    double error[NSNodeData::NUM_VARS]; // for MMS only

    /// --------------------------------------Setting Element Markers ------------------------------------------------------------///

    SetElementMarkerForEqs(inputData, nsEq, imga, octDA, distTree,
                           elementMarkers_withFalseIntercepted,
                           nodalFalseElement, domain, subDomain, &subDomainBoundary);

    /// -------------------------------------------------------------------------------------------------------------------------///
    IMGA *imga_TruePt_InterpPT = new IMGA(domain, IBM_METHOD::SBM);

//    OptimalSurrogateStruct optimalSurrogateData;
//    if (!imga->getGeometries().empty()) {
//        optimalSurrogateData =
//                PreProcessSBM(inputData, octDA, distTree, nodalFalseElement,
//                              domain, subDomain, &subDomainBoundary,
//                              imga, imga_TruePt_InterpPT, kd_trees, elementMarkers_withFalseIntercepted,
//                              &ti);
//    }
//
//    auto kd_tree_TrueGP = std::make_unique<my_kd_tree_t>(3 /*dim*/, optimalSurrogateData.TrueGPpos_PointCloud
//            , 10 /* max leaf */);
//    auto kd_tree_SurrogateGP = std::make_unique<my_kd_tree_t>(3 /*dim*/, optimalSurrogateData.SurrogateGPpos_PointCloud
//            , 10 /* max leaf */);

    std::unordered_map<std::string, uint32_t> queryToIndexMap;


    if (!inputData.IfScalingStudy) {
        /*Save CheckPoint at first*/
        if (octDA->isActive()) {
            std::vector<VecInfo> vecs;
            vecs.emplace_back(VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1));
            vecs.emplace_back(VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2));
            vecs.emplace_back(VecInfo(averaged_soln, NSNodeData::NS_DOF, NSNodeData::VEL_X));
            checkpointer.storeCheckpoint(octDA, &distTree, vecs, domain, &ti);
        }
    }

    ////--------------------------------------------------- Solve---------------------------------------------------////
    while (ti.getCurrentTime() < ti.getEndTime() - 1e-15) {
        PrintStatus("---------------------------------------------------");
        PrintStatus("t right now = ", ti.getCurrentTime());
        PrintStatus("Re right now = ", inputData.getReynolds(ti.getCurrentTime()));
        PrintStatus("---------------------------------------------------");

        if (inputData.CoefficientZeroForInitialGuess and ti.getTimeStepNumber() == 0) {
//            std::cout << "ti.getTimeStepNumber() = " << ti.getTimeStepNumber() << "\n";
            const double ramp_coe1 = 0.0;
            nsEq->equation()->setRampingCoefficient(&ramp_coe1);
            nsSolver->solve();
            const double ramp_coe2 = 1.0;
            nsEq->equation()->setRampingCoefficient(&ramp_coe2);
            VecCopy(nsSolver->getCurrentSolution(), prev1Solution);
            setInitialCondition(octDA, inputData, prev1Solution);
        }

        ti.increment();

        ti.print();
        nsSolver->solve();
        VecCopy(prev1Solution, prev2Solution);
        VecCopy(nsSolver->getCurrentSolution(), prev1Solution);
//        if (inputData.ifMMS) {
//            VecInfo v(prev1Solution, NSNodeData::NS_DOF, 0);
//            Analytic NSAnalytic(octDA, treePart, v, AnalyticalSolution, physDomain, ti.getCurrentTime());
//            NSAnalytic.getL2error(error);
//            NSAnalytic.getL2error();
//            if (not(rank)) {
//
//                for (int i = 0; i < NSNodeData::NUM_VARS; i++) {
//                    std::cout << error[i] << " ";
//                }
//                std::cout << "\n";
//            }
//        }

        ////------------------------------------ calculating time averaged velocity-----------------------------------------////


        VecGetArray(averaged_soln, &averaged_soln_array);
        VecGetArray(nsSolver->getCurrentSolution(), &soln_array);

        for (int i = 0; i < octDA->getLocalNodalSz() * NSNodeData::NS_DOF; i++) {

            averaged_soln_array[i] =
                    ((ti.getTimeStepNumber() - 1) * averaged_soln_array[i] + soln_array[i]) / ti.getTimeStepNumber();
        }
        VecRestoreArray(averaged_soln, &averaged_soln_array);
        VecRestoreArray(nsSolver->getCurrentSolution(), &soln_array);

        ////--------------------------------------------------- Writing---------------------------------------------------////

        if (ti.getTimeStepNumber() % inputData.checkpointFrequency == 0) {

            if (octDA->isActive()) {
                std::vector<VecInfo> vecs;
                vecs.emplace_back(VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1));
                vecs.emplace_back(VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2));
                vecs.emplace_back(VecInfo(averaged_soln, NSNodeData::NS_DOF, NSNodeData::VEL_X));
                checkpointer.storeCheckpoint(octDA, &distTree, vecs, domain, &ti);
            }
        }

        /// Output Cd after the time step is larger than CdStartTimeStep
//        if (ti.getTimeStepNumber() > inputData.CdStartTimeStep){
//            // close to PostProcessNu_Cd_BoundaryError_OptSurrogate
//            PostProcess_Cd_BoundaryError_OptSurrogate(inputData, octDA, distTree,
//                                                      nsSolver,
//                                                      domain, ti,
//                                                      imga, imga_TruePt_InterpPT,
//                                                      kd_tree_TrueGP.get(),
//                                                      kd_tree_SurrogateGP.get(),
//                                                      kd_trees,
//                                                      optimalSurrogateData,
//                                                      queryToIndexMap);
//        }

        if (ti.getTimeStepNumber() % ti.getCurrentOutputSpan() == 0) {
            PrintStatus("Output file");
            char folderName[PATH_MAX], fileName[PATH_MAX];
            std::snprintf(folderName, sizeof(folderName), "%s_%05d", "results", ti.getTimeStepNumber());
            std::snprintf(fileName, sizeof(fileName), "%s_%05d", "sol", ti.getTimeStepNumber());
            petscVectopvtu(octDA, distTree.getTreePartFiltered(), nsSolver->getCurrentSolution(), folderName, fileName,
                           varname,
                           domain, false, false,
                           NSNodeData::NS_DOF);

            std::snprintf(folderName, sizeof(folderName), "%s_%05d", "AvgResults", ti.getTimeStepNumber());
            std::snprintf(fileName, sizeof(fileName), "%s_%05d", "AvgSol", ti.getTimeStepNumber());
            petscVectopvtu(octDA, distTree.getTreePartFiltered(), averaged_soln, folderName, fileName, varname,
                           domain, false, false,
                           NSNodeData::NS_DOF);

        }

        int no_refine_return = 0;
        refineAndIntergrid(octDA, distTree, domain, inputData,
                           prev1Solution, prev2Solution, averaged_soln, subDomain, false,
                           inputData.meshDef.baseLevel, ti, no_refine_return);

        if (no_refine_return > 0) {
            elementMarkers_withFalseIntercepted.clear();
            VecDestroy(&nodalFalseElement);

            generateNewMarkers(octDA, distTree, domain, subDomain, &subDomainBoundary, imga,
                               elementMarkers_withFalseIntercepted, nodalFalseElement, &inputData);
            delete elementMarker_withoutFalseIntercept;
            elementMarker_withoutFalseIntercept = new Marker(octDA, distTree.getTreePartFiltered(), domain, imga,
                                                             MarkerType::GAUSS_POINT);
            dirichletNodes = GetIBMDirichletNodes(octDA, distTree, inputData, elementMarker_withoutFalseIntercept);

            delete nsEq;

            // update kd_tress again because we reset the nsEq
            nsEq = InitializeEquations(inputData, talyMesh, octDA, &subDomainBoundary, distTree,
                                       subDomain, ti, imga, boundaryDistanceData, queryToIndexMap_vector);

            SetElementMarkerForEqs(inputData, nsEq, imga, octDA, distTree,
                                   elementMarkers_withFalseIntercepted,
                                   nodalFalseElement, domain, subDomain, &subDomainBoundary);


            /// (Cheng-Hau) 2024/5/7: Previous dug is that we do not add below to distribute the new true GPs to new cell/elements.
            imga->initIMGAComputation(octDA, distTree.getTreePartFiltered());

            if (!imga->getGeometries().empty()) {
                delete imga_TruePt_InterpPT;
                imga_TruePt_InterpPT = new IMGA(domain, IBM_METHOD::SBM);
//                optimalSurrogateData.clear();
//                optimalSurrogateData =
//                        PreProcessSBM(inputData, octDA, distTree, nodalFalseElement,
//                                      domain, subDomain, &subDomainBoundary,
//                                      imga, imga_TruePt_InterpPT, kd_trees, elementMarkers_withFalseIntercepted,
//                                      &ti);
            }
//
//            kd_tree_SurrogateGP.reset();  // Deletes the object managed by kd_tree_SurrogateGP and sets the pointer to nullptr
//
//            kd_tree_SurrogateGP = std::make_unique<my_kd_tree_t>(3 /*dim*/, optimalSurrogateData.SurrogateGPpos_PointCloud
//                    , 10 /* max leaf */);


                updateSolver(&talyMesh, octDA, distTree, nsEq, domain, nsSolver,
                             &ti, &boundary, NSNodeData::NS_DOF, &inputData);

                resetDAVecs();


#ifdef IBM
                nsSolver->setIBMDirichletNodes(dirichletNodes);
#endif
                queryToIndexMap.clear();
                queryToIndexMap_vector.clear();
            }
        }

        /*Save CheckPoint at end*/
        if (!inputData.IfScalingStudy) {
            if (octDA->isActive()) {
                std::vector<VecInfo> vecs;
                vecs.emplace_back(VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1));
                vecs.emplace_back(VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2));
                vecs.emplace_back(VecInfo(averaged_soln, NSNodeData::NS_DOF, NSNodeData::VEL_X));
                checkpointer.storeCheckpoint(octDA, &distTree, vecs, domain, &ti);
            }
        }


        VecDestroy(&prev1Solution);
        VecDestroy(&prev2Solution);
        VecDestroy(&averaged_soln);
        delete nsEq;
        delete nsSolver;

        dendrite_finalize(octDA);


}