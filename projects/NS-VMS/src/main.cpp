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
    std::array<DENDRITE_UINT,DIM> PERIODS;
    for (int dim = 0; dim < DIM; dim++) {
        PERIODS[dim] = ((1u << m_uiMaxDepth) * inputData.meshDef.periodicScale[dim]);
    }
    periodic::PCoord<DENDRITE_UINT , DIM>::periods(PERIODS);


    ///---------------------------------------------- Time construct ------------------------------------------------////
    std::vector<DENDRITE_REAL> dt = inputData.dt;
    std::vector<DENDRITE_REAL> totalT = inputData.totalT;
    TimeInfo ti(0.0, dt, totalT);

    ////--------------------------------------------------- Vectors---------------------------------------------------////
    Vec prev1Solution, prev2Solution;


    ///------------------------------ Command line option to restart from a checkpoint -------------------------------////
    bool resume_from_checkpoint = false;
    {
        PetscBool resume = PETSC_FALSE;
        PetscOptionsGetBool (nullptr, nullptr, "-resume_from_checkpoint", &resume, nullptr);
        resume_from_checkpoint = (resume == PETSC_TRUE);
    }
    Checkpointer checkpointer (inputData.numberOfBackups, "CheckPoint");
    /// --------------------------------------------------------------------------------------------------------------////

    DomainInfo cubeDomain, physDomain;
    cubeDomain.min = inputData.meshDef.min;
    cubeDomain.max = inputData.meshDef.max;
    physDomain.min = inputData.physDomain.min;
    physDomain.max = inputData.physDomain.max;

    DomainExtents domain(cubeDomain, physDomain);
    SubDomain subDomain(domain);

    DistTREE distTree;
    DA *octDA = nullptr;


#if(DIM==2)
    double center[DIM]{0.5,0.5};
    VOXEL::Circle circle(center,0.15);
//    subDomain.addObject(circle);
#endif


#if(DIM==3)


//    double center[DIM]{0.5,0.5,0.5};
//    VOXEL::Sphere sphere(center,0.15);


//    const auto geomDef = [&](const double * coords) {
//        // Convert the physical coordinates to image coordinates using idata.physDomainMin, idata.physDomainMax
//        // and idata.image_size
//        double x = coords[0];
//        double y = coords[1];
//        double z = coords[2];
//
//        double x0 = x + 0.5;
//        double R = sqrt(pow(x,2) + pow(y,2));
//
//        if (R < 0.8 && z > - 1.5 && z < 1.5){
//                return ibm::Partition::OUT;
//            }
//            return ibm::Partition::IN;
//    };





    GEOMETRY::STL stl(inputData.stlFileName.c_str());
    Point<DIM> temp2(0.0);
    GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, inputData.stlRetainSide);



//    subDomain.addObject(geomDef);
    subDomain.addObject(stlGeo);



#endif

    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };


    if (!resume_from_checkpoint) {
        /// Create DA
        octDA = createSubDA(distTree,functionToRetain,inputData.meshDef.baseLevel,inputData.basisFunction,0.3);
//    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domain);
        performRefinement(octDA, distTree, domain, inputData.meshDef.refineLevelBoundary, subDomain);
//        performNoChangeRefinement(octDA, distTree, domain, 5, subDomain, inputData);
    }

    else{
        /// after loading vectors are at t+dt;
        TALYFEMLIB::PrintStatus ("Loading checkpoint");
        std::vector<VecInfo> vecs;
        checkpointer.loadFromCheckpoint (octDA, distTree, functionToRetain, vecs, domain, &ti, false);

        if (octDA->isActive()) {
            prev1Solution = vecs[0].v;
            prev2Solution = vecs[1].v;
        }
        TALYFEMLIB::PrintStatus ("Checkpoint Loaded");
    }



    const auto treePart = distTree.getTreePartFiltered();

    SubDomainBoundary boundary(&subDomain, octDA, domain);

    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "mesh", "mesh",domain);


    TALYFEMLIB::PrintStatus("total No of nodes in the mesh = ", octDA->getGlobalNodeSz());



    ////--------------------------------------------------- Equation---------------------------------------------------////

    NSBoundaryConditions nsBC(&inputData, &ti, &boundary);

    auto nsEq = new TalyEquation<NSEquation, NSNodeData>(octDA, treePart, domain, NSNodeData::NS_DOF, &ti, false,
                                                         nullptr, &inputData);

    NonlinearSolver *nsSolver = setNonLinearSolver(nsEq, octDA, distTree ,NSNodeData::NS_DOF, inputData.mfree);
    inputData.solverOptionsNS.apply_to_petsc_options("-");
    SNES m_snes = nsSolver->snes();
    SNESSetFromOptions(m_snes);
    if (inputData.ifMMS) {
        nsBC.setAnalyticalFunction(AnalyticalSolution);
    }

    nsSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
        Boundary b;
        nsBC.getMomentumBoundaryCondition(b, pos);
        return b;
    });


    if (!resume_from_checkpoint) {
        if (octDA->isActive()) {
            octDA->petscCreateVector(prev1Solution, false, false, NSNodeData::NS_DOF);
            octDA->petscCreateVector(prev2Solution, false, false, NSNodeData::NS_DOF);
            setInitialCondition(octDA, inputData, prev1Solution);
        }
    }

    nsEq->setVectors({VecInfo(PLACEHOLDER_GUESS, NSNodeData::NS_DOF, NSNodeData::VEL_X),
                      VecInfo(prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1),
                      VecInfo(prev2Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2)
                     });

    VecCopy(prev1Solution, nsSolver->getCurrentSolution());
#if(DIM == 3)
    static const char *varname[]{"u", "v", "w", "p"};
#endif
#if(DIM == 2)
    static const char *varname[]{"u", "v", "p"};
#endif
    petscVectopvtu(octDA, treePart, nsSolver->getCurrentSolution(), "initial", "ns_init", varname, domain, false, false,
                   NSNodeData::NS_DOF);

    double error[NSNodeData::NUM_VARS]; // for MMS only


    ////--------------------------------------------------- Solve---------------------------------------------------////
    while (ti.getCurrentTime() < ti.getEndTime() - 1e-15) {
        ti.increment();
        ti.print();
        nsSolver->solve();
        VecCopy(prev1Solution, prev2Solution);
        VecCopy(nsSolver->getCurrentSolution(), prev1Solution);
        if (inputData.ifMMS) {
            VecInfo v(prev1Solution, NSNodeData::NS_DOF, 0);
            Analytic NSAnalytic(octDA, treePart, v, AnalyticalSolution, physDomain, ti.getCurrentTime());
            NSAnalytic.getL2error(error);
            NSAnalytic.getL2error();
            if (not(rank)) {

                for (int i = 0; i < NSNodeData::NUM_VARS; i++) {
                    std::cout << error[i] << " ";
                }
                std::cout << "\n";
            }
        }

        ////--------------------------------------------------- Writing---------------------------------------------------////
        if (ti.getTimeStepNumber () % inputData.checkpointFrequency == 0) {

            if (octDA->isActive()) {
                std::vector<VecInfo> vecs;
                vecs.emplace_back (VecInfo (prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE1));
                vecs.emplace_back (VecInfo (prev1Solution, NSNodeData::NS_DOF, NSNodeData::VEL_X_PRE2));
                checkpointer.storeCheckpoint(octDA, &distTree, vecs, domain, &ti);
            }
        }

        if (ti.getTimeStepNumber() % inputData.OutputSpan == 0) {
            char folderName[PATH_MAX], fileName[PATH_MAX];
            std::snprintf (folderName, sizeof (folderName), "%s_%05d", "results", ti.getTimeStepNumber ());
            std::snprintf (fileName, sizeof (fileName), "%s_%05d", "sol", ti.getTimeStepNumber ());
            petscVectopvtu(octDA, distTree.getTreePartFiltered(), nsSolver->getCurrentSolution(), folderName, fileName, varname,
                           domain, false, false,
                           NSNodeData::NS_DOF);
        }

    }


    VecDestroy(&prev1Solution);
    VecDestroy(&prev2Solution);
    delete nsEq;
    delete nsSolver;

    dendrite_finalize(octDA);

}