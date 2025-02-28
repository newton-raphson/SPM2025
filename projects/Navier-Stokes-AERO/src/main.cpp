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
#include <NSPost.h>
#include "NSPost.h"

using namespace PETSc;



void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel, &inputData);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);

}

void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                       DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData){
    while(true) {
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, &inputData);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }

}

template<typename T>
void CreateMesh(NSInputData &inputData, T &geoObj, DENDRITE_UINT uniLevel, DENDRITE_UINT bdLevel, DENDRITE_UINT
eleOrder, std::string geoName){
    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = inputData.cubeDomainMin[dim];
        cubeDomain.max[dim] = inputData.cubeDomainMax[dim];
        physDomain.min[dim] = inputData.physDomainMin[dim];
        physDomain.max[dim] = inputData.physDomainMax[dim];
    }
    DomainExtents domainExtents(cubeDomain,physDomain);
    SubDomain subDomain(domainExtents);

    subDomain.addObject(geoObj);

    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DistTREE distTree;
    DA * octDA = createSubDA(distTree, functionToRetain, uniLevel, eleOrder, 0.3);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performRefinement(octDA, distTree, domainExtents, bdLevel, subDomain, inputData);
    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain, inputData);
    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), geoName.c_str(), geoName.c_str(), domainExtents);
    delete octDA;
}




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
        inputData.refine_lvl_uni = static_cast<DENDRITE_UINT>(std::atoi(argv[2]));
        inputData.mfree = static_cast<bool>(std::atoi(argv[3]));
    }

    TALYFEMLIB::PrintStatus("eleOrder ", inputData.basisFunction);
    TALYFEMLIB::PrintStatus("Level ", inputData.refine_lvl_uni);
    TALYFEMLIB::PrintStatus("Mfree ", inputData.mfree);
    TALYFEMLIB::PrintStatus("DIM =  ", DIM);


    ///---------------------------------------------- Time construct ------------------------------------------------////
    std::vector<DENDRITE_REAL> dt = inputData.dt;
    std::vector<DENDRITE_REAL> totalT = inputData.totalT;
    if (argc == 4) {
        dt = std::vector<DENDRITE_REAL>(1, 0.1);
        totalT = std::vector<DENDRITE_REAL>(1, 0.2);
    }
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




//    DomainInfo physDomain;
//    physDomain.min = inputData.meshDef.min;
//    physDomain.max = inputData.meshDef.max;
//
//
//    DomainExtents domain( physDomain);
//    SubDomain subDomain(domain);
////
//
//    DistTREE distTree;
//    DA *octDA = nullptr;





#if(DIM==2)
    double center[DIM]{0.0,0.0};
    VOXEL::Circle circle(center,8.0, RetainSide::IN);


    const auto geomDef = [&](const double * coords) {
        // Convert the physical coordinates to image coordinates using idata.physDomainMin, idata.physDomainMax
        // and idata.image_size
        double x = coords[0];
        double y = coords[1];
        double x0 = x + 0.5;

        double lower_p = - 0.6 * (0.2969 * sqrt(x0) - 0.1260 * x0 - 0.3516 * pow(x0, 2) + 0.2843*pow(x0, 3) - 0.1015*pow(x0, 4));

        double upper_p = 0.6 * (0.2969 * sqrt(x0) - 0.1260 * x0 - 0.3516 * pow(x0, 2) + 0.2843*pow(x0, 3) - 0.1015*pow(x0, 4));

        double R = sqrt(pow(x,2) + pow(y,2));

        if (R < 15.0){
            if (x > - 0.50 && x < 0.50 && y > lower_p && y < upper_p){
                return ibm::Partition::IN;
            }
            return ibm::Partition::OUT;
        }
        return ibm::Partition::IN;



//        if (x > 0 && x < 1){
//            if (y > lower_p && y < upper_p){
//                return ibm::Partition::IN;
//            }
//        }
//        return ibm::Partition::OUT;

    };


//    CreateMesh(inputData, circle, inputData.refine_lvl_uni, inputData.refine_lvl_bd, inputData.basisFunction, "circle");
//    subDomain.addObject(circle);
//    subDomain.addObject(geomDef);

#endif


#if(DIM==3)
//    double center[DIM]{0.5,0.5,0.5};
//    VOXEL::Sphere sphere(center,0.15);


    GEOMETRY::STL stl(inputData.stlFileName.c_str());
    Point<DIM> temp2(0.0);
    GEOMETRY::Geometry *stlGeo = new GEOMETRY::Geometry(&stl, temp2, inputData.stlRetainSide);



//    subDomain.addObject(sphere);
    subDomain.addObject(stlGeo);

#endif






    /// Physical dimensions.
    DomainInfo cubeDomain, physDomain;
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
        cubeDomain.min[dim] = inputData.cubeDomainMin[dim];
        cubeDomain.max[dim] = inputData.cubeDomainMax[dim];
        physDomain.min[dim] = inputData.physDomainMin[dim];
        physDomain.max[dim] = inputData.physDomainMax[dim];
    }
    DomainExtents domainExtents(cubeDomain,physDomain);
    SubDomain subDomain(domainExtents);

    subDomain.addObject(geomDef);
//    subDomain.addObject(circle);


    // Retain Function
    auto functionToRetain = [&](const double *physCoords,double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };

    /// Create DA
    DistTREE distTree;
    DA * octDA = createSubDA(distTree, functionToRetain, inputData.refine_lvl_uni, inputData.basisFunction, 0.3);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performRefinement(octDA, distTree, domainExtents, inputData.refine_lvl_bd, subDomain, inputData);
    performNoChangeRefinement(octDA, distTree, domainExtents, 5, subDomain, inputData);
    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "test", "test", domainExtents);


    SubDomainBoundary boundary(&subDomain, octDA, domainExtents);

    const auto treePart = distTree.getTreePartFiltered();




//    GenerateMeshInfo meshInfo(octDA,distTree.getTreePartFiltered(),domain,&boundary);
//    meshInfo.generate();
////    meshInfo.print("GMSH/nodes.txt","GMSH/elements.txt","GMSH/boundary.txt");
//    meshInfo.printGMSH("GMSH/mesh.msh");






/// -------------------------------------------------------------------------------------------------------------------

    NSBoundaryConditions nsBC(&inputData, &ti, &boundary);

    auto nsEq = new TalyEquation<NSEquation, NSNodeData>(octDA, treePart, domainExtents, NSNodeData::NS_DOF, &ti, false,
                                                         &boundary, &inputData);
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
    petscVectopvtu(octDA, treePart, nsSolver->getCurrentSolution(), "initial", "ns_init", varname, domainExtents, false, false,
                   NSNodeData::NS_DOF);

    double error[NSNodeData::NUM_VARS];
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
                checkpointer.storeCheckpoint(octDA, &distTree, vecs, domainExtents, &ti);
            }
        }

        if (ti.getTimeStepNumber() % inputData.OutputSpan == 0) {
            char folderName[PATH_MAX], fileName[PATH_MAX];
            std::snprintf (folderName, sizeof (folderName), "%s_%05d", "results", ti.getTimeStepNumber ());
            std::snprintf (fileName, sizeof (fileName), "%s_%05d", "sol", ti.getTimeStepNumber ());
            petscVectopvtu(octDA, distTree.getTreePartFiltered(), nsSolver->getCurrentSolution(), folderName, fileName, varname,
                           domainExtents, false, false,
                           NSNodeData::NS_DOF);

            double globalForce[DIM * 2]{};
            double projectedArea[DIM];

            VecInfo vec(nsSolver->getCurrentSolution(), NSNodeData::NS_DOF, 0);
            ForceCalc forceCalc(octDA, treePart, vec, domainExtents, inputData.Re, &boundary);
            forceCalc.getTotalForce(globalForce, projectedArea);





        }



    }

    petscVectopvtu(octDA, treePart, nsSolver->getCurrentSolution(), "solution","ns", varname, domainExtents, false, false,
                   NSNodeData::NS_DOF);
    /// Save vector file (only for regression test)
    if (inputData.dump_vec) {
        petscDumpFilesforRegressionTest(octDA, nsSolver->getCurrentSolution(), "solution_vec.vec");

    }


    VecDestroy(&prev1Solution);
    VecDestroy(&prev2Solution);
    delete nsEq;
    delete nsSolver;

    dendrite_finalize(octDA);

}