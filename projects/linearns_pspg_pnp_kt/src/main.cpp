
#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <NSPNPNodeData.hpp>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
#include <Checkpoint/Checkpointer.h>
//#include <PostProcessing/postProcessing.h>
#include <NSUtils.h>
#include <NSBoundaryConditions.h>
//#include <PressurePoissonEquation.h>
//#include <VelocityUpdateEquation.h>
#include <NSEquation.h>
#include <PNPEquation.h>
//#include <PNPEquation1.hpp>

using namespace PETSc;

int main(int argc, char *argv[]) {

	dendrite_init (argc, argv);
	int rank = TALYFEMLIB::GetMPIRank ();
//	if (argc < 3) {
//		if (not(rank))
//			std::cout << "Usage: " << argv[0]
//			          << " eleOrder level mfree"
//			          << std::endl;
//		return 0;
//	}

	NSPNPInputData inputData;

	if (!(inputData.ReadFromFile ())) {
		if (!rank) {
			throw TALYFEMLIB::TALYException () << "Can't read the config file \n";
		}
		MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
	}
	inputData.solverOptionsMomentum.apply_to_petsc_options("-momentum_");
	inputData.solverOptionsPNP.apply_to_petsc_options ("-pnp_");

//	NSParams nsParams (&inputData);

	TALYFEMLIB::PrintInfo ("Total number of processor = ", TALYFEMLIB::GetMPISize ());
	TALYFEMLIB::PrintInfo ("size of DendroInt ", sizeof (DendroIntL));
	TALYFEMLIB::PrintInfo ("size of PetscInt ", sizeof (PetscInt));


	TALYFEMLIB::PrintStatus ("eleOrder ", inputData.elemOrder);
	TALYFEMLIB::PrintStatus ("Level ", inputData.meshDef.baseLevel);
	TALYFEMLIB::PrintStatus ("Mfree ", inputData.ifMatrixFree);
	TALYFEMLIB::PrintStatus ("DIM =  ", DIM);






    ///================================ Things for periodic ================================================///
    std::array<DENDRITE_UINT,DIM> PERIODS;
    for (int dim = 0; dim < DIM; dim++) {
        PERIODS[dim] = ((1u << m_uiMaxDepth) * inputData.meshDef.periodicScale[dim]);
    }
    periodic::PCoord<DENDRITE_UINT , DIM>::periods(PERIODS);
    ///====================================================================================================///

    ///------------------------------ Command line option to restart from a checkpoint -------------------------------////
    bool resume_from_checkpoint = false;
    {
        PetscBool resume = PETSC_FALSE;
        PetscOptionsGetBool (nullptr, nullptr, "-resume_from_checkpoint", &resume, nullptr);
        resume_from_checkpoint = (resume == PETSC_TRUE);
    }
    Checkpointer checkpointer (inputData.numberOfBackups, "CheckPoint");
    /// --------------------------------------------------------------------------------------------------------------////
    ///-------------------------------------------------- Time construct-----------------------------------------------///
    std::vector<DENDRITE_REAL> dt(1, inputData.timeStep);
    std::vector<DENDRITE_REAL> totalT(1, inputData.TotalTime);
    TimeInfo ti(0.0, dt, totalT);
    /// --------------------------------------------------------------------------------------------------------------////

    ////--------------------------------------------------- Vectors---------------------------------------------------////
    /// vectors for NS
    Vec prev1VelocityPressureSolution, prev2VelocityPressureSolution, prev3VelocityPressureSolution;
    /// vectors for PNP
    Vec prev1PNPSolution, prev2PNPSolution, prev3PNPSolution;

    ///------------------------------------------Creation/loading of mesh----------------------------------------------///
    DomainInfo cubeDomain, physDomain;
    cubeDomain.min = inputData.meshDef.min;
    cubeDomain.max = inputData.meshDef.max;
    physDomain.min = inputData.physDomain.min;
    physDomain.max = inputData.physDomain.max;

    DomainExtents domain(cubeDomain, physDomain);
    SubDomain subDomain(domain);

    DistTREE distTree;
    DA *octDA = nullptr;

    /// carving subda
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };



    if (!resume_from_checkpoint) {
        octDA = createSubDA(distTree,functionToRetain,inputData.meshDef.baseLevel,inputData.elemOrder,0.3);


        /// perform initial refinement
        //performRefinement(octDA, distTree, domain, inputData.meshDef.refineLevelBoundary, subDomain, inputData);

    } else {/// Case with checkpointing
        /// after loading vectors are at t+dt;
        TALYFEMLIB::PrintStatus ("Loading checkpoint");
        std::vector<VecInfo> vecs;

        checkpointer.loadFromCheckpoint (octDA, distTree, functionToRetain, vecs, domain, &ti, false);

        if (octDA->isActive()) {
            prev1VelocityPressureSolution = vecs[0].v;
            prev2VelocityPressureSolution = vecs[1].v;
            prev3VelocityPressureSolution = vecs[2].v;

            prev1PNPSolution = vecs[3].v;
            prev2PNPSolution = vecs[4].v;
            prev3PNPSolution = vecs[5].v;

        }
        TALYFEMLIB::PrintStatus ("Checkpoint Loaded");
    }
    const auto treePart = distTree.getTreePartFiltered();
    TALYFEMLIB::PrintStatus ("total No of nodes in the mesh = ", octDA->getGlobalNodeSz ());
    //IO::writeBoundaryElements (octDA, treePart, "boundary", "subDA", domain);
    /// --------------------------------------------------------------------------------------------------------------////


	/// Boundary condition
	NSBoundaryConditions nsBC(&inputData, &ti);
	NSBoundaryConditions pnpBC(&inputData, &ti);


	/// Equation and solver setup
    auto nsEq = new TalyEquation<NSEquation, NSPNPNodeData>(octDA, treePart, domain, NSPNPNodeData::NS_DOF, &ti, false,
                                                         nullptr, &inputData);





	LinearSolver *nsSolver = setLinearSolver(nsEq, octDA, distTree, NSPNPNodeData::NS_DOF, inputData.ifMatrixFree);
	KSP momentumSolver_ksp = nsSolver->ksp();
	KSPSetOptionsPrefix(momentumSolver_ksp, "momentum_");
	KSPSetFromOptions(momentumSolver_ksp);

	/// Setting up the nonlinear PNP solver
	auto pnpEquation = new TalyEquation<PNPEquation, NSPNPNodeData> (octDA, treePart, domain, NSPNPNodeData::PNP_DOF, &ti, false,
                                                                      nullptr, &inputData);



       NonlinearSolver *pnpSolver = setNonLinearSolver(pnpEquation, octDA, distTree, NSPNPNodeData::PNP_DOF, inputData.ifMatrixFree);
       SNES pnpEquation_snes = pnpSolver->snes();
       SNESSetOptionsPrefix(pnpEquation_snes, "pnp_");
       SNESSetFromOptions(pnpEquation_snes);

       /// Analytical solution to BC for MMS
       if (inputData.ifMMS) {
           nsBC.setAnalyticalFunction (AnalyticalSolutionNS);
           pnpBC.setAnalyticalFunction (AnalyticalSolutionPNP);
       }

       /// Boundary condition setup for NS
       nsSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
               Boundary b;
               nsBC.getMomentumBoundaryCondition(b, pos);
               return b;
       });

       /// Boundary condition setup for PNP
       pnpSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV &pos, int nodeID) -> Boundary {
               Boundary b;
               pnpBC.getPNPBoundaryCondition (b, pos);
               return b;
       });


    if (!resume_from_checkpoint) {
        if (octDA->isActive()) {

            octDA->petscCreateVector(prev1VelocityPressureSolution, false, false, NSPNPNodeData::NS_DOF);
            octDA->petscCreateVector(prev2VelocityPressureSolution, false, false, NSPNPNodeData::NS_DOF);
            octDA->petscCreateVector(prev3VelocityPressureSolution, false, false, NSPNPNodeData::NS_DOF);
            octDA->petscCreateVector(prev1PNPSolution, false, false, NSPNPNodeData::PNP_DOF);
            octDA->petscCreateVector(prev2PNPSolution, false, false, NSPNPNodeData::PNP_DOF);
            octDA->petscCreateVector(prev3PNPSolution, false, false, NSPNPNodeData::PNP_DOF);

            /// Set initial conditions of all vectors
            VecSet(prev2VelocityPressureSolution, 0.0);
            VecSet(prev3VelocityPressureSolution, 0.0);
            VecSet(prev2PNPSolution, 0.0);
            VecSet(prev3PNPSolution, 0.0);

            setInitialConditionVelocity(octDA, inputData, prev1VelocityPressureSolution);
            //todo: Only implemented for MMS
            setInitialConditionPNP(octDA, inputData, prev1PNPSolution);
        }
    }

       /// Sync vectors: For pressure poisson and velocity update matrix does not need vector information
    const auto resetDAVecs = [&] {
        nsEq->setVectors({
                                 VecInfo(pnpSolver->getCurrentSolution(), NSPNPNodeData::PNP_DOF,
                                         NSPNPNodeData::POTENTIAL),
                                 VecInfo(prev1VelocityPressureSolution, NSPNPNodeData::NS_DOF,
                                         NSPNPNodeData::VEL_X_PRE1),
                                 VecInfo(prev2VelocityPressureSolution, NSPNPNodeData::NS_DOF,
                                         NSPNPNodeData::VEL_X_PRE2),
                                 VecInfo(prev3VelocityPressureSolution, NSPNPNodeData::NS_DOF,
                                         NSPNPNodeData::VEL_X_PRE3)
                         }, SYNC_TYPE::ALL);
        pnpEquation->setVectors({

                                        VecInfo(PLACEHOLDER_GUESS, NSPNPNodeData::PNP_DOF, NSPNPNodeData::POTENTIAL),

                                        VecInfo(prev1PNPSolution, NSPNPNodeData::PNP_DOF,
                                                NSPNPNodeData::POTENTIAL_PRE1),
                                        VecInfo(prev2PNPSolution, NSPNPNodeData::PNP_DOF,
                                                NSPNPNodeData::POTENTIAL_PRE2),
                                        VecInfo(prev3PNPSolution, NSPNPNodeData::PNP_DOF, NSPNPNodeData::POTENTIAL_PRE3),
                                        VecInfo(nsSolver->getCurrentSolution(), NSPNPNodeData::NS_DOF,
                                                NSPNPNodeData::VEL_X)

                                }, SYNC_TYPE::ALL);


        if (octDA->isActive()) {

            VecCopy(prev1VelocityPressureSolution, nsSolver->getCurrentSolution());
            VecCopy(prev1PNPSolution, pnpSolver->getCurrentSolution());
        }
    };
    ///Vectors are synched here
    resetDAVecs ();
   #if(DIM == 3)
       static const char *momentum_varname[]{"u", "v", "w", "p"};
       static const char *pnp_varname[]{"pot, c1, c2"};
   #endif
   #if(DIM == 2)
       static const char *momentum_varname[]{"u", "v", "p"};
       static const char *pnp_varname[]{"pot", "c1", "c2"};
   #endif
       //petscVectopvtu (octDA, treePart, nsSolver->getCurrentSolution (), "init", "momentum_init", momentum_varname, domain, false, false, NSPNPNodeData::NS_DOF);
       //petscVectopvtu (octDA, treePart, pnpSolver->getCurrentSolution (), "init","pnp_init", pnp_varname, domain, false,
        //               false, NSPNPNodeData::PNP_DOF);

       /// Setting first order for first step
       //nsParams.updateNSOrder(1);
       std::string fname = "Error.dat";
       DENDRITE_REAL nsError[NSPNPNodeData::NS_DOF],pnpError[NSPNPNodeData::PNP_DOF];
       std::ofstream fout(fname);
       fout.close();

       int config_nsOrder = inputData.nsOrder;
       int advectionExtrapolationOrder = inputData.advectionExtrapolationOrder;



//    while (ti.getCurrentTime() < ti.getEndTime ()) {
//        pnpSolver->solve();
//        VecCopy(prev2PNPSolution, prev3PNPSolution);
//        VecCopy(prev1PNPSolution, prev2PNPSolution);
//        VecCopy(pnpSolver->getCurrentSolution(), prev1PNPSolution);
//        ti.increment();
//        ti.print();
//
//    }printf("pnp solving finished");
//
//    ti.setCurrentTime(0.0);
//    ti.setTimeStepNumber(0);

       /// Time loop
       while (ti.getCurrentTime() < ti.getEndTime ()) {
           inputData.nsOrder = config_nsOrder;
           inputData.advectionExtrapolationOrder = advectionExtrapolationOrder;
           if (inputData.timescheme == NSPNPInputData::BDF && ti.getTimeStepNumber() == 0) {
               inputData.nsOrder = 1;
               inputData.advectionExtrapolationOrder = 1;
           }

           if ( ti.getCurrentTime() < inputData.restartTimeEof) {
               inputData.ifPurePNP = true;
           } else if ( ti.getCurrentTime() > inputData.restartTimeEof) {
               inputData.ifPurePNP = false;
           }





           ///------------------------------------------------- Solve-------------------------------------------------------///
           //VecSet(nsSolver->getCurrentSolution(),0.0);


           PrintLogStream(std::cout, "",
                          "/=============================== Solving momentum equations =================================/");
           nsSolver->solve();
           std::string tnumV = "momentum_";
           //std::string tnumP = "pressure_";
           //tnumV += std::to_string(ti.getTimeStepNumber());
           //tnumP += std::to_string(ti.getTimeStepNumber());
           //petscVectopvtu(octDA, nsSolver->getCurrentSolution(), tnumV.c_str(), momentum_varname, inputData
           //.physDomain, false,false, DIM + 1);
           //pressurePoissonSolver->solve();
           PrintLogStream(std::cout, "",
                          "/============================= Done solving momentum equations ==============================/");

           PrintLogStream(std::cout, "",
                          "/=========================== Solving Poisson-Nernst-Planck equations ========================/");
           pnpSolver->solve();

           PrintLogStream(std::cout, "",
                          "/========================= Done solving Poisson-Nernst-Planck equations =====================/");

           //petscVectopvtu (octDA, pressurePoissonSolver->getCurrentSolution(), tnumP.c_str(), pressure_varname, inputData
           //.physDomain, false,false,1);


           VecCopy(prev2VelocityPressureSolution, prev3VelocityPressureSolution);
           VecCopy(prev1VelocityPressureSolution, prev2VelocityPressureSolution);
           VecCopy(nsSolver->getCurrentSolution(), prev1VelocityPressureSolution);

           VecCopy(prev2PNPSolution, prev3PNPSolution);
           VecCopy(prev1PNPSolution, prev2PNPSolution);
           VecCopy(pnpSolver->getCurrentSolution(), prev1PNPSolution);

           //nsParams.updateNSOrder(inputData.nsOrder);
           ti.increment();
           ti.print();

           if (inputData.ifMMS) {
               {
                   VecInfo v (prev1VelocityPressureSolution, NSPNPNodeData::NS_DOF, NSPNPNodeData::VEL_X);
                   Analytic NSAnalytic(octDA, treePart, v, AnalyticalSolutionNS, domain, ti.getCurrentTime ());
                   TALYFEMLIB::PrintStatus ("Velocity and pressure error = ");
                   NSAnalytic.getL2error();
                   NSAnalytic.getL2error(nsError);
               }
               {
                   VecInfo v (prev1PNPSolution, NSPNPNodeData::PNP_DOF, NSPNPNodeData::VEL_X);
                   Analytic PNPAnalytic(octDA, treePart, v, AnalyticalSolutionPNP, domain, ti.getCurrentTime ());
                   TALYFEMLIB::PrintStatus ("PNP error = ");
                   PNPAnalytic.getL2error();
                   PNPAnalytic.getL2error(pnpError);
               }
               fout.open(fname,std::ios::app);
               if(not(rank)){
                   fout << ti.getCurrentTime () << " " << nsError[0] << " " << nsError[1] << " " << " " << nsError[2] << " " <<
                        pnpError[0] << " " << pnpError[1] << " " << pnpError[2] << "\n";
               }
               fout.close();
           }
           /// ----------------------------------------------file writing---------------------------------------------------///
           if (ti.getTimeStepNumber () % inputData.OutputSpan == 0) {
               /// writing solution for t+dt
               char folderName[PATH_MAX], fileNamePNP[PATH_MAX], fileNameNS[PATH_MAX];

               std::snprintf (folderName, sizeof (folderName), "%s_%05d", "results", ti.getTimeStepNumber ());
               std::snprintf (fileNameNS, sizeof (fileNameNS), "%s_%05d", "ns", ti.getTimeStepNumber ());
               std::snprintf (fileNamePNP, sizeof (fileNamePNP), "%s_%05d", "pnp", ti.getTimeStepNumber ());


               petscVectopvtu(octDA, treePart, prev1VelocityPressureSolution, folderName,fileNameNS, momentum_varname, domain,
                              false,
                              false, NSPNPNodeData::NS_DOF);
               petscVectopvtu(octDA, treePart, prev1PNPSolution, folderName, fileNamePNP, pnp_varname, domain, false, false,
                              NSPNPNodeData::PNP_DOF);

           }

       }


       delete nsEq;
       delete nsSolver;

       delete pnpEquation;
       delete pnpSolver;

       VecDestroy(&prev1PNPSolution);
       VecDestroy(&prev2PNPSolution);
       VecDestroy(&prev3PNPSolution);
       VecDestroy(&prev1VelocityPressureSolution);
       VecDestroy(&prev2VelocityPressureSolution);
       VecDestroy(&prev3VelocityPressureSolution);

	dendrite_finalize (octDA);

}