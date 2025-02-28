//
// Created by maksbh on 2/11/23.
//

#include "PETSc/PetscUtils.h"
#include "CHNodeData.h"
#include "CHEquation.h"
#include "CHInputData.h"
#include <DendriteUtils.h>
#include "CHInitialCondition.h"
#include "PETSc/IO/petscVTU.h"
#include "CHPostProcessing.h"
#include "CHUtils.h"
#include "MassComputation.h"
#include "MEquation.h"
using namespace PETSc;
int main(int argc, char *argv[]) {

    dendrite_init(argc, argv);
    int rank = TALYFEMLIB::GetMPIRank();
    int npes = TALYFEMLIB::GetMPISize();

    TALYFEMLIB::PrintStatus("DIM = ", DIM);

    CHInputData idata;
    if (!idata.ReadFromFile()) {  /// read from file named "config.txt"
        throw std::runtime_error("[ERR] Error reading input data, check the config file!");
    }
    if (!idata.CheckInputData()) {
        throw std::runtime_error("[ERR] Problem with input data, check the config file!");
    }
    idata.solverOptionsCH.apply_to_petsc_options ("-ch_");


    std::vector<double> dt(1,idata.dt);
    std::vector<double> totalT(1,idata.totalT);
    TimeInfo ti(0.0, dt,totalT);

    DistTREE distTree;

    DomainExtents domainExtents(idata.mesh_def.fullDADomain,idata.mesh_def.physDomain);
    SubDomain subDomain(domainExtents);

    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };
    DA *octDA = createSubDA(distTree,functionToRetain,idata.mesh_def.refine_lvl_interface, idata.elementOrder);

    MPI_Barrier(MPI_COMM_WORLD);
    auto chEquation = new TalyEquation<CHEquation, CHNodeData>(octDA,distTree.getTreePartFiltered(),domainExtents,CHNodeData::CH_DOF,&ti,false,
                                                               nullptr,&idata);
    auto chSolver = setNonLinearSolver(chEquation,octDA,distTree,CHNodeData::CH_DOF,false,false,false);

    Vec prevSolution;
    octDA->petscCreateVector(prevSolution,false,false,CHNodeData::CH_DOF);
    if(!rank) {
        std::ofstream fout("EnergyAndMass.txt");
        fout << "Time,MassAfterSolve,EnergyAferSolve,InterfacialEnergyAfterSolve,"
                "MassAfterRefine,EnergyAferRefine,InterfacialEnergyAfterRefine,"
                "MassAfterCoarsen,EnergyAferCoarsen,InterfacialEnergyAfterCoarsen\n";
        fout.close();
    }


    const auto resetDAVecs = [&] {
        chEquation->setVectors (
                {
                    VecInfo (PLACEHOLDER_GUESS, CHNodeData::CH_DOF, CHNodeData::PHI_DOF),
                    VecInfo(prevSolution, CHNodeData::CH_DOF, CHNodeData::CH_DOF, PLACEHOLDER_NONE)
                    }, SYNC_TYPE::ALL);

        if (octDA->isActive ()) {
            VecCopy (prevSolution, chSolver->getCurrentSolution ());
            SNES chEquation_snes = chSolver->snes ();
            SNESSetOptionsPrefix (chEquation_snes, "ch_");
            SNESSetFromOptions (chEquation_snes);
        }
    };


    chSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
        Boundary b;
        return b;
    });

    CHInitialCondition chInitialCondition(CHInitType::RANDOM);
    chInitialCondition.setInitialCondition(octDA,prevSolution);
    resetDAVecs();
    while (ti.getCurrentTime() < ti.getEndTime()){
        ti.print();

        chSolver->solve();
        auto snes = chSolver->snes();
        SNESConvergedReason reason;
        SNESGetConvergedReason(snes,&reason);
        if(reason < 0){
            throw std::runtime_error("Non linear solver diverged");
        }
        VecCopy(chSolver->getCurrentSolution(),prevSolution);



        if (ti.getTimeStepNumber () % idata.OutputInterval == 0) {
            /// writing solution for t+dt
            char folderName[PATH_MAX], fileNameCH[PATH_MAX];
            std::snprintf (folderName, sizeof (folderName), "%s_%05d", "results", ti.getTimeStepNumber ());
            std::snprintf (fileNameCH, sizeof (fileNameCH), "%s_%05d", "ch", ti.getTimeStepNumber ());
            petscVectopvtu (octDA, distTree.getTreePartFiltered(), chSolver->getCurrentSolution(),
                            folderName, fileNameCH,ch_varname, subDomain.domainExtents (),
                            false, false, DIM);

        }
        VecInfo vec(prevSolution,CHNodeData::CH_DOF,0);
        double globalQuantityAfterSolve[3];
        {
            CHPostProcessing chPostProcessing(octDA, distTree, domainExtents, vec);
            chPostProcessing.getGlobalQuantity(globalQuantityAfterSolve, octDA->getCommActive());
        }
        std::vector<double> dummy;
        bool didRefine = performRefinement<RefinementStage::REFINE>(octDA,distTree,vec,domainExtents,&idata,dummy,dummy,dummy,false);
        double globalQuantityAfterRefine[3];
        {
            CHPostProcessing chPostProcessing(octDA, distTree, domainExtents, vec);
            chPostProcessing.getGlobalQuantity(globalQuantityAfterRefine, octDA->getCommActive());
        }

        bool didCoarsen;
        if(idata.useMassConservingInterpolation) {
            std::vector<double> elementalMass,elementalMassPerChild,getCellCoarsened;
            MassCalculator calculator(octDA,distTree.getTreePartFiltered(),vec,domainExtents,elementalMass,elementalMassPerChild);
//            didCoarsen = performCoarsening<RefinementStage::COARSEN>(octDA,distTree, vec, domainExtents, &idata,elementalMass,
//                                                                    true);
            didCoarsen = performRefinement<RefinementStage::COARSEN>(octDA, distTree, vec, domainExtents, &idata,elementalMass,elementalMassPerChild,getCellCoarsened,
                                                                     true);


            if(didCoarsen) {
                auto massEquation = new TalyEquation<MEquation, CHNodeData>(octDA, distTree.getTreePartFiltered(),
                                                                            domainExtents, CHNodeData::CH_DOF,
                                                                            nullptr, false, nullptr);

                massEquation->equation()->assign(&elementalMass,&elementalMassPerChild,&getCellCoarsened);

                auto massSolver = setLinearSolver(massEquation, octDA, distTree, CHNodeData::CH_DOF, false, false,
                                                  false);
                KSP ksp_mass = massSolver->ksp();
                KSPSetOptionsPrefix(ksp_mass, "mass_");
                KSPSetFromOptions(ksp_mass);
                massSolver->solve();
                VecCopy(massSolver->getCurrentSolution(), vec.v);
                delete massSolver;
                delete massEquation;
            }

        }
        else{
            didCoarsen = performRefinement<RefinementStage::COARSEN>(octDA, distTree, vec, domainExtents, &idata,dummy,dummy,dummy,false);
        }
        double globalQuantityAfterCoarsen[3];
        {
            CHPostProcessing chPostProcessing(octDA, distTree, domainExtents, vec);
            chPostProcessing.getGlobalQuantity(globalQuantityAfterCoarsen, octDA->getCommActive());
        }
        if(!rank) {
            std::ofstream fout("EnergyAndMass.txt", std::ios::app);
            const auto default_precision {std::cout.precision()};

            fout <<  std::setprecision(15) << ti.getCurrentTime() << "," <<
            globalQuantityAfterSolve[0] << "," << globalQuantityAfterSolve[1] << "," << idata.epsilon*globalQuantityAfterSolve[2] << "," <<
            globalQuantityAfterRefine[0] << "," << globalQuantityAfterRefine[1] << "," << idata.epsilon*globalQuantityAfterRefine[2] << "," <<
            globalQuantityAfterCoarsen[0] << "," << globalQuantityAfterCoarsen[1] << "," << idata.epsilon*globalQuantityAfterCoarsen[2] << ","<<
            std::defaultfloat << "\n";
            std::setprecision(default_precision);
            fout.close();
        }
        prevSolution = vec.v;
        if(didRefine or didCoarsen) {
            updateEquationAndSolver(octDA, distTree, chEquation, domainExtents,chSolver,&ti,
                                    nullptr,CHNodeData::CH_DOF,&idata);
            resetDAVecs();
        }



        ti.increment();
    }

    delete chSolver;
    delete chEquation;
    VecDestroy(&prevSolution);
    dendrite_finalize(octDA);



}