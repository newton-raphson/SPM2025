//
// Created by maksbh on 1/27/23.
//

#include "DendriteUtils.h"
#include "PETSc/PetscUtils.h"
#include "PETSc/IO/petscVTU.h"
#include "MassComputation.h"
#include "Refinement.h"
#include "utils.h"
#include "MEquation.h"
#include "MNodeData.h"
#include "MassCalculatorTotal.h"
static const char *varName_ch[]{"ch"};

using namespace PETSc;
int main(int argc, char *argv[]) {
    dendrite_init(argc,argv);

    int rank = TALYFEMLIB::GetMPIRank();
    DENDRITE_UINT eleOrder = 1;
    DENDRITE_UINT level = 4;

    bool mfree = false;

    DomainInfo physDomain;
    physDomain.min.fill(0.0);
    physDomain.max.fill(1.0);
    DomainExtents domainExtents(physDomain);


    TALYFEMLIB::PrintInfo("Total number of processor = ", TALYFEMLIB::GetMPISize());
    TALYFEMLIB::PrintInfo("size of DendroInt ", sizeof(DendroIntL));
    TALYFEMLIB::PrintInfo("size of PetscInt ", sizeof(PetscInt));

    TALYFEMLIB::PrintStatus("eleOrder ", eleOrder);
    TALYFEMLIB::PrintStatus("Level ", level);
    TALYFEMLIB::PrintStatus("Mfree ", mfree);
    TALYFEMLIB::PrintStatus("DIM =  ", DIM);

    /// create DA
    SubDomain subDomain (domainExtents);
    std::function<ibm::Partition(const double *, double)> functionToRetain = [&](const double *octCoords, double scale) {
        return (subDomain.functionToRetain (octCoords, scale));
    };
    DistTREE  dTree;
    DA * octDA = createSubDA (dTree, functionToRetain, level, eleOrder);


    /// Problem size
    DENDRITE_UINT ndof = 1;

    std::function<void(const double *, double *)> initial_condition_ch = [&](const double *x, double *var) {
        var[0] = abs(cos(2*M_PI*x[0])) + 10;
//        var[0] = x[0] + x[1];
//          var[0] = 5;
    };

    Vec scalarValue;
    octDA->petscCreateVector(scalarValue,false,false,1);
    octDA->petscSetVectorByFunction(scalarValue,initial_condition_ch,false,false,1);

    petscVectopvtu(octDA,dTree.getTreePartFiltered(),scalarValue,"before","ch",varName_ch,domainExtents,false,false,1);

    VecInfo vec(scalarValue,1,0);

    std::vector<double> elemental_mass,element_mass_perChild;
    double originalMass;
    {

        TALYFEMLIB::PrintStatus("Num nodes = ", octDA->getGlobalNodeSz());
        MassCalculator massCalculator(octDA, dTree.getTreePartFiltered(), vec, domainExtents,elemental_mass,element_mass_perChild);
        MassCalculatorTotal massCalculatorTotal(octDA, dTree.getTreePartFiltered(), vec,domainExtents);
        TALYFEMLIB::PrintStatus("Total Mass = ", massCalculatorTotal.getTotalMass());
    }

    performRefinement(octDA, dTree,vec,domainExtents,subDomain,elemental_mass,element_mass_perChild);

    {
        MassCalculatorTotal massCalculatorTotal(octDA, dTree.getTreePartFiltered(), vec,domainExtents);
        TALYFEMLIB::PrintStatus("Mass by injection = ", massCalculatorTotal.getTotalMass());
//        std::vector<double> elementalMassByInjection;
//        TALYFEMLIB::PrintStatus("Num nodes = ", octDA->getGlobalNodeSz());
//        MassCalculator massCalculator(octDA, dTree.getTreePartFiltered(), vec, domainExtents,elementalMassByInjection,element_mass_perChild);
//        TALYFEMLIB::PrintStatus("Total Mass = ", massCalculator.getTotalMass());
////        auto difference = findDifference(elemental_mass,elementalMassByInjection,octDA->getCommActive());
//        TALYFEMLIB::PrintStatus("Max DIfference [Injection] = ", std::get<0>(difference), " TotalDifference = " , std::get<1>(difference));
//        TALYFEMLIB::PrintStatus("Global Mass difference = ", fabs(originalMass - massCalculator.getTotalMass()));


    }
    auto massEquation = new TalyEquation<MEquation,MNodeData> (octDA,dTree.getTreePartFiltered(),
                                                               subDomain.domainExtents(),1,
                                                               nullptr,false,nullptr
                                                               );
    massEquation->equation()->assign(&elemental_mass,&element_mass_perChild);

    auto massSolver = setLinearSolver(massEquation,octDA,dTree,1,false,false, true);

    massSolver->solve();



    Vec solverSolution = massSolver->getCurrentSolution();
    {
        vec.v = solverSolution;
        MassCalculatorTotal massCalculatorTotal(octDA, dTree.getTreePartFiltered(), vec,domainExtents);
        TALYFEMLIB::PrintStatus("Mass by projection = ", massCalculatorTotal.getTotalMass());
//        std::vector<double> elementalMassByProjection;
//
//
//        TALYFEMLIB::PrintStatus("Num nodes = ", octDA->getGlobalNodeSz());
//        MassCalculator massCalculator(octDA, dTree.getTreePartFiltered(), vec, domainExtents,elementalMassByProjection,element_mass_perChild);
//        TALYFEMLIB::PrintStatus("Total Mass = ", massCalculator.getTotalMass());
////        auto difference = findDifference(elemental_mass,elementalMassByProjection,octDA->getCommActive());
////        TALYFEMLIB::PrintStatus("Max DIfference [L2projection] = ", std::get<0>(difference), " TotalDifference = " , std::get<1>(difference));
////        TALYFEMLIB::PrintStatus("Global Mass difference = ", fabs(originalMass - massCalculator.getTotalMass()));

    }
    petscVectopvtu(octDA,dTree.getTreePartFiltered(),vec.v,"after","ch",varName_ch,domainExtents,false,false,1);



    VecDestroy(&vec.v);
    dendrite_finalize(octDA);



}