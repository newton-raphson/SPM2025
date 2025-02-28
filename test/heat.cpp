//
// Created by maksbh on 3/26/23.
//

#include "DendriteUtils.h"
#include "SubDA/SubDomain.h"
#include "HeatEquation/SSHTNodeData.h"
#include "HeatEquation/SSHTEquation.h"
#include "PETSc/VecInfo.h"
#include "IO/VTU.h"
#include "PETSc/IO/petscVTU.h"
#include "RefinementUtils.h"
#include "TalyEquation.h"
#include "PETSc/Solver/LinearSolver.h"
#include <PETSc/PetscUtils.h>
#include "Traversal/Analytic.h"
#include "Traversal/GenerateMesh.h"
#include "HeatEquation/Utils.h"

static const char * varname [] = {"nodes"};
using namespace PETSc;
int main(int argc, char *argv[]) {

    dendrite_init(argc,argv);
    DistTREE distTree;

    DomainInfo fullDADomain;
    DomainInfo physDomain;
    fullDADomain.min.fill(0.0);
    fullDADomain.max.fill(1.0);

    physDomain.min.fill(0.0);
    physDomain.max.fill(1.0);

    DomainExtents domainExtents(fullDADomain,physDomain);
    SubDomain subDomain(domainExtents);
    DENDRITE_REAL circleCenter[DIM]{0.5,0.5};
    VOXEL::Circle circle(circleCenter, 0.2, RetainSide::OUT);
    DENDRITE_REAL left[DIM]{0.4,0.4};
    DENDRITE_REAL right[DIM]{0.6,0.6};
//    VOXEL::Box box(left, right, RetainSide::OUT);
//    subDomain.addObject(box);


    int level = std::atoi(argv[1]);
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };
    DA *octDA = createSubDA(distTree,functionToRetain,level,1,0.3);

    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);

//    IO::writeBoundaryElements(octDA, distTree.getTreePartFiltered(), "test2", "test2", domainExtents);
    getnewDABoundaryRefine(octDA,distTree,domainExtents,1);

//    delete octDA;
//    DA *newDA = new DA(distTree,0,MPI_COMM_WORLD,1,100,0.1,DENDROMESH::NOHANGING_NODES);



    int ndof = 1;
    auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, distTree.getTreePartFiltered(), domainExtents);
    const int order = 1;
    sshtEq->setRelativeOrder(&order);

    LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA,distTree,  ndof, false);

    sshtSolver->setDirichletBoundaryCondition([&](const TALYFEMLIB::ZEROPTV pos, unsigned int nodeID) -> Boundary {
        Boundary b;
        static constexpr double eps = 1e-14;
        double x = pos.x();
        double y = pos.y();

        bool on_wall = (fabs(x - physDomain.min[0]) < eps) ||
                       (fabs(x - physDomain.max[0]) < eps) ||
                       (fabs(y - physDomain.min[1]) < eps) ||
                       (fabs(y - physDomain.max[1]) < eps);
#if (DIM >= 3)
        double z = pos.z();
    on_wall = on_wall || (fabs(z - physDomain.min[2]) < eps) || (fabs(z - physDomain.max[2]) < eps);
#endif
        if (on_wall) {
            b.addDirichlet(0, 0);
        }
        return b;
    });

    sshtSolver->solve();


    /// L2 error
    const auto analytic_sol = [](const TALYFEMLIB::ZEROPTV &pos, const DENDRITE_UINT dof, const DENDRITE_REAL time) {
#if (DIM == 2)
        return sin(M_PI * pos.x()) * sin(M_PI * pos.y());
//        return sin(2 * M_PI * pos.x()) * sin(2 * M_PI * pos.y());
//        return sin(4 * M_PI * pos.x()) * sin(4 * M_PI * pos.y());
#elif (DIM == 3)
        return sin(M_PI * pos.x()) * sin(M_PI * pos.y()) * sin(M_PI * pos.z());
#endif

    };
    {
        VecInfo v(sshtSolver->getCurrentSolution(), 1, 0);

        petscVectopvtu(octDA, distTree.getTreePartFiltered(), sshtSolver->getCurrentSolution(), "ssht", varname, domainExtents, false, false, ndof);


        Analytic sshtAnalytic(octDA, distTree.getTreePartFiltered(), v, analytic_sol, physDomain);

        sshtAnalytic.getL2error();
    }
    Vec localVec;
    octDA->petscCreateVector(localVec, false, false,ndof);
    std::function<void(const double *, double *)> initial_condition = [](const double *x, double *var) {
#if (DIM == 2)
        var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]);
//        var[0] = sin(2 * M_PI * x[0]) * sin(2 * M_PI * x[1]);
//        var[0] = sin(4 * M_PI * x[0]) * sin(4 * M_PI * x[1]);
#elif(DIM == 3)
        var[0] = sin(M_PI * x[0]) * sin(M_PI * x[1]) * sin (M_PI * x[2]);
#endif
    };
    octDA->petscSetVectorByFunction(localVec, initial_condition, false, false, ndof);
    PetscScalar * array, * soln_array;
    VecGetArray(localVec,&array);
    VecGetArray(sshtSolver->getCurrentSolution(),&soln_array);

    double sum_diff = 0;
    double max_diff = 0;
    for(int i = 0; i < octDA->getLocalNodalSz(); i++){
        sum_diff += std::abs(array[i] - soln_array[i]);
        max_diff = std::max(max_diff,std::abs(array[i] - soln_array[i]));
//        array[i] = std::abs((array[i] - soln_array[i])/(array[i]+1e-12));
        array[i] = std::abs((array[i] - soln_array[i]));
    }
    VecRestoreArray(localVec,&array);
    VecRestoreArray(sshtSolver->getCurrentSolution(),&soln_array);
    std::cout << octDA->getLocalNodalSz() <<"  "<< sum_diff/octDA->getLocalNodalSz() << " " << max_diff <<  "\n";
    {

        petscVectopvtu(octDA, distTree.getTreePartFiltered(), localVec, "error", varname, domainExtents, false, false, ndof);


    }


//    Utils Utils(newDA,distTree.getTreePartFiltered(),domainExtents);
//    Utils.WriteTrueGPToFile();  //


    SubDomainBoundary boundary(&subDomain,octDA,domainExtents);





//






    dendrite_finalize(octDA);

}
