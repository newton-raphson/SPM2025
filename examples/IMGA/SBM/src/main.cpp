//
// Created by maksbh on 7/11/20.
//
#include <vector>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
#include <IMGA/IMGA.h>
#include <IMGA/Marker.h>
#include <TalyEquation.h>
#include <sfcTreeLoop_matvec_io.h>
#include <PETSc/PetscUtils.h>
#include <PETSc/IO/petscVTU.h>
#include "utils.h"

#include "SSHTEquation.h"
#include "SSHTNodeData.h"
#include "IMGALoop.h"
#include "CheckSurface.h"
#include "DARefine.h"
#include "DACoarsen.h"
using namespace PETSc;


/**
 * @brief An example to test the background element of the processor.
 * In order for efficent implementation, Point should be used instead of ZEROPTV.
 */
int main(int argc, char *argv[]) {
    dendrite_init(argc, argv);
    int rank = TALYFEMLIB::GetMPIRank();

    DomainInfo domainInfoFullDomain,domainfoPhysDomain;
    domainInfoFullDomain.min.fill(0.0);
    domainInfoFullDomain.max.fill(1.0);

    domainfoPhysDomain.min.fill(0.0);
    domainfoPhysDomain.max.fill(1.0);

    DomainExtents domainExtents(domainInfoFullDomain,domainfoPhysDomain);
    const DENDRITE_UINT eleOrder = 1;
    const std::string filename = argv[1];
    const int baseLevel = std::atoi(argv[2]);
    const int refineLevel = std::atoi(argv[3]);
    bool doRefine = (bool) std::atoi(argv[4]);



#if(DIM == 2)
    SubDomain subDomain(domainExtents);
    double _center[DIM] {0.5,0.5};
    VOXEL::Circle circle(_center,0.5,RetainSide::IN);
    subDomain.addObject(circle);
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };
    DistTREE distTree;
    DA *octDA = createSubDA(distTree,functionToRetain,level,eleOrder);
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
    performRefinement(octDA,distTree,domainExtents,subDomain);
    IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"Boundary","boundary",domainExtents);
    Vec nodalFalseElement;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> elementMarkers;
    generateNewMarkers(octDA,distTree,domainExtents,subDomain,elementMarkers,nodalFalseElement);
    SubDomainBoundary subDomainBoundary(&subDomain,octDA,domainExtents);
    auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, distTree.getTreePartFiltered(), domainExtents,1, nullptr,
                                                               true,&subDomainBoundary,IBM_METHOD::SBM);

    sshtEq->setVectors({VecInfo(nodalFalseElement,1,NODE_ID)});
    IMGA imga(domainExtents,IBM_METHOD::SBM);
    sshtEq->assignIBMConstructs(&imga,elementMarkers.data(),NODE_ID);

    // Correct Cycles
    CheckSurface checkSurface(octDA,distTree.getTreePartFiltered(),{VecInfo(nodalFalseElement,1,NODE_ID)},domainExtents,elementMarkers,&subDomain);
    // This will correct the element markers
    checkSurface.correctCycles();

    LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA,distTree, 1, false);
    sshtSolver->solve();

    std::vector<TALYFEMLIB::ZEROPTV> positionOfpointToInterpolateOn, positionOnSurface;

    // TODO: Fill with Gauss points.

    IMGA imga1(domainExtents,IBM_METHOD::SBM);
    imga1.initIMGAComputation(octDA,distTree.getTreePartFiltered(),positionOfpointToInterpolateOn,positionOnSurface);

    VecInfo vec(sshtSolver->getCurrentSolution(),1,0);
    IMGALoop loop(octDA,&imga1,distTree.getTreePartFiltered(),vec,domainExtents);
    DENDRITE_REAL error;
    loop.computeBoundaryError(&error);

    delete sshtSolver;
    delete sshtEq;

#endif
#if(DIM == 3)


    SubDomain subDomain(domainExtents);
    GEOMETRY::STL stl(filename);
    Point<DIM> translate(0.0,0.0,0.0);
    GEOMETRY::Geometry geometry(&stl,translate,RetainSide::IN);
    subDomain.addObject(&geometry);
    std::function<ibm::Partition(const double *, double )> functionToRetain = [&](const double *physCoords, double physSize) {
        return (subDomain.functionToRetain(physCoords, physSize));
    };
    DistTREE distTree;
    DA *octDA = createSubDA(distTree,functionToRetain,baseLevel,eleOrder);
    subDomain.finalize(octDA,distTree.getTreePartFiltered(),domainExtents);
    performRefinement(octDA,distTree,domainExtents,subDomain);

    IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"Boundary","bnd",domainExtents);
//
    if(doRefine) {
        while (true) {
            DARefine refine(octDA, distTree.getTreePartFiltered(), domainExtents, refineLevel,false);

            DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateOutByIn);
            if (newDA == nullptr) {
                break;
            }
            std::swap(newDA, octDA);
            delete newDA;

            subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
        }

        while (true) {
            DARefine refine(octDA, distTree.getTreePartFiltered(), domainExtents, refineLevel+1,true);

            DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateOutByIn);
            if (newDA == nullptr) {
                break;
            }
            std::swap(newDA, octDA);
            delete newDA;

            subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
        }

    }
    {

        while (true) {
            DACoarse refine(octDA, distTree.getTreePartFiltered(), domainExtents,refineLevel, false);
            DA *newDA = refine.getRefineSubDA(distTree, 0.03, RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
            if (newDA == nullptr){
                break;
            }
            std::swap(newDA, octDA);
            delete newDA;

            subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainExtents);
        }
    }
//    std::cout << octDA->getLocalElementSz() << "\n";

    IO::writeBoundaryElements(octDA,distTree.getTreePartFiltered(),"BoundaryFinal","bnd",domainExtents);


    Vec nodalFalseElement;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> elementMarkers;
    generateNewMarkers(octDA,distTree,domainExtents,subDomain,elementMarkers,nodalFalseElement);
    SubDomainBoundary subDomainBoundary(&subDomain,octDA,domainExtents);
    auto sshtEq = new TalyEquation<SSHTEquation, SSHTNodeData>(octDA, distTree.getTreePartFiltered(), domainExtents,1, nullptr,
                                                               true,&subDomainBoundary,IBM_METHOD::SBM);

//    VecSet(nodalFalseElement,0);
    sshtEq->setVectors({VecInfo(nodalFalseElement,1,NODE_ID)});
    IMGA imga(domainExtents,IBM_METHOD::SBM);
    sshtEq->assignIBMConstructs(&imga,elementMarkers.data(),NODE_ID);
    LinearSolver *sshtSolver = setLinearSolver(sshtEq, octDA,distTree, 1, false);
    sshtSolver->solve();



#endif

    dendrite_finalize(octDA);
}