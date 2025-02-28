//
// Created by mehdi on 8/7/23.
//


#include "OctToPhysical.h"
#include <talyfem/grid/femelm.h>
#include "DataTypes.h"
#include "mathUtils.h"
#include "TalyDendroSync.h"
#include "Boundary/SubDomainBoundary.h"
#include "DataTypes.h"
#include "talyfem/grid/elem_types/elem2dbox.h"
#include "talyfem/grid/elem_types/elem2dtriangle.h"
#include "sfcTreeLoop_matvec.h"
#include "talyfem/utils/timers.h"
#include <Traversal/Traversal.h>
#include <mpi.h>
#include "Mesh/Mesh2D.h"

#ifndef DENDRITEKT_NEIGHBORSEARCH_H
#define DENDRITEKT_NEIGHBORSEARCH_H

#endif //DENDRITEKT_NEIGHBORSEARCH_H


class Utils : public MeshElements{


    static constexpr unsigned int numBits = intPow(3,DIM);
    void initiateOrder1();



    std::bitset<numBits> bits_;
    static constexpr int MAX_COORDS_PER_ELEMENT = 4;
    static constexpr int QUAD_COORDS_PER_ELEMENT = 4;
    static constexpr int TRI_COORDS_PER_ELEMENT = 3;
    static constexpr int MAX_COORDS_PER_FACE = 3;
    static constexpr int MIDDLE_NODE = 4;
    static constexpr int NUM_FACE = 4;


protected:
    enum GMSH_GRID:int{
        GMSH_TRIANGLE = 2,
        GMSH_QUAD = 3,
        GMSH_LINE = 1
    };

    std::vector<std::array<double,DIM>> m_vertexPosition;
    std::vector<std::tuple<GMSH_GRID,std::vector<PetscInt>>> gridElements;
    std::vector<std::tuple<GMSH_GRID,int,std::vector<PetscInt>>> boundary_elements;


    DA * m_octDA;
    OctToPhysical m_octToPhysical;
    const std::vector<TREENODE> & m_treePart;
    DENDRITE_REAL * m_coords;
    DENDRITE_REAL * target_coords;
    DENDRITE_REAL * neigh_coords;
    SubDomainBoundary * m_subDomainBoundary;


    TALYFEMLIB::ELEM *m_elem;
    TRAVERSAL_TYPE m_traversalType;
    DENDRITE_UINT m_level;
    DENDRITE_UINT m_level_neigh;
    bool m_BoundaryOctant = false;
    TALYFEMLIB::GRID m_grid;
    int m_childNum;

    DENDRITE_UINT nPe_;

    bool hasLocalGaussPoint = true;

    std::vector<TALYFEMLIB::ZEROPTV> DomainPtLocal;
    std::vector<TALYFEMLIB::ZEROPTV> DomainPtGlobal;
    std::vector<TALYFEMLIB::ZEROPTV> GaussPtLocal;
    std::vector<TALYFEMLIB::ZEROPTV> GaussPtGlobal;


    TalyDendroSync sync_;
    DENDRITE_UINT eleOrder_;
    DENDRITE_UINT elemID = 0;
    const int * relativeOrder_ = nullptr;

    inline int getRelativeOrder(int localElemNumber) const {
        int relativeOrder = 0;
        if(relativeOrder_) {
            relativeOrder = (relativeOrder_[localElemNumber]);
        }
        return relativeOrder;
    }

public:
    Utils(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents, SubDomainBoundary * subDomainBoundary);
    Utils(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents);
    void constructTalyGrid(const double *coords);
    void WriteTrueGPToFile();
    void checkNeighbors(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh,
                          bool &EastRefine, bool &WestRefine, bool &NorthRefine, bool &SouthRefine, bool &FrontRefine, bool &BackRefine);
    void checkNeighborsCoord(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh);
    bool CoordCheck(int target_idx, int neigh_idx);

};


//Utils::Utils(DA * octDA,  const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents,SubDomainBoundary * subDomainBoundary)
//        :m_octDA(octDA), m_octToPhysical(domainExtents), m_treePart(treePart),m_subDomainBoundary(subDomainBoundary){
//    target_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
//    neigh_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
//    nPe_ = m_octDA->getNumNodesPerElement();
//    eleOrder_ = m_octDA->getElementOrder();
//}
//
//
//Utils::Utils(DA * octDA,  const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents)
//        :m_octDA(octDA), m_octToPhysical(domainExtents), m_treePart(treePart){
//    target_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
//    neigh_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
//    nPe_ = m_octDA->getNumNodesPerElement();
//    eleOrder_ = m_octDA->getElementOrder();
//}





void Utils::constructTalyGrid(const double *coords){
    if(eleOrder_ == 1) {
        sync_.syncCoords<1>(coords, &m_grid);
        return;
    }
    if(eleOrder_ == 2) {
        sync_.syncCoords<2>(coords, &m_grid);
        return;
    }
    throw  TALYFEMLIB::TALYException() << "NotSupported in func" << " " << __func__ << "for eleOrder = " << eleOrder_ << "\n" ;

}




void Utils::WriteTrueGPToFile() {
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;

#ifdef BUILD_WITH_PETSC
    using ScalarT = PetscScalar;
    using IndexT = PetscInt;
#else
    using ScalarT = DendroScalar;
    using IndexT = long long unsigned;
#endif
    const size_t ghostedNodalSz = m_octDA->getTotalNodalSz();
    const ot::TreeNode<CoordT, DIM> *odaCoords = m_octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = m_octDA->getNodeLocalToGlobalMap();
    std::vector<std::tuple<std::array<int,2>,unsigned int>> boundaryFaceID;


    std::string fPrefix = "trueGP.vtk";FILE *fp = fopen(fPrefix.c_str(), "w");
    if (TALYFEMLIB::GetMPIRank() == 0) { // if:master cpu, then:print


#if (DIM == 2)
        fprintf(fp, "GP_x0,GP_y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "GP_x0,GP_y0,GP_z0\n");
#endif

    }

    ot::MatvecBase<DIM, RankI> treeloop(ghostedNodalSz, 1, m_octDA->getElementOrder(), odaCoords, &(*ghostedGlobalNodeId.cbegin()),
                                        &(*this->m_treePart.cbegin()),
                                        this->m_treePart.size(),
                                        *m_octDA->getTreePartFront(),
                                        *m_octDA->getTreePartBack());


    int counter = 0;
    while (!treeloop.isFinished()) {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
            const double *nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            m_BoundaryOctant = treeloop.subtreeInfo().isElementBoundary();
            m_level = treeloop.getCurrentSubtree().getLevel();
            m_childNum = treeloop.getCurrentSubtree().getMortonIndex();
            const int numNodes = treeloop.subtreeInfo().getNumNodesIn();

            std::cout << numNodes << std::endl;

            std::memcpy(m_coords,nodeCoordsFlat, sizeof(DENDRITE_REAL)*nPe_*DIM);
            m_octToPhysical.convertCoordsToPhys(m_coords,nPe_);
            constructTalyGrid(m_coords);
            TALYFEMLIB::FEMElm fe(&m_grid, TALYFEMLIB::BASIS_ALL);
            fe.refill(0,this->getRelativeOrder(elemID));
            while (fe.next_itg_pt()) {

#if (DIM == 2)
                fprintf(fp, "%.10e,%.10e\n",
                        fe.position().x(), fe.position().y());
#endif
#if (DIM == 3)
                fprintf(fp, "%.10e,%.10e,%.10e\n",
                        fe.position().x(), fe.position().y(), fe.position().z());
#endif


            }

            elemID++;
            treeloop.next();
        } else
            treeloop.step();
    }

}


