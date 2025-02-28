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
#include <mpi.h>

#ifndef DENDRITEKT_NEIGHBORSEARCH_H
#define DENDRITEKT_NEIGHBORSEARCH_H

#endif //DENDRITEKT_NEIGHBORSEARCH_H


class NodePrint {


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
    DENDRITE_REAL * target_coords;
    DENDRITE_REAL * neigh_coords;
    SubDomainBoundary * m_subDomainBoundary;


    TALYFEMLIB::ELEM *m_elem;
    TRAVERSAL_TYPE m_traversalType;
    DENDRITE_UINT m_level_target;
    DENDRITE_UINT m_level_neigh;
    bool m_BoundaryOctant = false;
    TALYFEMLIB::GRID m_grid;
    int m_childNum;

    DENDRITE_UINT nPe_;



public:
    NodePrint(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents, SubDomainBoundary * subDomainBoundary);
    void generate();
    void checkNeighbors(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh,
                          bool &EastRefine, bool &WestRefine, bool &NorthRefine, bool &SouthRefine, bool &FrontRefine, bool &BackRefine);
    void checkNeighborsCoord(const ot::RankI *target_nodeIdsFlat, const ot::RankI *neigh_nodeIdsFlat, DENDRITE_UINT m_level_target, DENDRITE_UINT m_level_neigh);
    bool CoordCheck(int target_idx, int neigh_idx);

};


NodePrint::NodePrint(DA * octDA,  const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents,SubDomainBoundary * subDomainBoundary)
        :m_octDA(octDA), m_octToPhysical(domainExtents), m_treePart(treePart),m_subDomainBoundary(subDomainBoundary){
    target_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
    neigh_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];
    nPe_ = m_octDA->getNumNodesPerElement();

}


void NodePrint::generate() {
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






    ot::MatvecBaseCoords<DIM> loop(ghostedNodalSz, m_octDA->getElementOrder(), false, 0, odaCoords,
                                   &(*m_treePart.cbegin()), m_treePart.size(), *m_octDA->getTreePartFront(),*m_octDA->getTreePartBack());



    int rank = TALYFEMLIB::GetMPIRank();
    std::string rankString = std::to_string(rank);
    std::string fPrefix = rankString + "_trueGP.csv";
    FILE *fp = fopen(fPrefix.c_str(), "w");
    fprintf(fp, "GP_x0,GP_y0,GP_z0\n");

    while (!loop.isFinished()) {
        if (loop.isPre() && loop.subtreeInfo().isLeaf()) {



            const double *target_nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();
//            const RankI *target_nodeIdsFlat = loop.subtreeInfo().readNodeValsIn();
            const bool isBoundaryOctant = loop.subtreeInfo().isElementBoundary();
            bits_ = loop.subtreeInfo().getLeafBitsetInfo();

            auto t = loop.subtreeInfo().getNumNonhangingNodes();
            m_level_target = loop.getCurrentSubtree().getLevel();
            std::memcpy(target_coords, target_nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe_ * DIM);

            m_octToPhysical.convertCoordsToPhys(target_coords,nPe_);

            for (int i = 0; i < nPe_; i++) {
                fprintf(fp, "%.10e,%.10e,%.10e\n", target_coords[DIM * i], target_coords[DIM * i + 1], target_coords[DIM * i + 2]);
            }



//            std::cout << "ID:" << counter <<  "  NumNonhangingNodes: " << t << std::endl;
//            std::cout << "x1:" << target_coords[6] <<  "  y1: " << target_coords[7] << std::endl;




            loop.next();

        } else {
            loop.step();
        }
    }







//
//
//
//
//
//
//    ot::MatvecBase<DIM, RankI> treeloop(ghostedNodalSz, 1, m_octDA->getElementOrder(), odaCoords, &(*ghostedGlobalNodeId.cbegin()),
//                                        &(*this->m_treePart.cbegin()),
//                                        this->m_treePart.size(),
//                                        *m_octDA->getTreePartFront(),
//                                        *m_octDA->getTreePartBack());
//
//
//    std::string fPrefix = "trueGP.vtk";
//    FILE *fp = fopen(fPrefix.c_str(), "w");
//    fprintf(fp, "GP_x0,GP_y0,GP_z0\n");
//
//
//
////    fprintf(fp, "%.10e,%.10e\n",
////           x, y);
//
//
////    MPITimer timer("loop object");
////    timer.Start();
////
////
////
////
////    timer.Stop();
////    timer.PrintTotalTimeSeconds();
//
//
//    int counter = 0;
//    while (!treeloop.isFinished()) {
//        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
//
//
//
//
//            const double *target_nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
//            const RankI *target_nodeIdsFlat = treeloop.subtreeInfo().readNodeValsIn();
//            const bool isBoundaryOctant = treeloop.subtreeInfo().isElementBoundary();
//            bits_ = treeloop.subtreeInfo().getLeafBitsetInfo();
//
//            auto t = treeloop.subtreeInfo().getNumNonhangingNodes();
//            m_level_target = treeloop.getCurrentSubtree().getLevel();
//            std::memcpy(target_coords, target_nodeCoordsFlat, sizeof(DENDRITE_REAL) * nPe_ * DIM);
//
//            for (int i = 0; i < nPe_; i++) {
//                fprintf(fp, "%.10e,%.10e,%.10e\n", target_coords[DIM * i], target_coords[DIM * i + 1], target_coords[DIM * i + 2]);
//            }
//
//
//
////            std::cout << "ID:" << counter <<  "  NumNonhangingNodes: " << t << std::endl;
////            std::cout << "x1:" << target_coords[6] <<  "  y1: " << target_coords[7] << std::endl;
//
//            MPITimer timer("loop object");
//            timer.Start();
//
//
//
//            timer.Stop();
////            timer.PrintTotalTimeSeconds();
//
//            counter ++;
//
//
//
//            treeloop.next();
//        } else
//            treeloop.step();
//    }



}







