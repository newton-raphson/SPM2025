//
// Created by maksbh on 3/26/23.
//

#pragma once

#include <talyfem/grid/femelm.h>
#include "DataTypes.h"
#include "mathUtils.h"
#include "TalyDendroSync.h"
#include "Boundary/SubDomainBoundary.h"

class MeshElements{


    static constexpr unsigned int numBits = intPow(3,DIM);
    void initiateOrder1();



    std::bitset<numBits> bits_;
    static constexpr int MAX_COORDS_PER_ELEMENT = 4;
    static constexpr int QUAD_COORDS_PER_ELEMENT = 4;
    static constexpr int TRI_COORDS_PER_ELEMENT = 3;
    static constexpr int MAX_COORDS_PER_FACE = 3;
    static constexpr int MIDDLE_NODE = 4;
    static constexpr int NUM_FACE = 4;




    template<int size>
    void createElemInfo(const std::array<int, size> &nodeBits, std::array<int, MAX_COORDS_PER_ELEMENT> &nodeCoords);

protected:
    const int FACE[NUM_FACE][MAX_COORDS_PER_FACE] = {{0,1,2},
                                                     {0,3,6},
                                                     {6,7,8},
                                                     {8,5,2}};

    TalyDendroSync m_sync;
    TALYFEMLIB::ELEM *m_triangle_elem;
    TALYFEMLIB::ELEM *m_quad_elem;

    TALYFEMLIB::GRID m_triangle_grid;
    TALYFEMLIB::GRID m_quad_grid;


    const int m_order;

    std::vector<std::tuple<ElementType_2D,std::array<int,MAX_COORDS_PER_ELEMENT>>> m_elemInfo;
//    std::array<int,numBits> m_bitsToNodeMap;



public:
    MeshElements(int order);

    int initElementType(const std::bitset<numBits> & val);

    ElementType_2D getFEMElemObject( const int id, const double * physicalCoords);

    void generateFaces(SubDomainBoundary * subDomainBoundary,const int id, const double * physCoords,
                       std::vector<std::tuple<std::array<int,2>,unsigned int>>& boundaryFaceID);

};


template<int size>
void MeshElements::createElemInfo(const std::array<int, size> &nodeBits, std::array<int, MAX_COORDS_PER_ELEMENT> &nodeCoords){
    static_assert(DIM == 2,"Cannot include this file");

    for(int i = 0; i < size; i++){
//        assert(this->m_bitsToNodeMap[nodeBits[i]] != -1);
        nodeCoords[i] = nodeBits[i];
    }
}