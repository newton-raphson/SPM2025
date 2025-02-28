//
// Created by maksbh on 3/26/23.
//

#pragma once

#include <talyfem/grid/femelm.h>
#include "DataTypes.h"
#include "mathUtils.h"
#include "TalyDendroSync.h"

#pragma once

class MeshElements{
    static constexpr unsigned int numBits = 27;
    void initiateOrder1();



    std::bitset<numBits> bits_;
    static constexpr int MAX_COORDS_PER_ELEMENT = 8;
    static constexpr int HEX_COORDS_PER_ELEMENT = 8;
    static constexpr int TET_COORDS_PER_ELEMENT = 4;
    static constexpr int PYRAMID_COORDS_PER_ELEMENT = 5;
    static constexpr int MAX_COORDS_PER_FACE = 9;
    static constexpr int MIDDLE_NODE = 13;
    static constexpr int NUM_FACE = 6;
    static constexpr int NUM_EDGE_PER_FACE = 4;


//    static constexpr int FACE[NUM_FACE][MAX_COORDS_PER_FACE] = {{0, 1, 2, 9, 10, 11, 18, 19, 20 },
//                                                                {0, 3, 6, 9, 12, 15, 18, 21, 24},
//                                                                {6, 7, 8, 15, 16, 17, 24, 25, 26},
//                                                                {8, 5, 2, 17, 14, 11, 26, 23, 20},
//                                                                {0, 1, 2, 3, 4, 5, 6, 7, 8},
//                                                                {18, 19, 20, 21, 22, 23, 24, 25, 26}};

    const int FACE[NUM_FACE][MAX_COORDS_PER_FACE] = {{0, 2, 18, 20, 1, 11, 19, 9, 10},
                                                     {0, 6, 18, 24, 3, 15, 21, 9, 12},
                                                     {6, 8, 24, 26, 7, 17, 25, 15, 16},
                                                     {2, 8, 20, 26, 5, 17, 23, 11, 14},
                                                     {0, 2, 6, 8, 1, 5, 7, 3, 4},
                                                     {18, 20, 24, 26, 19, 23, 25, 21, 22}};
    static constexpr int MAX_HANGING_PER_FACE = 5;
    static constexpr int MIN_NONHANGING_PER_FACE = 4;


    template<int size>
    void createElemInfo(const std::array<int, size> &nodeBits, std::array<int, MAX_COORDS_PER_ELEMENT> &nodeCoords);

//    https://docs.google.com/presentation/d/1MKc2jVD67TCAs2777ObxDblbhjuOVZrynmZh1p6Y07U/edit#slide=id.p
    enum CONFIGURATIONS:int{
        NO_HANGING_ELEMENT = 0,
        NO_HANGING_FACE = 1,
        ONE_NODE_HANGING = 2,
        TWO_NODE_HANGING_CONTIGUOUS = 3,
        TWO_NODE_HANGING_STRIDED = 4,
        THREE_NODE_HANGING = 5,
        FOUR_NODE_HANGING = 6,
        FACE_HANGING = 7,
        MAX_CONFIGURATION = 8
    };

protected:

    TalyDendroSync m_sync;

    TALYFEMLIB::ELEM *m_hex_elem;
    TALYFEMLIB::ELEM *m_pyramid_elem;
    TALYFEMLIB::ELEM *m_tet_elem;


    TALYFEMLIB::GRID m_hex_grid;
    TALYFEMLIB::GRID m_pyramid_grid;
    TALYFEMLIB::GRID m_tet_grid;




    const int m_order;

    std::vector<std::tuple<ElementType_3D,std::array<int,MAX_COORDS_PER_ELEMENT>>> m_elemInfo;
    std::array<CONFIGURATIONS,NUM_FACE> m_configurations;
    std::array<std::array<int,MAX_HANGING_PER_FACE>,NUM_FACE> hangingIDs;


public:
    MeshElements(int order);

    int initElementType(const std::bitset<numBits> &val);

    void generateConfigurations();

    virtual ~MeshElements();

    ElementType_3D getFEMElemObject( const int id, const double * physicalCoords);



};

template<int size>
void MeshElements::createElemInfo(const std::array<int, size> &nodeBits,
                                  std::array<int, MAX_COORDS_PER_ELEMENT> &nodeCoords) {
    for(int i = 0; i < size; i++){
//        assert(this->m_bitsToNodeMap[nodeBits[i]] != -1);
        nodeCoords[i] = nodeBits[i];
    }
}
