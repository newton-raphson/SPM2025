//
// Created by maksbh on 3/26/23.
//

#include "Mesh/Mesh3D.h"
#include "DataTypes.h"
#include "talyfem/grid/elem_types/elemPyramid.h"
#include "talyfem/grid/elem_types/elem3dtetrahedral.h"
#include "talyfem/grid/elem_types/elem3dhexahedral.h"

MeshElements::MeshElements(const int order)
        : m_order(order) {
    if (order != 1) {
        throw std::runtime_error("Only order = 1 currently supported");
    }
    if (order == 1) {
        this->initiateOrder1();
    }
}

void MeshElements::initiateOrder1() {
    // Tet element
    {
        static constexpr DENDRITE_UINT numNodes = 4;
        int node_id_array[numNodes];
        m_tet_grid.redimArrays(numNodes, 1);
        for (int i = 0; i < numNodes; i++) {
            m_tet_grid.node_array_[i] = new TALYFEMLIB::NODE();
            node_id_array[i] = i;
        }

        m_tet_elem = new TALYFEMLIB::ELEM3dTetrahedral();
        m_tet_grid.elm_array_[0] = m_tet_elem;


        m_tet_elem->redim(numNodes, node_id_array);
    }

    // Pyramid Element
    {
        static constexpr DENDRITE_UINT numNodes = 5;
        int node_id_array[numNodes];
        m_pyramid_grid.redimArrays(numNodes, 1);
        for (int i = 0; i < numNodes; i++) {
            m_pyramid_grid.node_array_[i] = new TALYFEMLIB::NODE();
            node_id_array[i] = i;
        }

        m_pyramid_elem = new TALYFEMLIB::ELEMPyramid();
        m_pyramid_grid.elm_array_[0] = m_pyramid_elem;

        m_pyramid_elem->redim(numNodes, node_id_array);
    }


    // Hex Element
    {
        static constexpr DENDRITE_UINT numNodes = 8;
        int node_id_array[numNodes];
        m_hex_grid.redimArrays(numNodes, 1);
        for (int i = 0; i < numNodes; i++) {
            m_hex_grid.node_array_[i] = new TALYFEMLIB::NODE();
            node_id_array[i] = i;
        }

        m_hex_elem = new TALYFEMLIB::ELEM3dHexahedral();
        m_hex_grid.elm_array_[0] = m_hex_elem;

        m_hex_elem->redim(numNodes, node_id_array);
    }
}

void MeshElements::generateConfigurations() {

    if (bits_.test(MIDDLE_NODE) == false) {
        assert(bits_.any() == false);
        std::fill(m_configurations.begin(),m_configurations.end(),CONFIGURATIONS::NO_HANGING_ELEMENT);
    }
    else{
        for(int face = 0; face < NUM_FACE; face++){
            int numHanging = 0;
            for(int coordsOnFace = MIN_NONHANGING_PER_FACE; coordsOnFace < MAX_COORDS_PER_FACE; coordsOnFace++){
                if(bits_.test(FACE[face][coordsOnFace])){
                    hangingIDs[face][numHanging] = coordsOnFace;
                    numHanging++;
                }
            }
            if(numHanging == 5){
                m_configurations[face] = CONFIGURATIONS::FACE_HANGING;
            }
            else if(numHanging == 0){
                m_configurations[face] = CONFIGURATIONS::NO_HANGING_FACE;
            }
            else if(numHanging == 1){
                assert(bits_.test(FACE[face][MAX_COORDS_PER_FACE - 1] == false));
                m_configurations[face] = CONFIGURATIONS::ONE_NODE_HANGING;
            }
            else if(numHanging == 2){
                assert(bits_.test(FACE[face][MAX_COORDS_PER_FACE - 1] == false));
                int val = hangingIDs[face][1] - hangingIDs[face][0];
                if((val == 1) or (val == 3)){
                    m_configurations[face] = CONFIGURATIONS::TWO_NODE_HANGING_CONTIGUOUS;
                }
                else{
                    assert(val == 2);
                    m_configurations[face] = CONFIGURATIONS::TWO_NODE_HANGING_STRIDED;
                }
            }
            else if(numHanging == 3){
                assert(bits_.test(FACE[face][MAX_COORDS_PER_FACE - 1] == false));
                m_configurations[face] = CONFIGURATIONS::THREE_NODE_HANGING;
            }
            else if(numHanging == 4){
                assert(bits_.test(FACE[face][MAX_COORDS_PER_FACE - 1] == false));
                m_configurations[face] = CONFIGURATIONS::FOUR_NODE_HANGING;
            }
        }
    }

}
int MeshElements::initElementType(const std::bitset<numBits> &val){
    bits_ = val;

    int numElements = 0;
    m_elemInfo.clear();

    if (bits_.test(MIDDLE_NODE) == false) {
        // No Hanging nodes in this element
        numElements += 1;
        m_elemInfo.reserve(numElements);
        const std::array<int, HEX_COORDS_PER_ELEMENT> nodeBits = {0, 2, 6, 8, 18, 20, 24, 26};
        std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
        this->createElemInfo<HEX_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
        m_elemInfo.push_back(std::make_tuple(ElementType_3D::HEX, nodeCoords));
    }
    else{
        this->generateConfigurations();
        for(int face = 0; face < NUM_FACE; face++){
            const auto & configuration = m_configurations[face];

            if(configuration == CONFIGURATIONS::NO_HANGING_FACE){
                static constexpr int NUM_PYRAMID = 1;
                numElements += NUM_PYRAMID;
                m_elemInfo.reserve(numElements);
                const std::array<int, PYRAMID_COORDS_PER_ELEMENT> nodeBits = {FACE[face][0],FACE[face][1],FACE[face][2],
                                                                              FACE[face][3],MIDDLE_NODE};
                std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
                this->createElemInfo<PYRAMID_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                m_elemInfo.push_back(std::make_tuple(ElementType_3D::PYRAMID, nodeCoords));
            }
            else if(configuration == CONFIGURATIONS::FACE_HANGING){
                static constexpr int NUM_PYRAMIDS = 4;
                static constexpr int PYRAMID_COORDS_ON_FACE = 4;
                numElements += NUM_PYRAMIDS;
                m_elemInfo.reserve(numElements);
                static constexpr int PYRAMID_CONFIGURATIONS[1][NUM_PYRAMIDS][PYRAMID_COORDS_ON_FACE] =
                        {
                                {{0, 4, 7, 8},{4, 1, 8, 5}, {7, 8, 2, 6},{8, 5, 6, 3}},
                        };
                for(int i = 0; i < NUM_PYRAMIDS; i++){
                    const std::array<int, PYRAMID_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][PYRAMID_CONFIGURATIONS[0][i][0]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][1]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][2]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][3]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<PYRAMID_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::PYRAMID, nodeCoords));
                }
            }
            // One node hanging
            else if(configuration == CONFIGURATIONS::ONE_NODE_HANGING){
                static constexpr int NUM_TETS = 3;
                static constexpr int TET_COORDS_ON_FACE = 3;
                numElements += NUM_TETS;
                m_elemInfo.reserve(numElements);
                static constexpr int TET_CONFIGURATIONS[NUM_EDGE_PER_FACE][NUM_TETS][TET_COORDS_ON_FACE] =
                        {
                            {{0, 4, 2},{4, 3, 2}, {4, 1, 3}},
                            {{0, 1, 5},{5, 3, 2}, {2, 0, 5}},
                            {{0, 2, 6},{0, 6, 1}, {1, 3, 6}},
                            {{0, 1, 7},{7, 1, 3}, {7, 2, 3}}
                        };
                int hangingNodeID = hangingIDs[face][0];
                for(int i = 0; i < NUM_TETS; i++){
                    const std::array<int, TET_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][TET_CONFIGURATIONS[hangingNodeID - 4][i][0]],
                             FACE[face][TET_CONFIGURATIONS[hangingNodeID - 4][i][1]],
                             FACE[face][TET_CONFIGURATIONS[hangingNodeID - 4][i][2]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TET_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::TET, nodeCoords));
                }
            }
            // Two hanging strided
            else if(configuration == CONFIGURATIONS::TWO_NODE_HANGING_STRIDED){
                static constexpr int NUM_PYRAMIDS = 2;
                static constexpr int PYRAMID_COORDS_ON_FACE = 4;
                numElements += NUM_PYRAMIDS;
                m_elemInfo.reserve(numElements);
                static constexpr int PYRAMID_CONFIGURATIONS[2][NUM_PYRAMIDS][PYRAMID_COORDS_ON_FACE] =
                        {
                                {{0, 4, 2, 6},{4, 1, 3, 6}},
                                {{0, 1, 7, 5},{7, 5, 2, 3}},
                        };
                int id = hangingIDs[face][0] - 4;
                for(int i = 0; i < NUM_PYRAMIDS; i++){
                    const std::array<int, PYRAMID_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][PYRAMID_CONFIGURATIONS[id][i][0]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][1]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][2]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][3]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<PYRAMID_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::PYRAMID, nodeCoords));
                }
            }
            else if(configuration == CONFIGURATIONS::TWO_NODE_HANGING_CONTIGUOUS){
                static constexpr int NUM_TETS = 4;
                static constexpr int TET_COORDS_ON_FACE = 3;
                numElements += NUM_TETS;
                m_elemInfo.reserve(numElements);
                static constexpr int TET_CONFIGURATIONS[NUM_EDGE_PER_FACE][NUM_TETS][TET_COORDS_ON_FACE] =
                        {
                                {{0, 4, 2},{4, 5, 2}, {5, 3, 2},{4, 1, 5}},
                                {{0, 1, 5},{0, 5, 6}, {5, 3, 6},{0, 6, 2}},
                                {{0, 1, 7},{1, 3, 6}, {6, 2, 7},{7, 1, 6}},
                                {{0, 4, 7},{4, 7, 3}, {4, 1, 3},{2, 7, 3}}
                        };
                int id = hangingIDs[face][0] - 4;
                if((hangingIDs[face][0] == 4) and (hangingIDs[face][1] == 7)){
                    id = 3;
                }
                for(int i = 0; i < NUM_TETS; i++){
                    const std::array<int, TET_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][TET_CONFIGURATIONS[id][i][0]],
                             FACE[face][TET_CONFIGURATIONS[id][i][1]],
                             FACE[face][TET_CONFIGURATIONS[id][i][2]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TET_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::TET, nodeCoords));
                }

            }
            else if(configuration == CONFIGURATIONS::THREE_NODE_HANGING){
                int id = -1;
                if( (hangingIDs[face][0] == 4) and (hangingIDs[face][1] == 5) and (hangingIDs[face][2] == 6)){
                    id = 0;
                }
                else if((hangingIDs[face][0] == 5) and (hangingIDs[face][1] == 6) and (hangingIDs[face][0] == 7)){
                    id = 1;
                }
                else if((hangingIDs[face][0] == 4) and (hangingIDs[face][1] == 5) and (hangingIDs[face][0] == 7)){
                    id = 2;
                }
                else if((hangingIDs[face][0] == 4) and (hangingIDs[face][1] == 6) and (hangingIDs[face][0] == 7)){
                    id = 3;
                }
                else{
                    throw std::runtime_error(" Unreachable statement");
                }

                static constexpr int NUM_TETS = 3;
                static constexpr int TET_COORDS_ON_FACE = 3;
                numElements += NUM_TETS;
                m_elemInfo.reserve(numElements);
                static constexpr int TET_CONFIGURATIONS[4][NUM_TETS][TET_COORDS_ON_FACE] =
                        {
                                {{4, 1, 5},{4, 5, 6}, {5, 3, 6}}, // (4,5,6)
                                {{5, 3, 6},{5, 6, 7}, {7, 6, 2}}, // (5,6,7)
                                {{0, 4, 7},{4, 5, 7}, {4, 1, 5}}, // (4,5,7)
                                {{0, 4, 7},{4, 6, 7}, {7, 6, 2}}, // (4,6,7)
                        };
                for(int i = 0; i < NUM_TETS; i++){
                    const std::array<int, TET_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][TET_CONFIGURATIONS[id][i][0]],
                             FACE[face][TET_CONFIGURATIONS[id][i][1]],
                             FACE[face][TET_CONFIGURATIONS[id][i][2]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TET_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::TET, nodeCoords));
                }
                static constexpr int NUM_PYRAMIDS = 1;
                static constexpr int PYRAMID_COORDS_ON_FACE = 4;
                numElements += NUM_PYRAMIDS;
                m_elemInfo.reserve(numElements);
                static constexpr int PYRAMID_CONFIGURATIONS[4][NUM_PYRAMIDS][PYRAMID_COORDS_ON_FACE] =
                        {
                                {{0, 4, 2, 6}},
                                {{0, 1, 5, 7}},
                                {{7, 5, 2, 3}},
                                {{4, 1, 6, 3}}
                        };

                for(int i = 0; i < NUM_PYRAMIDS; i++){
                    const std::array<int, PYRAMID_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][PYRAMID_CONFIGURATIONS[id][i][0]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][1]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][2]],
                             FACE[face][PYRAMID_CONFIGURATIONS[id][i][3]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<PYRAMID_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::PYRAMID, nodeCoords));
                }
            }
            else if(configuration == CONFIGURATIONS::FOUR_NODE_HANGING){
                static constexpr int NUM_TETS = 4;
                static constexpr int TET_COORDS_ON_FACE = 3;
                numElements += NUM_TETS;
                m_elemInfo.reserve(numElements);
                static constexpr int TET_CONFIGURATIONS[1][NUM_TETS][TET_COORDS_ON_FACE] =
                        {
                                {{0, 4, 7},{4, 1, 5}, {5, 3, 6},{7, 6, 2}},
                        };
                for(int i = 0; i < NUM_TETS; i++){
                    const std::array<int, TET_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][TET_CONFIGURATIONS[0][i][0]],
                             FACE[face][TET_CONFIGURATIONS[0][i][1]],
                             FACE[face][TET_CONFIGURATIONS[0][i][2]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TET_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::TET, nodeCoords));
                }

                static constexpr int NUM_PYRAMIDS = 1;
                static constexpr int PYRAMID_COORDS_ON_FACE = 4;
                numElements += NUM_PYRAMIDS;
                m_elemInfo.reserve(numElements);
                static constexpr int PYRAMID_CONFIGURATIONS[1][NUM_PYRAMIDS][PYRAMID_COORDS_ON_FACE] =
                        {
                                {{4, 5, 7, 6}}
                        };
                for(int i = 0; i < NUM_PYRAMIDS; i++){
                    const std::array<int, PYRAMID_COORDS_PER_ELEMENT> nodeBits =
                            {FACE[face][PYRAMID_CONFIGURATIONS[0][i][0]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][1]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][2]],
                             FACE[face][PYRAMID_CONFIGURATIONS[0][i][3]],
                             MIDDLE_NODE};

                    std::array<int, HEX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<PYRAMID_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_3D::PYRAMID, nodeCoords));
                }
            }
            else{
                throw std::runtime_error("No matching configuration on face found");
            }
        }
    }

    return numElements;
}

MeshElements::~MeshElements() {
    // Not needed as GRID also clears up the elements
}

ElementType_3D MeshElements::getFEMElemObject( const int id, const double * physicalCoords){
    const auto &elem = m_elemInfo[id];
    const auto &elemType = std::get<0>(elem);
    const auto &nodeCoords = std::get<1>(elem);

    switch (elemType) {

        case HEX:
            m_sync.syncCoords<1, static_cast<TALYFEMLIB::ElemType>(ElementType_3D::HEX)>(physicalCoords,
                                                                                              &m_hex_grid,
                                                                                              nodeCoords.data());
            break;

        case TET:
            m_sync.syncCoords<1, static_cast<TALYFEMLIB::ElemType>(ElementType_3D::TET)>(physicalCoords,
                                                                                          &m_tet_grid,
                                                                                          nodeCoords.data());
            break;

        case PYRAMID:
            m_sync.syncCoords<1, static_cast<TALYFEMLIB::ElemType>(ElementType_3D::PYRAMID)>(physicalCoords,
                                                                                         &m_pyramid_grid,
                                                                                         nodeCoords.data());
            break;
        default:
            throw std::runtime_error("Unreachable statement");

    }
    return elemType;
}