//
// Created by maksbh on 3/26/23.
//

#include "Mesh/Mesh2D.h"
#include "DataTypes.h"
#include "talyfem/grid/elem_types/elem2dbox.h"
#include "talyfem/grid/elem_types/elem2dtriangle.h"

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
    // Quad element
    {
        static constexpr DENDRITE_UINT numNodes = 4;
        int node_id_array[numNodes];
        m_quad_grid.redimArrays(numNodes, 1);
        for (int i = 0; i < numNodes; i++) {
            m_quad_grid.node_array_[i] = new TALYFEMLIB::NODE();
            node_id_array[i] = i;
        }

        m_quad_elem = new TALYFEMLIB::ELEM2dBox();
        m_quad_grid.elm_array_[0] = m_quad_elem;


        m_quad_elem->redim(numNodes, node_id_array);
    }
    // Tri Element

    {
        static constexpr DENDRITE_UINT numNodes = 3;
        int node_id_array[numNodes];
        m_triangle_grid.redimArrays(numNodes, 1);
        for (int i = 0; i < numNodes; i++) {
            m_triangle_grid.node_array_[i] = new TALYFEMLIB::NODE();
            node_id_array[i] = i;
        }

        m_triangle_elem = new TALYFEMLIB::ELEM2dTriangle();
        m_triangle_grid.elm_array_[0] = m_triangle_elem;

        m_triangle_elem->redim(numNodes, node_id_array);
    }


}

int MeshElements::initElementType(const std::bitset<numBits> &val) {
    bits_ = val;

//    m_bitsToNodeMap.fill(-1);
//    int counter = 0;
//    for (int i = 0; i < numBits; i++) {
//        if (bits_.test(i)) {
//            m_bitsToNodeMap[i] = counter;
//            counter++;
//        }
//
//    }
    int numElements = 0;
    m_elemInfo.clear();
    if (bits_.test(MIDDLE_NODE) == false) {
        // No Hanging nodes in this element
        numElements += 1;
        m_elemInfo.resize(numElements);
        const std::array<int, QUAD_COORDS_PER_ELEMENT> nodeBits = {0, 2, 6, 8};
        std::array<int, QUAD_COORDS_PER_ELEMENT> nodeCoords{};
        this->createElemInfo<QUAD_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
        m_elemInfo[0] = std::make_tuple(ElementType_2D::QUAD, nodeCoords);

        // No Hanging nodes in this element
//        numElements += 2;
//        m_elemInfo.resize(numElements);
//        {
//            const std::array<int, TRI_COORDS_PER_ELEMENT> nodeBits = {0, 2, 6};
//            std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
//            this->createElemInfo<TRI_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
//            m_elemInfo[0] = std::make_tuple(ElementType_2D::TRIANGLE, nodeCoords);
//        }
//        {
//            const std::array<int, TRI_COORDS_PER_ELEMENT> nodeBits = {2, 6, 8};
//            std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
//            this->createElemInfo<TRI_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
//            m_elemInfo[1] = std::make_tuple(ElementType_2D::TRIANGLE, nodeCoords);
//        }

    } else {

        // Hanging nodes is present in this element
        for (int face = 0; face < NUM_FACE; face++) {
            if (bits_.test(FACE[face][1])) {
                numElements += 2; // We split the face into 2 elements
                { // Element No. 1
                    const std::array<int, TRI_COORDS_PER_ELEMENT> nodeBits = {FACE[face][0], MIDDLE_NODE,
                                                                              FACE[face][1]};
                    std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TRI_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_2D::TRIANGLE, nodeCoords));
                }
                {
                    // Element No. 2
                    const std::array<int, TRI_COORDS_PER_ELEMENT> nodeBits = {FACE[face][1], MIDDLE_NODE,
                                                                              FACE[face][2]};
                    std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TRI_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_2D::TRIANGLE, nodeCoords));
                }
            } else {
                numElements += 1;
                // We create one triangle with middle node and coords on the faces
                {
                    // Element No. 2
                    const std::array<int, TRI_COORDS_PER_ELEMENT> nodeBits = {FACE[face][0], MIDDLE_NODE,
                                                                              FACE[face][2]};
                    std::array<int, MAX_COORDS_PER_ELEMENT> nodeCoords{};
                    this->createElemInfo<TRI_COORDS_PER_ELEMENT>(nodeBits, nodeCoords);
                    m_elemInfo.push_back(std::make_tuple(ElementType_2D::TRIANGLE, nodeCoords));
                }

            }

        }
    }
    return numElements;
}

ElementType_2D MeshElements::getFEMElemObject(const int id, const double *physicalCoords) {
    const auto &elem = m_elemInfo[id];
    const auto &elemType = std::get<0>(elem);
    const auto &nodeCoords = std::get<1>(elem);

    switch (elemType) {

        case TRIANGLE:
            m_sync.syncCoords<1, static_cast<TALYFEMLIB::ElemType>(ElementType_2D::TRIANGLE)>(physicalCoords,
                                                                                              &m_triangle_grid,
                                                                                              nodeCoords.data());
            break;

        case QUAD:
            m_sync.syncCoords<1, static_cast<TALYFEMLIB::ElemType>(ElementType_2D::QUAD)>(physicalCoords,
                                                                                              &m_quad_grid,
                                                                                              nodeCoords.data());
            break;

        default:
            throw std::runtime_error("Unreachable statement");

    }
    return elemType;
}


void MeshElements::generateFaces(SubDomainBoundary * subDomainBoundary, const int elementID, const double * physCoords,
                                 std::vector<std::tuple<std::array<int,2>,unsigned int>>& boundaryFaceID){

    const auto &elem = m_elemInfo[elementID];
    const auto &elemType = std::get<0>(elem);
    const auto &nodeCoords = std::get<1>(elem);
    switch (elemType) {

        case TRIANGLE: {
            static constexpr int NUMFACE = 3;
            static constexpr int NUM_COORDS_PER_FACE = 2;
            static constexpr int BOUNDARY_FACE[NUMFACE][NUM_COORDS_PER_FACE] = {{0, 1},
                                                                                {1, 2},
                                                                                {2, 0}};
            for (int numFace = 0; numFace < NUMFACE; numFace++) {
                unsigned int boundaryType;
                unsigned int boundaryID;
                int id0 = nodeCoords[BOUNDARY_FACE[numFace][0]];
                int id1 = nodeCoords[BOUNDARY_FACE[numFace][1]];
                double midCoord_X = 0.5*(physCoords[id0 * DIM + 0] +  physCoords[id1 * DIM + 0]);
                double midCoord_Y = 0.5*(physCoords[id0 * DIM + 1] +  physCoords[id1 * DIM + 1]);
                TALYFEMLIB::ZEROPTV pos(midCoord_X,midCoord_Y,0);
                subDomainBoundary->generateBoundaryFlags(pos,boundaryID);
                for(int i = 0; i < BoundaryTypes::VOXEL::MAX_BOUNDARY_TYPE; i++){
                    if(subDomainBoundary->checkBoundaryType(i)){
                        boundaryType = i;
                        std::array<int, 2> val = {id0, id1};
                        boundaryFaceID.push_back(std::make_tuple(val, boundaryType));
                        break;
                    }
                }

            }

            break;
        }

        case QUAD: {
            static constexpr int NUMFACE = 4;
            static constexpr int NUM_COORDS_PER_FACE = 2;
            static constexpr int BOUNDARY_FACE[NUMFACE][NUM_COORDS_PER_FACE] = {{0, 1},
                                                                                {1, 2},
                                                                                {2, 3},
                                                                                {3, 0}};
            for (int numFace = 0; numFace < NUMFACE; numFace++) {
                unsigned int boundaryType[NUM_COORDS_PER_FACE];
                unsigned int boundaryID[NUM_COORDS_PER_FACE];

                for (int numCoords = 0; numCoords < NUM_COORDS_PER_FACE; numCoords++) {
                    int id = BOUNDARY_FACE[numFace][numCoords];
                    TALYFEMLIB::ZEROPTV pos(physCoords[nodeCoords[id] * DIM], physCoords[nodeCoords[id] * DIM + 1],0);
                    subDomainBoundary->generateBoundaryFlags(pos, boundaryID[numCoords]);
                    for(int i = 0; i < BoundaryTypes::VOXEL::MAX_BOUNDARY_TYPE; i++){
                        if(subDomainBoundary->checkBoundaryType(i)){
                            boundaryType[numCoords] = i;
                        }
                    }
                }
                if ((boundaryType[0] == boundaryType[1]) and (boundaryID[0] == boundaryID[1])) {
                    int id0 = nodeCoords[BOUNDARY_FACE[numFace][0]];
                    int id1 = nodeCoords[BOUNDARY_FACE[numFace][1]];
                    std::array<int, 2> val = {id0, id1};
                    boundaryFaceID.push_back(std::make_tuple(val, boundaryType[0]));
                }
            }
            break;
        }
        default:
            throw std::runtime_error("Unreachable statement");

    }

}

