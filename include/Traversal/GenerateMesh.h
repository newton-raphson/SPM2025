//
// Created by maksbh on 4/4/23.
//

#ifndef DENDRITEKT_GENERATEMESH_H
#define DENDRITEKT_GENERATEMESH_H
#ifdef ENABLE_2D
#include "Mesh/Mesh2D.h"
#include "OctToPhysical.h"

class GenerateMeshInfo : public MeshElements {



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
    SubDomainBoundary * m_subDomainBoundary;


public:
    GenerateMeshInfo(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents, SubDomainBoundary * subDomainBoundary);
    void generate();
    void generateCoords(const double *octCoords) ;

    void print(const std::string & filename_vertex,const std::string & filename_element,const std::string & filename_boundary);
    void printGMSH(const std::string & GMSH_filename);

    template<int size, typename T>
    void copyIDS(const T * nodeIDsFlat, const int ndof, const int * nodeCoords, std::vector<PetscInt> & vals){
        vals.resize(size);
        for (int node = 0; node < size; node++) {
            vals[node] = nodeIDsFlat[nodeCoords[node]];

        }
    }

    void addFaces(const std::vector<std::tuple<std::array<int,2>,unsigned int>> & boundaryFaceID, const ot::RankI * nodeIDsFlat);

};

GenerateMeshInfo::GenerateMeshInfo(DA * octDA,  const std::vector<TREENODE> & treePart, const DomainExtents & domainExtents,SubDomainBoundary * subDomainBoundary)
: MeshElements(1),m_octDA(octDA), m_octToPhysical(domainExtents), m_treePart(treePart),m_subDomainBoundary(subDomainBoundary){
    m_coords = new DENDRITE_REAL[intPow(3, DIM) * DIM];

}

void GenerateMeshInfo::generate() {
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

    ot::MatvecBase<DIM, RankI> treeloop(ghostedNodalSz, 1, m_octDA->getElementOrder(), odaCoords, &(*ghostedGlobalNodeId.cbegin()),
                                        &(*this->m_treePart.cbegin()),
                                        this->m_treePart.size(),
                                        *m_octDA->getTreePartFront(),
                                        *m_octDA->getTreePartBack());


    while (!treeloop.isFinished()) {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
            const double *nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            const RankI *nodeIdsFlat = treeloop.subtreeInfo().readNodeValsIn();
            const bool isBoundaryOctant = treeloop.subtreeInfo().isElementBoundary();
            const auto bits = treeloop.subtreeInfo().getLeafBitsetInfo();

            int numElements = this->initElementType(bits);
            this->generateCoords(nodeCoordsFlat);
            const auto & m_elemInfo = this->m_elemInfo;
            for(int i = 0; i < numElements; i++){
                boundaryFaceID.clear();
                this->generateFaces(m_subDomainBoundary,i,m_coords,boundaryFaceID);
                this->addFaces(boundaryFaceID,nodeIdsFlat);
                const auto & elem = m_elemInfo[i];
                const auto &elemType = std::get<0>(elem);
                const auto &nodeCoords = std::get<1>(elem);
                if(elemType == ElementType_2D::QUAD){
                    std::vector<PetscInt> nodes;
                    this->copyIDS<4>(nodeIdsFlat,1,nodeCoords.data(),nodes);
                    gridElements.push_back(std::make_tuple(GMSH_GRID::GMSH_QUAD,nodes)); /// (element type, connectivity)
                }
                else if(elemType == ElementType_2D::TRIANGLE){
                    std::vector<PetscInt> nodes;
                    this->copyIDS<3>(nodeIdsFlat,1,nodeCoords.data(),nodes);
                    gridElements.push_back(std::make_tuple(GMSH_GRID::GMSH_TRIANGLE,nodes));
                }
                else{
                    throw std::runtime_error("Unreachable statement");
                }
            }
            treeloop.next(1);

        } else {
            treeloop.step(1);
        }

    }

    int eleOrder = m_octDA->getElementOrder();
    auto tnCoords = m_octDA->getTNCoords();
    std::vector<double> physcoords( DIM, 0 );
    m_vertexPosition.resize(m_octDA->getLocalNodalSz());
    for( int idx = 0; idx < m_octDA->getLocalNodalSz(); idx++ ) {
        ot::treeNode2Physical(*tnCoords, eleOrder + 1, &(*physcoords.begin()));
        m_vertexPosition[idx] = {physcoords[0],physcoords[1]};
        tnCoords++;
    }

}

void GenerateMeshInfo::generateCoords(const double *octCoords) {
//    static std::array<double, numBits * DIM> coords;
    double scalingVal[3] = {0.0, 0.5, 1.0};
#ifdef ENABLE_2D
    double DX[DIM];
    DX[0] = octCoords[6] - octCoords[0];
    DX[1] = octCoords[7] - octCoords[1];
    for(int dy = 0; dy < 3; dy++){
        for(int dx = 0; dx < 3; dx++){
            m_coords[(dy*3+dx)*DIM + 0] = octCoords[0] + scalingVal[dx]*DX[0];
            m_coords[(dy*3+dx)*DIM + 1] = octCoords[1] + scalingVal[dy]*DX[1];
        }
    }
#else
    double DX[DIM];
    DX[0] = octCoords[21] - octCoords[0];
    DX[1] = octCoords[22] - octCoords[1];
    DX[2] = octCoords[23] - octCoords[2];
    for (int dz = 0; dz < 3; dz++) {
        for (int dy = 0; dy < 3; dy++) {
            for (int dx = 0; dx < 3; dx++) {
                m_coords[(dz * 9 + dy * 3 + dx) * DIM + 0] = octCoords[0] + scalingVal[dx] * DX[0];
                m_coords[(dz * 9 + dy * 3 + dx) * DIM + 1] = octCoords[1] + scalingVal[dy] * DX[1];
                m_coords[(dz * 9 + dy * 3 + dx) * DIM + 2] = octCoords[2] + scalingVal[dz] * DX[2];
            }
        }
    }
#endif
//    int counter = 0;
//    for (int i = 0; i < m_bits.size(); i++) {
//        if (m_bits.test(i)) {
//            std::memcpy(&m_coords[counter * DIM], &coords[i * DIM], sizeof(double) * DIM);
//            counter++;
//        }
//    }

}

void GenerateMeshInfo::print(const std::string & filename_vertex,const std::string & filename_element,const std::string & filename_boundary){
    {
        std::ofstream fout(filename_vertex.c_str());
        fout << m_vertexPosition.size() << "\n";
        for (int i = 0; i < m_vertexPosition.size(); i++) {
            const auto & vertex = m_vertexPosition[i];
            fout << i+1 << " " << vertex[0] << " " << vertex[1] << "\n";
        }
        fout.close();
    }
    {
        std::ofstream fout(filename_element.c_str());
        fout << gridElements.size() << "\n";
        for (int i = 0; i < gridElements.size(); i++) {
            const auto &grid =  gridElements[i];
            const auto & conn = std::get<1>(grid);
            fout << i + 1 << " " <<  std::get<0>(grid) << " " << 15 << " " << 15 << " ";
            for(const auto & nodeIDs: conn){
                fout << nodeIDs + 1 << " ";
            }
            fout << "\n";
        }
        fout.close();
    }
    {
        std::ofstream fout(filename_boundary.c_str());
        fout << boundary_elements.size() << "\n";
        for (int i = 0; i < boundary_elements.size(); i++) {
            const auto & bnd = boundary_elements[i];
            const auto & elemType = std::get<0>(bnd);
            const auto & boundaryID = std::get<1>(bnd);
            const auto & boundaryNodes = std::get<2>(bnd);
            fout << i + 1 << " " << DIM << " "<<  elemType << " " << boundaryID + 1<< " " << boundaryID + 1 << " ";
            for(const auto & nodeIDs: boundaryNodes){
                fout << nodeIDs + 1 << " ";
            }
            fout << "\n";
        }
        fout.close();
    }
}

void GenerateMeshInfo::printGMSH(const std::string & gmsh_filename){
    std::ofstream fout(gmsh_filename.c_str());
    fout << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"  ;

    // Nodes
    fout << "$Nodes\n";
    fout << m_vertexPosition.size() << "\n";
    for (int i = 0; i < m_vertexPosition.size(); i++) {
        const auto &vertex = m_vertexPosition[i];
        fout << i + 1 << " " << vertex[0] << " " << vertex[1] << " 0 \n";
    }

    fout << "$EndNodes\n";

    // Elements
    fout << "$Elements\n";
    int totalNumElements = gridElements.size() + boundary_elements.size();
    fout << totalNumElements << "\n";
    // Boundary elements
    for (int i = 0; i < boundary_elements.size(); i++) {
        const auto & bnd = boundary_elements[i];
        const auto & elemType = std::get<0>(bnd);
        const auto & boundaryID = std::get<1>(bnd);
        const auto & boundaryNodes = std::get<2>(bnd);
        fout << i + 1 << " "  << elemType << " " << DIM << " " << boundaryID + 1<< " " << boundaryID + 1 << " ";
        for(const auto & nodeIDs: boundaryNodes){
            fout << nodeIDs + 1 << " ";
        }
        fout << "\n";
    }

    // Domain Elements
    for (int i = 0; i < gridElements.size(); i++) {
        const auto &grid =  gridElements[i];
        auto  conn = std::get<1>(grid); /// connectivity
        const auto & elem = std::get<0>(grid); /// element type
        fout << boundary_elements.size() + i + 1 << " " <<   std::get<0>(grid) << " " << DIM << " " << 15 << " " << 15 << " ";
        if(elem == GMSH_QUAD){
//          std::cout << "Inside\n";
          std::swap(conn[2],conn[3]);
        }
        for(const auto & nodeIDs: conn){
            fout << nodeIDs + 1 << " ";
        }
        fout << "\n";
    }

    fout << "$EndElements\n";
    fout.close();
}

void GenerateMeshInfo::addFaces(const std::vector<std::tuple<std::array<int,2>,unsigned int>> & boundaryFaceID, const ot::RankI * nodeIDsFlat){
    for (const auto & bID: boundaryFaceID){
        std::vector<PetscInt> nodeIDs(2);
        const auto & nodeCoords =  std::get<0>(bID);
        for(int i = 0; i < nodeCoords.size(); i++){
            nodeIDs[i] = nodeIDsFlat[nodeCoords[i]];
        }

        const auto & boundaryType =  std::get<1>(bID);
        boundary_elements.push_back(std::make_tuple(GMSH_GRID::GMSH_LINE,boundaryType,nodeIDs));
    }
}


#endif
#endif //DENDRITEKT_GENERATEMESH_H
