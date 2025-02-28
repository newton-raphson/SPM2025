//
// Created by maksbh on 12/5/22.
//

#ifndef DENDRITEKT_CHECKSURFACE_H
#define DENDRITEKT_CHECKSURFACE_H
#include "Traversal/Traversal.h"
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"

class CheckSurface:public Traversal{

    static constexpr  int numFaces = 2*DIM;
    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    std::vector<std::bitset<numFaces>> faceMarker_;
    SubDomainBoundary * subDomainBoundary_;

public:
    CheckSurface(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
                std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
                const SubDomain * subDomain);
    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    void correctCycles();
};

CheckSurface::CheckSurface(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
        std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, const SubDomain * subDomain)
        :eleMarkers_(eleMarkers), Traversal(octDA,treePart,v,domain),subdomain_(subDomain){

    std::bitset<numFaces> mark;
    mark.reset();

    faceMarker_.resize(octDA->getLocalElementSz());
    this->traverse();
}

void CheckSurface::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
#if (DIM == 3)
    static constexpr DENDRITE_UINT cornerMap[2][8]{{0,1,2,3,4,5,6,7},{0,2,6,8,18,20,24,26}};
  static constexpr DENDRITE_UINT numNodesPerFace = 4;
  static constexpr DENDRITE_UINT faceID[2*DIM][4]
  {
      {0,2,4,6}, // Left
      {1,3,5,7}, // Right
      {0,1,4,5}, // Bottom
      {2,3,6,7}, // Top
      {0,1,2,3}, // Back
      {4,5,6,7},  // Front
  };
#elif(DIM == 2)
    static constexpr DENDRITE_UINT cornerMap[2][numFaces]{{0,1,2,3},{0,2,6,8}};
    static constexpr DENDRITE_UINT numNodesPerFace = 2;
    static constexpr DENDRITE_UINT faceID[2*DIM][2]
            {
                    {2,0}, // Left
                    {1,3}, // Right
                    {0,1}, // Bottom
                    {3,2}, // Top
            };

#else
    throw std::runtime_error("Not implemented\n");
#endif
    const int eleOrder = m_octDA->getElementOrder();
    // We need to eliminate only intercepted elements.
    // Only Neighbors of intercepted cannot have cycles unless you are doing something stupid with very poor resolution.
    if(eleMarkers_[localElemID].test(ElementMarker::INTERCEPTED_ELEMENT)){
        std::vector<TALYFEMLIB::ZEROPTV> coords(m_octDA->getNumNodesPerElement());
        std::vector<bool> isBoundary(m_octDA->getNumNodesPerElement(),false);
        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            std::memcpy(coords[i].data(), &this->m_coords[i*DIM],sizeof(double)*DIM);
        }
        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            unsigned int id;
            subDomainBoundary_->generateBoundaryFlags(coords[i],id);
        }

        for(int i = 0; i < numFaces; i++) {
            bool isBoundaryFaceinterceptedElement = true;
            for(int j = 0; j < numNodesPerFace; j++) {
                unsigned int cornerNodesOnFace = faceID[i][j];
                unsigned int nodeID = cornerMap[eleOrder - 1][cornerNodesOnFace];
                isBoundaryFaceinterceptedElement = isBoundaryFaceinterceptedElement and isBoundary[nodeID];
            }

            bool SBMCheck = true;
            for (int j = 0; j < numNodesPerFace; j++) {
                SBMCheck = SBMCheck and ( (subdomain_->functionToRetainPhysNodes(&this->m_coords[cornerMap[eleOrder - 1][faceID[i][j]]*DIM]) == ibm::IN)// This node belong to inactive region
                                          or ((values[cornerMap[eleOrder - 1][faceID[i][j]]] == NodeMarker::SBM_FALSE_INTERCEPTED_NODES)));
            }

            bool isBoundaryFace = isBoundaryFaceinterceptedElement or SBMCheck;
            if(isBoundaryFace) {
                faceMarker_[localElemID].set(i,true);
            }
        }

    }
    localElemID++;


}

void CheckSurface::correctCycles(){
    int localElemSize = this->m_octDA->getLocalElementSz();
    for(int i = 0; i < localElemSize; i++){
        for(int nFace = 0; nFace < DIM; nFace++){
            if((faceMarker_[i].test(2*nFace)) and (faceMarker_[i].test(2*nFace + 1))){
                eleMarkers_[i].reset();
                eleMarkers_[i].set(ElementMarker::SBM_FALSE_INTERCEPTED);
            }
        }
    }
}





#endif //DENDRITEKT_CHECKSURFACE_H
