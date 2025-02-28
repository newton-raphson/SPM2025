//
// Created by maksbh on 11/6/22.
//

#pragma once
#include <Traversal/Traversal.h>
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
class BFS : public Traversal {

    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
public:
    BFS(DA * octDA, const std::vector<TREENODE> &treePart, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
        const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
};

BFS::BFS(DA * octDA, const std::vector<TREENODE> &treePart,std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
         const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain)
        : Traversal(octDA, treePart,v, domain),eleMarkers_(eleMarkers),subdomain_(subDomain){

    this->traverse();

}

void BFS::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
// Element marker for In Out Intercepted and False Intercepted must be set before hand
    int npe = m_octDA->getNumNodesPerElement();

    bool isNeighbor = std::any_of(values,values+npe,[](int i){return i == NodeMarker::SBM_FALSE_INTERCEPTED_NODES;});
    if(isNeighbor){
        eleMarkers_[localElemID].set(ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED, true);
    }



    localElemID++;
}


