//
// Created by maksbh on 11/6/22.
//

#pragma once
#include <Traversal/Traversal.h>
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
#include "NSInputData.h"

class BFS : public Traversal {

    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    NSInputData *idata_;
    const DENDRITE_UINT eleOrder_;
    SubDomainBoundary *m_subdomainBoundary;
    const IMGA *imga_;


public:
    BFS(DA * octDA, const std::vector<TREENODE> &treePart, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
        const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain, NSInputData *idata
            , SubDomainBoundary *subDomainBoundary, const IMGA *imga);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
};

BFS::BFS(DA * octDA, const std::vector<TREENODE> &treePart,std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
         const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain, NSInputData *idata, SubDomainBoundary *subDomainBoundary, const IMGA *imga)
        : Traversal(octDA, treePart,v, domain),eleMarkers_(eleMarkers),subdomain_(subDomain),idata_(idata),eleOrder_(octDA->getElementOrder()),m_subdomainBoundary(subDomainBoundary),imga_(imga)
{
    this->traverse();
}

void BFS::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
    int npe = m_octDA->getNumNodesPerElement();

    /// Element marker for In Out Intercepted and False Intercepted must be set before hand

    bool isNeighbor = std::any_of(values,values+npe,[](int i){return i == NodeMarker::SBM_FALSE_INTERCEPTED_NODES;});
    if(isNeighbor){
//        std::cout << "Element is neighbor" << std::endl;
        eleMarkers_[localElemID].set(ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED, true);
    }

    localElemID++;
}
