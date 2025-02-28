//
// Created by maksbh on 11/6/22.
//

#pragma once
#include <Traversal/Traversal.h>
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
#include "util.h"

class BFS : public Traversal {

    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    const LEInputData *idata_;

public:
    BFS(DA * octDA, const std::vector<TREENODE> &treePart, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
        const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain, const LEInputData *idata);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
};

BFS::BFS(DA * octDA, const std::vector<TREENODE> &treePart,std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
         const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain, const LEInputData *idata)
        : Traversal(octDA, treePart,v, domain),eleMarkers_(eleMarkers),subdomain_(subDomain),idata_(idata){
    std::bitset<ElementMarker::MAX_ELMENT_TYPE> init;
    init.reset();
    eleMarkers_.resize(octDA->getLocalElementSz(),init);
    this->traverse();

}

void BFS::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) {
    /// out means ACTIVE here
    int outPts = 0;

    int npe = m_octDA->getNumNodesPerElement();
    const auto coords = this->m_coords;
    for(int i = 0; i < npe; i++) {
        if (subdomain_->functionToRetainPhysNodes(&coords[DIM*i]) == ibm::Partition::OUT) {
            outPts++;
        }
    }

    if(outPts == npe){
        eleMarkers_[localElemID].set(ElementMarker::OUT_ELEMENT,true);
    }
    else if(outPts == 0){
        eleMarkers_[localElemID].set(ElementMarker::IN_ELEMENT, true);
    }
    else {
        eleMarkers_[localElemID].set(ElementMarker::INTERCEPTED_ELEMENT, true);
        int ActiveGP = 0;
        fe.refill(0, idata_->RelOrderCheckActive);
        while (fe.next_itg_pt()) {
            if (subdomain_->functionToRetainPhysNodes(fe.position().data()) == ibm::Partition::OUT) {
                ActiveGP++; // active GPs
            }
        }

        int InActiveGP = fe.n_itg_pts() - ActiveGP;
        //std::cout<<"fe.n_itg_pts() = " << fe.n_itg_pts() << "\n";
        //std::cout<<"ActiveGP = " << ActiveGP << "\n";

        //if ((double)ActiveGP/fe.n_itg_pts() < idata_->RatioGPSBM) {
        if ((1-(double)ActiveGP/fe.n_itg_pts()) > idata_->RatioGPSBM) {
            eleMarkers_[localElemID].set(ElementMarker::SBM_FALSE_INTERCEPTED, true); // Cut cell
        }
    }

    bool isNeighbor = std::any_of(values,values+npe,[](int i){return i == NodeMarker::SBM_FALSE_INTERCEPTED_NODES;});
    if(isNeighbor){
        eleMarkers_[localElemID].set(ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED, true);
    }

    localElemID++;
}

