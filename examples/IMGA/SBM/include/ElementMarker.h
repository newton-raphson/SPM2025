//
// Created by maksbh on 11/6/22.
//

#ifndef DENDRITEKT_ELEMENTMARKER_H
#define DENDRITEKT_ELEMENTMARKER_H
#include <Traversal/Traversal.h>
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
class SBMMarker : public Traversal {

    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    SubDomainBoundary * boundary_;
public:
    SBMMarker(DA * octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domain,
              std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
              const SubDomain * subDomain);
    void traverseOperation(TALYFEMLIB::FEMElm & fe) override;

};

SBMMarker::SBMMarker(DA * octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domain,
                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
                     const SubDomain * subDomain)
        : Traversal(octDA, treePart, domain),eleMarkers_(eleMarkers),subdomain_(subDomain){
    std::bitset<ElementMarker::MAX_ELMENT_TYPE> init;
    init.reset();
    eleMarkers_.resize(octDA->getLocalElementSz(),init);
    this->traverse();

}

void SBMMarker::traverseOperation(TALYFEMLIB::FEMElm & fe) {
    int outPts = 0;

    int npe = m_octDA->getNumNodesPerElement();
    const auto coords = this->m_coords;
    for(int i = 0; i < npe; i++){
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
        int outGP = 0;
        fe.refill(0,0);
        while (fe.next_itg_pt()){
            if (subdomain_->functionToRetainPhysNodes(fe.position().data()) == ibm::Partition::OUT) {
                outGP++;
            }

        }
        if(outGP < 2){
            eleMarkers_[localElemID].set(ElementMarker::SBM_FALSE_INTERCEPTED, true); // Cut cell
        }
    }




    localElemID++;
}


#endif //DENDRITEKT_ELEMENTMARKER_H
