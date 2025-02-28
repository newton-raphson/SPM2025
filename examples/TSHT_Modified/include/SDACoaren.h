//
// Created by dhruv on 6/8/23.
//

#ifndef DENDRITEKT_SDACOAREN_H
#define DENDRITEKT_SDACOAREN_H
#include <Traversal/Refinement.h>
class SDACoarsen : public Refinement {
    int count = -100;
    const int bLevel_;
public:
    SDACoarsen(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const int blevel);
    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

SDACoarsen::SDACoarsen(DA *da, const std::vector<TREENODE> &treePart, const DomainExtents &domainInfo, const int blevel)
        :Refinement(da,treePart,domainInfo),bLevel_(blevel){
    this->traverse();
}

ot::OCT_FLAGS::Refine
SDACoarsen::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
    if(this->m_BoundaryOctant and (bLevel_ < this->m_level)){
//        count++;
        return ot::OCT_FLAGS::Refine::OCT_COARSEN;
    }
    else{
//        count++;
        return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
    }
}
#endif //DENDRITEKT_SDACOAREN_H
