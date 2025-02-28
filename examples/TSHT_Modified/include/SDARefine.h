//
// Created by maksbh on 10/29/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>
class SDARefine : public Refinement {
    int count = -100;
    const int bLevel_;
public:
    SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const int blevel);
    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

SDARefine::SDARefine(DA *da, const std::vector<TREENODE> &treePart, const DomainExtents &domainInfo, const int blevel)
:Refinement(da,treePart,domainInfo),bLevel_(blevel){
    this->traverse();
}

ot::OCT_FLAGS::Refine
SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
    if(this->m_BoundaryOctant and (bLevel_ > this->m_level)){
//        count++;
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    else{
//        count++;
        return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
    }
}
#endif //DENDRITEKT_SDAREFINE_H
