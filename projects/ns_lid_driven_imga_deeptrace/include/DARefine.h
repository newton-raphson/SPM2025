//
// Created by maksbh on 3/31/21.
//

#ifndef SBM_DAREFINE_H
#define SBM_DAREFINE_H

#include <Traversal/Refinement.h>
#include <Boundary/SubDomainBoundary.h>
class DARefine: public Refinement{
    const int refineLevel_;
    const bool onlyBoundary_;
public:
    DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,const int refineLevel,const bool onlyBoundary);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

DARefine::DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,const int refineLevel,const bool onlyBoundary)
        :Refinement(da,treePart,domainInfo),refineLevel_(refineLevel),onlyBoundary_(onlyBoundary){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine DARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
    if(onlyBoundary_){
        if(this->m_BoundaryOctant and this->m_level < refineLevel_) {
            return ot::OCT_FLAGS::Refine::OCT_REFINE;
        }
        else{
            return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
        }
    }
    else{
        if(this->m_level < refineLevel_) {
            return ot::OCT_FLAGS::Refine::OCT_REFINE;
        }
        else{
            return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
        }
    }

}

#endif //SBM_DAREFINE_H
