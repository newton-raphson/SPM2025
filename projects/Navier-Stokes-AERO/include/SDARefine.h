//
// Created by maksbh on 11/8/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>

class SDARefine: public Refinement{
  const DENDRITE_UINT maxLevel_;
public:
    NSInputData *inputData_;

    SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, DENDRITE_UINT maxLevel,NSInputData *inputData);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

SDARefine::SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const DENDRITE_UINT maxLevel,NSInputData *inputData)
:Refinement(da,treePart,domainInfo),maxLevel_(maxLevel),inputData_(inputData){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
    const double eps = 1e-13;
    unsigned int levelForRect = this->m_level;
    unsigned int levelForBd = this->m_level;

    for (const auto &p : coords)
    {
        if (p.y() - inputData_->rectLowerLeft[1] > eps and p.y() - inputData_->rectUpperRight[1] < eps and
            p.x() - inputData_->rectLowerLeft[0] > eps and p.x() - inputData_->rectUpperRight[0] < eps and
        this->m_level < inputData_->rect_lvl)
        {
            levelForRect = inputData_->rect_lvl;
            break;
        }

        if (p.y() - inputData_->rectLowerLeft[1] > eps and p.y() - inputData_->rectUpperRight[1] < eps and
            p.x() - inputData_->rectLowerLeft[0] > eps and p.x() - inputData_->rectUpperRight[0] < eps and
            this->m_level < inputData_->refine_lvl_bd) {

            levelForBd = inputData_->refine_lvl_bd;
            break;
        }
    }


    if(((this->m_BoundaryOctant) and (this->m_level < levelForBd)) or
            (this->m_level < levelForRect)){
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }


    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}

//ot::OCT_FLAGS::Refine SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
//    if((this->m_BoundaryOctant) and (this->m_level < maxLevel_)){
//        return ot::OCT_FLAGS::Refine::OCT_REFINE;
//    }
//    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
//}







#endif //DENDRITEKT_SDAREFINE_H
