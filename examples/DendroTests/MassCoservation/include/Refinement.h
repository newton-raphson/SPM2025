//
// Created by maksbh on 1/27/23.
//
#pragma once
#include <Traversal/Refinement.h>
#include "DataTypes.h"


class MRefine: public Refinement{
    std::vector<double> & coarseFlag_;
    int elemID = 0;
public:
   MRefine(DA *da, const std::vector<TREENODE> & treenode, const DomainExtents & domainExtents, std::vector<double> & coarseFlag);

   ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords) override;

};

MRefine::MRefine(DA *da, const std::vector<TREENODE> & treenode, const DomainExtents & domainExtents, std::vector<double> & coarseFlag)
:Refinement(da,treenode,domainExtents),coarseFlag_(coarseFlag){
    coarseFlag.resize(da->getLocalElementSz(),1.0);
    this->initRefinement();
}

ot::OCT_FLAGS::Refine MRefine::getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords){
    ot::OCT_FLAGS::Refine flag;
    if(this->m_BoundaryOctant) {
        flag = ot::OCT_FLAGS::OCT_COARSEN;
    }
    else{
        flag = ot::OCT_FLAGS::OCT_COARSEN;
    }
    elemID++;
    return flag;
}



