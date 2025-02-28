//
// Created by maksbh on 2/11/23.
//

#ifndef DENDRITEKT_CHREFINE_H
#define DENDRITEKT_CHREFINE_H

#include "Traversal/Refinement.h"
#include "CHInputData.h"

template<RefinementStage stage>
class CHRefine:public Refinement{
    const DENDRITE_UINT ndof_;
    const CHInputData * chInputData_;
    const bool coarseningTest_ = false;
public:
    CHRefine(DA * octDA, DistTREE & dtree, VecInfo & vec, DomainExtents domainExtents, const CHInputData * chInputData, const bool coarseningTest);
    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords, const PetscScalar * values) override;

};
template<RefinementStage stage>
CHRefine<stage>::CHRefine(DA * octDA, DistTREE & dtree, VecInfo & vec, DomainExtents domainExtents, const CHInputData * chInputData,const bool coarseningTest)
:Refinement(octDA,dtree.getTreePartFiltered(),vec,domainExtents,false),ndof_(vec.ndof),chInputData_(chInputData),coarseningTest_(coarseningTest){
    this->traverse();

}
template<RefinementStage stage>
ot::OCT_FLAGS::Refine CHRefine<stage>::getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords, const PetscScalar * values){
    if(coarseningTest_){
        return ot::OCT_FLAGS::OCT_COARSEN;
    }
    const int nPe = m_octDA->getNumNodesPerElement();
    std::vector<double> chValues(nPe);
    for(int i = 0; i < nPe; i++){
        chValues[i] = values[i*ndof_ + 0];
    }
    int reqRefineLevel = chInputData_->mesh_def.refine_lvl_base;
    if(std::any_of(chValues.begin(),chValues.end(),[&](double phi){
        return (fabs(phi) < chInputData_->mesh_def.refine_tol);
    }) ){
        reqRefineLevel = chInputData_->mesh_def.refine_lvl_interface;
    }
    if(stage == RefinementStage::REFINE){
        if(reqRefineLevel > this->m_level){
            return ot::OCT_FLAGS::OCT_REFINE;
        }
        else{
            return ot::OCT_FLAGS::OCT_NO_CHANGE;
        };
    }
    if(stage == RefinementStage::COARSEN){
        if(reqRefineLevel < this->m_level){
            return ot::OCT_FLAGS::OCT_COARSEN;
        }
        else{
            return ot::OCT_FLAGS::OCT_NO_CHANGE;
        }
    }
}
#endif //DENDRITEKT_CHREFINE_H
