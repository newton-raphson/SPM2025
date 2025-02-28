//
// Created by maksbh on 8/12/21.
//

#ifndef DENDRITEKT_PREFINE_H
#define DENDRITEKT_PREFINE_H
#include <Traversal/Refinement.h>
#include <SubDA/Voxel.h>

class PRefine : public Refinement {


public:
  PRefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;


};

PRefine::PRefine(DA *da, const std::vector<TREENODE> & treePart,const DomainExtents &domainInfo)
  : Refinement(da, treePart,domainInfo) {

  this->traverse();
}

ot::OCT_FLAGS::Refine
PRefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {
  if(this->m_BoundaryOctant){
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  }
  else{
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }
}

#endif //DENDRITEKT_PREFINE_H
