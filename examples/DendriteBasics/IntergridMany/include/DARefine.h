//
// Created by maksbh on 2/28/22.
//

#ifndef DENDRITEKT_DAREFINE_H
#define DENDRITEKT_DAREFINE_H
#include <Traversal/Refinement.h>

class DARefine: public Refinement{

public:


  DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo);

  int getNumLevelRefine(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

DARefine::DARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo)
  :Refinement(da,treePart,domainInfo,true){
  this->initRefinement();
}

int DARefine::getNumLevelRefine(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
  using namespace TALYFEMLIB;
  if(this->m_BoundaryOctant){
    return 2;
  }
  else{
    return 0;
  }
}
#endif //DENDRITEKT_DAREFINE_H
