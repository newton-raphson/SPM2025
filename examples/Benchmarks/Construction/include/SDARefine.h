//
// Created by maksbh on 6/18/21.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>
#include <random>
#include <BenchDataTypes.h>
class Refine: public Refinement{
  const DENDRITE_UINT boundaryLevel;
  const DENDRITE_UINT baseLevel;
  const bool allRefine_ = false;
  const DENDRITE_REAL fractionRefine_;
  std::mt19937_64 rng;
  std::uniform_real_distribution<double>  unif;
  const CaseType caseType_;
  const SubDomain * subDomain_;
public:
  Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const CaseType caseType, const SubDomain * subDomain,const DENDRITE_UINT blevel, const DENDRITE_UINT baseLevel, const DENDRITE_REAL fracRefine, const bool allRefine = false);

  ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

Refine::Refine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo,const CaseType caseType, const SubDomain * subDomain,const DENDRITE_UINT blevel, const DENDRITE_UINT baselevel, const DENDRITE_REAL fracRefine, const bool allRefine)
  :Refinement(da,treePart,domainInfo),boundaryLevel(blevel),baseLevel(baselevel),allRefine_(allRefine),caseType_(caseType),fractionRefine_(fracRefine),unif(0,1),subDomain_(subDomain){
  uint64_t timeSeed = 731 * da->getRankActive();
  std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
  rng.seed(ss);
  unif.reset();
  this->initRefinement();
}

ot::OCT_FLAGS::Refine Refine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){

  if(caseType_ == CaseType::CHANNEL) {
    if (allRefine_) {
      return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    if (unif(rng) < fractionRefine_) {
      return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }
  else{
    const auto & spheres = subDomain_->getVoxelSpheres();
    for(const auto & sphere:spheres){
      const auto&center= sphere.center;
      for(const auto  & coord:coords){
        double dist = 0;
        for(int d = 0; d < DIM; d++){
          dist +=(coord[d] - center[d])*(coord[d] - center[d]);
        }
        dist= sqrt(dist);
        if((dist < 0.9) and (this->m_level < boundaryLevel)){
          return ot::OCT_FLAGS::Refine::OCT_REFINE;
        }
      }
    }
    return ot::OCT_FLAGS::OCT_NO_CHANGE;
  }
}

#endif //DENDRITEKT_SDAREFINE_H
