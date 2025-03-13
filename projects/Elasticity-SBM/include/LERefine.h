#pragma once

#include <Traversal/Refinement.h>
#include <Boundary/SubDomainBoundary.h>
#include "LEInputData.h"

class LERefine : public Refinement {
  LEInputData *inputData_;
  const DomainExtents &domainExtents_;
  SubDomainBoundary *subDomainBoundary_;

 public:
  LERefine(DA *octDA,
             const std::vector<TREENODE> &treePart,
             const DomainExtents &domainExtents,
             LEInputData *inputData,
             SubDomainBoundary *subDomainBoundary1);

  virtual ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

  ~LERefine() {}

};

LERefine::LERefine(DA *octDA,
                       const std::vector<TREENODE> &treePart,
                       const DomainExtents &domainExtents,
                       LEInputData *inputData,
                       SubDomainBoundary *subDomainBoundary)
    : Refinement(octDA, treePart, domainExtents), inputData_(inputData), domainExtents_(domainExtents), subDomainBoundary_(subDomainBoundary) {
  this->traverse();
}

ot::OCT_FLAGS::Refine LERefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) {

  const DomainInfo &physDomain = domainExtents_.physicalDADomain;


  /// refine walls (maximum to the refine_h level)
  const double eps = 1e-13;
  unsigned int baselevel = inputData_->mesh_def.refine_lvl_base;


  DENDRITE_UINT id = -1;
  bool isObject = false;
  if (this->m_BoundaryOctant) {
    for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
      subDomainBoundary_->generateBoundaryFlags(coords[i], id);
      if (subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or subDomainBoundary_->checkBoundaryType(BoundaryTypes::FUNCTION)) {
          isObject = true;
        break;
      }
    }
  }

  unsigned int maxlevelForCarvedOutGeom = inputData_->mesh_def.refine_lvl_base;
  if (isObject) {
    maxlevelForCarvedOutGeom = std::max(inputData_->ibm_geom_def.at(id).refine_lvl, maxlevelForCarvedOutGeom);
  }

//  std::cout << "maxlevelForCarvedOutGeom = " << maxlevelForCarvedOutGeom << "\n";

  // the regional refine should not refine inside the geometry
//  if (insideGeo) {
//    maxlevelForRegion = inputData_->mesh_def.refine_l;
//  }
//  if (outsideRetain) {
//    levelForWall = inputData_->mesh_def.refine_l;
//  }
    if(this->m_level < baselevel or this->m_level < maxlevelForCarvedOutGeom){
        return ot::OCT_FLAGS::Refine::OCT_REFINE;
    }
    else {
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }

}
