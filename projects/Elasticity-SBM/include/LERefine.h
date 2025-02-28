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
  const unsigned int nodeNo = m_octDA->getNumNodesPerElement();
  const DENDRITE_UINT currLevel = this->m_level;

  /// refine walls (maximum to the refine_h level)
  const double eps = 1e-13;
  unsigned int levelForWall = inputData_->mesh_def.refine_lvl_base;
  const auto &bc = inputData_->boundary_def;
  if (inputData_->mesh_def.refine_walls) {
    if (bc[BoundaryDef::X_MINUS].ifRefine) {
      for (const auto &p : coords) {
        if (p.x() - inputData_->mesh_def.physDomain.min[0] < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
    if (bc[BoundaryDef::Y_MINUS].ifRefine) {
      for (const auto &p : coords) {
        if (p.y() - inputData_->mesh_def.physDomain.min[1] < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
    if (bc[BoundaryDef::X_PLUS].ifRefine) {
      for (const auto &p : coords) {
        if (inputData_->mesh_def.physDomain.max[0] - p.x() < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
    if (bc[BoundaryDef::Y_PLUS].ifRefine) {
      for (const auto &p : coords) {
        if (inputData_->mesh_def.physDomain.max[1] - p.y() < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
#if (DIM == 3)
    if (bc[BoundaryDef::Z_MINUS].ifRefine) {
      for (const auto &p : coords) {
        if (p.z() - inputData_->mesh_def.physDomain.min[2] < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
    if (bc[BoundaryDef::Z_PLUS].ifRefine) {
      for (const auto &p : coords) {
        if (inputData_->mesh_def.physDomain.max[2] - p.z() < eps and currLevel < inputData_->mesh_def.refine_lvl_channel_wall) {
          levelForWall = inputData_->mesh_def.refine_lvl_channel_wall;
          break;
        }
      }
    }
#endif
  }
  /// extra wall refine
  DENDRITE_UINT levelForWallExtra = inputData_->mesh_def.refine_lvl_base;

  /// region refine
  auto &rf = inputData_->region_refine;
  unsigned int maxlevelForRegion = inputData_->mesh_def.refine_lvl_base;
  std::vector<unsigned int> levelForRegions;
  levelForRegions.resize(rf.size());
  for (int i = 0; i < rf.size(); i++) {
    auto &r = rf[i];
    if (!r.forRetain) {
      for (const auto &p : coords) {
        if (r.in_region(p)) {
          levelForRegions[i] = r.refine_region_lvl;
          break;
        } else {
          levelForRegions[i] = inputData_->mesh_def.refine_lvl_base;
        }
      }
    }
  }
  if (!levelForRegions.empty()) {
    maxlevelForRegion = *max_element(levelForRegions.begin(), levelForRegions.end());
  }
  // for region with geometry, refine the boundries.
  std::vector<unsigned int> levelForRegionsGeomBoundary;
  std::vector<int> countOfInPointsEachRegionGeom;
  levelForRegionsGeomBoundary.resize(rf.size());
  countOfInPointsEachRegionGeom.resize(rf.size());
  for (unsigned int i = 0; i < rf.size(); i++) {
    auto &r = rf[i];
    if (!r.forRetain and r.GetRefineType() == RegionalRefine::MESHOBJECT) {
      for (const auto &p : coords) {
        if (r.in_region(p)) {
          countOfInPointsEachRegionGeom[i]++;
        }
      }
    }
  }

  for (unsigned int i = 0; i < rf.size(); i++) {
    if (not(countOfInPointsEachRegionGeom[i] == nodeNo || countOfInPointsEachRegionGeom[i] == 0)) {
      levelForRegionsGeomBoundary[i] = (rf[i].refine_region_lvl_boundary);
    }
  }
  levelForRegionsGeomBoundary.push_back(maxlevelForRegion);
  if (!levelForRegionsGeomBoundary.empty()) {
    maxlevelForRegion = *max_element(levelForRegionsGeomBoundary.begin(), levelForRegionsGeomBoundary.end());
  }


  /// region retain (for complete outside elements outside retain region, we don't want to refine those)
  bool outsideRetain = false;
  for (auto &r : rf) {
    if (r.forRetain) {
      bool all_out = true;
      for (const auto &p : coords) {
        all_out = all_out and (r.out_retain(p) == RegionalRefine::OUTSIDE);
      }
      outsideRetain = outsideRetain or all_out;
    }
  }

  DENDRITE_UINT id = -1;
  bool isObject = false;
  if (this->m_BoundaryOctant) {
    for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
      subDomainBoundary_->generateBoundaryFlags(coords[i], id);
      if (subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
          subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE)) {
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
  if (currLevel < maxlevelForCarvedOutGeom or
      currLevel < maxlevelForRegion or
          currLevel < levelForWall or
          currLevel < maxlevelForCarvedOutGeom or
          currLevel < levelForWallExtra) {
    return ot::OCT_FLAGS::Refine::OCT_REFINE;
  } else {
    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
  }

}
