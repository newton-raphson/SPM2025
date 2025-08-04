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
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::FUNCTION)) {
                isObject = true;
                break;
            }
        }
    }
    bool offset_region = false;
    bool inside_radial_circle = false;
    double offset_value = 0.05;
    if (inputData_->SbmGeo == LEInputData::PLANT) {
        //////////////////////////// this is for offset //////////////////////////////////

//    for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++)
//    {
////      double dx = (coords[i].x() - 0.6) / 0.3;
////      double dy = (coords[i].y() - 0.6) / 0.3;
////      double ellipse_phi = (std::sqrt(dx*dx + dy*dy) - 1.0) * 0.3;
////
////      double circle_phi = std::sqrt(coords[i].x()*coords[i].x()+ coords[i].y()*coords[i].y()) - 0.9;
//
//      double signed_distance = inputData_->plantGeometry.signedDistance(coords[i].x(), coords[i].y());
//
//      if (signed_distance > -inputData_->plantGeometry.offset_threshold)
//      {
//        offset_region = true;
//        break;
//      }
//    }
        offset_region = inputData_->plantGeometry.isInOffsetRegion(coords);
        inside_radial_circle = inputData_->plantGeometry.isInsideRadialCircle(coords);

//    /////////////////////////////////////////////////////
//    using Circle = std::pair<double, double>;
//
//    std::default_random_engine rng(200);  // Fixed seed for reproducibility
//    std::uniform_real_distribution<double> angle_jitter(-0.15, 0.15);     // jitter in radians
//    std::uniform_real_distribution<double> radius_jitter(-0.04, 0.04);    // jitter in radius
//
//    double x_origin = 0.0;
//    double y_origin = 0.0;
//    double radius_max = 0.8;
//    int num_rings = 8;
//    int num_sectors = 8;
//    double radius_tolerance = 0.02;  // actual radius of each small circle
//    double angular_step = 2.0 * M_PI / num_sectors;
//    double base_radius = radius_max / num_rings;
//
//    std::vector<Circle> circle_centers;
//
//    for (int r = 1; r <= num_rings; ++r) {
//      double radial_distance = r * base_radius;
//      for (int a = 0; a < num_sectors; ++a) {
//        double jittered_angle = a * angular_step + angle_jitter(rng);
//        double jittered_radius = radial_distance + radius_jitter(rng);
//
//        double x = x_origin + jittered_radius * std::cos(jittered_angle);
//        double y = y_origin + jittered_radius * std::sin(jittered_angle);
//
//        // Ensure no overlap with existing circles
//        bool overlaps = false;
//        for (const auto& c : circle_centers) {
//          double dx = x - c.first;
//          double dy = y - c.second;
//          if (std::sqrt(dx * dx + dy * dy) < 2.1 * radius_tolerance) {
//            overlaps = true;
//            break;
//          }
//        }
//
//        if (!overlaps) {
//          circle_centers.emplace_back(x, y);
//        }
//      }
//    }
//
//    // -----------------------
//    // Check if element nodes fall in any circle
//    // -----------------------
//
//
//    for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
//      for (const auto& center : circle_centers) {
//        double dx = coords[i].x() - center.first;
//        double dy = coords[i].y() - center.second;
//        double dist = std::sqrt(dx * dx + dy * dy);
//
//        if (dist <= radius_tolerance) {
//          inside_radial_circle = true;
//          break;
//        }
//      }
//
//      if (inside_radial_circle)
//        break;
//    }
//  }
        // double dx = (x - 0.6) / 0.3;
        // double dy = (y - 0.6) / 0.3;
        // double ellipse_phi = (std::sqrt(dx*dx + dy*dy) - 1.0) * 0.3;
        //
        // double circle_phi = std::sqrt(x*x + y*y) - 0.9;
        //
        // double signed_distance = std::max(circle_phi, -ellipse_phi);


        unsigned int maxlevelForCarvedOutGeom = inputData_->mesh_def.refine_lvl_base;
        if (isObject) {
            maxlevelForCarvedOutGeom = std::max(inputData_->ibm_geom_def.at(id).refine_lvl, maxlevelForCarvedOutGeom);
        }
        if (inside_radial_circle or offset_region) {
            maxlevelForCarvedOutGeom = inputData_->mesh_def.refine_lvl_channel_wall;
        }

//  std::cout << "maxlevelForCarvedOutGeom = " << maxlevelForCarvedOutGeom << "\n";

        // the regional refine should not refine inside the geometry
//  if (insideGeo) {
//    maxlevelForRegion = inputData_->mesh_def.refine_l;
//  }
//  if (outsideRetain) {
//    levelForWall = inputData_->mesh_def.refine_l;
//  }

        if (this->m_level < baselevel or this->m_level < maxlevelForCarvedOutGeom) {
            return ot::OCT_FLAGS::Refine::OCT_REFINE;
        } else {
            return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
        }

    }
}

