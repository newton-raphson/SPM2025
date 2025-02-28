#pragma once

#include <Traversal/Refinement.h>
#include <Boundary/SubDomainBoundary.h>

#include <utility>
#include "NSInputData.h"
#include "SBMCalc.h"

class NSRefine : public Refinement
{
    NSInputData *inputData_;
    const DomainExtents &domainExtents_;
    SubDomainBoundary *subDomainBoundary_;
    bool doCoarsen_ = false;
    const IMGA *imga_;
    const my_kd_tree_t *kd_tree_;
    const TimeInfo ti_;

public:
    NSRefine(DA *octDA,
               const std::vector<TREENODE> &treePart,
               const DomainExtents &domainExtents,
               NSInputData *inputData,
               SubDomainBoundary *subDomainBoundary,
               const IMGA *imga,
               const my_kd_tree_t *kd_tree,
             const TimeInfo& ti,
             bool doCoarsen = false);

    NSRefine(DA *octDA,
               const std::vector<TREENODE> &treePart,
               const DomainExtents &domainExtents,
               NSInputData *inputData,
               SubDomainBoundary *subDomainBoundary,
               const TimeInfo& ti,
               bool doCoarsen = false);
    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

    ~NSRefine() = default;
};
/**
 * @brief constructor
 * @param octDA
 * @param treePart tree partition
 * @param domainExtents domain information
 * @param inputData input data
 * @param subDomainBoundary boundary information (generateBoundaryFlags)
 */
NSRefine::NSRefine(DA *octDA,
                       const std::vector<TREENODE> &treePart,
                       const DomainExtents &domainExtents,
                       NSInputData *inputData,
                       SubDomainBoundary *subDomainBoundary,
                       const IMGA *imga,
                       const my_kd_tree_t *kd_tree,
                   const TimeInfo& ti,
                   bool doCoarsen)
        : Refinement(octDA, treePart, domainExtents), inputData_(inputData), domainExtents_(domainExtents),
          subDomainBoundary_(subDomainBoundary), doCoarsen_(doCoarsen),imga_(imga),kd_tree_(kd_tree), ti_(ti)
{
    this->traverse();
}
NSRefine::NSRefine(DA *octDA,
                       const std::vector<TREENODE> &treePart,
                       const DomainExtents &domainExtents,
                       NSInputData *inputData,
                       SubDomainBoundary *subDomainBoundary,
                       const TimeInfo& ti,
                       bool doCoarsen)
        : Refinement(octDA, treePart, domainExtents), inputData_(inputData), domainExtents_(domainExtents),
          subDomainBoundary_(subDomainBoundary), doCoarsen_(doCoarsen) , ti_(ti)
{
    this->traverse();
}



/**
 * @brief The refinement flags per element by element. If nothing is returned, its set to NO_CHANGE.
 * @param fe the element
 * @param coords vector of coords
 * @return the flags for each elements
 */
ot::OCT_FLAGS::Refine NSRefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords)
{

    const DomainInfo &physDomain = domainExtents_.physicalDADomain;
    const unsigned int nodeNo = m_octDA->getNumNodesPerElement();
    const DENDRITE_UINT currLevel = this->m_level;

    /// refine walls (maximum to the refine_h level)
    const double eps = 1e-13;
    unsigned int levelForWall = inputData_->meshDef.baseLevel;
    const auto &bc = inputData_->boundary_def;
//    std::cout<<"inputData_->meshDef.refineLevelBoundary = "<<inputData_->meshDef.refineLevelBoundary<<"\n";
    if (inputData_->meshDef.refineLevelBoundary)
    {
            for (const auto &p : coords)
            {
                if (p.x() - inputData_->physDomain.min[0] < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
                {
                    levelForWall = inputData_->meshDef.refineLevelBoundary;
                    break;
                }
                if (p.y() - inputData_->physDomain.min[1] < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
                {
                    levelForWall = inputData_->meshDef.refineLevelBoundary;
                    break;
                }
                if (inputData_->physDomain.max[0] - p.x() < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
                {
                    levelForWall = inputData_->meshDef.refineLevelBoundary;
                    break;
                }
                if (inputData_->physDomain.max[1] - p.y() < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
                {
                    levelForWall = inputData_->meshDef.refineLevelBoundary;
                    break;
                }
#if (DIM == 3)
        if (p.z() - inputData_->physDomain.min[2] < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
        {
          levelForWall = inputData_->meshDef.refineLevelBoundary;
          break;
        }

        if (inputData_->physDomain.max[2] - p.z() < eps and currLevel < inputData_->meshDef.refineLevelBoundary)
        {
          levelForWall = inputData_->meshDef.refineLevelBoundary;
          break;
        }
#endif

            }

    }


    /// region refine
    auto &rf = inputData_->region_refine;
    unsigned int maxlevelForRegion = inputData_->meshDef.baseLevel;
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
                    levelForRegions[i] = inputData_->meshDef.baseLevel;
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
    for (unsigned int i = 0; i < rf.size(); i++)
    {
        auto &r = rf[i];
        if (!r.forRetain and (r.GetRefineType() == RegionalRefine::MESHOBJECT))
        {
            for (const auto &p : coords)
            {
                if (r.in_region(p))
                {
                    countOfInPointsEachRegionGeom[i]++;
                }
            }
        }
    }

    for (unsigned int i = 0; i < rf.size(); i++)
    {
        if (not(countOfInPointsEachRegionGeom[i] == nodeNo || countOfInPointsEachRegionGeom[i] == 0))
        {
            levelForRegionsGeomBoundary[i] = (rf[i].refine_region_lvl_boundary);
        }
    }
    levelForRegionsGeomBoundary.push_back(maxlevelForRegion);
    if (!levelForRegionsGeomBoundary.empty())
    {
        maxlevelForRegion = *max_element(levelForRegionsGeomBoundary.begin(), levelForRegionsGeomBoundary.end());
    }

    /// region retain (for complete outside elements outside retain region, we don't want to refine those)
    bool outsideRetain = false;
    for (auto &r : rf)
    {
        if (r.forRetain)
        {
            bool all_out = true;
            for (const auto &p : coords)
            {
                all_out = all_out and (r.out_retain(p) == RegionalRefine::OUTSIDE);
            }
            outsideRetain = outsideRetain or all_out;
        }
    }

    DENDRITE_UINT id = -1;
    bool isObject = false;
    if (this->m_BoundaryOctant)
    {
        for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++)
        {
            subDomainBoundary_->generateBoundaryFlags(coords[i], id);
            if (subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::FUNCTION) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE))
            {
//                std::cout<<"isObject = true\n";
                isObject = true;
                break;
            }
        }
    }
    unsigned int maxlevelForCurvedOutGeom = inputData_->meshDef.baseLevel;
    if (isObject)
    {
        maxlevelForCurvedOutGeom = std::max(inputData_->carved_out_geoms_def.at(id).refine_lvl, maxlevelForCurvedOutGeom);
//        std::cout<<"maxlevelForCurvedOutGeom = "<<maxlevelForCurvedOutGeom<<"\n";
    }

    int refine_lvl_shift = 0;

//    refine_lvl_shift = inputData_->lvl[0];


//    std::cout << "refine_lvl_shift = " << refine_lvl_shift << "\n";


    //  if (outsideRetain) {
    //    levelForWall = inputData_->meshDef.refine_l;
    //  }
    if (!doCoarsen_)
    {
        if (currLevel < maxlevelForRegion + refine_lvl_shift or
            currLevel < levelForWall or
            currLevel < maxlevelForCurvedOutGeom + refine_lvl_shift)
        {
            return ot::OCT_FLAGS::Refine::OCT_REFINE;
        }
        else
        {
            return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
        }
    }
    else
    {

        if (currLevel > maxlevelForRegion + refine_lvl_shift and
            currLevel > levelForWall and
            currLevel > maxlevelForCurvedOutGeom + refine_lvl_shift)
        {
            return ot::OCT_FLAGS::Refine::OCT_COARSEN;
        }
    }

    return ot::OCT_FLAGS::Refine::OCT_NO_CHANGE;
}
