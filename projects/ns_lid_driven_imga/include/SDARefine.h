//
// Created by maksbh on 11/8/20.
//

#ifndef DENDRITEKT_SDAREFINE_H
#define DENDRITEKT_SDAREFINE_H
#include <Traversal/Refinement.h>

class SDARefine: public Refinement{
    const DENDRITE_UINT maxLevel_;
    const TimeInfo ti_;
public:
    NSInputData *inputData_;
    SubDomainBoundary *subDomainBoundary_;

    SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, DENDRITE_UINT maxLevel,NSInputData *inputData,
              SubDomainBoundary *subDomainBoundary,
              const TimeInfo& ti);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

    ~SDARefine() = default;

};

SDARefine::SDARefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo, const DENDRITE_UINT maxLevel,NSInputData *inputData,
                     SubDomainBoundary *subDomainBoundary,
                     const TimeInfo& ti)
    :Refinement(da, treePart, domainInfo), maxLevel_(maxLevel), inputData_(inputData), ti_(ti),subDomainBoundary_(subDomainBoundary)
{
    this->traverse();
}

ot::OCT_FLAGS::Refine SDARefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){
    const double eps = 1e-13;
    unsigned int BaseLevel;
    unsigned int levelForBd = this->m_level;

    for (const auto &p : coords)
    {
        if (inputData_->lvl.empty()) {
//            std::cout << "empty\n";
//            if ((not(p.y() - inputData_->rectLowerLeft[1] > eps and p.y() - inputData_->rectUpperRight[1] < eps and
//                     p.x() - inputData_->rectLowerLeft[0] > eps and p.x() - inputData_->rectUpperRight[0] < eps)) and
//                this->m_level < inputData_->rect_lvl) {
//                levelForRect = inputData_->rect_lvl;
//                break;
//            }
            BaseLevel = maxLevel_;
        } else {

//            std::cout << "not empty\n";

            int refine_lvl = inputData_->lvl[0];

            for (int idx = inputData_->totalT.size() - 1; idx > 0; idx--) {
                if (ti_.getCurrentTime() >= inputData_->totalT.at(idx - 1)) {
//                    std::cout << "inputData_->totalT.at(idx - 1) = " << inputData_->totalT.at(idx - 1) << "\n";
                    refine_lvl = inputData_->lvl[idx];
                    break;
                }
            }

            BaseLevel = refine_lvl;

//            std::cout << "ti_.getCurrentTime() = " << ti_.getCurrentTime() << "\n";
//            std::cout << "refine_lvl = " << refine_lvl << "\n";
//            if ((not(p.y() - inputData_->rectLowerLeft[1] > eps and p.y() - inputData_->rectUpperRight[1] < eps and
//                     p.x() - inputData_->rectLowerLeft[0] > eps and p.x() - inputData_->rectUpperRight[0] < eps)) and
//                this->m_level < refine_lvl) {
//                levelForRect = refine_lvl;
//                break;
//            }
        }


//        std::cout << "this->m_level = " << this->m_level << "\n";
//        std::cout << "levelForRect = " << levelForRect << "\n";

//        if (p.y() - inputData_->rectLowerLeft[1] > eps and p.y() - inputData_->rectUpperRight[1] < eps and
//            p.x() - inputData_->rectLowerLeft[0] > eps and p.x() - inputData_->rectUpperRight[0] < eps and
//            this->m_level < inputData_->meshDef.refineLevelBoundary) {
//
//            levelForBd = inputData_->meshDef.refineLevelBoundary;
//            break;
//        }
    }

    DENDRITE_UINT id = -1;
    bool isObject = false;
    if (this->m_BoundaryOctant)
    {
        for (DENDRITE_UINT i = 0; i < m_octDA->getNumNodesPerElement(); i++)
        {
            subDomainBoundary_->generateBoundaryFlags(coords[i], id);
            if (subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or
                subDomainBoundary_->checkBoundaryType(BoundaryTypes::FUNCTION))
            {
                isObject = true;
                break;
            }
        }
    }
    unsigned int maxlevelForCurvedOutGeom = BaseLevel;
    if (isObject)
    {
//        std::cout << "here\n";
//        std::cout << "id = " << id << "\n";
//        std::cout << "inputData_->carved_out_geoms_def.at(id).refine_lvl = " << inputData_->carved_out_geoms_def.at(id).refine_lvl << "\n";
        maxlevelForCurvedOutGeom = std::max(inputData_->carved_out_geoms_def.at(id).refine_lvl, maxlevelForCurvedOutGeom);
    }

//    std::cout << "maxlevelForCurvedOutGeom = " << maxlevelForCurvedOutGeom << "\n";

    if(this->m_level < BaseLevel or this->m_level < maxlevelForCurvedOutGeom){
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
