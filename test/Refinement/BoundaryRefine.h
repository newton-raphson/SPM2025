//
// Created by maksbh on 3/26/23.
//

#ifndef DENDRITEKT_BOUNDARYREFINE_H
#define DENDRITEKT_BOUNDARYREFINE_H
#include <Traversal/Refinement.h>

class BoundaryRefine: public Refinement{


public:


    BoundaryRefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo);

    ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords) override;

};

BoundaryRefine::BoundaryRefine(DA *da,  const std::vector<TREENODE> & treePart, const DomainExtents &domainInfo)
        :Refinement(da,treePart,domainInfo){
    this->initRefinement();
}

ot::OCT_FLAGS::Refine BoundaryRefine::getRefineFlags(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &coords){

    DENDRITE_REAL Left[DIM]{0.4,0.4};
    DENDRITE_REAL Right[DIM]{0.6,0.6};
    DENDRITE_UINT rect_lvl = 6;


    const double eps = 1e-13;
    unsigned int levelForRect = this->m_level;
    unsigned int levelForBd = this->m_level;

    for (const auto &p : coords)
    {
        if (p.y() - Left[1] > eps and p.y() - Right[1] < eps and
            p.x() - Left[0] > eps and p.x() - Right[0] < eps and
            this->m_level < rect_lvl)
        {
            levelForRect = rect_lvl;
            break;
        }

//        if (p.y() - Left[1] > eps and p.y() - Right[1] < eps and
//            p.x() - Left[0] > eps and p.x() - Right[0] < eps and
//            this->m_level < inputData_->refine_lvl_bd) {
//
//            levelForBd = refine_lvl_bd;
//            break;
//        }
    }




   if(this->m_BoundaryOctant or  (this->m_level < levelForRect)){
       return ot::OCT_FLAGS::OCT_REFINE;
   }

   if(
     (FEQUALS(coords[0].x(),0.0)) and (FEQUALS(coords[0].y(),0)) or
     (FEQUALS(coords[1].x(),1.0)) and (FEQUALS(coords[1].y(),0)) or
     (FEQUALS(coords[2].x(),0.0)) and (FEQUALS(coords[2].y(),1)) or
     (FEQUALS(coords[3].x(),1.0)) and (FEQUALS(coords[3].y(),1))
     )
   {
     //return ot::OCT_FLAGS::OCT_REFINE;
   }
   return ot::OCT_FLAGS::OCT_NO_CHANGE;
}
#endif //DENDRITEKT_BOUNDARYREFINE_H
