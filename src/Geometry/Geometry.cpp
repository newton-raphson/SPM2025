//
// Created by maksbh on 8/19/20.
//

#include <Geometry/Geometry.h>

namespace GEOMETRY{
#if(DIM == 3)
Geometry::Geometry(const STL *_stl, const Point<DIM> & point,const RetainSide _retainSide):
m_stl(_stl),geomType_(GeomType::STL_3D),retainSide(_retainSide){
  m_translation.clear();
  m_translation.push_back(point);
}
#endif
#if(DIM == 2)
Geometry::Geometry(const MSH *_msh, const Point<DIM> & point,const RetainSide _retainSide)
:m_msh(_msh),geomType_(GeomType::MSH_2D),retainSide(_retainSide){
  m_translation.clear();
  m_translation.push_back(point);
}
#endif

void Geometry::addTranslations(const Point<DIM> & point) {
  m_translation.push_back(point);
}


bool Geometry::ifInside(const DENDRITE_REAL * position) const{
#ifndef DNDEBUG
  if(m_translation.empty()){
    throw std::runtime_error("Need to translate the Geometry. If not translating, still need to pass the point with 0's");
  }
  assert(m_translation.size() != 0);
#endif
  DENDRITE_REAL _position[3];
  bool isInside = false;
  for(const auto & pt : m_translation) {
    // 1. Translate the points
    for(int dim = 0; dim < DIM; dim++){
      _position[dim] = position[dim] - pt.x(dim);
    }

    // Mark In Out: If any of the translated point marks it inside. Return inside.
#if(DIM == 3)
    if (geomType_ == GeomType::STL_3D) {
      if(m_stl->ifInside(_position)){
        isInside = true;
      }
    }
#endif
#if(DIM == 2)
    if (geomType_ == GeomType::MSH_2D) {
      if(m_msh->ifInside(_position)){
        isInside = true;
      }
    }
#endif

  }

    isInside = (retainSide==RetainSide::OUT) ? isInside:not(isInside);
    return isInside;



}

}