//
// Created by maksbh on 8/19/20.
//

#ifndef DENDRITEKT_GEOMETRY_H
#define DENDRITEKT_GEOMETRY_H

#include <Geometry/STL.h>
#include <Geometry/MSH.h>

namespace GEOMETRY {

enum GeomType:unsigned short {
  STL_3D = 0,
  MSH_2D = 1
};

/**
 * @brief class to contain the GEOMETRY information.
 * This class is designed to accumulate different classes of same objects which will have
 * same boundary conditions.
 */


class Geometry {
  GeomType geomType_;
  const RetainSide retainSide;
 protected:
#if (DIM == 3)
  const STL * m_stl = nullptr;
#endif
#if (DIM == 2)
  const MSH * m_msh = nullptr;
#endif
  std::vector<Point<DIM>> m_translation;
 public:
#if (DIM == 3)
 Geometry(const STL * _stl,const Point<DIM> & point, const RetainSide _retainSide = RetainSide::OUT);
#endif
#if (DIM == 2)
  Geometry(const MSH * _msh,const Point<DIM> & point,const RetainSide _retainSide = RetainSide::OUT);
#endif
  void addTranslations(const Point<DIM> & point);

  bool ifInside(const DENDRITE_REAL * point) const;

#if (DIM == 3)
  inline const InOutTest getInOutTestType() const {
    return m_stl->getInOutTest();
  }
#endif
#if (DIM == 2)
  inline const InOutTest2D getInOutTestType() const {
    return m_msh->getInOutTest();
  }
#endif



#if (DIM == 3)
  inline const STL * getSTL() const{
    return m_stl;
  }
#endif
#if (DIM == 2)
  inline const MSH * getMSH() const{
    return m_msh;
  }
#endif

  inline const std::vector<Point<DIM>> & getTranslations() const{
    return m_translation;
  }
  inline const RetainSide & getRetainSide() const{
    return retainSide;
  }
};
}
#endif //DENDRITEKT_GEOMETRY_H
