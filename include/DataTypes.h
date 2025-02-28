//
// Created by maksbh on 9/16/19.
//

#ifndef DENDRITEKT_DATATYPES_H
#define DENDRITEKT_DATATYPES_H

#include <oda.h>
#include <talyfem/grid/zeroptv.h>
#include <array>
#include<string.h>
#include "talyfem/grid/elem_common.h"

#ifdef ENABLE_4D
#define DIM 4
#endif

#ifdef ENABLE_2D
#define DIM 2
static constexpr bool SCALINGX[4]{false, true, false, true};
static constexpr bool SCALINGY[4]{false, false, true, true};

enum ElementType_2D:int{
    TRIANGLE = TALYFEMLIB::kElem2dTriangle,
    QUAD = TALYFEMLIB::kElem2dBox
};

#endif

#ifdef ENABLE_3D
#define DIM 3
  static constexpr bool SCALINGX[8]{false,true,false,true,false,true,false,true};
  static constexpr bool SCALINGY[8]{false,false,true,true,false,false,true,true};
  static constexpr bool SCALINGZ[8]{false,false,false,false,true,true,true,true};

  enum ElementType_3D:int{
    TET = TALYFEMLIB::kElem3dTetrahedral,
    HEX = TALYFEMLIB::kElem3dHexahedral,
    PYRAMID = TALYFEMLIB::kElemPyramid
};

#endif
#include "point.h"

typedef unsigned int DENDRITE_UINT;
typedef double DENDRITE_REAL;
typedef ot::DA<DIM> DA;

typedef ot::TreeNode<DENDRITE_UINT, DIM> TREENODE;
typedef ot::DistTree<DENDRITE_UINT, DIM> DistTREE;

#define FEQUALS(x, y) fabs((x) - (y)) < 1E-10 ? true : false

typedef std::function<double(const TALYFEMLIB::ZEROPTV & pos, const DENDRITE_UINT dof, const DENDRITE_REAL time)> AnalyticFunction;

typedef std::function<ibm::Partition(const double *elemPhysCoords, double elemPhysSize)> PhysicalDomainDecider;

typedef std::function<ibm::Partition(const std::vector<Point<DIM>> & coords)> FunctionToRetain;

typedef std::function<ibm::Partition (const double * position)> InOutClassifier;


struct DomainInfo{
  std::array<DENDRITE_REAL,DIM> min; /// Minimum of the domain
  std::array<DENDRITE_REAL,DIM> max; /// maximum of the domain
};

struct DomainExtents{
  DomainInfo  & fullDADomain;
  DomainInfo  & physicalDADomain;
  DomainExtents( DomainInfo & fullDA, DomainInfo & physDomain)
  :fullDADomain(fullDA),physicalDADomain(physDomain){
  }
  DomainExtents(DomainInfo & fullDA)
      :fullDADomain(fullDA),physicalDADomain(fullDA){
  }

};

enum RetainSide:bool{
    OUT = false,
    IN = true
};

template <const char ** enumVarname,int maxEnums>
static int convertEnumToStrings(const char * stringName){
  for(int i = 0; i < maxEnums; i++){
    if(strcmp(stringName,enumVarname[i]) == 0){
      return i;
    }
  }
  throw std::logic_error("String did not match not found");
}

extern std::array<DENDRITE_UINT,DIM> PERIODS;
#endif //DENDRITEKT_DATATYPES_H
