#pragma once

#include "LEInputData.h"
#include "LENodeData.h"
#include <cmath>
/**
* Class calculates Appropriate boundary conditions for LE objects
* Usage:
* BoundaryConditions bc(b);
* bc.setBC;
*/

class LEBCSetup
{
private:
  LEInputData *input_data_;
  SubDomainBoundary *boundaries_;
  void returnNormalTractionBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnDisplacementBothSideBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnFixedAtWallBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnFixedWallBottomForceBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnHalfBeamBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnCsvTractionBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
  void returnCarvedOutBoundary(PETSc::Boundary &b, const ZEROPTV &pos);

public:
  LEBCSetup(SubDomainBoundary *boundary, LEInputData *inputData)
      : input_data_(inputData), boundaries_(boundary)
  {
    PrintInfo("Setting up boundary conditions");
  }

  /**
   * Method to setup Boundary conditions for different cases of HT equations.
   * @param b Boundary object to setup dirichlet or other conditions
   * @param position position
   * @param boundary_def input_data_ boundary definition vector
   */
  void setBoundaryConditions(PETSc::Boundary &b, const ZEROPTV &position, BCTYPE type)
  {
    if (type == e100)
    {
      b.addDirichlet(0,position.x());
      b.addDirichlet(1,0);
      return;
    }
    else if (type == e010)
    {
      b.addDirichlet(0,0);
      b.addDirichlet(1,position.y());
      return;
    }
    else if (type == e001)
    {
      b.addDirichlet(0,0.5*position.y());
      b.addDirichlet(1,0.5*position.x());
      return;

    }
    throw std::runtime_error("Unsupported boundary condition type");

  }
};



