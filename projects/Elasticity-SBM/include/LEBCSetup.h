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
  void setBoundaryConditions(PETSc::Boundary &b, const ZEROPTV &position)
  {

    if (input_data_->bccaseType == BCCaseType::NORMAL_TRACTION)
    {
      returnNormalTractionBoundary(b, position);
    }

    else if (input_data_->bccaseType == BCCaseType::DISPLACEMENT_BOTH_SIDE)
    {
      returnDisplacementBothSideBoundary(b, position);
    }

    else if (input_data_->bccaseType == BCCaseType::FIXED_AT_WALL)
    {
      returnFixedAtWallBoundary(b, position);
    }

    else if (input_data_->bccaseType == BCCaseType::HALF_BEAM)
    {
      returnHalfBeamBoundary(b, position);
    }

    else if (input_data_->bccaseType == BCCaseType::BOTTOM_FORCE and input_data_->SbmGeo == LEInputData::SBMGeo::NONE)
    {
        returnFixedWallBottomForceBoundary(b, position);
    }

    else if (input_data_->bccaseType == BCCaseType::CSV_FORCE)
    {
        returnCsvTractionBoundary(b, position);
    }

      returnCarvedOutBoundary(b, position);
  };

};

void LEBCSetup::returnDisplacementBothSideBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos)
{
  static const double eps = 1e-15;

  bool x_minus_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.min[0]) < eps;
  bool x_max_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.max[0]) < eps;
  if (x_max_wall)
  {
    b.addDirichlet(0, input_data_->DisplacementBothSide.displacement); // x-dir displacement
    b.addDirichlet(1, 0.0);
  }
  if (x_minus_wall)
  {
    b.addDirichlet(0, -input_data_->DisplacementBothSide.displacement); // x-dir displacement
    b.addDirichlet(1, 0.0);
  }
}

void LEBCSetup::returnFixedAtWallBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos)
{
  static const double eps = 1e-15;

  bool x_minus_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.min[0]) < eps;
  if (x_minus_wall)
  {
    b.addDirichlet(0, 0.0); // x-dir displacement
    b.addDirichlet(1, 0.0); // y-dir displacement
#if (DIM == 3)
    b.addDirichlet(2, 0.0); // z-dir displacement
#endif
  }
}

void LEBCSetup::returnNormalTractionBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos)
{
  static const double eps = 1e-15;

  bool x_minus_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.min[0]) < eps;
  bool y_minus_wall = fabs(pos.y() - input_data_->mesh_def.physDomain.min[1]) < eps;
  if (x_minus_wall)
  {
    b.addDirichlet(0, 0.0); // x-dir displacement
  }
  if (y_minus_wall)
  {
    b.addDirichlet(1, 0.0); // y-dir displacement
  }

#if (DIM == 3)
  bool z_minus_wall = fabs(pos.z() - input_data_->mesh_def.physDomain.min[2]) < eps;
  if (z_minus_wall)
  {
    b.addDirichlet(2, 0.0); // z-dir displacement
  }
#endif
}

void LEBCSetup::returnHalfBeamBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos)
{
  static const double eps = 1e-15;

  bool x_minus_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.min[0]) < eps;
  bool x_max_wall = fabs(pos.x() - input_data_->mesh_def.physDomain.max[0]) < eps;
  if (x_max_wall)
  {
    b.addDirichlet(0, 0.0); 
  }
  if (x_minus_wall)
  {
    b.addDirichlet(1, 0.0);
#if (DIM == 3)
    b.addDirichlet(2, 0.0); 
#endif    
  }
}

void LEBCSetup::returnFixedWallBottomForceBoundary(PETSc::Boundary &b, const ZEROPTV &pos) {

    static const double eps = 1e-5;

//    std::cout <<"---------------------------------\n";
//    std::cout <<" input_data_->mesh_def.fullDADomain.max[0] = " << input_data_->mesh_def.fullDADomain.max[0] << "\n";
//    std::cout <<" input_data_->mesh_def.physDomain.max[0] = " << input_data_->mesh_def.physDomain.max[0] << "\n";
//    std::cout << "input_data_->mesh_def.channel_max[0] = " << input_data_->mesh_def.channel_max[0] << "\n";
//    std::cout <<"---------------------------------\n";

    bool x_max_wall = fabs(pos.x() - input_data_->DomainMax(0)) < eps;
    if (x_max_wall) {
        b.addDirichlet(0, 0.0); // x-dir displacement
        b.addDirichlet(1, 0.0); // y-dir displacement
#if (DIM == 3)
        b.addDirichlet(2, 0.0); // z-dir displacement
#endif
    }
}


void LEBCSetup::returnCsvTractionBoundary(PETSc::Boundary &b, const ZEROPTV &pos) {

    static const double eps = 1e-15;

    for (int dim = 0; dim<DIM;dim++) {
        if (input_data_->CsvForce.fixside/2 == dim) {
            bool wall = (input_data_->CsvForce.fixside%2)?fabs(pos(dim) - input_data_->DomainMax(dim)) < eps:
                    fabs(pos(dim) - input_data_->DomainMin(dim)) < eps;
            if (wall) {
#ifndef NDEBUG
#endif
                std::cout << "setting Dirichlet BC\n";

                b.addDirichlet(0, 0.0); // x-dir displacement
                b.addDirichlet(1, 0.0); // y-dir displacement
#if (DIM == 3)
                b.addDirichlet(2, 0.0); // z-dir displacement
#endif
            }
        }
    }
}

void LEBCSetup::returnCarvedOutBoundary(PETSc::Boundary &b, const ZEROPTV &pos) {

    bool DirichletHaveSet = false;
    double DirichletBCValue[DIM];

    double x = pos.x();
    double y= pos.y();

    DENDRITE_UINT objectID = -1;
    boundaries_->generateBoundaryFlags(pos, objectID);
    //std::cout << "NSHTNodeData::TEMPERATURE - NSHTNodeData::NS_DOF = " << NSHTNodeData::TEMPERATURE - NSHTNodeData::NS_DOF <<"\n";
    if (boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::SPHERE) or
        boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::CIRCLE) or
        boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::BOX) or
        boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY)) {
        const auto &carved_geo = input_data_->ibm_geom_def.at(objectID);

//        if (carved_geo.bc_type_D[0] == IBMGeomDef::STRONG_Dirichlet) {
//
//            switch (input_data_->SbmGeo) {
//                case LEInputData::SBMGeo::CIRCLE: {
//
//
//                    double rad = 0.5;
//                    double radius_gp = sqrt((x - carved_geo.InitialDisplacement.x()) * (x - carved_geo.InitialDisplacement.x())
//                            + (y - carved_geo.InitialDisplacement.y()) * (y - carved_geo.InitialDisplacement.y()));
//                    double x_proj = (rad * (x - carved_geo.InitialDisplacement.x()) / radius_gp + carved_geo.InitialDisplacement.x());
//                    double y_proj = (rad * (y - carved_geo.InitialDisplacement.y()) / radius_gp + carved_geo.InitialDisplacement.y());
//                    DirichletBCValue[0] = sin(M_PI*x_proj)*cos(M_PI*y_proj)/10.0;
//                    DirichletBCValue[1] = cos(M_PI*x_proj)*sin(M_PI*y_proj)/10.0;
//
//                    DirichletHaveSet = true;
//                    break;
//                }
//            }
//
//            if (DirichletHaveSet) {
//
//                for (int dim = 0; dim < DIM; dim++) {
//                    b.addDirichlet(dim, DirichletBCValue[dim]);
//                }
//            } else {
//                b.addDirichlet(0, carved_geo.getBC_D(0)[0]); // x-dir displacement
//                b.addDirichlet(1, carved_geo.getBC_D(1)[0]); // y-dir displacement
//#if (DIM == 3)
//                b.addDirichlet(2, carved_geo.getBC_D(2)[0]); // z-dir displacement
//#endif
//            }
//        }

    }
}


