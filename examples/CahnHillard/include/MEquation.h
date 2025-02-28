/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //



#pragma once
#include <talyfem/fem/cequation.h>
#include "CHNodeData.h"

class MEquation : public TALYFEMLIB::CEquation<CHNodeData> {
    const std::vector<double> * cellMass_;
    const std::vector<double> * cellMassPerChild_;
    const std::vector<double> * cellCoarsened_;
    const int numItg = 1u << DIM;
 public:

  MEquation()
      : TALYFEMLIB::CEquation<CHNodeData>(false, TALYFEMLIB::kAssembleGaussPoints)
              {}


  // never called, but required for PPEquation to be instantiable (pure virtual in CEquation)
  virtual void Solve(double dt, double t) { assert(false); }
  void assign(const std::vector<double> * cellVector, const std::vector<double> * cellVectorPerChild, const std::vector<double> * cellCoarsen){
      cellMass_ = cellVector;
      cellMassPerChild_ = cellVectorPerChild;
      cellCoarsened_ = cellCoarsen;
  }

  bool isCoarsenedCell(const int localElemID)const {
      if(FEQUALS(cellCoarsened_->at(localElemID),1)){
          return false;
    }
      else{
          return true;
      }
  }
  inline double getCellPhi(const int localElemID, int currentItg) const{
      double val = cellMass_->at(
                  localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::PHI_DOF));
//      if(not isCoarsenedCell(localElemID)) {
//          val = cellMass_->at(
//                  localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::PHI_DOF));
//      }
//      else{
//
//          val = cellMassPerChild_->at(
//                  localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::PHI_DOF));
//          assert(val!=0);
//      }

      return val;
  }

  inline double getCellMu(const int localElemID, int currentItg) const{
      double val = cellMass_->at(
              localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::MU_DOF));
//      if(not isCoarsenedCell(localElemID)) {
//          val = cellMass_->at(
//                  localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::MU_DOF));
//      }
//      else{
//          val = cellMassPerChild_->at(
//                  localElemID * numItg * CHNodeData::CH_DOF + currentItg * CHNodeData::CH_DOF + (CHNodeData::MU_DOF));
//      }

      return val;
  }
  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,double * h) {
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();



    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.
    const int localElemID = fe.elem()->elm_id();
    if(isCoarsenedCell(localElemID)){
        for (int a = 0; a < n_basis_functions; a++) {
            for (int b = 0; b < n_basis_functions; b++) {
                double M = fe.N(a) * fe.N(b) * fe.jacc();
                Ae(a*CHNodeData::CH_DOF + CHNodeData::PHI_DOF,b*CHNodeData::CH_DOF + CHNodeData::PHI_DOF) += M;
                Ae(a*CHNodeData::CH_DOF + CHNodeData::MU_DOF, b*CHNodeData::CH_DOF + CHNodeData::MU_DOF) += M;
            }
        }
    }
    else{
        for (int a = 0; a < n_basis_functions; a++) {
            for (int b = 0; b < n_basis_functions; b++) {
                double M = fe.N(a) * fe.N(b) * detJxW;
                Ae(a*CHNodeData::CH_DOF + CHNodeData::PHI_DOF,b*CHNodeData::CH_DOF + CHNodeData::PHI_DOF) += M;
                Ae(a*CHNodeData::CH_DOF + CHNodeData::MU_DOF, b*CHNodeData::CH_DOF + CHNodeData::MU_DOF) += M;
            }
        }
    }


  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
      const int localElemID = fe.elem()->elm_id();
      const double cellPhi = this->getCellPhi(localElemID,fe.cur_itg_pt_num());
      const double cellMu =  this->getCellMu(localElemID,fe.cur_itg_pt_num());

      const int n_basis_functions = fe.nbf();
      if(isCoarsenedCell(localElemID)) {
          for (int a = 0; a < n_basis_functions; a++) {
              be(a * CHNodeData::CH_DOF + CHNodeData::PHI_DOF) += fe.N(a) * cellPhi * fe.jacc();
              be(a * CHNodeData::CH_DOF + CHNodeData::MU_DOF) += fe.N(a) * cellMu * fe.jacc();
          }
      }
      else{
          for (int a = 0; a < n_basis_functions; a++) {
              be(a * CHNodeData::CH_DOF + CHNodeData::PHI_DOF) += fe.N(a) * cellPhi * fe.detJxW();
              be(a * CHNodeData::CH_DOF + CHNodeData::MU_DOF) += fe.N(a) * cellMu * fe.detJxW();
          }
      }
  }

};
