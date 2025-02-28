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
#include "MNodeData.h"

class MEquation : public TALYFEMLIB::CEquation<MNodeData> {
    const std::vector<double> * cellMass_;
    const std::vector<double> * cellMassPerChild_;
 public:

  MEquation()
      : TALYFEMLIB::CEquation<MNodeData>(false, TALYFEMLIB::kAssembleGaussPoints)
              {}


  // never called, but required for PPEquation to be instantiable (pure virtual in CEquation)
  virtual void Solve(double dt, double t) { assert(false); }
  void assign(const std::vector<double> * cellMassVector,const std::vector<double> * cellMassPerChild ){
      cellMass_ = cellMassVector;
      cellMassPerChild_ = cellMassPerChild;
  }

  inline double getCellMass(const int localElemID, int currItg) const{
      static constexpr int numItg = 1u << DIM;

//      std::cout << localElemID << " " << cellMass_->at(localElemID) << "\n";
//      return cellMass_->at(localElemID*numItg + currItg);
      return cellMassPerChild_->at(localElemID*numItg + currItg);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,double * h) {
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();


    // in order to assemble the gauss point, we loop over each pair of basis
    // functions and calculate the individual contributions to the 'Ae' matrix
    // and 'be' vector.


    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        Ae(a,b) += fe.N(a) * fe.N(b) * fe.detJxW();
      }
    }

  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
      const double cellMass =  this->getCellMass(fe.elem()->elm_id(), fe.cur_itg_pt_num());
//      for(int a = 0; a < be.size(); a++){
//          be(a) = cellMass;
//      }
      const int n_basis_functions = fe.nbf();
////      if(fe.cur_itg_pt_num() == 0) {
          for (int a = 0; a < n_basis_functions; a++) {
              be(a) += fe.N(a) * cellMass * fe.detJxW();
//              be(a) = cellMass;
          }
//      }
//      if(fe.cur_itg_pt_num() == 3){
//          for(int i = 0; i < be.size();i++){
//              std::cout << be(i) << " ";
//          }
//          std::cout << "\n";
//      }
  }

};
