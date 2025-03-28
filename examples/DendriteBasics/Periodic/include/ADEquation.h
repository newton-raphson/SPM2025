#pragma once

#include <talyfem/fem/cequation.h>
#include "ADNodeData.h"
#include <DataTypes.h>
#include <talyfem/talyfem.h>

class ADEquation : public TALYFEMLIB::CEquation<ADNodeData> {
  double diffusion = 1E-2;


 public:
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae) {

    using namespace TALYFEMLIB;
    const int n_dimensions = DIM;

    // # of basis functions
    const int n_basis_functions = fe.nbf();
    TezduyarUpwindFE SUPG;
    TALYFEMLIB::ZEROPTV netVelocity(1.0,0.0,0.0);
    SUPG.calcSUPGWeightingFunction(fe, netVelocity, diffusion);
     // (determinant of J) cross W
    const double detJxW = fe.detJxW();
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double M = (fe.N(a) + SUPG.SUPG(a))*fe.N(b)*detJxW/dt_;
        double N = 0;
        double A = 0;
        for (int k = 0; k < n_dimensions; k++) {
          A += (fe.N(a) + SUPG.SUPG(a))*fe.dN(b,k)*netVelocity(k)*detJxW;
          N += diffusion*fe.dN(a, k) * fe.dN(b, k)*detJxW;
        }
        Ae(a, b) += A + N + M;
      }
    }



  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be) {
    using namespace TALYFEMLIB;
    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W
    const double detJxW = fe.detJxW();

    const ZEROPTV p = fe.position();
    double force = calc_d2u_at(p);
    TezduyarUpwindFE SUPG;
    TALYFEMLIB::ZEROPTV netVelocity(1.0,0.0,0.0);
    SUPG.calcSUPGWeightingFunction(fe, netVelocity, diffusion);

    double C = this->p_data_->valueFEM(fe,0);
    for (int a = 0; a < n_basis_functions; a++) {
      be(a) += fe.N(a)*C*detJxW/dt_;
    }

  }



 protected:
  double calc_d2u_at(const TALYFEMLIB::ZEROPTV &pt) const {
    return 0;
  }
};
