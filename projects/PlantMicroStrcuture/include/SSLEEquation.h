#pragma once

#include <cmath>
#include <talyfem/fem/cequation.h>
#include <talyfem/talyfem.h>
#include "LEInputData.h"
#include "LENodeData.h"

class SSLEEquation : public TALYFEMLIB::CEquation<LENodeData> {
public:
  explicit SSLEEquation(LEInputData *idata)
      : TALYFEMLIB::CEquation<LENodeData>(false, TALYFEMLIB::kAssembleGaussPoints),
        idata_(idata) {}

  void Solve(double, double) override {}

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe,
                     TALYFEMLIB::ZEROMATRIX<double> &Ae,
                     const double *h) override {
    (void)h;
    const int nbf = fe.nbf();
    const int nd = DIM;
    const int strain_sz = 3 * (DIM - 1);
    const double detJxW = fe.detJxW();

    for (int a = 0; a < nbf; ++a) {
      for (int b = 0; b < nbf; ++b) {
        for (int i = 0; i < nd; ++i) {
          for (int k = 0; k < nd; ++k) {
            double entry = 0.0;
            for (int alpha = 0; alpha < strain_sz; ++alpha) {
              const double Ba = strainComponent(fe, a, i, alpha);
              if (Ba == 0.0) {
                continue;
              }
              for (int beta = 0; beta < strain_sz; ++beta) {
                const double Bb = strainComponent(fe, b, k, beta);
                if (Bb == 0.0) {
                  continue;
                }
                entry += Ba * idata_->Cmatrix[alpha][beta] * Bb;
              }
            }
            Ae(nd * a + i, nd * b + k) += entry * detJxW;
          }
        }
      }
    }
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe,
                     TALYFEMLIB::ZEROARRAY<double> &be,
                     const double *h) override {
    (void)h;
    ZEROPTV bodyForce{0.0, 0.0, 0.0};
    bool forceSet = false;
    CalcForce(fe.position(), bodyForce, forceSet);
    if (!forceSet) {
      return;
    }

    const int nbf = fe.nbf();
    const int nd = DIM;
    const double detJxW = fe.detJxW();
    for (int a = 0; a < nbf; ++a) {
      for (int d = 0; d < nd; ++d) {
        be(nd * a + d) += fe.N(a) * bodyForce[d] * detJxW;
      }
    }
  }

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe,
                          const unsigned int side_idx,
                          const unsigned int id,
                          TALYFEMLIB::ZeroMatrix<double> &Ae,
                          const double *h) override {
    (void)fe;
    (void)side_idx;
    (void)id;
    (void)Ae;
    (void)h;
  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe,
                          const unsigned int side_idx,
                          const unsigned int id,
                          TALYFEMLIB::ZEROARRAY<double> &be,
                          const double *h) override {
    (void)fe;
    (void)side_idx;
    (void)id;
    (void)be;
    (void)h;
  }

private:
  LEInputData *idata_;

  static double pi() {
    return 3.14159265358979323846;
  }

  double strainComponent(const TALYFEMLIB::FEMElm &fe,
                         int node,
                         int comp,
                         int strainIdx) const {
#if (DIM == 2)
    switch (strainIdx) {
    case 0:
      return (comp == 0) ? fe.dN(node, 0) : 0.0;
    case 1:
      return (comp == 1) ? fe.dN(node, 1) : 0.0;
    case 2:
      return (comp == 0) ? fe.dN(node, 1) : fe.dN(node, 0);
    default:
      return 0.0;
    }
#else
    switch (strainIdx) {
    case 0:
      return (comp == 0) ? fe.dN(node, 0) : 0.0;
    case 1:
      return (comp == 1) ? fe.dN(node, 1) : 0.0;
    case 2:
      return (comp == 2) ? fe.dN(node, 2) : 0.0;
    case 3:
      return (comp == 1) ? fe.dN(node, 2) : (comp == 2 ? fe.dN(node, 1) : 0.0);
    case 4:
      return (comp == 0) ? fe.dN(node, 2) : (comp == 2 ? fe.dN(node, 0) : 0.0);
    case 5:
      return (comp == 0) ? fe.dN(node, 1) : (comp == 1 ? fe.dN(node, 0) : 0.0);
    default:
      return 0.0;
    }
#endif
  }

  void CalcForce(const TALYFEMLIB::ZEROPTV &p,
                 ZEROPTV &bodyForce,
                 bool &forceSet) const {
    forceSet = false;

    switch (idata_->SbmGeo) {
    case LEInputData::SBMGeo::PLANT: {
      forceSet = true;
      const double x = p.x();
      const double y = p.y();
      const double pi_v = pi();
      const double E = idata_->planeStrain.young;
      const double nu = idata_->planeStrain.poisson;
      bodyForce.x() = -E * (pi_v * pi_v * (3 * nu + 2) * std::sin(pi_v * x) * std::sin(pi_v * y) +
                            pi_v * pi_v * nu * std::cos(pi_v * x) * std::cos(pi_v * y)) /
                      (100 * (1 + nu) * (1 - 2 * nu));
      bodyForce.y() = -E * (pi_v * pi_v * (3 * nu + 2) * std::cos(pi_v * x) * std::cos(pi_v * y) +
                            pi_v * pi_v * nu * std::sin(pi_v * x) * std::sin(pi_v * y)) /
                      (100 * (1 + nu) * (1 - 2 * nu));
      bodyForce.z() = 0.0;
      break;
    }
    case LEInputData::SBMGeo::SPHERE: {
      forceSet = true;
      const double pi_v = pi();
      const double lam = idata_->lame.lamda;
      const double mu = idata_->lame.mu;
      const double x = p.x();
      const double y = p.y();
      const double z = p.z();
      const double pi2 = pi_v * pi_v;
      bodyForce.x() = -0.1 * lam * pi2 * std::sin(pi_v * x) * std::sin(pi_v * y) * std::sin(pi_v * z) +
                      0.2 * mu * pi2 * std::sin(pi_v * x) * std::sin(pi_v * y) * std::sin(pi_v * z);
      bodyForce.y() = 0.1 * lam * pi2 * std::sin(pi_v * z) * std::cos(pi_v * x) * std::cos(pi_v * y) +
                      0.4 * mu * pi2 * std::sin(pi_v * z) * std::cos(pi_v * x) * std::cos(pi_v * y);
      bodyForce.z() = 0.1 * lam * pi2 * std::sin(pi_v * y) * std::cos(pi_v * x) * std::cos(pi_v * z) +
                      0.4 * mu * pi2 * std::sin(pi_v * y) * std::cos(pi_v * x) * std::cos(pi_v * z);
      break;
    }
    case LEInputData::SBMGeo::BUNNY:
      bodyForce = {0.0, 0.0, -4900.0};
      forceSet = true;
      break;
    case LEInputData::SBMGeo::EIFFEL:
      bodyForce = {0.0, 0.0, -2000.0};
      forceSet = true;
      break;
    default:
      break;
    }
  }
};
