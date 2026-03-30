#pragma once

#include <cmath>
#include <talyfem/fem/cequation.h>
#include <talyfem/talyfem.h>
#include "LEInputData.h"
#include "LENodeData.h"


#include "MaterialMap2D.h"



class SSLEEquation : public TALYFEMLIB::CEquation<LENodeData> {
public:
  explicit SSLEEquation(LEInputData *idata)
      : TALYFEMLIB::CEquation<LENodeData>(false, TALYFEMLIB::kAssembleGaussPoints),
        idata_(idata)
  {

    MaterialMap2D::Domain2D dom;
    dom.x0 = 0;
    dom.x1 = 1;
    dom.y0 = 0;
    dom.y1 = 1;

    // --- Phase properties ---
    MaterialMap2D::PhaseProps props;
    props.E_f = idata_->planeFiberProp.FiberE;
    props.nu_f = idata_->planeFiberProp.Fibermu;
    props.E_m = idata_->planeFiberProp.MatrixE;
    props.nu_m = idata_->planeFiberProp.Matrixmu;

    props.threshold01 = idata_->planeFiberProp.hardThresholdValue;
    props.mode = (!idata_->planeFiberProp.hardThreshold)
                   ? MaterialMap2D::MixingMode::SMOOTH_MIXTURE
                   : MaterialMap2D::MixingMode::HARD_THRESHOLD;

    // --- Load the PNG once ---
    // If it's binary/grayscale you can set forceGray=true (faster + simpler).
    const std::string pngPath = idata_->planeFiberProp.Image_Path;
    auto img = MaterialMap2D::Image::LoadPNG(pngPath, /*forceGray=*/true);

    // flipY: true if your image origin is top-left and y in physics increases upward
    const bool flipY = true;

    materialMap_ = MaterialMap2D::Map(dom, props, std::move(img), flipY);
  }

  MaterialMap2D::MatProps queryMaterialAt(double x, double y) const
  {
    return materialMap_.query(x, y);
  }

  void buildLocalCmatrix(const TALYFEMLIB::ZEROPTV &pos,
                         double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)]) const
  {
    MaterialMap2D::MatProps matProp = materialMap_.query(pos.x(), pos.y());
#if (DIM == 2)
    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = matProp.E;
      double poisson = matProp.nu;
      Cmatrix[0][0] = young / (1 - pow(poisson, 2));
      Cmatrix[0][1] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[1][1] = young / (1 - pow(poisson, 2));
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (2 * (1 + poisson));
    }
#endif

    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = matProp.E;
      double poisson = matProp.nu;
#if (DIM == 2)
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
#endif
#if (DIM == 3)
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }

    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
#if (DIM == 2)
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
#endif
#if (DIM == 3)
      double young = mu * (3 * lamda + 2 * mu) / (mu + lamda);
      double poisson = lamda / 2 / (lamda + mu);

      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }
  }

  void Solve(double, double) override {}

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe,
                     TALYFEMLIB::ZEROMATRIX<double> &Ae,
                     const double *h) override {
    (void)h;
    const int nbf = fe.nbf();
    const int nd = DIM;
    const int strain_sz = 3 * (DIM - 1);
    const double detJxW = fe.detJxW();

    double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];

    CalcCmatrixFiber(fe, Cmatrix);

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
                entry += Ba * Cmatrix[alpha][beta] * Bb;
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
  MaterialMap2D::Map materialMap_;

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

    // switch (idata_->SbmGeo) {
    // case LEInputData::SBMGeo::PLANT: {
    //   forceSet = true;
    //   const double x = p.x();
    //   const double y = p.y();
    //   const double pi_v = pi();
    //   const double E = idata_->planeStrain.young;
    //   const double nu = idata_->planeStrain.poisson;
    //   bodyForce.x() = -E * (pi_v * pi_v * (3 * nu + 2) * std::sin(pi_v * x) * std::sin(pi_v * y) +
    //                         pi_v * pi_v * nu * std::cos(pi_v * x) * std::cos(pi_v * y)) /
    //                   (100 * (1 + nu) * (1 - 2 * nu));
    //   bodyForce.y() = -E * (pi_v * pi_v * (3 * nu + 2) * std::cos(pi_v * x) * std::cos(pi_v * y) +
    //                         pi_v * pi_v * nu * std::sin(pi_v * x) * std::sin(pi_v * y)) /
    //                   (100 * (1 + nu) * (1 - 2 * nu));
    //   bodyForce.z() = 0.0;
    //   break;
    // }
    // case LEInputData::SBMGeo::SPHERE: {
    //   forceSet = true;
    //   const double pi_v = pi();
    //   const double lam = idata_->lame.lamda;
    //   const double mu = idata_->lame.mu;
    //   const double x = p.x();
    //   const double y = p.y();
    //   const double z = p.z();
    //   const double pi2 = pi_v * pi_v;
    //   bodyForce.x() = -0.1 * lam * pi2 * std::sin(pi_v * x) * std::sin(pi_v * y) * std::sin(pi_v * z) +
    //                   0.2 * mu * pi2 * std::sin(pi_v * x) * std::sin(pi_v * y) * std::sin(pi_v * z);
    //   bodyForce.y() = 0.1 * lam * pi2 * std::sin(pi_v * z) * std::cos(pi_v * x) * std::cos(pi_v * y) +
    //                   0.4 * mu * pi2 * std::sin(pi_v * z) * std::cos(pi_v * x) * std::cos(pi_v * y);
    //   bodyForce.z() = 0.1 * lam * pi2 * std::sin(pi_v * y) * std::cos(pi_v * x) * std::cos(pi_v * z) +
    //                   0.4 * mu * pi2 * std::sin(pi_v * y) * std::cos(pi_v * x) * std::cos(pi_v * z);
    //   break;
    // }
    // case LEInputData::SBMGeo::BUNNY:
    //   bodyForce = {0.0, 0.0, -4900.0};
    //   forceSet = true;
    //   break;
    // case LEInputData::SBMGeo::EIFFEL:
    //   bodyForce = {0.0, 0.0, -2000.0};
    //   forceSet = true;
    //   break;
    // default:
    //   break;
    // }
  }

  void CalcCmatrixFiber(const TALYFEMLIB::FEMElm &fe,double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)] ) const
  {
    buildLocalCmatrix(fe.position(), Cmatrix);
  }

    void CalcCmatrix(double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)])
  {
    /*
     * 3D do not have plane stress case
     */
#if (DIM == 2)
    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      // C for plane stress
      Cmatrix[0][0] = young / (1 - pow(poisson, 2));
      Cmatrix[0][1] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[1][1] = young / (1 - pow(poisson, 2));
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (2 * (1 + poisson));
    }
#endif

    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      // C for plane strain
#if (DIM == 2)
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
#endif
#if (DIM == 3)
      /*
       * this formulation is from FEM book page 241
       */
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      // [previous bug here]
      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      // C for lame parameters
#if (DIM == 2)
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
#endif
#if (DIM == 3)
      double young = mu * (3 * lamda + 2 * mu) / (mu + lamda);
      double poisson = lamda / 2 / (lamda + mu);

      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][3] = 0;
      Cmatrix[0][4] = 0;
      Cmatrix[0][5] = 0;

      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][3] = 0;
      Cmatrix[1][4] = 0;
      Cmatrix[1][5] = 0;

      Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[2][3] = 0;
      Cmatrix[2][4] = 0;
      Cmatrix[2][5] = 0;

      Cmatrix[3][0] = 0;
      Cmatrix[3][1] = 0;
      Cmatrix[3][2] = 0;
      Cmatrix[3][3] = young / (1 + poisson) / 2;
      Cmatrix[3][4] = 0;
      Cmatrix[3][5] = 0;

      // [previous bug here]
      Cmatrix[4][0] = 0;
      Cmatrix[4][1] = 0;
      Cmatrix[4][2] = 0;
      Cmatrix[4][3] = 0;
      Cmatrix[4][4] = young / (1 + poisson) / 2;
      Cmatrix[4][5] = 0;

      Cmatrix[5][0] = 0;
      Cmatrix[5][1] = 0;
      Cmatrix[5][2] = 0;
      Cmatrix[5][3] = 0;
      Cmatrix[5][4] = 0;
      Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
    }
  }
};
