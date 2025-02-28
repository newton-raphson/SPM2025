#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"

class NSEquation : public TALYFEMLIB::CEquation<NSNodeData> {
  TALYFEMLIB::ZEROPTV forcing_;
  TALYFEMLIB::ZEROPTV forcingPre_;
  DENDRITE_REAL theta_ = 1.0;
  NSInputData *idata_;

protected:

    struct VMSParams {
        double tauM;
        double tauC;

        VMSParams() : tauM(1.0), tauC(1.0) {}
    };


 public:
  explicit NSEquation(NSInputData *idata)
      : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
    idata_ = idata;
    if (idata->timeStepping == NSInputData::TIME_STEPPING::CRANK_NICHOLSON) {
      theta_ = 0.5;
    }
  }
  void Solve(double dt, double t) override {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override {
    assert(false);
  }


  void calcAe_vms(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const std::vector<double> &bdf) {
    using namespace TALYFEMLIB;
    double Re = idata_->Re;
    double Ci_f = 36;
    const int nsd = DIM;
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double dt = dt_;
    const double lambda = - 2.0 / 3.0;

    if (idata_->ifMMS) {
      calcForcing(forcing_, fe.position(), t_);
    }

      VMSParams vms_params = calc_tau(fe);


      /// Calculate the tau_M based on the projections calculated, Equation 63 in
      /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
      /// in the paper it is set to 1.
      const double tauM = vms_params.tauM;
      /// Calculate continuity residual based on
      const double tauC = vms_params.tauC;

//    MPITimer timer("total Ae time");
//    timer.Start();
//
//
      const double temp_coeff = 1.0;
      const double adv_coeff = 1.0;
      const double press_coeff = 1.0;
      const double diff_coeff = 1.0 / Re;
      const double force_coeff = 1.0;


      /// time_scheme
      double force[nsd];

      for (int i = 0; i < nsd; i++) {
          force[i] = forcing_(i);
      }


      /// Calculate the tau_M based on the projections calculated, Equation 63 in
      /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
      /// in the paper it is set to 1.
//    const double tauM = vms_params.tauM;
      /// Calculate continuity residual based on
//    const double tauC = vms_params.tauC;

      /// field variables
      /** Get u and u_pre from the NodeData
       * NOTE: ValueFEM is defined on NodeData object, therefore it follows
       * indices defined in NodeData subclass.
       * NOTE: Every iteration the solution is copied into NodeData, we use these
       * solved fields from the NodeData object to calculate the be (NS residual)
       * of current
       * iteration
       */


      ////////////////////////////// precalculate velocity terms to simplify the expressions


      double u[nsd];
      double u_pre[nsd];
      double u_pre_pre[nsd];

      double res_M[nsd];
      for (int i = 0; i < nsd; i++) {
          u[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
          u_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
          u_pre_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);

      }

      /// Define velocity gradient tensor
      double du[nsd][nsd];
      double dp[nsd];
      double ui_ujj[nsd];
      double uj_uij[nsd];
      for (int i = 0; i < nsd; i++) {
          dp[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
          for (int j = 0; j < nsd; j++) {
              du[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
          }
      }


      double d2u_ijj[nsd];
      double d2u_jij[nsd];
      /// loop for three directions of velocity (MAX NSD is no of velocity
      /// directions in the problem)
      /// PS: This assumes the first three degrees are always velocity, it is
      /// wise to keep it that way
      for (int dof = 0; dof < nsd; dof++) {
          /// Summing over three directions of velocity as it is laplacian
          d2u_ijj[dof] = 0;
          d2u_jij[dof] = 0;
          for (int dir = 0; dir < nsd; dir++) {
              d2u_ijj[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dof, dir, dir);
              d2u_jij[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dir, dof, dir);
          }
      }


      double ukk;
      for (int i = 0; i < nsd; i++) {
          ui_ujj[i] = 0;
          uj_uij[i] = 0;
          ukk = 0;
          for (int j = 0; j < nsd; j++) {
              ui_ujj[i] += u[i] * du[j][j];
              uj_uij[i] += u[j] * du[i][j];
              ukk += du[j][j];
          }
          res_M[i] = temp_coeff * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt + adv_coeff * (ui_ujj[i] + uj_uij[i]) +
                     press_coeff * dp[i] - diff_coeff * (d2u_ijj[i] + d2u_jij[i]) - force_coeff * force[i];

          /// since d2u_jij = d2u_kik we don't need to calculate it again
          if (idata_->SecondViscosity){
              res_M[i] -= diff_coeff * lambda * d2u_jij[i];
          }
      }


      ////////////////////////////////// end of precalculations




      ////////////////////////////////// calculate the derivative of the momentunm residual


      double dres_du_M_Diag[nsd][nbf];
      double dres_du_M_offDiag[nsd][nsd][nbf];

      for (int a = 0; a < nbf; a++) {

          for (int i = 0; i < nsd; i++) {

              dres_du_M_Diag[i][a] = temp_coeff * (fe.N(a) / dt);

              for (int j = 0; j < nsd; j++) {
                  dres_du_M_Diag[i][a] += adv_coeff * ((du[j][j] * fe.N(a)) + (u[j] * fe.dN(a, j)))  -
                                       (diff_coeff * fe.d2N(a, j, j));



                  dres_du_M_offDiag[i][j][a] =
                          adv_coeff * ((u[i] * fe.dN(a, j)) + (du[i][j] * fe.N(a))) -
                          (diff_coeff * fe.d2N(a, i, j));



              } // end of j loop

              if (idata_->SecondViscosity){
                  dres_du_M_Diag[i][a] -= diff_coeff * lambda * fe.d2N(a, i, i);
              }

          } // end of i loop
      } // end of a loop


      ////////////////////////////////// end of calculation of the derivative of the momentunm residual

//    printf("inside the integrands\n");


      for (int a = 0; a < nbf; a++) {
          for (int b = 0; b < nbf; b++) {
              for (int i = 0; i < nsd; i++) {


                  /// first part of momentum equation
                  /// Temporal term -- Term No. in the Document= 1
                  Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += temp_coeff * fe.N(a) * (fe.N(b) / dt) * detJxW; ///ok!



                  for (int j = 0; j < nsd; j++) {


                      /// coarseAdvection 1 -- Term No. in the Document= 2
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * (u[j] * fe.N(b)) * detJxW; ///ok!

                      /// coarseAdvection 2 -- Term No. in the Document= 3
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * (u[i] * fe.N(b)) * detJxW; ///ok!

                      /// fineAdvection 1 -- Term No. in the Document= 4
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * (tauM * res_M[j] * fe.N(b)) * detJxW; ///ok!

                      /// fineAdvection 2 -- Term No. in the Document= 5
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * u[i] * tauM * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) * detJxW; ///ok!


                      /// fineAdvection 3a  -- Term No. in the Document= 6
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * tauM * dres_du_M_Diag[i][b] * u[j] * detJxW; ///ok!

                      /// fineAdvection 3b  -- Term No. in the Document= 7
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * tauM * dres_du_M_offDiag[i][j][b] * u[j] * detJxW; ///ok!


                      /// fineAdvection 4 -- Term No. in the Document= 8
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * (tauM * res_M[i] * fe.N(b)) * detJxW; ///ok!



                      /// fineAdvection 5a -- Term No. in the Document= 9
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * tauM * tauM * dres_du_M_Diag[i][b] * res_M[j] * detJxW; ///ok!


                      /// fineAdvection 5b -- Term No. in the Document= 10
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * tauM * tauM * dres_du_M_offDiag[i][j][b] * res_M[j] * detJxW; ///ok!


                      /// fineAdvection 6  -- Term No. in the Document= 11
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * tauM * tauM * res_M[i] * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) *  detJxW; ///ok!

                      /// finePressure 1 -- Term No. in the Document= 12
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += press_coeff * fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW; ///ok!

                      /// coarseDiffusion 1 -- Term No. in the Document= 13
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * fe.dN(a, j) * fe.dN(b, j) * detJxW; ///ok!

                      /// coarseDiffusion 2-- Term No. in the Document= 14
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.dN(a, j) * fe.dN(b, i) * detJxW; ///ok!


                      if (idata_->DiffFineTerms) {
                          /// fineDiffusion 1  -- Term No. in the Document= 15
                          Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.d2N(a, j, i) * tauM * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) *  detJxW; ///ok!

                          /// second part of momentum equation
                          /// fineDiffusion 3-- Term No. in the Document= 23
                          Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += diff_coeff * fe.d2N(a, i, j) * tauM * fe.dN(b, j) * detJxW;

                      }


                      /// second part of momentum equation
                      /// fineAdvection 7 -- Term No. in the Document= 18
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a, j) * u[i] * tauM * fe.dN(b, j) * detJxW; ///ok!


                      /// fineAdvection 8 -- Term No. in the Document= 19
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a, j) * u[j] * tauM * fe.dN(b, i) * detJxW; ///ok!

                      /// fineAdvection 9 -- Term No. in the Document= 20
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= adv_coeff * fe.dN(a, j) * tauM * tauM * fe.dN(b, i) *
                                                                    res_M[j] * detJxW; ///ok!

                      /// fineAdvection 10 -- Term No. in the Document= 21
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= adv_coeff * fe.dN(a, j) * tauM * tauM * fe.dN(b, j) *
                                                                    res_M[i] * detJxW; ///ok!


                      /// fineContinuityVelocity-a -- Term No. in the Document= 26
                      Ae((nsd + 1) * a + nsd, (nsd + 1) * b + j) += fe.dN(a, i) * tauM * dres_du_M_offDiag[i][j][b] * detJxW; ///ok!



                  } // end of j loop


                  /// coarsePressure 1 -- Term No. in the Document= 22
                  Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= press_coeff * fe.dN(a, i) * fe.N(b) * detJxW; ///ok!

                  /// first part of continuity equation
                  /// coarseContinuityVelocity -- Term No. in the Document= 25
                  Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.N(a) * fe.dN(b, i) * detJxW; ///ok!


                  /// fineContinuityVelocity-b -- Term No. in the Document= 27
                  Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * dres_du_M_Diag[i][b] * detJxW; ///ok!


                  /// second part of continuity equation
                  /// fineContinuityPressure -- Term No. in the Document= 28
                  Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) += fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW; ///ok!



                  if (idata_->SecondViscosity) {

                      /// coarseDiffusion 3-- Term No. in the Document= 16
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * lambda * fe.dN(a, i) * fe.dN(b, i) * detJxW; ///ok!

                      /// fineDiffusion 2-- Term No. in the Document= 17
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * lambda * fe.d2N(a, i, i) * tauM * dres_du_M_Diag[i][b] * detJxW; ///ok!

                      /// second part of continuity equation
                      /// fineDiffusion 4-- Term No. in the Document= 24
                      Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += diff_coeff * lambda * fe.d2N(a, i, i) * tauM * fe.dN(b, i) * detJxW; ///ok!

                  }





              } // end of i loop
          } // end of b loop
      } // end of a loop

//
//    timer.Stop();
//    timer.PrintTotalTimeSeconds();





//    MPITimer timer("total Ae time");
//    timer.Start();



//    timer.Stop();
//    timer.PrintTotalTimeSeconds();


  }

  void calcbe_vms(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const std::vector<double> &bdf) {
    using namespace TALYFEMLIB;
      double Re = idata_->Re;
      double Ci_f = 36;
      const int nsd = DIM;
      const int nbf = fe.nbf();
      const double detJxW = fe.detJxW();
      const double dt = dt_;
      const double lambda = - 2.0 / 3.0;

      if (idata_->ifMMS) {
          calcForcing(forcing_, fe.position(), t_);
      }

      VMSParams vms_params = calc_tau(fe);

      /// Calculate the tau_M based on the projections calculated, Equation 63 in
      /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
      /// in the paper it is set to 1.
      const double tauM = vms_params.tauM;
      /// Calculate continuity residual based on
      const double tauC = vms_params.tauC;

//    MPITimer timer("total Ae time");
//    timer.Start();
//
//
      const double temp_coeff = 1.0;
      const double adv_coeff = 1.0;
      const double press_coeff = 1.0;
      const double diff_coeff = 1.0 / Re;
      const double force_coeff = 1.0;


      /// time_scheme
      double force[nsd];

      for (int i = 0; i < nsd; i++) {
          force[i] = forcing_(i);
      }

      /// Calculate the tau_M based on the projections calculated, Equation 63 in
      /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
      /// in the paper it is set to 1.
//    const double tauM = vms_params.tauM;
      /// Calculate continuity residual based on
//    const double tauC = vms_params.tauC;

      /// field variables
      /** Get u and u_pre from the NodeData
       * NOTE: ValueFEM is defined on NodeData object, therefore it follows
       * indices defined in NodeData subclass.
       * NOTE: Every iteration the solution is copied into NodeData, we use these
       * solved fields from the NodeData object to calculate the be (NS residual)
       * of current
       * iteration
       */


      ////////////////////////////// precalculate velocity terms to simplify the expressions


      double u[nsd];
      double u_pre[nsd];
      double u_pre_pre[nsd];
      double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
      for (int i = 0; i < nsd; i++) {
          u[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
          u_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
          u_pre_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);

      }

      /// Define velocity gradient tensor
      double du[nsd][nsd];
      double dp[nsd];
      double ui_ujj[nsd];
      double uj_uij[nsd];
      for (int i = 0; i < nsd; i++) {
          dp[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
          for (int j = 0; j < nsd; j++) {
              du[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
          }
      }


      double d2u_ijj[nsd];
      double d2u_jij[nsd];
      /// loop for three directions of velocity (MAX NSD is no of velocity
      /// directions in the problem)
      /// PS: This assumes the first three degrees are always velocity, it is
      /// wise to keep it that way
      for (int dof = 0; dof < nsd; dof++) {
          /// Summing over three directions of velocity as it is laplacian
          d2u_ijj[dof] = 0;
          d2u_jij[dof] = 0;
          for (int dir = 0; dir < nsd; dir++) {
              d2u_ijj[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dof, dir, dir);
              d2u_jij[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dir, dof, dir);

          }
      }






      double res_M[nsd];
      double res_C = 0;
      double ukk;
      for (int i = 0; i < nsd; i++) {
          ui_ujj[i] = 0;
          uj_uij[i] = 0;
          ukk = 0;
          for (int j = 0; j < nsd; j++) {
              ui_ujj[i] += u[i] * du[j][j];
              uj_uij[i] += u[j] * du[i][j];
              ukk += du[j][j];
          }
          res_M[i] = temp_coeff * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt + adv_coeff * (ui_ujj[i] + uj_uij[i]) +
                     press_coeff * dp[i] - diff_coeff * (d2u_ijj[i] + d2u_jij[i]) - force_coeff * force[i];
          if (idata_->SecondViscosity){
              res_M[i] -= diff_coeff * lambda * d2u_jij[i];
          }
          res_C += du[i][i];
      }


      ////////////////////////////////// end of precalculations



      double R_NS[nsd];
      double R_C;
      for (int a = 0; a < nbf; a++) {
          R_C = 0;
          for (int i = 0; i < nsd; i++) {
              R_NS[i] = 0;
              for (int j = 0; j < nsd; j++) {
                  /// expression:  R_NS[i] += advection[i][j] + diffusion[i][j]
                  R_NS[i] += adv_coeff * ((-fe.dN(a, j) * u[i] * u[j]) + (fe.dN(a, j) * u[i] * tauM * res_M[j]) +
                                          (fe.dN(a, j) * tauM * res_M[i] * u[j]) - (fe.dN(a, j) * tauM * res_M[i] * tauM * res_M[j])) +
                             diff_coeff * (fe.dN(a, j) * (du[i][j] + du[j][i]));
                  if (idata_->DiffFineTerms){
                      R_NS[i] += diff_coeff *  (fe.d2N(a, j, i) * tauM * res_M[j]);
                  }

                  if (idata_->SecondViscosity) {
                      R_NS[i] += diff_coeff * lambda * (fe.dN(a, i) * du[j][j] + fe.d2N(a, i, j) * tauM * res_M[j]);
                  }

              } // end of j loop
              /// expression:  R_NS[i] += temporal[i] - pressure[i] - force[i]
              R_NS[i] += temp_coeff * (fe.N(a) * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt) +
                         press_coeff * ((-fe.dN(a, i) * p) + (fe.dN(a, i) * tauC * res_C)) -
                         force_coeff * fe.N(a) * force[i];
              R_C += (fe.N(a) * du[i][i]) + (fe.dN(a, i) * tauM * res_M[i]);

              be((nsd + 1) * a + i) +=  R_NS[i] * detJxW;

          } // end of i loop
          be((nsd + 1) * a + nsd) += R_C * detJxW;
      } // end of a loop

  }


    VMSParams calc_tau(const TALYFEMLIB::FEMElm &fe) const {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;
        const double Ci_f = 36;
        const int nsd = DIM;

        VMSParams params;

        ZEROPTV u;
        for (int i = 0; i < nsd; i++) {
            u(i) = this->p_data_->valueFEM(fe, i);
        }

        ZeroMatrix<double> ksiX;
        ksiX.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                ksiX(i, j) = fe.cof(j, i) / fe.jacc();
            }
        }

        ZeroMatrix<double> Ge;
        Ge.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                Ge(i, j) = 0.0;
                for (int k = 0; k < nsd; k++)
                    Ge(i, j) += ksiX(k, i) * ksiX(k, j);
            }
        }

        double u_Gu = 0.0;
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                u_Gu += u(i) * Ge(i, j) * u(j);
            }
        }

        double G_G_u = 0.0;
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
            }
        }

        params.tauM = 1.0 / sqrt(4.0 / (this->dt_ * this->dt_) + u_Gu + G_G_u);

        ZEROARRAY<double> ge;
        ge.redim(nsd);
        for (int i = 0; i < nsd; i++) {
            ge(i) = 0.0;
            for (int j = 0; j < nsd; j++)
                ge(i) += ksiX(j, i);
        }

        double g_g = 0.0;
        for (int i = 0; i < nsd; i++) {
            g_g += ge(i) * ge(i);
        }

        params.tauC = 1.0 / (params.tauM * g_g);

        return params;
    }



    void ReRamping (double t) {
        if (t <= idata_->RampingTime) {
            idata_->Re = ((idata_->RampingRe - 10) / idata_->RampingTime) * t + 10; /*ramping from Re = 10 */
        } else {
            idata_->Re = idata_->Re_tmp;
        }
    }


  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {

      if (idata_->DoReRamping) {
          ReRamping(t_);
      }
      std::vector<double> bdf;
    if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
            t_ < 1.5 * dt_) {
        bdf = {1.0, -1.0, 0.0};
    } else if (idata_->timeStepping == NSInputData::BDF2) {
        bdf = {1.5, -2.0, 0.5};
    }
      calcAe_vms(fe, Ae, bdf);
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
      if (idata_->DoReRamping) {
          ReRamping(t_);
      }
      std::vector<double> bdf;
    if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
        t_ < 1.5 * dt_) {
        bdf = {1.0, -1.0, 0.0};
    } else if (idata_->timeStepping == NSInputData::BDF2) {
        bdf = {1.5, -2.0, 0.5};
    }
      calcbe_vms(fe, be, bdf);
    }



  void calcForcing(TALYFEMLIB::ZEROPTV &forcing, const TALYFEMLIB::ZEROPTV &location, const double time) {
    double Re(idata_->Re);
    forcing.x() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI)
        + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI) * 2.0
        + ((M_PI * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) *
            2.0) / Re
        + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
            (cos(location.y() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0)
        + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
            (cos(location.x() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0);

    forcing.y() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI) * (-2.0)
        + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI)
        - ((M_PI * M_PI) * cos(location.x() * M_PI) * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) *
            2.0) / Re
        + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
            (cos(location.x() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0)
        + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
            (cos(location.y() * M_PI)
                * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0);
    forcing.z() = 0;

  }

};

#endif //DENDRITEKT_NSEQUATION_H
