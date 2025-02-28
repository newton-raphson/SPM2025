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
        double force_pre[nsd];
        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
            force_pre[i] = forcing_(i);
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


        double u_pre1[nsd];
        double u_pre2[nsd];
        double u_pre3[nsd];


        for (int i = 0; i < nsd; i++) {
            u_pre1[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre2[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
            u_pre3[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE3 + i);

        }

        std::vector<double> u_extrapolated = calc_linear_advection(u_pre1, u_pre2, idata_->VelocityExtrapolationOrder);
        VMSParams vms_params = calc_tau(fe, u_extrapolated);


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        /// Define velocity gradient tensor
        double du_pre1[nsd][nsd];
        double du_pre2[nsd][nsd];
        double dp_pre1[nsd];
        double dp_pre2[nsd];
        for (int i = 0; i < nsd; i++) {
            dp_pre1[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
            dp_pre2[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
            for (int j = 0; j < nsd; j++) {
                du_pre1[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
                du_pre2[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + i, j);
            }
        }


        double d2u_ijj_pre1[nsd];
        double d2u_jij_pre1[nsd];
        double d2u_ijj_pre2[nsd];
        double d2u_jij_pre2[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj_pre1[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dof, dir, dir);
                d2u_jij_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dir, dof, dir);
                d2u_ijj_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dof, dir, dir);
                d2u_jij_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dir, dof, dir);

            }
        }

        double u_full_pre1[nsd];
        double u_full_pre2[nsd];
        double res_M_pre1[nsd], res_M_pre2[nsd];
        double ui_ujj, uj_uij, ui_ujj_pre, uj_uij_pre;
        double ukk, ukk_pre;
        for (int i = 0; i < nsd; i++) {
            ui_ujj = 0;
            uj_uij = 0;
            ui_ujj_pre = 0;
            uj_uij_pre = 0;
            ukk = 0;
            ukk_pre = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj += u_pre1[i] * du_pre1[j][j];
                uj_uij += u_pre1[j] * du_pre1[i][j];
                ui_ujj_pre += u_pre2[i] * du_pre2[j][j];
                uj_uij_pre += u_pre2[j] * du_pre2[i][j];
                ukk += du_pre1[j][j];
                ukk_pre += du_pre2[j][j];
            }
            res_M_pre1[i] = temp_coeff * (bdf[0] * u_pre1[i] + bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt + adv_coeff * (ui_ujj + 0.5 * uj_uij) +
                            press_coeff * dp_pre1[i] - diff_coeff * (d2u_ijj_pre1[i] + d2u_jij_pre1[i]) - force_coeff * force[i];

            /// for temporal term we consider backward euler scheme to avoid storing one more degree of freedom
            res_M_pre2[i] = temp_coeff * (u_pre2[i] - u_pre3[i]) / dt + adv_coeff * (ui_ujj_pre + 0.5 * uj_uij_pre) +
                            press_coeff * dp_pre2[i] - diff_coeff * (d2u_ijj_pre2[i] + d2u_jij_pre2[i]) - force_coeff * force_pre[i];

            if (idata_->SecondViscosity) {
                res_M_pre1[i] -= diff_coeff * lambda * ukk;
                res_M_pre2[i] -= diff_coeff * lambda * ukk_pre;
            }


            u_full_pre1[i] = u_pre1[i] - tauM * res_M_pre1[i];
            u_full_pre2[i] = u_pre2[i] - tauM * res_M_pre2[i];

        }


        std::vector<double> u_star = calc_linear_advection(u_full_pre1, u_full_pre2, idata_->VelocityExtrapolationOrder);

//        if (u[0] != 0)
//        printf ("u_star = %f %f\n", u[0], Re);




        ////////////////////////////////// end of precalculations




        ////////////////////////////////// calculate the derivative of the momentunm residual


        double res_M_Diag[nsd][nbf];
        double res_M_offDiag[nsd][nsd][nbf];

        for (int a = 0; a < nbf; a++) {
            for (int i = 0; i < nsd; i++) {
                res_M_Diag[i][a] = temp_coeff * (bdf[0] * fe.N(a) / dt);

                for (int j = 0; j < nsd; j++) {
                    res_M_Diag[i][a] += adv_coeff * (u_star[j] *  fe.dN(a, j) + 0.5 * du_pre1[j][j] * fe.N(a)) - diff_coeff * fe.d2N(a, i, j);
                    res_M_offDiag[i][j][a] = - diff_coeff *  fe.d2N(a, j, i);

                    if (idata_->SecondViscosity) {
                        res_M_offDiag[i][j][a] -= diff_coeff * lambda * fe.d2N(a, j, j);
                    }
                } // end of j loop
            } // end of i loop
        } // end of a loop


        ////////////////////////////////// end of calculation of the momentunm residual


//    printf("inside the integrands\n");


        for (int a = 0; a < nbf; a++) {
            for (int b = 0; b < nbf; b++) {
                for (int i = 0; i < nsd; i++) {


                    /// first part of momentum equation
                    /// Temporal term -- Term No. in the Document= 1
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += temp_coeff * fe.N(a) * (fe.N(b) / dt) * detJxW; ///ok!



                    for (int j = 0; j < nsd; j++) {


                        /// coarseAdvection 1 -- Term No. in the Document= 2
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * (u_star[j] * fe.N(b)) * detJxW; ///ok!



                        /// fineAdvection 1a -- Term No. in the Document= 3
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * res_M_Diag[i][b] * detJxW; ///ok!


                        /// fineAdvection 1b -- Term No. in the Document= 4
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * res_M_offDiag[i][j][b] * detJxW; ///ok!


                        /// fineAdvection 1c -- Term No. in the Document= 5
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * fe.dN(b, i) * detJxW; ///ok!


                        /// coarseAdvection 2 -- Term No. in the Document= 6
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * fe.N(b) * detJxW; ///ok!


                        /// fineAdvection 2a -- Term No. in the Document= 7
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * res_M_Diag[i][b] * detJxW; /// ok!


                        /// fineAdvection 2b -- Term No. in the Document= 8
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * res_M_offDiag[i][j][b] * detJxW; /// ok!


                        /// fineAdvection 2c -- Term No. in the Document= 9
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * fe.dN(b, i) * detJxW; /// ok!


                        /// coarseDiffusion 1 -- Term No. in the Document= 12
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * fe.dN(a,j) * fe.dN(b, j) * detJxW; ///ok!

                        /// coarseDiffusion 2 -- Term No. in the Document= 13
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.dN(a,j) * fe.dN(b, i) * detJxW; ///ok!

                        if (idata_->SecondViscosity) {
                            /// coarseDiffusion 3 -- Term No. in the Document= 14
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * lambda * fe.dN(a,j) * fe.dN(b, j) * detJxW; ///ok!

                            /// fineDiffusion 2a -- Term No. in the Document= 17
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                    diff_coeff * lambda * fe.d2N(a, j, i) * tauM * (res_M_Diag[j][b] + res_M_offDiag[j][j][b]) *
                                    detJxW; ///ok!

                            /// fineDiffusion 2b -- Term No. in the Document= 18
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                                    diff_coeff * lambda * fe.d2N(a, j, i) * tauM * fe.dN(b, i) * detJxW; ///ok!
                        }



                        if (idata_->DiffFineTerms) {
                            /// fineDiffusion 1a -- Term No. in the Document= 15
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                    diff_coeff * fe.d2N(a, j, i) * tauM * (res_M_Diag[j][b] + res_M_offDiag[j][j][b]) *
                                    detJxW; ///ok!

                            /// fineDiffusion 1b -- Term No. in the Document= 16
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                                    diff_coeff * fe.d2N(a, j, i) * tauM * fe.dN(b, i) * detJxW; ///ok!

                        }

                        /// fineContinuity 1b -- Term No. in the Document= 21
                        Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * res_M_offDiag[i][j][b] * detJxW; ///ok!



                    } // end of j loop



                    /// coarsePressure 1 -- Term No. in the Document= 10
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= press_coeff * fe.dN(a, i) * fe.N(b) * detJxW; ///ok!

                    /// finePressure 1 -- Term No. in the Document= 11
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += press_coeff * fe.dN(a, i) * tauC * fe.dN(b, i) * detJxW; ///ok!


                    /// coarseContinuity 1 -- Term No. in the Document= 19
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.N(a) * fe.dN(b, i) * detJxW; ///ok!


                    /// fineContinuity 1a -- Term No. in the Document= 20
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * res_M_Diag[i][b] * detJxW; ///ok!


                    /// fineContinuity 1c -- Term No. in the Document= 22
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) += fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW; ///ok!


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
        double force_pre[nsd];
        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
            force_pre[i] = forcing_(i);
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


        double u_pre1[nsd];
        double u_pre2[nsd];
        double u_pre3[nsd];


        for (int i = 0; i < nsd; i++) {
            u_pre1[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre2[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
            u_pre3[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE3 + i);

        }

        std::vector<double> u_extrapolated = calc_linear_advection(u_pre1, u_pre2, idata_->VelocityExtrapolationOrder);
        VMSParams vms_params = calc_tau(fe, u_extrapolated);


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        /// Define velocity gradient tensor
        double du_pre1[nsd][nsd];
        double du_pre2[nsd][nsd];
        double dp_pre1[nsd];
        double dp_pre2[nsd];
        for (int i = 0; i < nsd; i++) {
            dp_pre1[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
            dp_pre2[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
            for (int j = 0; j < nsd; j++) {
                du_pre1[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
                du_pre2[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + i, j);
            }
        }


        double d2u_ijj_pre1[nsd];
        double d2u_jij_pre1[nsd];
        double d2u_ijj_pre2[nsd];
        double d2u_jij_pre2[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj_pre1[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dof, dir, dir);
                d2u_jij_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dir, dof, dir);
                d2u_ijj_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dof, dir, dir);
                d2u_jij_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dir, dof, dir);

            }
        }

        double u_full_pre1[nsd];
        double u_full_pre2[nsd];
        double res_M_pre1[nsd], res_M_pre2[nsd];
        double ui_ujj, uj_uij, ui_ujj_pre, uj_uij_pre;
        double ukk, ukk_pre;
        for (int i = 0; i < nsd; i++) {
            ui_ujj = 0;
            uj_uij = 0;
            ui_ujj_pre = 0;
            uj_uij_pre = 0;
            ukk = 0;
            ukk_pre = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj += u_pre1[i] * du_pre1[j][j];
                uj_uij += u_pre1[j] * du_pre1[i][j];
                ui_ujj_pre += u_pre2[i] * du_pre2[j][j];
                uj_uij_pre += u_pre2[j] * du_pre2[i][j];
                ukk += du_pre1[j][j];
                ukk_pre += du_pre2[j][j];
            }
            res_M_pre1[i] = temp_coeff * (bdf[0] * u_pre1[i] + bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt + adv_coeff * (ui_ujj + 0.5 * uj_uij) +
                       press_coeff * dp_pre1[i] - diff_coeff * (d2u_ijj_pre1[i] + d2u_jij_pre1[i]) - force_coeff * force[i];

            /// for temporal term we consider backward euler scheme to avoid storing one more degree of freedom
            res_M_pre2[i] = temp_coeff * (u_pre2[i] - u_pre3[i]) / dt + adv_coeff * (ui_ujj_pre + 0.5 * uj_uij_pre) +
                           press_coeff * dp_pre2[i] - diff_coeff * (d2u_ijj_pre2[i] + d2u_jij_pre2[i]) - force_coeff * force_pre[i];

            if (idata_->SecondViscosity) {
                res_M_pre1[i] -= diff_coeff * lambda * ukk;
                res_M_pre2[i] -= diff_coeff * lambda * ukk_pre;
            }

            u_full_pre1[i] = u_pre1[i] - tauM * res_M_pre1[i];
            u_full_pre2[i] = u_pre2[i] - tauM * res_M_pre2[i];

        }


        std::vector<double> u_star = calc_linear_advection(u_full_pre1, u_full_pre2, idata_->VelocityExtrapolationOrder);



        ////////////////////////////////// end of precalculations




        double res_M_RHS[nsd];
        for (int i = 0; i < nsd; i++) {
            res_M_RHS[i] = temp_coeff * (bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt - force_coeff * force[i];
        }


        for (int a = 0; a < nbf; a++) {
            for (int i = 0; i < nsd; i++) {
                for (int j = 0; j < nsd; j++) {

                    /// fineAdvection RHS 1 -- Term No. in the Document= 1
                    be((nsd + 1) * a + i) -= fe.dN(a, j) * adv_coeff * u_star[i] * tauM * res_M_RHS[i] * detJxW; ///ok!

                    /// fineAdvection RHS 2 -- Term No. in the Document= 2
                    be((nsd + 1) * a + i) -= fe.N(a) * adv_coeff * tauM * 0.5 * du_pre1[j][j] * res_M_RHS[i] * detJxW; ///ok!

                    if (idata_->DiffFineTerms) {
                        /// fineDiffusion RHS 1 -- Term No. in the Document= 3
                        be((nsd + 1) * a + i) -= diff_coeff * fe.d2N(a, j, i) * tauM * res_M_RHS[j] * detJxW; ///ok!
                    }

                    if (idata_->SecondViscosity) {

                        /// fineDiffusion RHS 2 -- Term No. in the Document= 4
                        be((nsd + 1) * a + i) -= diff_coeff * lambda * fe.d2N(a, j, i) * tauM * res_M_RHS[j] * detJxW; ///ok!
                    }

                    /// force -- Term No. in the Document= 5
                    be((nsd + 1) * a + i) += force_coeff * fe.N(a) * force[i] * detJxW; ///ok!


                } // end of j loop

                /// fineContinuity RHS 1 -- Term No. in the Document= 6
                be((nsd + 1) * a + nsd) -= fe.dN(a, i) * tauM * res_M_RHS[i] * detJxW; ///ok!

            } // end of i loop

        } // end of a loop





    }


    VMSParams calc_tau(const TALYFEMLIB::FEMElm &fe, std::vector<double> u) const {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;
        const double Ci_f = 36;
        const int nsd = DIM;

        VMSParams params;

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
                u_Gu += u[i] * Ge(i, j) * u[j];
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


    std::vector<double> calc_linear_advection(double u[DIM], double u_prev[DIM], int option) {

        std::vector<double> u_star(DIM);

        switch (option) {
            case 1:
                for (int i = 0; i < DIM; i++) {
                    u_star[i] = u[i];
                }
                break;
            case 2:
                for (int i = 0; i < DIM; i++) {
                    u_star[i] = 2.0 * u[i] - u_prev[i];
                }
                break;
        }
        return u_star;
    }



    void ReRamping (double t) {
        if (t <= idata_->RampingTime) {
            idata_->Re = ((idata_->RampingRe - 10) / idata_->RampingTime) * t + 10;
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
