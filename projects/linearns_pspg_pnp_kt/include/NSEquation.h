//
// Created by makrand on 9/5/17.
//

#pragma once

#include <talyfem/talyfem.h>
#include <link.h>
#include "NSPNPInputData.h"
#include "NSPNPNodeData.hpp"
#include "NSPNPExpr.hpp"
#include <fstream>

#define CHECK_INF_NAN

class NSEquation : public TALYFEMLIB::CEquation<NSPNPNodeData> {
 protected:
	struct VMSParams {
			double tauM;
			double tauC;

			VMSParams()
					: tauM (1.0), tauC (1.0) {}
	};

 private:
//	GridField<NSPNPNodeData> *p_data_;
	NSPNPInputData *input_data_;

	// this flag is responsible for switching the contribution of dp_prev, especially in PP and VU
	double dp_prev_flag;
	// only applicable to the fine scale of the previous step in theta method
	double res_prev_flag;
	// whether fine scale contribution goes into PP and VU
	double stab_in_ppe_vu;
	// flag to switch on/off rotational form terms
	double rot_form_flag;
	// Rotational form: flag for terms on the boundary
	double rot_form_boundary_flag;
	// Theta method coefficient
	double theta;
	// BDF coefficient vector
	std::vector<double> bdf_c;
	// pressure extrapolation order
	int pExtrapOrder;
	// pressure extrapolation coefficients
	std::vector<double> p_extrap_c;


	/// Right now this only pertains to the linearized set of integrands.
	/// Coefficients to handle all non-dimensionalizations for NS
	/// ndcf := non dimensional coefficients
	/// Check constructor for assignment
	double ndcf_time_;
	double ndcf_conv_;
	double ndcf_diff_;
	double ndcf_pres_;
	double ndcf_pnp_coupling_;

	//// PNP Stuff
	DENDRITE_REAL pi = M_PI;
	std::vector<int> z;
	std::vector<DENDRITE_REAL> A_;
	std::vector<DENDRITE_REAL> B_;
	std::vector<DENDRITE_REAL> C_;
	std::vector<DENDRITE_REAL> D_;
	DENDRITE_REAL E_;
	DENDRITE_UINT nSpcs_ = 2;

 public:
#pragma mark Constructor
	/// Class constructor
    explicit NSEquation(NSPNPInputData *inputData)
			: TALYFEMLIB::CEquation<NSPNPNodeData>(false, TALYFEMLIB::kAssembleGaussPoints),
			        input_data_ (inputData), bdf_c (3, 0.0), p_extrap_c (4, 0.0) {
		setNonDimensionalCoefficients();

		/// PNP setup
		z.resize(2);
		z[0] = 1;
		z[1] = -1;
	}

    void Solve(double dt, double t) override {
        assert(false);
    }

    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }

#pragma mark Internal-functions

	void setNonDimensionalCoefficients() {
		ndcf_time_ = input_data_->ndcf_time_;
		ndcf_conv_ = input_data_->ndcf_conv_;
		ndcf_diff_ = input_data_->ndcf_diff_;
		ndcf_pres_ = input_data_->ndcf_pres_;
		ndcf_pnp_coupling_ = input_data_->ndcf_pnp_coupling_;

	}

	void setIntegrandsFlags() {
		if (input_data_->nsOrder == 1) {
			dp_prev_flag = 0.0;
		} else if (input_data_->nsOrder == 2) {
			dp_prev_flag = 1.0;
		}

		//dp_prev_flag = 1.0;
		res_prev_flag = 0.0;
		stab_in_ppe_vu = 1.0;

		rot_form_flag = 0.0;
		if (input_data_->ifUseRotationalForm) {
			rot_form_flag = 1.0;
		}

		rot_form_boundary_flag = 1.0;

		this->pExtrapOrder = input_data_->pExtrapOrder;

		if (input_data_->nsOrder == 1) {
			pExtrapOrder = 0;
			rot_form_flag = 0.0;
			theta = 1.0;
			bdf_c[0] = 1.0;
			bdf_c[1] = -1.0;
			bdf_c[2] = 0.0;
		} else if (input_data_->nsOrder == 2) {
			theta = 0.5;
			bdf_c[0] = 1.5;
			bdf_c[1] = -2.0;
			bdf_c[2] = 0.5;
		}

		if (pExtrapOrder == 0) {
			p_extrap_c[0] = 0.0;///< for t_current = t_n
			p_extrap_c[1] = 0.0;///< for t_n-1
			p_extrap_c[2] = 0.0;///< for t_n-2
			p_extrap_c[3] = 0.0;///< for t_n-3
		} else if (pExtrapOrder == 1) {
			p_extrap_c[0] = 0.0;///< for t_current = t_n
			p_extrap_c[1] = 1.0;///< for t_n-1
			p_extrap_c[2] = 0.0;///< for t_n-2
			p_extrap_c[3] = 0.0;///< for t_n-3
		} else if (pExtrapOrder == 2) {
			p_extrap_c[0] = 0.0;///< for t_current = t_n
			p_extrap_c[1] = 2.0;///< for t_n-1
			p_extrap_c[2] = -1.0;///< for t_n-2
			p_extrap_c[3] = 0.0;///< for t_n-3
		}
	}

	/** This function returns the linear extrapolated velocity
	 * calculated from velocities from previous steps
	 * @param vel_pre_list This vector should contain
	 * velocities at t_n, t_(n-1), t_(n-2) and so on in that order
	 * @return
	 */
    ZEROPTV calc_linear_advection(std::vector<ZEROPTV> &vel_pre_list) {
        int extrapolation_order = input_data_->advectionExtrapolationOrder;

        ZEROPTV u_star;
        for (int i = 0; i < DIM; i++) {
            if (extrapolation_order == 1) {
                u_star (i) = vel_pre_list[0] (i);
            } else if (extrapolation_order == 2) {
                u_star (i) = 2.0 * vel_pre_list[0] (i) - vel_pre_list[1] (i);
            }
        }
        return u_star;
    }

#pragma mark Integrands-interface-functions NS

	/**
	 * Returns the linear system (LHS part) of the integrands
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @param freeEnergy: freeEnergy object to calculate free energy on the fly
	 */
	void getIntegrandsNSAe(const FEMElm &fe, ZeroMatrix<double> &Ae, double dt, double t) {
		setIntegrandsFlags ();
		if (input_data_->timescheme == NSPNPInputData::THETA) {
			if (input_data_->ifLinearNS){
				Integrands_NS_Momentum_ThetaAe(fe, Ae, dt, t);
			} else {
				throw TALYException () << "Not implemented yet";
			}
		} else if (input_data_->timescheme == NSPNPInputData::BDF) {
			if (input_data_->ifLinearNS){
				Integrands_NS_Momentum_PSPG_BDF12Ae(fe, Ae, dt, t);
			} else {
				Integrands_NS_Momentum_BDF12Ae(fe, Ae, dt, t);
			}
		} else {
			throw TALYException () << "Time scheme not implemented. Should be either \"theta\" or \"bdf\"";
		}
	}


	/**
	 * Returns the linear system (RHS part) of the integrands
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @param freeEnergy: freeEnergy object to calculate free energy on the fly
	 */
	void getIntegrandsNSbe(const FEMElm &fe, ZEROARRAY<double> &be, double dt, double t) {
		setIntegrandsFlags ();
		if (input_data_->timescheme == NSPNPInputData::THETA) {
			if (input_data_->ifLinearNS) {
				Integrands_NS_Momentum_Thetabe(fe, be, dt, t);
			} else {
				throw TALYException () << "Not implemented yet";
			}
		} else if (input_data_->timescheme == NSPNPInputData::BDF) {
			if (input_data_->ifLinearNS) {
				Integrands_NS_Momentum_PSPG_BDF12be(fe, be, dt, t);
			} else {
				Integrands_NS_Momentum_BDF12be(fe, be, dt, t);
			}
		} else {
			throw TALYException () << "Time scheme not implemented. Should be either \"theta\" or \"bdf\"";
		}
	}


	/** Function to return the integrands for Pressure Poisson (Just Ae)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 */
	void getIntegrandsPPAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		setIntegrandsFlags ();
		NSIntegrandsPressurePoissonAe(fe, Ae, dt, t);
	}

	/** Function to return the integrands for Pressure Poisson (Just be)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 */
	void getIntegrandsPPbe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {
		setIntegrandsFlags ();
		NSIntegrandsPressurePoissonbe(fe, be, dt, t);
	}

	/** Function to return the integrands for boundary terms for Pressure Poisson (only contributions to RHS)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be matrix to be filled in the method
	 */
	void getIntegrands4sidePPbe(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be, const double dt, const double t) {
		setIntegrandsFlags ();
		NSIntegrands4sidePressurePoissonbe (fe, sideInd, be, dt, t);
	}

	/** Function to return the integrands for Velocity update (Just Ae)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 */
	 void getIntegrandsVelocityUpdateAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		setIntegrandsFlags ();
		NSIntegrandsVelocityUpdateAe (fe, Ae, dt, t);
	}

	/** Function to return the integrands for Velocity update (Just be)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be matrix to be filled in the method
	 */
	void getIntegrandsVelocityUpdatebe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {
		setIntegrandsFlags ();
		NSIntegrandsVelocityUpdatebe(fe, be, dt, t);
	}

#pragma mark tau-calculation

	/**
	 * This is a function to calculate tau at the current velocity, useful for non-linear BDF integrands
	 * @param fe finite element object
	 * @param vel velocity for which tau needs to be calculated
	 * @param Re Reynolds no.
	 * @param Ci_f Coefficient to control diffusive parameter in tau
	 * @param dt timestep
	 * @return
	 */
	virtual VMSParams calc_tau(const FEMElm &fe, double Re, double Ci_f, double dt) const {
		const double Coe_diff = 1.0 / Re;
		const int nsd = DIM;

		VMSParams params;

		ZEROPTV u;
		for (int i = 0; i < nsd; i++) {
			u (i) = this->p_data_->valueFEM (fe, i);
		}

		ZeroMatrix<double> ksiX;
		ksiX.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				ksiX (i, j) = fe.cof (j, i) / fe.jacc ();
			}
		}

		ZeroMatrix<double> Ge;
		Ge.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				Ge (i, j) = 0.0;
				for (int k = 0; k < nsd; k++)
					Ge (i, j) += ksiX (k, i) * ksiX (k, j);
			}
		}

		double u_Gu = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				u_Gu += u (i) * Ge (i, j) * u (j);
			}
		}

		double G_G_u = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				G_G_u += Ci_f * Coe_diff * Coe_diff * Ge (i, j) * Ge (i, j);
			}
		}

		params.tauM = 1.0 / sqrt (4.0 / (dt * dt) + u_Gu + G_G_u);

		ZEROARRAY<double> ge;
		ge.redim (nsd);
		for (int i = 0; i < nsd; i++) {
			ge (i) = 0.0;
			for (int j = 0; j < nsd; j++)
				ge (i) += ksiX (j, i);
		}

		double g_g = 0.0;
		for (int i = 0; i < nsd; i++) {
			g_g += ge (i) * ge (i);
		}

		params.tauC = 1.0 / (params.tauM * g_g);

		return params;
	}

	/**
	 * This is a function to calculate tau for given velocity, useful for linearized (semi-Implicit NS), as tau need to
	 * be calculated at different time levels in those integrands
	 * @param fe finite element object
	 * @param vel velocity for which tau needs to be calculated
	 * @param Re Reynolds no.
	 * @param Ci_f Coefficient to control diffusive parameter in tau
	 * @param dt timestep
	 * @return
	 */
	VMSParams calc_tau_vel(const FEMElm &fe, ZEROPTV &vel, double Re, double Ci_f, double dt) const {
		const double Coe_diff = 1.0 / Re;
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
				u_Gu += vel(i) * Ge(i, j) * vel(j);
			}
		}

		double G_G_u = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
			}
		}

		double denom_squared = (ndcf_time_ * ndcf_time_) * (4.0 / (dt * dt))
		                       + (ndcf_conv_ * ndcf_conv_) * u_Gu
		                       + (ndcf_diff_ * ndcf_diff_) * G_G_u;

		params.tauM = 1.0 / sqrt(denom_squared);

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

	virtual VMSParams calc_tau_surf(const FEMElm &fe, double Re, double Ci_f, double dt) const {

		const double Coe_diff = 1.0 / Re;
		const int nsd = DIM;

		VMSParams params;

		ZEROPTV u;
		for (int i = 0; i < nsd; i++) {
			u (i) = this->p_data_->valueFEM (fe, i);
		}

		double hy = 2.0 * fe.jacc ();
		double hx = 4.0 * fe.volume_jacc () / hy;

		ZeroMatrix<double> ksiX;
		ksiX.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				ksiX (i, j) = 0.0; //1.0 / fe.jacc(); //fe.cof(j, i) / fe.jacc();
				/*std::cout << "(i, j) = (" << i << ", " << j <<"), ksiX = "
                << ksiX(i,j) << ", J = " << fe.jacc()
                << ", rank = " << GetMPIRank()
                << "\n";*/
				//PrintInfo("(i, j) = (", i, ", ", j, "), ksiX = ", ksiX(i,j));
			}
		}

		ksiX (0, 0) = 2.0 / hx;
		ksiX (1, 1) = 2.0 / hy;

		ZeroMatrix<double> Ge;
		Ge.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				Ge (i, j) = 0.0;
				for (int k = 0; k < nsd; k++)
					Ge (i, j) += ksiX (k, i) * ksiX (k, j);
			}
		}

		double u_Gu = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				u_Gu += u (i) * Ge (i, j) * u (j);
			}
		}

		double G_G_u = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				G_G_u += Ci_f * Coe_diff * Coe_diff * Ge (i, j) * Ge (i, j);
			}
		}

		params.tauM = 1.0 / sqrt (4.0 / (dt * dt) + u_Gu + G_G_u);

		ZEROARRAY<double> ge;
		ge.redim (nsd);
		for (int i = 0; i < nsd; i++) {
			ge (i) = 0.0;
			for (int j = 0; j < nsd; j++)
				ge (i) += ksiX (j, i);
		}

		double g_g = 0.0;
		for (int i = 0; i < nsd; i++) {
			g_g += ge (i) * ge (i);
		}

		params.tauC = 1.0 / (params.tauM * g_g);

		return params;
	}

#pragma mark NS-integrands

	/** This is a variational multiscale solver for Navier-Stokes momentum
	 * equations with the BDF1 and BDF2 methods. Refer to Guermond et. al (2006),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void Integrands_NS_Momentum_BDF12Ae(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt

		double Re = input_data_->Re;           ///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1.0 / Re;
		double Ci_f = input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		/// field variables
		/** Get u and u_pre from the NodeData
		 * NOTE: ValueFEM is defined on NodeData object, therefore it follows indices defined in NodeData subclass.
		 * NOTE: Every iteration the solution is copied into NodeData, we use these solved fields from the NodeData
		 * object to calculate the be (NS residual)
		 * of current
		 * iteration
		 */
		ZEROPTV u, u_pre1, u_pre2;
		for (int i = 0; i < nsd; i++) {
			u (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X + i);
			u_pre1 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre2 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE2 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du;
		du.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
			}
		}

		/** Calculate the laplacian of velocity. This is required for course scale residual of Navier-Stokes for
		 * diffusion term
		 */
		ZEROPTV d2u;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += p_data_->value2DerivativeFEM (fe, dof, dir, dir);
			}
		}

		/// Pressure fields setup
		double p       = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE);
		double p_prev1 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_prev2 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE2);
		double pStar = p_extrap_c[0] * p + p_extrap_c[1] * p_prev1 + p_extrap_c[2] * p_prev2;

		/// Get gradient of pressure at previous timestep
		ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE     , i);
			dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
			dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE2, i);
			dpStar (i) = p_extrap_c[0] * dp (i) + p_extrap_c[1] * dpPrev1(i) + p_extrap_c[2] * dpPrev2(i);
		}


		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x () = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y () = 0;
		} else {
			gravity.y () = -1.0 / Fr;
		}
		gravity.z () = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of NS
		ZEROPTV forcingNS;
		if (input_data_->ifMMS) {
			forcingNS = calcForcingManufacturedSolutionNS (fe, t_);
		} else {
			forcingNS = {0.0, 0.0, 0.0};
		}

		VMSParams vms_params = calc_tau (fe, Re, Ci_f, dt);
		const double tauM = vms_params.tauM;

		// NS terms
		/** We define the convection of Navier Stokes here from
		 *  using inner product of gradient tensor and fields we acquired above.
		 */
		ZEROPTV convec;
		ZEROPTV diffusion;
		for (int i = 0; i < nsd; i++) {
			convec(i) = 0.0;
			diffusion(i) = 0.0;
			for (int j = 0; j < nsd; j++)
				convec(i) += du (i, j) * u (j);
				diffusion(i) += Coe_diff * d2u (i);
		}

		/** Construct the Navier Stokes strong form residual (used for calculating fine scale velocitys)
		 * Here diffusion is not present as its contribution to stabilizers for
		 * linear basis functions is zero
		 * Diffusion term is added for higher order basis functions
		 * Equation 61 of Bazilevs et al. (2007)
		 */
		ZEROPTV NS;
		for (int i = 0; i < nsd; i++) {
			NS (i) = (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_
			         + convec (i) + dpStar (i) - diffusion (i) - gravity (i) - forcingNS (i);
		}

		/** Calculate continuity residual for PSPG stabilizer
		 * This residual is a essential for calculating the PSPG stabilizer
		 * Equation 62 of Bazilevs et al. (2007)
		 */
		double cont = 0.0;
		for (int i = 0; i < nsd; i++) {
			cont += du (i, i);
		}

		/** Calculating the Elemental matrices requires loop over basis functions
		 * We loop over basis functions and calculate the Jacobian.
		 * The equation we are solving is \vect{J} \cdot \delta u = -E
		 * Here \vect{J} is the Jacobian matrix operating on \delta u and E is the residual.  We will be using stabilized
		 * forms of Jacobian matrix and residual
		 */
		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Jacobian operator--------------------///
			// NS terms

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52
			 * in Bazilevs et al. (2007) PS: these are not calculated using the SUPG class, and this is an approximation of
			 * (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u (i);
			}

			/// Actual fine scale contribution to cross terms from tauM(inverse
			/// estimate)
			double crossTermFineScalePart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermFineScalePart += fe.dN (a, i) * tauM * NS (i);
			}

			for (int b = 0; b < n_basis_functions; b++) {

				/// Convection term
				double conv = 0.0;
				for (int i = 0; i < nsd; i++) {
					conv += fe.dN (b, i) * u (i);
				}

				/// Adding terms to the Jacobian matrix.
				for (int i = 0; i < nsd; i++) {

					/// Transient terms and t
					/// he convection term and part of the stress
					/// tensor in the diagonals
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							(fe.N (a) * (bdf_c[0] * fe.N (b) / dt_ + conv)) * detJxW;

					for (int j = 0; j < nsd; j++) {
						/// This term calculates (w, (delta_u . grad{u_n}))
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								fe.N (a) * du (i, j) * fe.N (b) * detJxW;
						/** This term calculates (grad{w},grad{delta_u}), goes only in diagonals PS: In this case we are using
						 * the diffusion form not the stress tensor form of the momentun equation
						 */
						Ae ((ndof) * a + i, (ndof) * b + i) +=
								Coe_diff * fe.dN (a, j) * fe.dN (b, j) * detJxW;
					}
				}

				/** crossTerm 1: This term is written as, (w, grad{u u'}), which is
				 * weakened to (grad{w}, u\innerproduct u') which is then linearised.
				 * Here u' is fine scale velocity and u is resolved velocity. The fine scale velocity is approximated
				 * as -tau*Res{NS} (Equation 58 bazilevs et al. (2007))
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}
					/** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from  Bazilevs et al. (2007)) you get two
					 * terms, one goes in diagonal and other in the non-diagonal parts of the matrix
					 * PS: crossTermVelocityPart is (u. grad{w} => u*dN(a))
					 */
					/// Diagonal part
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							crossTermVelocityPart * tauM * (bdf_c[0] * fe.N (b) / dt_ + conv - diff_J) *
							detJxW;
					for (int j = 0; j < nsd; j++) {
						/// Off diagonal part
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								crossTermVelocityPart * tauM * du (i, j) * fe.N (b) * detJxW;
						/** this term is essentially, (grad{w}*tauM*Residual(NS), delta u) Equation 57 Bazilevs et al. (2007)
						 * PS: u' = (tauM*residual(NS))
						 */
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								tauM * fe.dN (a, j) * NS (i) * fe.N (b) * detJxW;
					}
				}

				/** Crossterm2:This term can be mathematically written as,
				 * (w, u'.grad{u}). In this term we do not further weaken it as done in
				 * cross term 1
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}

					/** This term represents the contribution of, (w, tau \delta{u} . grad{u}) term to the cross term, here
					 * \delta{} is the linearisation operator (basically delta{u'} represents the change in u').  This is the
					 * first term which arise due to linearisation of the nonliner term of the Navier-Stokes residual when u' is
					 * substitute for in (w,\delta{u'}.grad{u}). PS: u' = - tauM*res{NS}
					 */
					for (int j = 0; j < nsd; j++) {
						/// k is the dummy index
						for (int k = 0; k < nsd; k++) {
							Ae ((ndof) * a + i, (ndof) * b + j) +=
									-du (i, k) * fe.N (a) * tauM * fe.N (b) * du (k, j) * detJxW;
						}
					}

					for (int j = 0; j < nsd; j++) {
						/// This term represents the contribution of, (w, \delta{u} . grad{u\delta{u'}})
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								-du (i, j) * fe.N (a) * tauM * (bdf_c[0] * fe.N (b) / dt_ + conv - diff_J) *
								detJxW;
						/// This term is the contribution of (w, u'.grad{delta{u}}) to the Jacobian operator.
						Ae ((ndof) * a + i, (ndof) * b + i) +=
								-fe.N (a) * tauM * NS (j) * fe.dN (b, j) * detJxW;
					}
				}

				/// The Reynolds Stress term: (w, u'.grad{u'}), we subsitute u' as -tau*Res(NS) and expand.
				for (int i = 0; i < nsd; i++) {
					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}

					/** Three terms arising from (w,\delta(u').grad{u'})
					 * u' has to be expanded to -tau*res{NS}: when the linearisation acts on the res{NS} it gives three terms.
					 * First term which goes in the diagonals of the matrix is
					 * (-grad{w}* tauM*res(NS), tauM*(d_t{\delta{u}} + conv -diff_J))
					 * PS: grad{w}*tau*res(NS) is corssTermFineScalePart
					 */
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							-crossTermFineScalePart * tauM * (bdf_c[0] * fe.N (b) / dt_ + conv - diff_J) *
							detJxW;

					/** Second term from (w,\delta(u').grad{u'}) which goes to the off
					 * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
					 */
					for (int j = 0; j < nsd; j++) {
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								-crossTermFineScalePart * tauM * du (i, j) * fe.N (b) * detJxW;
					}

					for (int j = 0; j < nsd; j++)
						for (int k = 0; k < nsd; k++)
							Ae ((ndof) * a + i, (ndof) * b + j) +=
									-tauM * NS (i) * fe.dN (a, k) * tauM * fe.N (b) * du (k, j) *
									detJxW;

					/// Just as above terms which arise when(w,(u').grad{\delta{u'}}) is expanded is given below.
					for (int j = 0; j < nsd; j++) {
						Ae ((ndof) * a + i, (ndof) * b + j) +=
								-tauM * NS (i) * fe.dN (a, j) * tauM *
								(bdf_c[0] * fe.N (b) / dt_ + conv - diff_J) * detJxW;
					}
				}
			}
			///--------------------Done with Jacobian operator Matrix----------------///
		}
	}

	/** This is a variational multiscale solver for Navier-Stokes momentum
	 * equations with the BDF1 and BDF2 methods. Refer to Guermond et. al (2006),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return void
	 */
	void Integrands_NS_Momentum_BDF12be(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt

		double Re = input_data_->Re;           ///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1.0 / Re;
		double Ci_f = input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		/// field variables
		/** Get u and u_pre from the NodeData
		 * NOTE: ValueFEM is defined on NodeData object, therefore it follows indices defined in NodeData subclass.
		 * NOTE: Every iteration the solution is copied into NodeData, we use these solved fields from the NodeData
		 * object to calculate the be (NS residual)
		 * of current
		 * iteration
		 */
		ZEROPTV u, u_pre1, u_pre2;
		for (int i = 0; i < nsd; i++) {
			u (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X + i);
			u_pre1 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre2 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE2 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du;
		du.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
			}
		}

		/** Calculate the laplacian of velocity. This is required for course scale residual of Navier-Stokes for
		 * diffusion term
		 */
		ZEROPTV d2u;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += p_data_->value2DerivativeFEM (fe, dof, dir, dir);
			}
		}

		/// Pressure fields setup
		double p       = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE);
		double p_prev1 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_prev2 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE2);
		double pStar = p_extrap_c[0] * p + p_extrap_c[1] * p_prev1 + p_extrap_c[2] * p_prev2;

		/// Get gradient of pressure at previous timestep
		ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE     , i);
			dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
			dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE2, i);
			dpStar (i) = p_extrap_c[0] * dp (i) + p_extrap_c[1] * dpPrev1(i) + p_extrap_c[2] * dpPrev2(i);
		}

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x () = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y () = 0;
		} else {
			gravity.y () = -1.0 / Fr;
		}
		gravity.z () = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of NS
		ZEROPTV forcingNS;
		if (input_data_->ifMMS) {
			forcingNS = calcForcingManufacturedSolutionNS (fe, t_);
		} else {
			forcingNS = {0.0, 0.0, 0.0};
		}

		VMSParams vms_params = calc_tau (fe, Re, Ci_f, dt);
		const double tauM = vms_params.tauM;

		// NS terms
		/** We define the convection of Navier Stokes here from
		 *  using inner product of gradient tensor and fields we acquired above.
		 */
		ZEROPTV convec;
		ZEROPTV diffusion;
		for (int i = 0; i < nsd; i++) {
			convec (i) = 0.0;
			diffusion (i) = 0.0;
			for (int j = 0; j < nsd; j++)
				convec (i) += du (i, j) * u (j);
				diffusion (i) += Coe_diff * d2u (i);
		}

		/** Construct the Navier Stokes strong form residual (used for calculating fine scale velocitys)
		 * Here diffusion is not present as its contribution to stabilizers for
		 * linear basis functions is zero
		 * Diffusion term is added for higher order basis functions
		 * Equation 61 of Bazilevs et al. (2007)
		 */
		ZEROPTV NS;
		for (int i = 0; i < nsd; i++) {
			NS (i) = (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_
			         + convec (i) + dpStar (i) - diffusion (i) - gravity (i) - forcingNS (i);
		}

		for (int a = 0; a < n_basis_functions; a++) {
			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52
			 * in Bazilevs et al. (2007) PS: these are not calculated using the SUPG class, and this is an approximation of
			 * (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u (i);
			}

			/** Constructing the RHS vector of the residual of the equation. This involves just defining the terms involved
			 * in the equation without any linearisation, as this is calculated from the existing data on the NodeData
			 * structure
			 */
			/// All the terms involved in the equation
			ZEROPTV cross1, cross2, reystress, pressurep, normal_NS, M;

			/// Diffusion term
			ZEROPTV diff;
			for (int i = 0; i < nsd; i++) {
				diff (i) = 0.0;
				for (int j = 0; j < nsd; j++) {
					diff (i) += Coe_diff * fe.dN (a, j) * (du (i, j));
				}
			}

			/// Terms of Normal Navier Stokes without any additional VMS terms
			for (int i = 0; i < nsd; i++) {
				normal_NS (i) = fe.N (a) * (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_ * detJxW
				                + fe.N (a) * (convec (i)) * detJxW
				                + (-fe.dN (a, i) * pStar) * detJxW
				                + (diff (i)) * detJxW
				                - (fe.N (a) * gravity (i) * detJxW)
				                - (fe.N (a) * forcingNS (i) * detJxW);

				/// first cross term
				cross1 (i) = (tauM * crossTermVelocityPart * NS (i)) * detJxW;

				cross2 (i) = 0.0;
				reystress (i) = 0.0;

				for (int j = 0; j < nsd; j++) {
					/// Second cross term
					cross2 (i) += fe.N (a) * (du (i, j) * tauM * NS (j)) * detJxW;

					/// The Reynolds stress term
					reystress (i) += fe.dN (a, j) * (tauM * NS (i) * tauM * NS (j)) * detJxW;
				}

				/// Add all the terms to the vector
				M (i) = normal_NS (i) + cross1 (i) - cross2 (i) - reystress (i);
			}

			/// Adding all the calculated terms to the vector
			for (int i = 0; i < nsd; i++) {
				be ((ndof) * a + i) += M (i);
			}
		}
		///-------------------Done with populating be vector-----------------------///
	}



	/** This is a variational multiscale solver for Navier-Stokes momentum
	 * equations with the Theta methods (i.e., BE and CN). Refer to John Volker's book (2016),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return void
	 */
	void Integrands_NS_Momentum_ThetaAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM;
		double Re = input_data_->Re;
		double Coe_diff = 1. / Re;

		double Ci_f = this->input_data_->Ci_f;

		const int nsd = DIM;
		const int n_basis_functions = fe.nbf ();
		const double detJxW = fe.detJxW ();
		/// Calculate the tau_M based on the projections calculated, Equation 63 in
		/// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
		/// in the paper it is set to 1.
		VMSParams vms_params = calc_tau (fe, Re, Ci_f, dt);
		const double tauM = vms_params.tauM;
		/// Calculate continuity residual based on
		const double tauC = vms_params.tauC;

		//unsigned int pre_vars_offset = NSPNPNodeData::NUM_VARS + NSPNPNodeData::MANUFAC_SOL_VARS;
		//int solenoidalOffset = pre_vars_offset + NSPNPNodeData::NS_DOF;

		/// field variables
		/** Get u and u_pre from the NodeData
		 * NOTE: ValueFEM is defined on NodeData object, therefore it follows indices defined in NodeData subclass.
		 * NOTE: Every iteration the solution is copied into NodeData, we use these solved fields from the NodeData
		 * object to calculate the be (NS residual) of current iteration
		 */
		ZEROPTV u, u_pre;
		for (int i = 0; i < nsd; i++) {
			u(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X + i);
			u_pre(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du, du_pre;
		du.redim (nsd, nsd);
		du_pre.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
				du_pre (i, j) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
			}
		}

		/** Calculate the laplacian of velocity
		 * This is required for course scale residual of Navier-Stokes for diffusion term
		 */
		ZEROPTV d2u, d2u_pre;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += this->p_data_->value2DerivativeFEM (fe, dof, dir, dir);
				d2u_pre (dof) += this->p_data_->value2DerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + dof, dir, dir);
			}
		}

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test
		ZEROPTV forcing, forcing_prev;
		if (input_data_->ifMMS) {
			forcing = calcForcingManufacturedSolutionNS (fe, t + dt);
			forcing_prev = calcForcingManufacturedSolutionNS (fe, t);
		}

		double p = p_data_->valueFEM (fe, NSPNPNodeData::PRESSURE);
		double p_pre = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		// double theta = p_data_->valueFEM(fe, nsd + 1);

		/// Get gradient of pressure
		ZEROPTV dp, dp_pre;
		for (int i = 0; i < nsd; i++) {
			dp (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::PRESSURE, i);
			dp_pre (i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
		}

		// NS terms
		ZEROPTV convec, convec_pre;
		ZEROPTV diffusion, diffusion_pre;
		for (int i = 0; i < nsd; i++) {
			convec (i) = 0.0;
			diffusion (i) = 0.0;

			convec_pre (i) = 0.0;
			diffusion_pre (i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec (i) += du (i, j) * u(j);
				convec_pre (i) += du_pre (i, j) * u_pre(j);
			}
			diffusion (i) += Coe_diff * d2u (i);
			diffusion_pre (i) += Coe_diff * d2u_pre (i);
		}

		/** Construct the Navier Stokes equation without diffusion. Here diffusion is not present as its contribution to
		 * stabilizers for linear basis functions is zero. Diffusion term is added for higher order basis functions
		 * Equation 61 of Bazilevs et al. (2007)
		 */
		ZEROPTV NS, NS_pre;
		for (int i = 0; i < nsd; i++) {
			NS (i) = (u(i) - u_pre(i)) / dt + convec (i) + dp_prev_flag * dp_pre (i) - diffusion (i) - forcing (i);
			NS_pre (i) = (u (i) - u_pre (i)) / dt + convec_pre (i) + dp_pre (i) - diffusion_pre (i) - forcing_prev (i);
		}

		/** Calculate continuity residual for PSPG stabilizer
         * This residual is a essential for calculating the PSPG stabilizer
         * Equation 62 of Bazilevs et al. (2007)
         */
		double cont = 0.0, cont_pre = 0.0;
		for (int i = 0; i < nsd; i++) {
			cont += du (i, i);
			cont_pre += du_pre (i, i);
		}

		/** Calculating the Elemental matrices requires loop over basis functions. We loop over basis functions and
		 * calculate the Jacobian and the Residual. The equation we are solving is \vect{J} \cdot \delta u = -E
		 * Here \vect{J} is the Jacobian matrix operating on \delta u and E is the residual.  We will be using stabilized
		 * forms of Jacobian matrix and residual
		 */
		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Jacobian operator--------------------///

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an
			 * approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0, crossTermVelocityPart_pre = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u (i);
				crossTermVelocityPart_pre += fe.dN (a, i) * u_pre (i);
			}

			/// Actual fine scale contribution to cross terms from tauM(inverse
			/// estimate)
			double crossTermFineScalePart = 0.0, crossTermFineScalePart_pre = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermFineScalePart += fe.dN (a, i) * tauM * NS (i);
				crossTermFineScalePart_pre += fe.dN (a, i) * tauM * NS_pre (i);
			}

			for (int b = 0; b < n_basis_functions; b++) {

				/// Convection term
				double conv = 0.0;
				for (int i = 0; i < nsd; i++) {
					conv += fe.dN (b, i) * u (i);
				}

				/// Adding terms to the Jacobian matrix.
				for (int i = 0; i < nsd; i++) {

					/// Transient terms and the convection term and part of the stress
					/// tensor in the diagonals
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							(fe.N (a) * (fe.N (b) / dt + theta * conv)) * detJxW;

					for (int j = 0; j < nsd; j++) {
						/// This term calculates (w, (delta_u . grad{u_n}))
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (fe.N (a) * du (i, j) * fe.N (b) * detJxW);
						/** This term calculates (grad{w},grad{delta_u}), goes only in diagonals PS: In this case we are using
						 * the diffusion form not the stress tensor form of the momentun equation
						 */
						/// 1/Re (\partial_j w_i, \partial v_i)
						/// linearised : Coe_diff (\partial_j w_i, \partial_j \delta(v_i))
						/// COeff_diff * fe.dN(a,j) * fe.dN(b,j) *detJxW summed over j for each i in diagonals
						Ae ((ndof) * a + i, (ndof) * b + i) += theta *
						                                       (Coe_diff * fe.dN (a, j) * fe.dN (b, j) * detJxW);
					}
				}

				/** crossTerm 1: This term is written as, (w, grad{u u'}), which is weakened to (grad{w}, u\innerproduct u')
				 * which is then linearised. Here u' is fine scale velocity and u is resolved velocity. The fine scale
				 * velocity is approximated as -tau*Res{NS} (Equation 58 bazilevs et al. (2007))
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}
					/** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from Bazilevs et al. (2007)) you get two
					 * terms, one goes in diagonal and other in the non-diagonal parts of the matrix
					 * PS: crossTermVelocityPart is (u. grad{w} => u*dN(a))
					 */
					/// Diagonal part
					Ae ((ndof) * a + i, (ndof) * b + i) += theta *
					                                       (crossTermVelocityPart * tauM *
					                                        (fe.N (b) / dt + conv - diff_J) *
					                                        detJxW) + res_prev_flag * (1.0 - theta) *
					                                                  (crossTermVelocityPart_pre * tauM *
					                                                   (fe.N (b) / dt) * detJxW);
					for (int j = 0; j < nsd; j++) {
						/// Off diagonal part
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (crossTermVelocityPart * tauM * du (i, j) * fe.N (b) *
						                                        detJxW);
						/** this term is essentially, (grad{w}*tauM*Residual(NS), delta u) Equation 57 Bazilevs et al. (2007)
						 * PS: u' = (tauM*residual(NS))
						 */
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (tauM * fe.dN (a, j) * NS (i) * fe.N (b) * detJxW);
					}
				}

				/** Crossterm2:This term can be mathematically written as, (w, u'.grad{u}). In this term we do not further
				 * weaken it as done in cross term 1
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}

					/** This term represents the contribution of, (w, tau \delta{u} . grad{u}) term to the cross term, here
					 * \delta{} is the linearisation operator (basically delta{u'} represents the change in u').  This is the
					 * first term which arise due to linearisation of the nonliner term of the Navier-Stokes residual when u'
					 * is substitute for in (w,\delta{u'}.grad{u}). PS: u' = - tauM*res{NS}
					 */
					for (int j = 0; j < nsd; j++) {
						/// k is the dummy index
						for (int k = 0; k < nsd; k++) {
							Ae ((ndof) * a + i, (ndof) * b + j) += theta *
							                                       (-du (i, k) * fe.N (a) * tauM * fe.N (b) *
							                                        du (k, j) * detJxW);
						}
					}

					for (int j = 0; j < nsd; j++) {
						/// This term represents the contribution of, (w, \delta{u} . grad{u\delta{u'}})
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (-du (i, j) * fe.N (a) * tauM *
						                                        (fe.N (b) / dt + conv - diff_J) *
						                                        detJxW) + res_prev_flag * (1.0 - theta) *
						                                                  (-du_pre (i, j) * fe.N (a) * tauM *
						                                                   (fe.N (b) / dt) * detJxW);
						/// This term is the contribution of (w, u'.grad{delta{u}}) to the Jacobian operator.
						Ae ((ndof) * a + i, (ndof) * b + i) += theta *
						                                       (-fe.N (a) * tauM * NS (j) * fe.dN (b, j) * detJxW);

						/// Pressure part of the term which arise from (w,\delta{u'}
						/// .grad{u}), always the last diagonal term
					}
				}

				/// The Reynolds Stress term: (w, u'.grad{u'}), we subsitute u' as -tau*Res(NS) and expand.
				for (int i = 0; i < nsd; i++) {
					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}

					/** Three terms arising from (w,\delta(u').grad{u'}) u' has to be expanded to -tau*res{NS}: when the
					 * linearisation acts on the res{NS} it gives three terms. First term which goes in the diagonals of the
					 * matrix is (-grad{w}* tauM*res(NS), tauM*(d_t{\delta{u}} + conv -diff_J)) PS: grad{w}*tau*res(NS) is
					 * corssTermFineScalePart
					 */
					Ae ((ndof) * a + i, (ndof) * b + i) += theta *
					                                       (-crossTermFineScalePart * tauM *
					                                        (fe.N (b) / dt + conv - diff_J) *
					                                        detJxW) + res_prev_flag * (1.0 - theta) *
					                                                  (-crossTermFineScalePart_pre * tauM *
					                                                   (fe.N (b) / dt) * detJxW);

					/** Second term from (w,\delta(u').grad{u'}) which goes to the off diagonals is: (-grad{w}* tauM*res(NS),
					 * tauM*(\delta{u}. grad{u}))
					 */
					for (int j = 0; j < nsd; j++) {
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (-crossTermFineScalePart * tauM * du (i, j) *
						                                        fe.N (b) * detJxW);
					}

					/** Third term from (w,\delta(u').grad{u'}) which is the contribution of pressure and goes in the last
					 * diagonal term
					 */
					for (int j = 0; j < nsd; j++) {
						for (int k = 0; k < nsd; k++) {
							Ae ((ndof) * a + i, (ndof) * b + j) += theta *
							                                       (-tauM * NS (i) * fe.dN (a, k) * tauM * fe.N (b) *
							                                        du (k, j) *
							                                        detJxW);
						}
					}

					/// Just as above terms which arise when(w,(u').grad{\delta{u'}}) is expanded is given below.
					for (int j = 0; j < nsd; j++) {
						Ae ((ndof) * a + i, (ndof) * b + j) += theta *
						                                       (-tauM * NS (i) * fe.dN (a, j) * tauM *
						                                        (fe.N (b) / dt + conv - diff_J) * detJxW) +
						                                       res_prev_flag * (1.0 - theta) *
						                                       (-tauM * NS_pre (i) * fe.dN (a, j) * tauM *
						                                        (fe.N (b) / dt) * detJxW);
					}
				}
			}
		}
		/// Done with assembling elemental Jacobian
	}

	/** This is a variational multiscale solver for Navier-Stokes momentum equations with the Theta methods (i.e., BE
	 * and CN). Refer to John Volker's book (2016),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return void
	 */
	void Integrands_NS_Momentum_Thetabe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM;
		double Re = input_data_->Re;
		double Coe_diff = 1. / Re;

		double Ci_f = this->input_data_->Ci_f;

		const int nsd = DIM;
		const int n_basis_functions = fe.nbf ();
		const double detJxW = fe.detJxW ();
		/// Calculate the tau_M based on the projections calculated, Equation 63 in
		/// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
		/// in the paper it is set to 1.
		VMSParams vms_params = calc_tau (fe, Re, Ci_f, dt);
		const double tauM = vms_params.tauM;
		/// Calculate continuity residual based on
		const double tauC = vms_params.tauC;

		//unsigned int pre_vars_offset = NSPNPNodeData::NUM_VARS + NSPNPNodeData::MANUFAC_SOL_VARS;
		//int solenoidalOffset = pre_vars_offset + NSPNPNodeData::NS_DOF;

		/// field variables
		/** Get u and u_pre from the NodeData
		 * NOTE: ValueFEM is defined on NodeData object, therefore it follows indices defined in NodeData subclass.
		 * NOTE: Every iteration the solution is copied into NodeData, we use these solved fields from the NodeData
		 * object to calculate the be (NS residual) of current iteration
		 */
		ZEROPTV u, u_pre;
		for (int i = 0; i < nsd; i++) {
			u(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X + i);
			u_pre(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du, du_pre;
		du.redim (nsd, nsd);
		du_pre.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
				du_pre (i, j) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
			}
		}

		/** Calculate the laplacian of velocity
		 * This is required for course scale residual of Navier-Stokes for diffusion term
		 */
		ZEROPTV d2u, d2u_pre;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += this->p_data_->value2DerivativeFEM (fe, dof, dir, dir);
				d2u_pre (dof) += this->p_data_->value2DerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + dof, dir, dir);
			}
		}

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test
		ZEROPTV forcing, forcing_prev;
		if (input_data_->ifMMS) {
			forcing = calcForcingManufacturedSolutionNS (fe, t + dt);
			forcing_prev = calcForcingManufacturedSolutionNS (fe, t);
		}

		double p = p_data_->valueFEM (fe, NSPNPNodeData::PRESSURE);
		double p_pre = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		// double theta = p_data_->valueFEM(fe, nsd + 1);

		/// Get gradient of pressure
		ZEROPTV dp, dp_pre;
		for (int i = 0; i < nsd; i++) {
			dp (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::PRESSURE, i);
			dp_pre (i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
		}

		// NS terms
		ZEROPTV convec, convec_pre;
		ZEROPTV diffusion, diffusion_pre;
		for (int i = 0; i < nsd; i++) {
			convec (i) = 0.0;
			diffusion (i) = 0.0;

			convec_pre (i) = 0.0;
			diffusion_pre (i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec (i) += du (i, j) * u(j);
				convec_pre (i) += du_pre (i, j) * u_pre(j);
			}
			diffusion (i) += Coe_diff * d2u (i);
			diffusion_pre (i) += Coe_diff * d2u_pre (i);
		}

		/** Construct the Navier Stokes equation without diffusion. Here diffusion is not present as its contribution to
		 * stabilizers for linear basis functions is zero. Diffusion term is added for higher order basis functions
		 * Equation 61 of Bazilevs et al. (2007)
		 */
		ZEROPTV NS, NS_pre;
		for (int i = 0; i < nsd; i++) {
			NS (i) = (u(i) - u_pre(i)) / dt + convec (i) + dp_prev_flag * dp_pre (i) - diffusion (i) - forcing (i);
			NS_pre (i) = (u (i) - u_pre (i)) / dt + convec_pre (i) + dp_pre (i) - diffusion_pre (i) - forcing_prev (i);
		}

		/** Calculate continuity residual for PSPG stabilizer
         * This residual is a essential for calculating the PSPG stabilizer
         * Equation 62 of Bazilevs et al. (2007)
         */
		double cont = 0.0, cont_pre = 0.0;
		for (int i = 0; i < nsd; i++) {
			cont += du (i, i);
			cont_pre += du_pre (i, i);
		}

		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Jacobian operator--------------------///

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an
			 * approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0, crossTermVelocityPart_pre = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u (i);
				crossTermVelocityPart_pre += fe.dN (a, i) * u_pre (i);
			}

			/// Actual fine scale contribution to cross terms from tauM(inverse
			/// estimate)
			double crossTermFineScalePart = 0.0, crossTermFineScalePart_pre = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermFineScalePart += fe.dN (a, i) * tauM * NS (i);
				crossTermFineScalePart_pre += fe.dN (a, i) * tauM * NS_pre (i);
			}

			/** Constructing the RHS vector of the residual of the equation. This involves just defining the terms involved
			 * in the equation without any linearisation, as this is calculated from the existing data on the NodeData
			 * structure
			 */
			/// All the terms involved in the equation
			ZEROPTV time_der;
			ZEROPTV cross1, cross2, reystress, pressurep, diff, press, normal_NS, NScontrib;
			ZEROPTV cross1_pre, cross2_pre, reystress_pre, pressurep_pre, diff_pre, press_pre, normal_NS_pre, NScontrib_pre;

			/// Diffusion term
			for (int i = 0; i < nsd; i++) {
				diff (i) = 0.0;
				diff_pre (i) = 0.0;
				for (int j = 0; j < nsd; j++) {
					//diff(i) += Coe_diff * fe.dN(a, j) * (du(i, j) + du(j, i));
					diff (i) += Coe_diff * fe.dN (a, j) * du (i, j);
					//diff_pre(i) += Coe_diff * fe.dN(a, j) * (du_pre(i, j) + du_pre(j, i));
					diff_pre (i) += Coe_diff * fe.dN (a, j) * du_pre (i, j);
				}
			}

			/// Terms of Normal Navier Stokes without any additional VMS terms
			for (int i = 0; i < nsd; i++) {
				time_der (i) = fe.N (a) * (u (i) - u_pre (i)) / dt * detJxW;
				//press(i) = (-fe.dN(a, i) * p) * detJxW;
				press_pre (i) = (-fe.dN (a, i) * p_pre) * detJxW;
				normal_NS (i) = fe.N (a) * (convec (i)) * detJxW
				                + (diff (i)) * detJxW
				                - fe.N (a) * forcing (i) * detJxW;
				normal_NS_pre (i) =
						fe.N (a) * convec_pre (i) * detJxW
						+ diff_pre (i) * detJxW
						- fe.N (a) * forcing_prev (i) * detJxW;

				/// first cross term
				cross1 (i) = (tauM * crossTermVelocityPart * NS (i)) * detJxW;
				cross1_pre (i) = (tauM * crossTermVelocityPart_pre * NS_pre (i)) * detJxW;

				cross2 (i) = 0.0, cross2_pre (i) = 0.0;
				reystress (i) = 0.0, reystress_pre (i) = 0.0;

				for (int j = 0; j < nsd; j++) {
					/// Second cross term
					cross2 (i) += fe.N (a) * (du (i, j) * tauM * NS (j)) * detJxW;
					cross2_pre (i) += fe.N (a) * (du_pre (i, j) * tauM * NS_pre (j)) * detJxW;

					/// The Reynolds stress term
					reystress (i) += fe.dN (a, j) * (tauM * NS (i) * tauM * NS (j)) * detJxW;
					reystress_pre (i) += fe.dN (a, j) * (tauM * NS_pre (i) * tauM * NS_pre (j)) * detJxW;
				}

				/// Add all the terms to the vector
				NScontrib (i) = time_der (i) +
				                theta * (normal_NS (i) + cross1 (i) - cross2 (i) - reystress (i));
				NScontrib_pre (i) = (1.0 - theta) *
				                    (normal_NS_pre (i) + res_prev_flag * cross1_pre (i) - res_prev_flag * cross2_pre (i) -
				                     res_prev_flag * reystress_pre (i));
			}

			/// Adding all the calculated terms to the vector
			for (int i = 0; i < nsd; i++) {
				be ((ndof) * a + i) += NScontrib (i) + NScontrib_pre (i) + press_pre (i) * dp_prev_flag;
			}
		}
	}

	#pragma mark NS-integrands-linearized
    /// linearized VMS NS
	/** This is a variational multiscale solver for Navier-Stokes momentum (semi implicit linear solver) equations with
	 * the BDF1 and BDF2 methods. Refer to Guermond et. al (2006),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void Integrands_NS_Momentum_BDF12_semiImplLinearAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt,
	                                                   const double t) {
		const int ndof = DIM;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt
		int degOfFreedom = DIM; ///equal to no of vel components degrees of freedom

		double Re = input_data_->Re;///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1 / Re;

		//todo this is temperorary till a proper way to set these coeff is added
		ndcf_time_ = 1.0;
		ndcf_conv_ = 1.0;
		ndcf_diff_ = 1.0;
		ndcf_pres_ = 1.0;


		double Ci_f = this->input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf(); // # of basis functions
		double detJxW = fe.detJxW();

		/// FIELD VARIABLES
		ZEROPTV u_pre_1, u_pre_2, u_pre_3;
		for (int i = 0; i < nsd; i++) {
			u_pre_1(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre_2(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
			u_pre_3(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE3 + i);
		}


		/// LINEARIZED CONVECTION

		/// u_star is the approximate linear advection at this time step
		std::vector<ZEROPTV> vel_pre_list = {u_pre_1, u_pre_2};
		ZEROPTV u_star = calc_linear_advection(vel_pre_list);

		VMSParams vms_params_current = calc_tau_vel(fe, u_star, Re, Ci_f, dt);
		double tauM = vms_params_current.tauM;

		/// Calculating the Elemental matrices requires loop over basis functions
		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the LHS Elemental matrix--------------------///

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN(a, i) * u_star(i);
			}

			for (int b = 0; b < n_basis_functions; b++) {

				/// Convection term
				double conv = 0.0;
				for (int i = 0; i < nsd; i++) {
					conv += fe.dN(b, i) * u_star(i);
				}

				/// Adding terms to the Elemental matrix.
				for (int i = 0; i < nsd; i++) {
					// Time derivative
					Ae((ndof) * a + i, (ndof) * b + i) +=
							(fe.N(a) * (ndcf_time_ * bdf_c[0] * fe.N(b) / dt_)) * detJxW;
					//validateMatrix(Ae);

					// Coarse scale convection term
					Ae((ndof) * a + i, (ndof) * b + i) += (fe.N(a) * (ndcf_conv_ * conv)) * detJxW;
					//validateMatrix(Ae);

					// Coarse scale (only) diffusion term
					for (int j = 0; j < nsd; j++) {
						Ae((ndof) * a + i, (ndof) * b + i) +=
								Coe_diff * fe.dN(a, j) * (ndcf_diff_ * fe.dN(b, j)) * detJxW;
					}
				}


				/** crossTerm 1: This term is written as, (w, grad{u u'}), which is
				 * weakened to (grad{w}, u\innerproduct u') which is then linearised.
				 * Here u' is fine scale velocity
				 * and u is resolved velocity. The fine scale velocity is approximated
				 * as -tau*Res{NS} (Equation 58 bazilevs et al. (2007)),
				 *
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N(b, j, j);
					}
					/** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from
					 * Bazilevs et al. (2007)) you get two terms, one goes
					 * in diagonal and other in the non-diagonal parts of the matrix
					 * PS: crossTermVelocityPart is (u. grad{w} => u*dN(a))
					 */
					/// Diagonal part
					Ae ((ndof) * a + i, (ndof) * b + i) += crossTermVelocityPart * tauM *
							((ndcf_time_ * bdf_c[0] * fe.N (b) / dt_) + (ndcf_conv_ * conv) - (ndcf_diff_ * diff_J)) * detJxW;
					//validateMatrix(Ae);
				}
			}
			///--------------------Done with Elemental operator Matrix----------------///
		}
	}

	/** This is a variational multiscale solver for Navier-Stokes momentum (semi implicit linear solver) equations with
	 * the BDF1 and BDF2 methods. Refer to Guermond et. al (2006) (just be),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector (LHS) to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void Integrands_NS_Momentum_BDF12_semiImplLinearbe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt,
	                                                   const double t) {
		const int ndof = DIM;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt
		int degOfFreedom = DIM; ///equal to no of vel components degrees of freedom

		double Re = input_data_->Re;///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1 / Re;

		//todo this is temperorary till a proper way to set these coeff is added
		ndcf_time_ = 1.0;
		ndcf_conv_ = 1.0;
		ndcf_diff_ = 1.0;
		ndcf_pres_ = 1.0;


		double Ci_f = this->input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf(); // # of basis functions
		double detJxW = fe.detJxW();

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x() = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y() = 0;
		} else {
			gravity.y() = -1.0 / Fr;
		}
		gravity.z() = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcing_current;
		if (input_data_->ifMMS) {
			forcing_current = calcForcingManufacturedSolutionNS(fe, t_);
		} else {
			forcing_current = {0.0, 0.0, 0.0};
		}

		/// FIELD VARIABLES
		ZEROPTV u_pre_1, u_pre_2, u_pre_3;
		for (int i = 0; i < nsd; i++) {
			u_pre_1(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre_2(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
			u_pre_3(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE3 + i);
		}

		/// Pressure values at previous steps
		double p_pre_1_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_pre_2_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE2);
		double p_pre_3_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE3);
		/// p_extrap_c[0] is not there because there is no contribution from current timestep. pressure is completely
		/// explicit
		double pStar = p_extrap_c[1] * p_pre_1_c + p_extrap_c[2] * p_pre_2_c + p_extrap_c[3] * p_pre_3_c;

		/// Derivatives of pressure at previous steps
		ZEROPTV dp_pre_1_c, dp_pre_2_c, dp_pre_3_c, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp_pre_1_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE1), i);
			dp_pre_2_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE2), i);
			dp_pre_3_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE3), i);
			dpStar(i) = p_extrap_c[1] * dp_pre_1_c(i) + p_extrap_c[2] * dp_pre_2_c(i) + p_extrap_c[3] * dp_pre_3_c(i);
		}

		/// CALCULATED VARIABLES (USING FIELD VALUES)

		/// LINEARIZED CONVECTION

		/// u_star is the approximate linear advection at this time step
		std::vector<ZEROPTV> vel_pre_list = {u_pre_1, u_pre_2};
		ZEROPTV u_star = calc_linear_advection(vel_pre_list);
		VMSParams vms_params_current = calc_tau_vel(fe, u_star, Re, Ci_f, dt);
		double tauM = vms_params_current.tauM;

		/// Calculating the Elemental matrices requires loop over basis functions
		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the LHS Elemental matrix--------------------///

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN(a, i) * u_star(i);
			}

			/// Constructing the RHS vector

			/// All the terms involved in the equation
			ZEROPTV convFine, normal_NS, momentum;

			/// Terms of Normal Navier Stokes without any additional VMS terms
			for (int i = 0; i < nsd; i++) {
				normal_NS(i) =
						-fe.N(a) * (ndcf_time_ * (bdf_c[1] * u_pre_1(i) + bdf_c[2] * u_pre_2(i)) / dt_) * detJxW
						+ (fe.dN(a, i) * pStar) * detJxW
						+ (fe.N(a) * gravity(i) * detJxW)
						+ (fe.N(a) * forcing_current(i) * detJxW);

				/// first cross term
				convFine(i) =
						- (tauM * crossTermVelocityPart * (ndcf_time_ * (bdf_c[1] * u_pre_1(i) + bdf_c[2] * u_pre_2(i)) / dt_)) * detJxW
						- (tauM * crossTermVelocityPart * dpStar(i)) * detJxW
						+ (tauM * crossTermVelocityPart * (forcing_current(i))) * detJxW;

				/// Add all the terms to the vector
				momentum(i) = normal_NS(i) + convFine(i);
			}

			/// Adding all the calculated terms to the vector
			for (int i = 0; i < nsd; i++) {
				be((ndof) * a + i) += momentum(i);
			}
		}
		///-------------------Done with populating be vector-----------------------///
	}

	/** This is a linear-PSPG variational multiscale solver for Navier-Stokes momentum equations with the BDF1 and BDF2
	 * methods. Refer to Guermond et. al (2006),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void Integrands_NS_Momentum_PSPG_BDF12Ae(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM + 1;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt

		double Re = input_data_->Re;           ///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1 / Re;
		double Ci_f = this->input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf(); // # of basis functions
		double detJxW = fe.detJxW();

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x() = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y() = 0;
		} else {
			gravity.y() = -1.0 / Fr;
		}
		gravity.z() = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcing_current, forcing_pre_1, forcing_pre_2;
		if (input_data_->ifMMS) {
			forcing_pre_1 = calcForcingManufacturedSolutionNS(fe, t_ - dt_);
			if (t_ >= 2 * dt_) {
				forcing_pre_2 = calcForcingManufacturedSolutionNS(fe, t_ - 2 * dt_);
			}
		} else {
			forcing_pre_1 = {0.0, 0.0, 0.0};
			forcing_pre_2 = {0.0, 0.0, 0.0};
		}

		// FIELD VARIABLES
		ZEROPTV u_pre_1_c, u_pre_2_c, u_pre_3_c;
		for (int i = 0; i < nsd; i++) {
			u_pre_1_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre_2_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
			u_pre_3_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE3 + i);
		}

		// Velocity gradient tensors at previous steps
		ZeroMatrix<double> du_pre_1_c, du_pre_2_c;
		du_pre_1_c.redim(nsd, nsd);
		du_pre_2_c.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du_pre_1_c(i, j) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
				du_pre_2_c(i, j) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i, j);
			}
		}

		// Laplacian of velocities at previous steps
		ZEROPTV d2u_pre_1_c, d2u_pre_2_c;
		for (int dof = 0; dof < nsd; dof++) {
			for (int dir = 0; dir < nsd; dir++) {
				d2u_pre_1_c(dof) += p_data_->value2DerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE1 + dof, dir, dir);
				d2u_pre_2_c(dof) += p_data_->value2DerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE2 + dof, dir, dir);
			}
		}

		// Pressure values at previous steps
		double p_pre_1_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_pre_2_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_pre_3_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);

		// Derivatives of pressure at previous steps
		ZEROPTV dp_pre_1_c, dp_pre_2_c;
		for (int i = 0; i < nsd; i++) {
			dp_pre_1_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE1), i);
			dp_pre_2_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE2), i);
		}

		// CALCULATED VARIABLES (USING FIELD VALUES)
		// Time derivatives at previous steps
		ZEROPTV t_der_pre_1, t_der_pre_2;
		for (int i = 0; i < nsd; i++) {
			t_der_pre_1(i) = (bdf_c[0] * u_pre_1_c(i) + bdf_c[1] * u_pre_2_c(i) + bdf_c[2] * u_pre_3_c(i)) / dt_;
			t_der_pre_2(i) = (u_pre_2_c(i) - u_pre_3_c(i)) / dt_; // backward Euler approximation at t_(n-1)
		}

		// Convection & Diffusion terms at previous steps
		ZEROPTV convec_pre_1, convec_pre_2;
		ZEROPTV diffusion_pre_1, diffusion_pre_2;
		for (int i = 0; i < nsd; i++) {
			convec_pre_1(i) = 0.0;
			convec_pre_2(i) = 0.0;
			diffusion_pre_1(i) = 0.0;
			diffusion_pre_2(i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec_pre_1(i) += du_pre_1_c(i, j) * u_pre_1_c(j);
				convec_pre_2(i) += du_pre_2_c(i, j) * u_pre_2_c(j);
			}
			diffusion_pre_1(i) += Coe_diff * d2u_pre_1_c(i);
			diffusion_pre_2(i) += Coe_diff * d2u_pre_2_c(i);
		}

		///------------------------------------Coupling with PNP terms--------------------------------------------------///
		///Get gradient of electric potential
		ZEROPTV gradPotential, gradPotentialPrev1, gradPotentialPrev2;
		for (int i = 0; i < nsd; i++) {
			gradPotential(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL, i);
			gradPotentialPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL_PRE1, i);
			gradPotentialPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL_PRE2, i);
		}

		/// Charge density
		double rho_e = 0.0, rho_e_prev1, rho_e_prev2;
		for (int species = 0; species < nSpcs_; species++)  {
			double speciesCon = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1 + species);
			double speciesConPrev1 = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE1 + species);
			double speciesConPrev2 = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE2 + species);
			rho_e += z[species] * speciesCon;
			rho_e_prev1 += z[species] * speciesConPrev1;
			rho_e_prev2 += z[species] * speciesConPrev2;
		}

		/// body force of the NS equation. This couples NS equation with PNP equation.
		ZEROPTV pnpCoupling(0.0, 0.0, 0.0),
						pnpCouplingPrev1(0.0, 0.0, 0.0),
						pnpCouplingPrev2(0.0, 0.0, 0.0);
		if (!input_data_->ifPurePNP){
			for (int i = 0; i < nsd; i++){
				pnpCoupling(i) = gradPotential(i) * rho_e * -1.0;
				pnpCouplingPrev1(i) = gradPotentialPrev1(i) * rho_e_prev1 * -1.0;
				pnpCouplingPrev2(i) = gradPotentialPrev2(i) * rho_e_prev2 * -1.0;
			}
		}
		///--------------------------------------------------------------------------------------------------------------///

		// Residuals at previous steps
		ZEROPTV res_pre_1, res_pre_2;
		for (int i = 0; i < nsd; i++) {
			res_pre_1(i) = ndcf_time_ * t_der_pre_1(i) + ndcf_conv_ * convec_pre_1(i) - ndcf_diff_ * diffusion_pre_1(i) +
					ndcf_pres_ * dp_pre_1_c(i) - forcing_pre_1(i) - ndcf_pnp_coupling_ * pnpCouplingPrev1(i);
			res_pre_2(i) = ndcf_time_ * t_der_pre_2(i) + ndcf_conv_ * convec_pre_2(i) - ndcf_diff_ * diffusion_pre_2(i) +
					ndcf_pres_ * dp_pre_2_c(i) - forcing_pre_2(i) - ndcf_pnp_coupling_ * pnpCouplingPrev2(i);
		}

		// Calculate total velocities u = (u^c + u^f) at previous steps
		VMSParams vms_params1 = calc_tau_vel(fe, u_pre_1_c, Re, Ci_f, dt);
		VMSParams vms_params2 = calc_tau_vel(fe, u_pre_2_c, Re, Ci_f, dt);
		ZEROPTV u_pre_1, u_pre_2;
		for (int i = 0; i < nsd; i++) {
			u_pre_1(i) = u_pre_1_c(i) - vms_params1.tauM * res_pre_1(i);
			u_pre_2(i) = u_pre_2_c(i) - vms_params2.tauM * res_pre_2(i);
		}


		/// LINEARIZED CONVECTION
		// u_star is the approximate linear advection at this time step
		std::vector<ZEROPTV> vel_pre_list = {u_pre_1, u_pre_2};
		ZEROPTV u_star = calc_linear_advection(vel_pre_list);
		// and u_star_c is the coarse scale counterpart of u_star
		std::vector<ZEROPTV> vel_pre_c_list = {u_pre_1_c, u_pre_2_c};
		ZEROPTV u_star_c = calc_linear_advection(vel_pre_c_list);
		VMSParams vms_params_current = calc_tau_vel(fe, u_star_c, Re, Ci_f, dt);
		double tauM = vms_params_current.tauM;
		double tauC = vms_params_current.tauC;

		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Jacobian operator--------------------///
			// NS terms

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u_star (i);
			}

			for (int b = 0; b < n_basis_functions; b++) {

				/// Convection term
				double conv = 0.0;
				for (int i = 0; i < nsd; i++) {
					conv += fe.dN (b, i) * u_star (i);
				}

				/// Adding terms to the Jacobian matrix.
				for (int i = 0; i < nsd; i++) {
					// Time derivative
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							(fe.N (a) * (ndcf_time_ * bdf_c[0] * fe.N (b) / dt_)) * detJxW;
					//validateMatrix(Ae);

					// Coarse scale convection term
					Ae ((ndof) * a + i, (ndof) * b + i) += (fe.N (a) * (ndcf_conv_ * conv)) * detJxW;
					//validateMatrix(Ae);

					// Coarse scale (only) diffusion term
					for (int j = 0; j < nsd; j++) {
						Ae ((ndof) * a + i, (ndof) * b + i) +=
								Coe_diff * fe.dN (a, j) * (ndcf_diff_ * fe.dN (b, j)) * detJxW;
					}
					//validateMatrix(Ae);

					Ae (ndof * a + i, ndof * b + nsd) +=
							-fe.dN (a, i) * (ndcf_pres_ * fe.N (b)) * detJxW;
					//validateMatrix(Ae);

				}


				/** crossTerm 1: This term is written as, (w, grad{u u'}), which is weakened to (grad{w}, u\innerproduct u')
				 * which is then linearised. Here u' is fine scale velocity and u is resolved velocity. The fine scale
				 * velocity is approximated as -tau*Res{NS} (Equation 58 bazilevs et al. (2007)),
				 */
				for (int i = 0; i < nsd; i++) {

					/// Contribution of laplacian of velocity(diffusion) to the diagonal
					double diff_J = 0;
					for (int j = 0; j < nsd; j++) {
						diff_J += Coe_diff * fe.d2N (b, j, j);
					}
					/** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from Bazilevs et al. (2007)) you get two
					 * terms, one goes in diagonal and other in the non-diagonal parts of the matrix PS: crossTermVelocityPart
					 * is (u. grad{w} => u*dN(a))
					 */
					/// Diagonal part
					Ae ((ndof) * a + i, (ndof) * b + i) +=
							crossTermVelocityPart * tauM *
							((ndcf_time_ * bdf_c[0] * fe.N (b) / dt_)
							 + (ndcf_conv_ * conv)
							 - (ndcf_diff_ * diff_J)) * detJxW;
					//validateMatrix(Ae);

					Ae (ndof * a + i, ndof * b + nsd) +=
							crossTermVelocityPart * tauM * (ndcf_pres_ * fe.dN (b, i)) * detJxW;
					//validateMatrix(Ae);

				}

				// Grad-div equivalent
				/** This term represents the fine scale pressure given by
				 * (grad{w}, \delta{tauC*res{Cont}}), where Cont is the continuity
				 * equation. See third term in equation (71) in Bazilev et al. (2007)
				 */
				for (int i = 0; i < nsd; i++) {
					for (int j = 0; j < nsd; j++) {
						Ae (ndof * a + i, ndof * b + j) +=
								fe.dN (a, i) * tauC * fe.dN (b, j) * detJxW;
						//validateMatrix(Ae);
					}
				}

				// Divergence eqution below
				/// pspg term from VMS theory: See second term in equation (71) from Bazilevs et al. (2007)
				for (int j = 0; j < nsd; j++) {
					double diff_J = 0;
					for (int k = 0; k < nsd; k++) {
						diff_J += Coe_diff * fe.d2N (b, k, k);
					}
					Ae (ndof * a + nsd, ndof * b + j) +=
							///(q, \partial_j(u_j^{c,n+1}))
							fe.N (a) * fe.dN (b, j) * detJxW +
							///(\partial_j(q), tau_M * Res(u_j^{c,n+1},p^{c,n+1}))
							/// (\partial_j(q), tau_M * (time_der + (u^*_k)\partial_j(u_i^{c,n+1}) - diffusion))
							fe.dN (a, j) * tauM
							* ((ndcf_time_ * bdf_c[0] * fe.N (b) / dt) + (ndcf_conv_ * conv) - (ndcf_diff_ * diff_J)) * detJxW;
					//validateMatrix(Ae);
					/// (\partial_i(q), tau_M * \partial_i p^{c,n+1})
					Ae (ndof * a + nsd, ndof * b + nsd) +=
							fe.dN (a, j) * tauM * (ndcf_pres_ * fe.dN (b, j)) * detJxW;
					//validateMatrix(Ae);
				}
			}
			///--------------------------------------Done with Elemental operator Matrix-----------------------------------///
		}
//        for (int a = 0; a < Ae.nx_; a++) {
//            for (int b = 0; b < Ae.ny_; b++) {
//
//                std::cout << Ae (a, b) << " ";
//            }
//            std::cout << "\n";
//
//        }



	}

	/** This is a linear-PSPG variational multiscale solver for Navier-Stokes momentum equations with the BDF1 and BDF2
	 * methods. Refer to Guermond et. al (2006),
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void Integrands_NS_Momentum_PSPG_BDF12be(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {
		const DENDRITE_UINT ndof = DIM + 1;
		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt



		double Re = input_data_->Re;           ///<Reynolds number
		double Fr = input_data_->Fr;

		/// double Re = this->input_data_->integralReNo;
		double Coe_diff = 1 / Re;
		double Ci_f = this->input_data_->Ci_f;
		int nsd = DIM;
		int n_basis_functions = fe.nbf(); // # of basis functions
		double detJxW = fe.detJxW();

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x() = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y() = 0;
		} else {
			gravity.y() = -1.0 / Fr;
		}
		gravity.z() = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcing_current, forcing_pre_1, forcing_pre_2;
		if (input_data_->ifMMS) {
			forcing_current = calcForcingManufacturedSolutionNS(fe, t_);
			forcing_pre_1 = calcForcingManufacturedSolutionNS(fe, t_ - dt_);
			if (t_ >= 2 * dt_) {
				forcing_pre_2 = calcForcingManufacturedSolutionNS(fe, t_ - 2 * dt_);
			}
		} else {
			forcing_pre_1 = {0.0, 0.0, 0.0};
			forcing_pre_2 = {0.0, 0.0, 0.0};
		}

		// FIELD VARIABLES
		ZEROPTV u_pre_1_c, u_pre_2_c, u_pre_3_c;
		for (int i = 0; i < nsd; i++) {
			u_pre_1_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre_2_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
			u_pre_3_c(i) = p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE3 + i);
		}

		// Velocity gradient tensors at previous steps
		ZeroMatrix<double> du_pre_1_c, du_pre_2_c;
		du_pre_1_c.redim(nsd, nsd);
		du_pre_2_c.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du_pre_1_c(i, j) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
				du_pre_2_c(i, j) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i, j);
			}
		}

		// Laplacian of velocities at previous steps
		ZEROPTV d2u_pre_1_c, d2u_pre_2_c;
		for (int dof = 0; dof < nsd; dof++) {
			for (int dir = 0; dir < nsd; dir++) {
				d2u_pre_1_c(dof) += p_data_->value2DerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE1 + dof, dir, dir);
				d2u_pre_2_c(dof) += p_data_->value2DerivativeFEM(fe, NSPNPNodeData::VEL_X_PRE2 + dof, dir, dir);
			}
		}

		// Pressure values at previous steps
		double p_pre_1_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_pre_2_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_pre_3_c = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);

		// Derivatives of pressure at previous steps
		ZEROPTV dp_pre_1_c, dp_pre_2_c;
		for (int i = 0; i < nsd; i++) {
			dp_pre_1_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE1), i);
			dp_pre_2_c(i) = this->p_data_->valueDerivativeFEM(fe, (NSPNPNodeData::PRESSURE_PRE2), i);
		}

		// CALCULATED VARIABLES (USING FIELD VALUES)
		// Time derivatives at previous steps
		ZEROPTV t_der_pre_1, t_der_pre_2;
		for (int i = 0; i < nsd; i++) {
			t_der_pre_1(i) = (bdf_c[0] * u_pre_1_c(i) + bdf_c[1] * u_pre_2_c(i) + bdf_c[2] * u_pre_3_c(i)) / dt_;
			t_der_pre_2(i) = (u_pre_2_c(i) - u_pre_3_c(i)) / dt_; // backward Euler approximation at t_(n-1)
		}

		// Convection & Diffusion terms at previous steps
		ZEROPTV convec_pre_1, convec_pre_2;
		ZEROPTV diffusion_pre_1, diffusion_pre_2;
		for (int i = 0; i < nsd; i++) {
			convec_pre_1(i) = 0.0;
			convec_pre_2(i) = 0.0;
			diffusion_pre_1(i) = 0.0;
			diffusion_pre_2(i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec_pre_1(i) += du_pre_1_c(i, j) * u_pre_1_c(j);
				convec_pre_2(i) += du_pre_2_c(i, j) * u_pre_2_c(j);
			}
			diffusion_pre_1(i) += Coe_diff * d2u_pre_1_c(i);
			diffusion_pre_2(i) += Coe_diff * d2u_pre_2_c(i);
		}

		///------------------------------------Coupling with PNP terms--------------------------------------------------///
		///Get gradient of electric potential
		ZEROPTV gradPotential, gradPotentialPrev1, gradPotentialPrev2;
		for (int i = 0; i < nsd; i++) {
			gradPotential(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL, i);
			gradPotentialPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL_PRE1, i);
			gradPotentialPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL_PRE2, i);
		}

		/// Charge density
		double rho_e = 0.0, rho_e_prev1, rho_e_prev2;
		for (int species = 0; species < nSpcs_; species++)  {
			double speciesCon = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1 + species);
			double speciesConPrev1 = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE1 + species);
			double speciesConPrev2 = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE2 + species);
			rho_e += z[species] * speciesCon;
			rho_e_prev1 += z[species] * speciesConPrev1;
			rho_e_prev2 += z[species] * speciesConPrev2;
//            if (speciesCon!=0)
//            printf("Integrands_NS_Momentum_PSPG_BDF12be %f\n", speciesCon);
		}
//        printf("phi=   %f\n", this->p_data_->valueFEM(fe, NSPNPNodeData::POTENTIAL));
//        printf("phi=   %f , %f\n", this->p_data_->valueFEM(fe, NSPNPNodeData::POTENTIAL), this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X));

		/// body force of the NS equation. This couples NS equation with PNP equation.
		ZEROPTV pnpCoupling(0.0, 0.0, 0.0),
				pnpCouplingPrev1(0.0, 0.0, 0.0),
				pnpCouplingPrev2(0.0, 0.0, 0.0);
		if (!input_data_->ifPurePNP){
			for (int i = 0; i < nsd; i++){
				pnpCoupling(i) = gradPotential(i) * rho_e * -1.0;
				pnpCouplingPrev1(i) = gradPotentialPrev1(i) * rho_e_prev1 * -1.0;
				pnpCouplingPrev2(i) = gradPotentialPrev2(i) * rho_e_prev2 * -1.0;

			}
		}


		///--------------------------------------------------------------------------------------------------------------///


		// Residuals at previous steps
		ZEROPTV res_pre_1, res_pre_2;
		for (int i = 0; i < nsd; i++) {
			res_pre_1(i) = ndcf_time_ * t_der_pre_1(i) + ndcf_conv_ * convec_pre_1(i) - ndcf_diff_ * diffusion_pre_1(i) +
			               ndcf_pres_ * dp_pre_1_c(i) - forcing_pre_1(i) - ndcf_pnp_coupling_ * pnpCouplingPrev1(i);
			res_pre_2(i) = ndcf_time_ * t_der_pre_2(i) + ndcf_conv_ * convec_pre_2(i) - ndcf_diff_ * diffusion_pre_2(i) +
			               ndcf_pres_ * dp_pre_2_c(i) - forcing_pre_2(i) - ndcf_pnp_coupling_ * pnpCouplingPrev2(i);
		}

		// Calculate total velocities u = (u^c + u^f) at previous steps
		VMSParams vms_params1 = calc_tau_vel(fe, u_pre_1_c, Re, Ci_f, dt);
		VMSParams vms_params2 = calc_tau_vel(fe, u_pre_2_c, Re, Ci_f, dt);
		ZEROPTV u_pre_1, u_pre_2;
		for (int i = 0; i < nsd; i++) {
			u_pre_1(i) = u_pre_1_c(i) - vms_params1.tauM * res_pre_1(i);
			u_pre_2(i) = u_pre_2_c(i) - vms_params2.tauM * res_pre_2(i);
		}

		// LINEARIZED CONVECTION
		// u_star is the approximate linear advection at this time step
		std::vector<ZEROPTV> vel_pre_list = {u_pre_1, u_pre_2};
		ZEROPTV u_star = calc_linear_advection(vel_pre_list);
		// and u_star_c is the coarse scale counterpart of u_star
		std::vector<ZEROPTV> vel_pre_c_list = {u_pre_1_c, u_pre_2_c};
		ZEROPTV u_star_c = calc_linear_advection(vel_pre_c_list);
		VMSParams vms_params_current = calc_tau_vel(fe, u_star_c, Re, Ci_f, dt);
		double tauM = vms_params_current.tauM;
		double tauC = vms_params_current.tauC;

		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Jacobian operator--------------------///
			// NS terms

			/** Course scale terms for the cross terms (terms 4 and 5 in equation 52 in Bazilevs et al. (2007)
			 * PS: these are not calculated using the SUPG class, and this is an approximation of (u, grad{u})
			 */
			double crossTermVelocityPart = 0.0;
			for (int i = 0; i < nsd; i++) {
				crossTermVelocityPart += fe.dN (a, i) * u_star (i);
			}
			/// All the terms involved in the equation
			ZEROPTV convFine, normal_NS, momentum;

			/// Terms of Normal Navier Stokes without any additional VMS terms
			for (int i = 0; i < nsd; i++) {
				normal_NS (i) =
						- fe.N (a) * (ndcf_time_ * (bdf_c[1] * u_pre_1_c (i) + bdf_c[2] * u_pre_2_c (i)) / dt_) * detJxW
						+ (fe.N (a) * gravity (i) * detJxW)
						+ (fe.N (a) * forcing_current (i) * detJxW)
						- (fe.N (a) * ndcf_pnp_coupling_ * pnpCoupling(i) * detJxW)
						;

				/// first cross term
				convFine (i) =
						-(tauM * crossTermVelocityPart * (ndcf_time_ * (bdf_c[1] * u_pre_1_c (i) + bdf_c[2] * u_pre_2_c (i)) / dt_))
						* detJxW
						+ (tauM * crossTermVelocityPart * (forcing_current (i))) * detJxW
						- (tauM * crossTermVelocityPart * ndcf_pnp_coupling_ * pnpCoupling(i) * detJxW)
						;

				/// Add all the terms to the vector
				momentum (i) = normal_NS (i) + convFine (i);
			}

			double pspg = 0.0;
			for (int i = 0; i < nsd; i++) {
				/// Fine scale contribution to pressure
				pspg +=
						-(fe.dN (a, i) * tauM * (ndcf_time_ * (bdf_c[1] * u_pre_1_c (i) + bdf_c[2] * u_pre_2_c (i)) / dt_)) * detJxW
						+ (fe.dN (a, i) * tauM * forcing_current (i)) * detJxW;
			}

			/// Adding all the calculated terms to the vector
			for (int i = 0; i < nsd; i++) {
				be ((ndof) * a + i) += momentum (i);
			}
			be (ndof * a + nsd) += pspg;
		}
		///------------------------------------------Done with populating be vector--------------------------------------///
	}

	#pragma mark Projection-integrands
	/** Helmholtz decomposition projection step. It updates pressure
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void NSIntegrandsPressurePoissonAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		int degOfFreedom = 1;

		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		double time_coeff = (input_data_->timescheme == NSPNPInputData::BDF) ? (1.0 / bdf_c[0]) : theta;

		double Coeff_pressure_poisson = 1.0;

		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Calculating the Elemental operator------------///
			for (int b = 0; b < n_basis_functions; b++) {
				for (int i = 0; i < nsd; i++) {
					Ae ((degOfFreedom * a) + 0, (degOfFreedom * b) + 0) +=
							///(\partial_i(w), \partial_i(p))
							(Coeff_pressure_poisson * fe.dN (a, i) * fe.dN (b, i)) * detJxW;
				}
			}
		}

	}

	/** Helmholtz decomposition projection step. It updates pressure
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void NSIntegrandsPressurePoissonbe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {

		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt
		double Re = this->input_data_->Re;           ///<Reynolds number
		double Fr = this->input_data_->Fr;
		int degOfFreedom = 1;

		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		double Coe_diff = 1.0 / Re;
		double Ci_f = this->input_data_->Ci_f;


		/// field variables
		ZEROPTV u, u_pre1, u_pre2;
		for (int i = 0; i < nsd; i++) {
			u(i)      = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X + i);
			u_pre1(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre2(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
		}

		VMSParams vms_params = calc_tau_vel(fe, u, Re, Ci_f, dt);
		const double tauM = vms_params.tauM;

		/// Define velocity gradient tensor
		ZeroMatrix<double> du;
		du.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
			}
		}


		ZeroMatrix<double> duPrev;
		duPrev.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				duPrev (i, j) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
			}
		}

		/// Pressure fields setup
		double p       = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE);
		double p_prev1 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_prev2 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE2);
//		double pStar = p_extrap_c[0] * p + p_extrap_c[1] * p_prev1 + p_extrap_c[2] * p_prev2;
        assert(FEQUALS(p,p_prev1));
		/// Get gradient of pressure at previous timestep
		ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE, i);
			dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
			dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE2, i);
//			dpStar (i) = p_extrap_c[0] * dp(i) + p_extrap_c[1] * dpPrev1(i) + p_extrap_c[2] * dpPrev2(i);
		}

		ZEROPTV d2u;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += p_data_->value2DerivativeFEM (fe, dof, dir, dir);
			}
		}

		/// Work in progress on rotational version:
		/// the below term is used for the  rotational form
		/// it is \grad(\divergence(u)) which is a vector of size nsd
		ZEROPTV grad_of_div;
		for (int idof = 0; idof < nsd; idof++) {
			for (int i = 0; i < nsd; i++) {
				grad_of_div (idof) += p_data_->value2DerivativeFEM (fe, i, idof, i);
			}
		}

		ZEROPTV convec;
		ZEROPTV diffusion;
		for (int i = 0; i < nsd; i++) {
			convec(i) = 0.0;
			diffusion(i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec(i) += du(i, j) * u(j);
			}
			diffusion(i) += Coe_diff * d2u(i);
		}

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x () = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y () = 0;
		} else {
			gravity.y () = -1.0 / Fr;
		}
		gravity.z () = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcingNS;
		if (input_data_->ifMMS) {
			forcingNS = calcForcingManufacturedSolutionNS (fe, t_);
		} else {
			forcingNS = {0.0, 0.0, 0.0};
		}

		ZEROPTV timeder;
		for (int i = 0; i < nsd; i++) {
			if (input_data_->timescheme == NSPNPInputData::THETA) {
				timeder (i) = (u (i) - u_pre1 (i)) / dt_;
			} else if (input_data_->timescheme == NSPNPInputData::BDF) {
				timeder (i) = (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_;
			}
		}

		ZEROPTV NS;
		for (int i = 0; i < nsd; i++) {
			NS (i) = timeder(i) + convec(i) + dpStar(i) - diffusion(i) - gravity(i) - forcingNS(i);
		}
		double time_coeff = (input_data_->timescheme == NSPNPInputData::BDF) ? (1.0 / bdf_c[0]) : theta;

		double Coeff_pressure = time_coeff * dt_;
		double Coeff_vel_update = 1.0;
		double Coeff_velPred = (1.0 / Coeff_pressure);
		double Coeff_pressure_poisson = 1.0;

		for (int a = 0; a < n_basis_functions; a++) {
			///--------------Elemental vector-------------///
			double pressurePoissonRHS1 = 0.0;
			double pressurePoissonRHS2 = 0.0;
			double pressurePoissonGradDiv = 0.0;
			for (int i = 0; i < nsd; i++) {
				pressurePoissonRHS1 += -(Coeff_velPred * fe.N (a) * du(i, i));
				//pressurePoissonRHS1 += (Coeff_velPred * fe.dN(a, i) * u(i));
				pressurePoissonRHS1 += -Coeff_velPred * fe.dN (a, i) * tauM * NS (i);
				pressurePoissonRHS2 += (Coeff_pressure_poisson * fe.dN (a, i) * dpStar (i));
				pressurePoissonGradDiv += -Coeff_pressure_poisson * Coe_diff * fe.dN (a, i) * grad_of_div (i);
			}
			/// Contribution of this term goes to last element of the vector
			double contContrib = (pressurePoissonRHS1 * detJxW) + (pressurePoissonRHS2 * detJxW)
			                     + (rot_form_flag * pressurePoissonGradDiv * detJxW);

			/// Adding all the calculated terms to the vector
			be ((degOfFreedom * a)) += contContrib;
		}
	}

	/** Helmholtz decomposition projection step. Boundary terms for elemental vector
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void NSIntegrands4sidePressurePoissonbe(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be, const double dt,
			const double t) {
		enum {
				LEFT = 1,
				RIGHT = 2,
				BOTTOM = 3,
				TOP = 4,
				CIRCLE = 7,
		};

		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt
		int degOfFreedom = 1; ///4 degrees of freedom

		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		ZEROPTV position = fe.position ();
		double detJxW = fe.detJxW ();

		double Fr = input_data_->Fr;
		double Re = input_data_->Re;

		double Coe_diff = 1.0 / Re;
		double Ci_f = this->input_data_->Ci_f;

		VMSParams vms_params = calc_tau_surf (fe, Re, Ci_f, dt);
		const double tauM = stab_in_ppe_vu * vms_params.tauM;

		/// field variables
		ZEROPTV u, u_pre1, u_pre2;
		for (int i = 0; i < nsd; i++) {
			u(i)      = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X + i);
			u_pre1(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre2(i) = this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X_PRE2 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du;
		du.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
			}
		}


		ZeroMatrix<double> duPrev;
		duPrev.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				duPrev (i, j) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
			}
		}

		/// Pressure fields setup
		double p       = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE);
		double p_prev1 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_prev2 = p_data_->valueFEM(fe, NSPNPNodeData::PRESSURE_PRE2);
		double pStar = p_extrap_c[0] * p + p_extrap_c[1] * p_prev1 + p_extrap_c[2] * p_prev2;

		/// Get gradient of pressure at previous timestep
		ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp(i)      = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE, i);
			dpPrev1(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE1, i);
			dpPrev2(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::PRESSURE_PRE2, i);
			dpStar (i) = p_extrap_c[0] * dp(i) + p_extrap_c[1] * dpPrev1(i) + p_extrap_c[2] * dpPrev2(i);
		}

		ZEROPTV d2u;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += p_data_->value2DerivativeFEM (fe, dof, dir, dir);
			}
		}

		ZEROPTV convec;
		ZEROPTV diffusion;
		for (int i = 0; i < nsd; i++) {
			convec (i) = 0.0;
			diffusion (i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec (i) += du (i, j) * u (j);
			}
			diffusion (i) += Coe_diff * d2u (i);
		}

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x () = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y () = 0;
		} else {
			gravity.y () = -1.0 / Fr;
		}
		gravity.z () = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcingNS;
		if (input_data_->ifMMS) {
			forcingNS = calcForcingManufacturedSolutionNS (fe, t_);
		} else {
			forcingNS = {0.0, 0.0, 0.0};
		}

		ZEROPTV timeder;
		for (int i = 0; i < nsd; i++) {
			if (input_data_->timescheme == NSPNPInputData::THETA) {
				timeder (i) = (u (i) - u_pre1 (i)) / dt_;
			} else if (input_data_->timescheme == NSPNPInputData::BDF) {
				timeder (i) = (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_;
			}
		}

		ZEROPTV NS;
		for (int i = 0; i < nsd; i++) {
			NS (i) = timeder (i)
			         + convec (i) + dpStar (i) - diffusion (i) - gravity (i) - forcingNS (i);
		}

		ZEROPTV normal = fe.surface ()->normal ();

		double u_dot_n = 0.0;
		for (int i = 0; i < nsd; i++) {
			u_dot_n += -tauM * NS (i) * normal (i);
		}

		double time_coeff = (input_data_->timescheme == NSPNPInputData::BDF) ? (1.0 / bdf_c[0]) : theta;

		double Coeff_pressure = time_coeff * dt_;
		double Coeff_vel_update = 1.0;
		double Coeff_velPred = (1.0 / Coeff_pressure);
		double Coeff_pressure_poisson = 1.0; //(1.0 / rhoMix);

		for (int a = 0; a < n_basis_functions; a++) {
			be(degOfFreedom * a) += -Coeff_velPred * fe.N (a) * u_dot_n * detJxW;
		}
	}

	#pragma mark Velocity update-integrands
	/** Helmholtz decomposition projection step. It updates pressure and velocity
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void NSIntegrandsVelocityUpdateAe(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t) {
		double dt_ = dt;
		double t_ = t + dt_;  ///< current time is elapsed + dt
		double Re = input_data_->Re;///<Reynolds number
		double Fr = input_data_->Fr;
		int degOfFreedom = DIM; /// degrees of freedom are the same as no of velocity components

		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		double time_coeff = (input_data_->timescheme == NSPNPInputData::BDF) ? (1.0 / bdf_c[0]) : theta;

		double Coeff_pressure = time_coeff * dt_;
		double Coeff_vel_update = 1.0;
		double Coeff_velPred = (1.0 / Coeff_pressure);
		double Coeff_pressure_poisson = 1.0;

		for (int a = 0; a < n_basis_functions; a++) {
			///---------------Filling the elemental operator ------------///
			for (int b = 0; b < n_basis_functions; b++) {
				/// All the terms which does not contain any fine scale corrections for the first three rows of the elemental
				/// matrix
				for (int i = 0; i < degOfFreedom; i++) {
					Ae ((degOfFreedom * a) + i, (degOfFreedom * b) + i) +=
							///(w_i, v_i)
							(Coeff_vel_update * fe.N (a) * fe.N (b)) * detJxW;
				}
			}
		}
	}

	/** Helmholtz decomposition velocity update step. It updates velocity to the solenoidal variation
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be vector to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @return
	 */
	void NSIntegrandsVelocityUpdatebe(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t) {

		double dt_ = dt;
		double t_ = t + dt_; ///< current time is elapsed + dt
		double Re = input_data_->Re; ///<Reynolds number
		double Fr = input_data_->Fr;
		int degOfFreedom = DIM; ///4 degrees of freedom

		int nsd = DIM;
		int n_basis_functions = fe.nbf (); // # of basis functions
		double detJxW = fe.detJxW ();

		double Coe_diff = 1.0 / Re;
		double Ci_f = this->input_data_->Ci_f;

		VMSParams vms_params = calc_tau (fe, Re, Ci_f, dt);
		const double tauM = stab_in_ppe_vu * vms_params.tauM;

		/// field variables
		ZEROPTV u, u_pre1, u_pre2;
		for (int i = 0; i < nsd; i++) {
			u (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X + i);
			u_pre1 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i);
			u_pre2 (i) = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X_PRE2 + i);
		}

		/// Define velocity gradient tensor
		ZeroMatrix<double> du;
		du.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				du (i, j) = this->p_data_->valueDerivativeFEM (fe, i, j);
			}
		}


		ZeroMatrix<double> duPrev;
		duPrev.redim (nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				duPrev (i, j) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::VEL_X_PRE1 + i, j);
			}
		}

		/// Pressure fields setup
		double p = p_data_->valueFEM (fe, NSPNPNodeData::PRESSURE);
		double p_prev1 = p_data_->valueFEM (fe, NSPNPNodeData::PRESSURE_PRE1);
		double p_prev2 = p_data_->valueFEM (fe, NSPNPNodeData::PRESSURE_PRE2);
		double pStar = p_extrap_c[0] * p + p_extrap_c[1] * p_prev1 + p_extrap_c[2] * p_prev2;

		/// Get gradient of pressure at previous timestep
		ZEROPTV dp, dpPrev1, dpPrev2, dpStar;
		for (int i = 0; i < nsd; i++) {
			dp (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::PRESSURE, i);
			dpPrev1 (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::PRESSURE_PRE1, i);
			dpPrev2 (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::PRESSURE_PRE2, i);
			dpStar (i) = p_extrap_c[0] * dp(i) + p_extrap_c[1] * dpPrev1(i) + p_extrap_c[2] * dpPrev2(i);
		}

		double contAvg = 0.0;
		for (int i = 0; i < nsd; i++) {
			contAvg += 0.5 * (du (i, i) + duPrev (i, i));
		}

		ZEROPTV d2u;
		/// loop for three directions of velocity (MAX NSD is no of velocity
		/// directions in the problem)
		/// PS: This assumes the first three degrees are always velocity, it is
		/// wise to keep it that way
		for (int dof = 0; dof < nsd; dof++) {
			/// Summing over three directions of velocity as it is laplacian
			for (int dir = 0; dir < nsd; dir++) {
				d2u (dof) += p_data_->value2DerivativeFEM (fe, dof, dir, dir);
			}
		}

		// the below term is used for the  rotational form
		// it is \grad(\divergence(u)) which is a vector of size nsd
		ZEROPTV grad_of_div;
		for (int idof = 0; idof < nsd; idof++) {
			for (int i = 0; i < nsd; i++) {
				grad_of_div(idof) += p_data_->value2DerivativeFEM (fe, i, idof, i);
			}
		}

		ZEROPTV convec;
		ZEROPTV diffusion;
		for (int i = 0; i < nsd; i++) {
			convec (i) = 0.0;
			diffusion (i) = 0.0;
			for (int j = 0; j < nsd; j++) {
				convec (i) += du (i, j) * u (j);
			}
			diffusion (i) += Coe_diff * d2u (i);
		}

		/// Contribution of gravity in terms of Froude number
		ZEROPTV gravity;
		gravity.x () = 0.0;
		/// if Fr is specified to be zero in the config file turn gravity off
		if (Fr < 1e-10) {
			gravity.y () = 0;
		} else {
			gravity.y () = -1.0 / Fr;
		}
		gravity.z () = 0.0;

		/// External forcing: In the current implementation this term is non-zero
		/// only in the case of Manufactured solutions test of CHNS
		ZEROPTV forcingNS;
		if (input_data_->ifMMS) {
			forcingNS = calcForcingManufacturedSolutionNS (fe, t_);
		} else {
			forcingNS = {0.0, 0.0, 0.0};
		}

		ZEROPTV timeder;
		for (int i = 0; i < nsd; i++) {
			if (input_data_->timescheme == NSPNPInputData::THETA) {
				timeder (i) = (u (i) - u_pre1 (i)) / dt_;
			} else if (input_data_->timescheme == NSPNPInputData::BDF) {
				timeder (i) = (bdf_c[0] * u (i) + bdf_c[1] * u_pre1 (i) + bdf_c[2] * u_pre2 (i)) / dt_;
			}
		}

		ZEROPTV NS;
		for (int i = 0; i < nsd; i++) {
			NS (i) = timeder(i) + convec(i) + dpStar(i) - diffusion(i) - gravity(i) - forcingNS(i);
		}

		double time_coeff = (input_data_->timescheme == NSPNPInputData::BDF) ? (1.0 / bdf_c[0]) : theta;

		double Coeff_pressure = time_coeff * dt_;
		double Coeff_vel_update = 1.0;
		double Coeff_velPred = (1.0 / Coeff_pressure);
		double Coeff_pressure_poisson = 1.0;

		for (int a = 0; a < n_basis_functions; a++) {
			///--------------Filling the elemental vector-------------///
			for (int i = 0; i < degOfFreedom; i++) {
				be ((degOfFreedom * a) + i) +=
						/// (w_i, v_i) for each i
						(Coeff_vel_update * fe.N (a) * u (i)) * detJxW
						- Coeff_vel_update * fe.N (a) * tauM * NS (i) * detJxW
						/// (w_i, \delta t/2 \partial_i(p_prev) ) for each i
						- (Coeff_pressure * fe.N (a) * (dp (i) - dpStar (i))) * detJxW
						- (rot_form_flag * Coeff_pressure * Coe_diff * fe.N (a) * grad_of_div (i)) * detJxW;
			}
		}
	}

#pragma mark PNP Integrands

	/**
	 * References for VMS
	 * 1. Simplified VMS with charge neutrality assumption for NS-PNP
	 *    - G. Bauer et al. / Comput. Methods Appl. Mech. Engrg. 223–224 (2012) 199–210
	 * 2. Whole VMS terms for PNP
	 *    - Shadid, John Nicolas, et al. Simulation of neutron radiation damage in silicon semiconductor devices, 2007.
	 * 3. VMS for Navier-Stokes
	 *    - Bazilevs et al. / Comput. Methods Appl. Mech. Engrg. 197 (2007) 173–201
	 */


#pragma mark Manufactured solutions helper NS

	/**Function to calculate the external forcing in manufactured solutions case for Navier Stokes
	 * @param fe FEMElm object for positions
	 * @param t current time
	 * @return forcing due to manufactured solutions
	 */
	ZEROPTV calcForcingManufacturedSolutionNS(const FEMElm &fe, const double &t) {
		ZEROPTV forcingNS;
		ZEROPTV location = fe.position ();

		forcingNS.x () = forcingNSXdir (location, t);
		forcingNS.y () = forcingNSYdir (location, t);
		forcingNS.z () = 0.0;
		return forcingNS;
	}


	/// Function which calculate forcing for manufactured solutions case for
	/// Navier-Stokes in X direction
	double forcingNSXdir(const ZEROPTV &location, const double &t) {
		double time = t;

		/// Get non dimensional numbers
		double Re = this->input_data_->Re;
		double We = 1.0;
		double Fr = this->input_data_->Fr; /// Froude number

		double forcingX = 0.0;
		double x = location.x(), y = location.y(), z = location.z();

		if (input_data_->ifPurePNP){
			forcingX =
					ndcf_conv_ *
					((M_PI * M_PI * M_PI) * cos (x * M_PI) * pow (sin (x * M_PI), 3.0) * pow (sin (y * M_PI * 2.0), 2.0) *
					 pow (sin (t), 2.0) * 2.0 -
					 (M_PI * M_PI * M_PI) * cos (y * M_PI * 2.0) * pow (sin (x * M_PI), 2.0) * sin (x * M_PI * 2.0) *
					 pow (sin (y * M_PI), 2.0) * pow (sin (t), 2.0) * 2.0) - (ndcf_diff_ * ((M_PI * M_PI * M_PI) *
					                                                                        pow (cos (x * M_PI), 2.0) *
					                                                                        sin (y * M_PI * 2.0) * sin (t) *
					                                                                        2.0 - (M_PI * M_PI * M_PI) *
					                                                                              pow (sin (x * M_PI), 2.0) *
					                                                                              sin (y * M_PI * 2.0) *
					                                                                              sin (t) * 6.0)) / Re -
					ndcf_pres_ * M_PI * sin (x * M_PI) * sin (y * M_PI) * sin (t) +
					ndcf_time_ * M_PI * pow (sin (x * M_PI), 2.0) * sin (y * M_PI * 2.0) * cos (t);
		} else {
			double t2 = sin (t);
			double t3 = pi * x;
			double t4 = pi * y;
			double t5 = pi * pi * pi;
			double t6 = t2 * t2;
			double t7 = t3 * 2.0;
			double t8 = t4 * 2.0;
			double t9 = cos (t3);
			double t10 = sin (t3);
			double t11 = sin (t4);
			double t12 = cos (t8);
			double t13 = sin (t7);
			double t14 = sin (t8);
			double t15 = t10 * t10;
			forcingX = ndcf_conv_ * (t5 * t6 * t9 * (t10 * t10 * t10) * (t14 * t14) * 2.0
			                          - t5 * t6 * (t11 * t11) * t12 * t13 * t15 * 2.0)
			            + (ndcf_diff_ * (t2 * t5 * t14 * t15 * 6.0 - t2 * t5 * (t9 * t9) * t14 * 2.0)) / Re
			            - ndcf_pres_ * pi * t2 * t10 * t11 + ndcf_time_ * pi * t14 * t15 * cos (t)
			            + ndcf_pnp_coupling_ * pi * t2 * t13 * t14 * (t2 * t14 * cos (t7) - t2 * t12 * t13) * 2.0;
		}
		return forcingX;
	}

	/// Function which calculate forcing for manufactured solutions case for
	/// Navier-Stokes in Y direction
	double forcingNSYdir(const ZEROPTV &location, const double &t) {
		double time = t;

		/// Get non dimensional numbers
		double Re = this->input_data_->Re;
		double We = 1.0;
		double Fr = this->input_data_->Fr; /// Froude number

		// Forcing in y
		double forcingY = 0.0;
		double x = location.x(), y = location.y(), z = location.z();

		if (input_data_->ifPurePNP){
			forcingY = ndcf_conv_ *
			           ((M_PI * M_PI * M_PI) * cos (y * M_PI) * pow (sin (x * M_PI * 2.0), 2.0) * pow (sin (y * M_PI), 3.0) *
			            pow (sin (t), 2.0) * 2.0 -
			            (M_PI * M_PI * M_PI) * cos (x * M_PI * 2.0) * pow (sin (x * M_PI), 2.0) * pow (sin (y * M_PI), 2.0) *
			            sin (y * M_PI * 2.0) * pow (sin (t), 2.0) * 2.0) + (ndcf_diff_ *
			                                                                ((M_PI * M_PI * M_PI) * pow (cos (y * M_PI), 2.0)
			                                                                 *
			                                                                 sin (x * M_PI * 2.0) * sin (t) * 2.0 -
			                                                                 (M_PI * M_PI * M_PI) * sin (x * M_PI * 2.0) *
			                                                                 pow (sin (y * M_PI), 2.0) * sin (t) * 6.0)) / Re
			           +
			           ndcf_pres_ * M_PI * cos (x * M_PI) * cos (y * M_PI) * sin (t) -
			           ndcf_time_ * M_PI * sin (x * M_PI * 2.0) * pow (sin (y * M_PI), 2.0) * cos (t);
		} else {
			double t2 = sin (t);
			double t3 = pi * x;
			double t4 = pi * y;
			double t5 = pi * pi * pi;
			double t6 = t2 * t2;
			double t7 = t3 * 2.0;
			double t8 = t4 * 2.0;
			double t9 = cos (t4);
			double t10 = sin (t4);
			double t11 = cos (t7);
			double t12 = cos (t8);
			double t13 = sin (t7);
			double t14 = sin (t8);
			double t15 = t10 * t10;
			forcingY = ndcf_conv_ * (t5 * t6 * t9 * (t10 * t10 * t10) * (t13 * t13) * 2.0
			                          - t5 * t6 * t11 * t14 * t15 * pow (sin (t3), 2.0) * 2.0)
			            - (ndcf_diff_ * (t2 * t5 * t13 * t15 * 6.0 - t2 * t5 * (t9 * t9) * t13 * 2.0)) / Re
			            + ndcf_pres_ * pi * t2 * t9 * cos (t3) - ndcf_time_ * pi * t13 * t15 * cos (t)
			            - ndcf_pnp_coupling_ * pi * t2 * t11 * t12 * (t2 * t11 * t14 - t2 * t12 * t13) * 2.0;
		}
		return forcingY;
	}

#pragma mark NS Integrands

    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {
        setIntegrandsFlags ();
        if (input_data_->timescheme == NSPNPInputData::THETA) {
            if (input_data_->ifLinearNS){
                Integrands_NS_Momentum_ThetaAe(fe, Ae, dt_, t_);
            } else {
                throw TALYException () << "Not implemented yet";
            }
        } else if (input_data_->timescheme == NSPNPInputData::BDF) {
            if (input_data_->ifLinearNS){
                Integrands_NS_Momentum_PSPG_BDF12Ae(fe, Ae, dt_, t_);
            } else {
                Integrands_NS_Momentum_BDF12Ae(fe, Ae, dt_, t_);
            }
        } else {
            throw TALYException () << "Time scheme not implemented. Should be either \"theta\" or \"bdf\"";
        }
    }

    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {

        setIntegrandsFlags ();
        if (input_data_->timescheme == NSPNPInputData::THETA) {
            if (input_data_->ifLinearNS) {
                Integrands_NS_Momentum_Thetabe(fe, be, dt_, t_);
            } else {
                throw TALYException () << "Not implemented yet";
            }
        } else if (input_data_->timescheme == NSPNPInputData::BDF) {
            if (input_data_->ifLinearNS) {
                Integrands_NS_Momentum_PSPG_BDF12be(fe, be, dt_, t_);
            } else {
                Integrands_NS_Momentum_BDF12be(fe, be, dt_, t_);
            }
        } else {
            throw TALYException () << "Time scheme not implemented. Should be either \"theta\" or \"bdf\"";
        }
    }










};
