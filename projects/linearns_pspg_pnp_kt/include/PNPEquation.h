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

class PNPEquation : public TALYFEMLIB::CEquation<NSPNPNodeData> {
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
    explicit PNPEquation(NSPNPInputData *inputData)
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
		ndcf_time_ = 1.0;
		ndcf_conv_ = 1.0;
		ndcf_diff_ = 1.0;
		ndcf_pres_ = 1.0;
		ndcf_pnp_coupling_ = 1.0;

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


#pragma mark Integrands-interface-functions PNP

	/**
	 * Returns the linear system (LHS part) of the integrands
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param Ae Empty Ae matrix to be filled in the method
	 * @param dt timestep
	 * @param t current time
	 * @param freeEnergy: freeEnergy object to calculate free energy on the fly
	 */
	void getIntegrandsPNPAe(const FEMElm &fe, ZeroMatrix<double> &Ae, double dt, double t) {
		setIntegrandsFlags ();
		if (input_data_->timescheme == NSPNPInputData::THETA) {
			throw TALYException () << "Not implemented yet";
		} else if (input_data_->timescheme == NSPNPInputData::BDF) {
			PNPIntegrands_BDF12Ae(fe, Ae, dt, t);
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
	void getIntegrandsPNPbe(const FEMElm &fe, ZEROARRAY<double> &be, double dt, double t) {
		setIntegrandsFlags ();
		if (input_data_->timescheme == NSPNPInputData::THETA) {
			throw TALYException () << "Not implemented yet";
		} else if (input_data_->timescheme == NSPNPInputData::BDF) {
			PNPIntegrands_BDF12be(fe, be, dt, t);
		} else {
			throw TALYException () << "Time scheme not implemented. Should be either \"theta\" or \"bdf\"";
		}
	}


	/** Function to return the integrands for boundary terms for Pressure Poisson (only contributions to RHS)
	 * @param fe FEMElm object for looping over elements and gauss points
	 * @param be Empty be matrix to be filled in the method
	 */
	void getIntegrands4sidePNPAe(const FEMElm &fe, int sideInd, ZeroMatrix<double> &Ae, const double dt, const double t) {
		setIntegrandsFlags ();
		PNPIntegrands4sideAe(fe, sideInd, Ae, dt, t);
	}

	/** Function to return the integrands for boundary terms for Pressure Poisson (only contributions to RHS)
 * @param fe FEMElm object for looping over elements and gauss points
 * @param be Empty be matrix to be filled in the method
 */
	void getIntegrands4sidePNPbe(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be, const double dt, const double t) {
		setIntegrandsFlags ();
		PNPIntegrands4sidebe(fe, sideInd, be, dt, t);
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
	void PNPIntegrands_BDF12Ae(const FEMElm &fe, ZeroMatrix<double> &Ae, const double dt, const double t){

		load_param(); /// set parameters
		const double detJxW = fe.detJxW();  /// (determinant of J) cross W
		const DENDRITE_UINT nsd = DIM;

		const DENDRITE_UINT ndof_ = NSPNPNodeData::PNP_DOF;
		DENDRITE_REAL dt_ = dt;
		DENDRITE_REAL t_ = t + dt_; ///< current time is elapsed + dt
		DENDRITE_UINT n_basis_functions = fe.nbf(); // # of basis functions

		/// load fluid velocity
		const double u = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
		const double v = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity

//        printf("phi=   %f , %f\n", this->p_data_->valueFEM(fe, NSPNPNodeData::POTENTIAL), this->p_data_->valueFEM(fe, NSPNPNodeData::VEL_X));

		double w = 0.0;
		#if (DIM == 3)
		w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
		#endif
		ZEROPTV U(u,v,w);  /// non-dimensional velocity
		/// For other than NS_PNP and NS_P_NP eq, set velocity to zero
		if (input_data_->ifPurePNP) {
			U = U * .0;   /// fluid velocity is zero for uncoupled equation
		}

		/// load values from previous time step and previous guess
		double Ci_g[nSpcs_];   /// Concentration guess for Newton-Rapson method
		double Ci_pre1[nSpcs_];   /// Concentration of previous time step (n1)
		double Ci_pre2[nSpcs_];   /// Concentration of previous time step (n0)
		for (int i = 0; i < nSpcs_; i++) {
			Ci_g[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i);
			Ci_pre1[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE1 + i);
			Ci_pre2[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE2 + i);
		}

		double test = this->p_data_->value2DerivativeFEM (fe, 0, 1, 1);
		/// load gradient of conc. and phi from previous time step
		ZEROPTV dCi_c[nSpcs_]; /// grad {Ci_g}
		std::vector<double> ddCi_c (nSpcs_); /// grad^2 {Ci_g}
		for (int sp = 0; sp < nSpcs_; sp++) {
			ddCi_c[sp] = 0.;
		}
		ZEROPTV dPhi_g;   /// grad {phi_g}
		double ddPhi_g=0.0;  /// grad.{dPhi_g}
		for (int k = 0; k < nsd; k++) {
			dPhi_g(k) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL, k);
			ddPhi_g += this->p_data_->value2DerivativeFEM(fe, NSPNPNodeData::POTENTIAL, k, k);
			for (int i = 0; i < nSpcs_; i++) {
				dCi_c[i](k) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i, k);
				ddCi_c[i]+=this->p_data_->value2DerivativeFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i, k, k);
			}
		}

		/// vec{u].grad{ck}
		std::vector<double> UdotdCk (nSpcs_);
		for (int sp = 0; sp < nSpcs_; sp++) {
			UdotdCk[sp] = 0.;
			for (int k = 0; k < nsd; k++) {
				UdotdCk[sp] += U[k] * dCi_c[sp] (k);
			}
		}

		/// dCdPhi
		std::vector<double> dCdPhi (nSpcs_);
		for (int sp = 0; sp < nSpcs_; sp++) {
			dCdPhi[sp] = 0.;
			for (int k = 0; k < nsd; k++) {
				dCdPhi[sp] += dCi_c[sp] (k) * dPhi_g (k);
			}
		}

		/// set forcing term for MMS
		double forcing[3];
		for (int i = 0; i < ndof_; i++)  /// forcing term is zero, if mms is not used
			forcing[i] = 0.0;
		if (input_data_->ifMMS) { /// if mms is being used
			forcing[0] = Force_Poisson (t_, fe.position (), E_);
			forcing[1] = PNPCpForce (t_, fe.position (), A_[0], B_[0], C_[0], D_[0]);
			forcing[2] = PNPCnForce (t_, fe.position (), A_[1], B_[1], C_[1], D_[1]);
		}

		/// get VMS parameter
		std::vector<double> tauK (nSpcs_);
		tauK = calc_tauK (fe, U, dPhi_g, ddPhi_g);

		////================================ Starting the Elemental assembly ============================================///
		for (int a = 0; a < n_basis_functions; a++) {  /// loop over Gaussian points a

			double dNaDotU = 0.; ///grad{w}.U
			double dNaDotPhi = 0.; ///grad{Na}.grad{phi}

			for (int dir = 0; dir < nsd; dir++) {
				dNaDotU+= fe.dN(a, dir) * U(dir);
				dNaDotPhi+=fe.dN(a, dir) * dPhi_g(dir);
			}

			for (int b = 0; b < n_basis_functions; b++) {  /// loop over Gaussian points b

				/// Mass with SUPG basis function N*. This zero when solving steady state eqn
				double Mstar = fe.N(a) * fe.N(b); /// mass matrix
				double K = .0;      /// Stiffness matrix
				double Conv = .0;  /// Convection matrix
				double MI = .0;     /// migration matrix

				/// for VMS Res terms
				double UdNb = .0; /// U.grad{Nb}
				double ddNb = .0; /// grad.grad{Nb}
				double dNbdPhi = .0; /// grad{Nb}.grad{Phi}
				std::vector<double> dCdNb(nSpcs_); /// grad{ci_0}.grad{Nb}
				for (int x : dCdNb) {
					dCdNb.push_back(0.0);
				}

				for (int dir = 0; dir < fe.nsd (); dir++) {
					K += fe.dN (a, dir) * fe.dN (b, dir);
					Conv += fe.N (a) * U (dir) * fe.dN (b, dir);
					MI += fe.dN (a, dir) * dPhi_g (dir) * fe.N (b);
					UdNb += U (dir) * fe.dN (b, dir);
					ddNb += fe.d2N (b, dir, dir);
					dNbdPhi += fe.dN (b, dir) * dPhi_g (dir);



					for (int sp = 0; sp < nSpcs_; sp++) {
						dCdNb[sp] += dCi_c[sp] (dir) * fe.dN (b, dir);
					}
				}

                //printf("E_ = %d\n", z[1]);

				/// looping over sub-matrix j X k, j corresponds to Potential and species, k corresponds to unknowns
				for (int j = 0; j < ndof_; j++)  {
					for (int k = 0; k < ndof_; k++)  {

						/**
						 * first row of sub matrix is for Poisson equation
						 * E_*( grad{w}, grad{ delta {phi} } ) - sigma ( w, zi * delta {ci} ) = -R + f
						 */
						if (j == 0) {
							if (k == 0) /// at sub matrix(0,0), 1st term of the Poisson equation is defined
								Ae(a * ndof_,b * ndof_) += E_ * K * detJxW;
							else  /// rest of terms in Poisson equation is defined here
								Ae(a * ndof_,b * ndof_ + k) += - z[k - 1] * Mstar * detJxW;
						}

							/**
							 * rest of rows are for NP equations for multiple species
							 *  Galerkin
							 *    A * coef_BDF[0] * 1/dt * ( w, d delta{Ci_n2} ) + B * ( w, vec{u}.grad{delta{ci}} ) + C * ( grad{w}, grad{delta{ci}} )
							 *    + D * ( grad{w}, zi*delta{ci}*grad{phi_k} ) +D * ( grad{w}, zi*Ci_k*grad{delta{phi}} )
							 *  VMS
							 *    - B * ( grad{w}*vec{u}, delta{ci'} ) + D * ( grad{w}*zi*grad{phi_k}, delta{ci'} )
							 *    + D * ( grad{w}*zi*ci'_k,  grad{delta{phi}} )
							 *  RHS
							 *    = -(w,R_np)
							 *    = -A*(w,dc_i/dt) - B*(w,vec{u}.grad{c_i,k}) - C*(grad{w},grad{c_i,k}) - D*(grad{w},zi*c_i,k * grad{phi,k})
							 *        + B*(grad{w},vec{u}*c'_i,k) - D*(grad{w},zi*c'_i,k*grad{phi_k})
							 *  where,
							 *    ci' = -tau * Res_k(delta{ci})
							 *    ci'_k = -tau * Res_k(ci_k)
							 */
						else { //checkcheck
							if (k == 0) {
								/** First column in matrix (Ae) corresponds to phi in vector (xe)
								 *    Galerkin
								 *        D * ( grad{w}, zi*Ci_g*grad{delta{phi}} )
								 */
								Ae(a * ndof_ + j, b * ndof_) += D_[j-1] * z[j - 1] * (Ci_g[j - 1] * K) * detJxW;

								/** VMS
								 *      D_ * ( grad{w}, zi*C'_k * grad{delta{phi}} )
								 *      where
								 *      C'_k = -tauK * Res(ci_k)
								 */
								Ae (a * ndof_ + j, b * ndof_) += -D_[j - 1] * z[j - 1] * K * tauK[j - 1]
								                                 /// Res(c_k)
								                                 /// A*1/dt*( 3/2*c_n+2 - 2*c_n+1 + 1/2*c_n
								                                 * (1. / dt * A_[j - 1] *
								                                 (bdf_c[0] * Ci_g[j - 1] + bdf_c[1] * Ci_pre1[j - 1]
								                                 + bdf_c[2] * Ci_pre2[j - 1])
								                                 /// B*vec{u}.grad{c_k} - C_*grad^2{c_k}
								                                 + B_[j - 1] * UdotdCk[j - 1] - C_[j - 1] * ddCi_c[j - 1]
								                                 /// -D*z_i*grad{c_k}.grad{phi_k} -D*z_i*c_k.grad^2{phi_k}
								                                 - D_[j - 1] * z[j - 1] * dCdPhi[j - 1]
								                                 - D_[j - 1] * z[j - 1] * Ci_g[j - 1] * ddPhi_g) * detJxW;
							}

								/** Diagonal terms in matrix (Ae) corresponds to ci in vector (xe)
								 *  Galerkin
								 *    A * coef_BDF[0] * 1/dt * ( w, d delta{Ci_n2} ) + B * ( w, vec{u}.grad{delta{ci}} ) + C * ( grad{w}, grad{delta{ci}} )
								 *    + D * ( grad{w}, zi*delta{ci}*grad{phi_g} )
								 *  VMS
								 *    - B * ( grad{w}.vec{u}, delta{ci'} ) + ( grad{w}*zi*grad.grad{phi}, ci' )
								 */
							else if (k == j)  {
								/// Galerkin
								Ae(a * ndof_ + j, b * ndof_ + k) += bdf_c[0]*(A_[j-1] * Mstar)/dt * detJxW;
								Ae(a * ndof_ + j, b * ndof_ + k) += (B_[j-1] *  Conv + C_[j-1] * K + D_[j-1] * z[j - 1] * MI) * detJxW;

								/** VMS convection
								 *      - B * ( grad{w}.vec{u}, delta{ci'} )
								 *  where delta{ci'} =
								 *          -tauK * delta{Res(ci)}
								 */
								Ae (a * ndof_ + j, b * ndof_ + k) +=
																											(B_[j - 1] * dNaDotU * tauK[j - 1]
								                                      /// delta{ Res(ci} }
								                                      * (bdf_c[0] * (A_[j - 1] * fe.N (b)) / dt
								                                      ///
								                                      + B_[j - 1] * UdNb   ///
								                                      // transient and convection
								                                      /// - C_ * grad^2{ delta{ci} } - D_*zi*grad{delta{ci}}.grad{phi_k}
								                                      - C_[j - 1] * ddNb - D_[j - 1] * z[j - 1] * dNbdPhi
								                                      /// -D_*zi*delta{ci}.grad^2{phi_k}
								                                      - D_[j - 1] * z[j - 1] * fe.N (b) * ddPhi_g) ) * detJxW;
								/// Add remaining terms in Res(ci) (contribution from phi)
								Ae (a * ndof_ + j, b * ndof_) +=
																									(B_[j - 1] * dNaDotU * tauK[j - 1]
								                                  /// -D_*zi*grad{ck}*grad{delta{phi}} - D_*zi*ck*grad^2{delta{phi}}
								                                  * (-D_[j - 1] * z[j - 1] * dCdNb[j - 1]
								                                  - D_[j - 1] * z[j - 1] * Ci_g[j - 1] * ddNb)) * detJxW;

								/// VMS 2 ( zi*grad{w}.grad{phi}, delta{ci'} )
								Ae (a * ndof_ + j, b * ndof_ + k) +=
																											-(D_[j - 1] * z[j - 1] * dNaDotPhi * tauK[j - 1]
								                                       /// delta{ Res(ci} }
								                                       * (bdf_c[0] * (A_[j - 1] * fe.N (b)) / dt
								                                       + B_[j - 1] * UdNb
								                                       /// transient and convection
								                                       /// - C_* grad^2{ delta{ci} } - D_*zi*grad{delta{ci}}.grad{phi_k}
								                                       - C_[j - 1] * ddNb - D_[j - 1] * z[j - 1] * dNbdPhi
								                                       /// -D_*zi*delta{ci}.grad^2{phi_k}
								                                       - D_[j - 1] * z[j - 1] * fe.N (b) * ddPhi_g) ) * detJxW;
								/// Add remaining terms in Res(ci) (contribution from phi)
								Ae (a * ndof_ + j, b * ndof_) += -(D_[j - 1] * z[j - 1] * dNaDotPhi * tauK[j - 1]
								                                   /// -D_*zi*grad{ck}*grad{delta{phi}} - D_*zi*ck*grad^2{delta{phi}}
								                                   * (-D_[j - 1] * z[j - 1] * dCdNb[j - 1]
								                                      - D_[j - 1] * z[j - 1] * Ci_g[j - 1] * ddNb)) * detJxW;
							}
						}
					}
				}
			}
			/// ==========================================Done with elemental assembly ====================================///
		}
	}
	/**
	 * References for VMS
	 * 1. Simplified VMS with charge neutrality assumption for NS-PNP
	 * - G. Bauer et al. / Comput. Methods Appl. Mech. Engrg. 223–224 (2012) 199–210
	 * 2. Whole VMS terms for PNP
	 * - Shadid, John Nicolas, et al. Simulation of neutron radiation damage in silicon semiconductor devices, 2007.
	 * 3. VMS for Navier-Stokes
	 * - Bazilevs et al. / Comput. Methods Appl. Mech. Engrg. 197 (2007) 173–201
	 */
	void PNPIntegrands_BDF12be(const FEMElm &fe, ZEROARRAY<double> &be, const double dt, const double t){

		load_param(); /// set parameters
		const double detJxW = fe.detJxW();  /// (determinant of J) cross W
		const DENDRITE_UINT nsd = DIM;

		const DENDRITE_UINT ndof_ = NSPNPNodeData::PNP_DOF;
		DENDRITE_REAL dt_ = dt;
		DENDRITE_REAL t_ = t + dt_; ///< current time is elapsed + dt
		DENDRITE_UINT n_basis_functions = fe.nbf(); // # of basis functions

		/// load fluid velocity
		const double u = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
		const double v = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
		double w = 0.0;
		#if (DIM == 3)
		w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
		#endif
		ZEROPTV U(u,v,w);  /// non-dimensional velocity
		/// For other than NS_PNP and NS_P_NP eq, set velocity to zero
		if (input_data_->ifPurePNP) {
			U = U * .0;   /// fluid velocity is zero for uncoupled equation
		}

		/// load values from previous time step and previous guess
		double Ci_g[nSpcs_];   /// Concentration guess for Newton-Rapson method
		double Ci_pre1[nSpcs_];   /// Concentration of previous time step (n1)
		double Ci_pre2[nSpcs_];   /// Concentration of previous time step (n0)
		for (int i = 0; i < nSpcs_; i++) {
			Ci_g[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i);
			Ci_pre1[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE1 + i);
			Ci_pre2[i] = this->p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1_PRE2 + i);
		}

		double test = this->p_data_->value2DerivativeFEM (fe, 0, 1, 1);
		/// load gradient of conc. and phi from previous time step
		ZEROPTV dCi_c[nSpcs_]; /// grad {Ci_g}
		std::vector<double> ddCi_c (nSpcs_); /// grad^2 {Ci_g}
		for (int sp = 0; sp < nSpcs_; sp++) {
			ddCi_c[sp] = 0.;
		}
		ZEROPTV dPhi_g;   /// grad {phi_g}
		double ddPhi_g=0.0;  /// grad.{dPhi_g}
		for (int k = 0; k < nsd; k++) {
			dPhi_g(k) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL, k);
			ddPhi_g += this->p_data_->value2DerivativeFEM(fe, NSPNPNodeData::POTENTIAL, k, k);
			for (int i = 0; i < nSpcs_; i++) {
				dCi_c[i](k) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i, k);
				ddCi_c[i]+=this->p_data_->value2DerivativeFEM(fe, NSPNPNodeData::CONCENTRATION_1 + i, k, k);
			}
		}

		/// vec{u].grad{ck}
		std::vector<double> UdotdCk (nSpcs_);
		for (int sp = 0; sp < nSpcs_; sp++) {
			UdotdCk[sp] = 0.;
			for (int k = 0; k < nsd; k++) {
				UdotdCk[sp] += U[k] * dCi_c[sp] (k);
			}
		}

		/// dCdPhi
		std::vector<double> dCdPhi (nSpcs_);
		for (int sp = 0; sp < nSpcs_; sp++) {
			dCdPhi[sp] = 0.;
			for (int k = 0; k < nsd; k++) {
				dCdPhi[sp] += dCi_c[sp] (k) * dPhi_g (k);
			}
		}

		/// set forcing term for MMS
		double forcing[3];
		for (int i = 0; i < ndof_; i++)  /// forcing term is zero, if mms is not used
			forcing[i] = 0.0;
		if (input_data_->ifMMS) { /// if mms is being used
			forcing[0] = Force_Poisson (t_, fe.position (), E_);
			forcing[1] = PNPCpForce (t_, fe.position (), A_[0], B_[0], C_[0], D_[0]);
			forcing[2] = PNPCnForce (t_, fe.position (), A_[1], B_[1], C_[1], D_[1]);
		}

		/// get VMS parameter
		std::vector<double> tauK (nSpcs_);
		tauK = calc_tauK (fe, U, dPhi_g, ddPhi_g);


		for (int a = 0; a < n_basis_functions; a++) {  /// loop over Gaussian points a

			double dNaDotU = 0.; ///grad{w}.U
			double dNaDotPhi = 0.; ///grad{Na}.grad{phi}

			for (int dir = 0; dir < nsd; dir++) {
				dNaDotU+= fe.dN(a, dir) * U(dir);
				dNaDotPhi+=fe.dN(a, dir) * dPhi_g(dir);
			}

			/// =====================================Beginning vector assembly=============================================///
			double udotGradCig[ndof_-1];       /// vec{u}.grad{Ci_g}
			double GradNadotGradCig[ndof_-1];  /// grad{Na}.grad{Ci_g}
			for (int species = 0; species < ndof_-1; species++)  {
				udotGradCig[species] = .0;
				GradNadotGradCig[species] = .0;
			}

			double GradNadotGradPhig = .0;      /// grad{Na}.grad{phi_g}

			for (int dir = 0; dir < nsd; dir++) {
				for (int species = 0; species < ndof_-1; species++)  {
					udotGradCig[species] += U(dir) * dCi_c[species](dir);
					GradNadotGradCig[species] += fe.dN(a,dir) * dCi_c[species](dir);
				}
				GradNadotGradPhig += fe.dN(a,dir) * dPhi_g(dir);
			}


			std::vector<double> dCg_dPhig(nSpcs_);   /// grad{cg}.grad{phi_g}
			for (int x:dCg_dPhig) {
				dCg_dPhig.push_back(0.0);
			}
			for (int dir = 0; dir < nsd; dir++) {
				for (int sp = 0; sp < nSpcs_; sp++)  {
					dCg_dPhig[sp] += dCi_c[sp](dir)*dPhi_g(dir);
				}
			}

			std::vector<double> Cg_ddPhi(nSpcs_);    /// cg * grad.grad{phi}
			for (int x:Cg_ddPhi) {
				Cg_ddPhi.push_back(0.0);
			}
			for (int sp = 0; sp < nSpcs_; sp++)  {
				Cg_ddPhi[sp] += Ci_g[sp]*ddPhi_g;
			}


			double rho_e = 0.; /// electric charge density
			DENDRITE_UINT numSpecies = 2;
			for (int species = 0; species < numSpecies; species++){
				rho_e += z[species] * Ci_g[species];
			}

			/// Residual of Poisson equation and forcing term
			be(a * ndof_) += (E_ * GradNadotGradPhig
			                  - fe.N(a) * rho_e
			                  - fe.N(a) * forcing[0]) * detJxW;

			/// Residual of Nernst-Planck equation and forcing term
			for (int j = 1; j < ndof_; j++)  {
				/// Galerkin
				be(a * ndof_ + j) += (A_[j - 1] * (bdf_c[0] * fe.N(a) * Ci_g[j-1]
				                                   + bdf_c[1] * fe.N(a) * Ci_pre1[j - 1]
				                                   + bdf_c[2] * fe.N(a) * Ci_pre2[j - 1]) / dt) * detJxW;
				be(a * ndof_ + j) += (B_[j-1]*fe.N(a)*udotGradCig[j-1]
				                      + C_[j-1]*GradNadotGradCig[j-1]
				                      + D_[j-1] * z[j - 1] * Ci_g[j - 1] * GradNadotGradPhig
				                      - fe.N(a)*forcing[j] ) * detJxW;

				/// VMS 1 - B * ( grad{w}.vec{u}, c'_{i,k} ) where ci' = -tauK * Res(ci)
				be(a * ndof_ + j) += (B_[j - 1] * dNaDotU * tauK[j - 1] *
				                      (A_[j - 1] * (1.0/dt) *
				                       (bdf_c[0] * Ci_g[j - 1] + bdf_c[1] * Ci_pre1[j - 1] + bdf_c[2] * Ci_pre2[j - 1])
				                       + B_[j - 1] * udotGradCig[j - 1] - C_[j - 1] * ddCi_c[j - 1]
				                       - D_[j - 1] * z[j - 1] * dCdPhi[j - 1]
				                       - D_[j - 1] * z[j - 1] * Cg_ddPhi[j - 1] )) * detJxW;

				/// VMS 2 ( z*grad{w}.grad{phi}, ci'_k )
				be(a * ndof_ + j) += -(D_[j - 1] * GradNadotGradPhig * z[j - 1] * tauK[j - 1]
				                     * ( A_[j - 1] * (1.0/dt) *
				                         (bdf_c[0] * Ci_g[j - 1] + bdf_c[1] * Ci_pre1[j - 1] + bdf_c[2] * Ci_pre2[j - 1])
				                         + B_[j - 1] * udotGradCig[j - 1] - C_[j - 1] * ddCi_c[j - 1]
				                         - D_[j - 1] * z[j - 1] * dCdPhi[j - 1]
				                         - D_[j - 1] * z[j - 1] * Cg_ddPhi[j - 1]) ) * detJxW;
			}
		}
		///=======================================Done assmebling elemental vector=======================================///
	}

	//todo
	void PNPIntegrands4sideAe(const FEMElm &fe, int sideInd, ZeroMatrix<double> &Ae, const double dt, const double t) {
		DENDRITE_REAL dt_ = dt;
		DENDRITE_REAL t_ = t + dt_;
		///==============================================================================================================///
		///============================================ Weak BC (Phi) ===================================================///
		///==============================================================================================================///
		if (input_data_->numWeakBCPhi != 0) {  /// loop for weak BC
			for (int num = 0; num < input_data_->numWeakBCPhi; num++) {  /// loop over weakBC index
				if (input_data_->weakBCPhiIdx[num] == sideInd) { /// if current node is on weakBC index
					//TODO set following variable as class members

					double detJxW = fe.detJxW ();
					int nsd = DIM; /// Number of space dimensions
					int nbf = fe.nbf (); /// Number of basis function
					int ndof = NSPNPNodeData::PNP_DOF;

					const int numSpecies = 2;  /// number of species
					double coeff_E = input_data_->E_; /// non-dimensional Debye length (lambda / L)

					const ZEROPTV &normal = fe.surface ()->normal ();   /// normal to surface

					double gamma = input_data_->gamma_NP;  /// coeff for adjoint term. either 1 or -1
					double Cb = input_data_->Cb_Poisson;   /// coeff of penalty term

					double h_b = BdElmSize (fe, normal);  /// boundary element size

					/// Define expression variable to link functions from 'config.txt' with actual value
					NSPNPExpr Expr;
					Expr.set_variable (t_, fe);
					double phiBC = input_data_->weakBCPhi[num].value (); /// phi boundary condition


					double c0[numSpecies]; /// previous concentration guess
					ZEROPTV gradC0[numSpecies]; /// grad{c0}
					double gradCiDotN[numSpecies]; /// grad{c0}.vec{n}
					for (int species = 0; species < numSpecies; species++) {
						c0[species] = this->p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species);
						gradCiDotN[species] = 0.0;
						for (int i = 0; i < nsd; i++) {
							gradC0[species] (i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::CONCENTRATION_1 + species, i);
							gradCiDotN[species] += gradC0[species] (i) * normal (i);
						}
					}

					/// load fluid velocity
					const double u = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
					const double v = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
					double w = 0.0;
					#if (DIM == 3)
						w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
					#endif
					ZEROPTV U(u,v,w);  /// non-dimensional velocity


					/// load gradient values from previous time step
					ZEROPTV dPhi_guess;
					for (int i = 0; i < nsd; i++)
						dPhi_guess(i) = this->p_data_->valueDerivativeFEM(fe, NSPNPNodeData::POTENTIAL, i);
					double phi0 = this->p_data_->valueFEM(fe, NSPNPNodeData::POTENTIAL);

					/// SUPG calculation
					ZEROPTV conv = U - dPhi_guess;
					TezduyarUpwindFE adv;
					adv.calcSUPGWeightingFunction(fe, conv, 1.0); /// calculate coefficient tau


					///=============================== Beginning Elemental assembly ===========================================///
					for (int a = 0; a < fe.nbf (); a++) {


						double delta = adv.SUPG (a);
						double NaStar = fe.N (a) + delta;
						double gradNaDotN = .0;
						for (int k = 0; k < fe.nsd (); k++)
							gradNaDotN += fe.dN (a, k) * normal (k);

						double gradPhi0 = .0; /// grad {phi0}
						double gradPhi0DotN = 0.0;  /// grad{phi0} . vec{n}
						for (int i = 0; i < nsd; i++) {
							gradPhi0 += p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
							gradPhi0DotN += dPhi_guess (i) * normal (i);
						}

						for (int b = 0; b < fe.nbf (); b++) {

							double gradNbDotN = .0;
							for (int k = 0; k < nsd; k++){
								gradNbDotN += fe.dN (b, k) * normal (k);
							}

							/// looping over sub-matrix j X k
							for (int j = 0; j < ndof; j++) {
								for (int k = 0; k < ndof; k++) {
									/// first row of sub-matrix is for Poisson equation
									if (j == 0) { /// for Poisson equation
										if (k == 0) { /// BC terms for Poisson equation
											/// consistency terms
											Ae (a * ndof, b * ndof) += -coeff_E * NaStar * gradNbDotN * detJxW;
											/// weak bc
											Ae (a * ndof, b * ndof) += -gamma * coeff_E * gradNaDotN * fe.N (b) * detJxW;
											/// penalty term
											Ae (a * ndof, b * ndof) += fabs (coeff_E) * NaStar * Cb / h_b * fe.N (b) * detJxW;
										}
									} else {  /// for Nernst-Planck equation
										double zj = z[j - 1];
										double cj = p_data_->valueFEM(fe, NSPNPNodeData::CONCENTRATION_1 + (j - 1));
										if (k == 0) { /// BC terms
											/*Ae(a*ndof+j, b*ndof) += - coeff_D[j-1] * NaStar * zj * cj * gradNbDotN * detJxW; /// consistency term
											*//*Ae(a*ndof+j, b*ndof) += - gamma * coeff_D[j-1] * zj * (NaStar * gradCiDotN[j-1]
                                                                   + cj * gradNaDotN) * fe.N(b) * fe.detJxW();*//* /// weak bc
                      Ae(a*ndof+j, b*ndof) += - gamma *  coeff_D[j-1]* zj * cj * gradNaDotN * fe.N(b) * detJxW; /// weak bc
                      Ae(a*ndof+j, b*ndof) += NaStar * Cb * fabs(zj*cj*coeff_D[j-1])/h_b * fe.N(b) * detJxW;*/  /// penalty term
										} else if (k == j) { /// diagonal terms, derivative with respect to each species
											/// consistency term
											//Ae(a*n_dof_+j, b*n_dof_+k) += - NaStar * gradNbDotN * fe.detJxW();
											//Ae(a*n_dof_+j, b*n_dof_+k) += - fe.N(a) * zj * gradPhi0DotN * fe.N(b) * fe.detJxW();
											/// weak bc
											//Ae(a*n_dof_+j, b*n_dof_+k) += - gamma * gradNaDotN * fe.N(b) * fe.detJxW();
											//Ae(a*n_dof_+j, b*n_dof_+k) += - gamma * zj * gradPhiGuess * gradNaDotN * fe.N(b) * fe.detJxW();
											/// penalty term
											//Ae(a*n_dof_+j, b*n_dof_+k) += NaStar * Cnp/h * fe.N(b) * fe.detJxW();
											//Ae(a*n_dof_+j, b*n_dof_+k) += NaStar * Cnp1 * fabs(zj*gradPhiGuess)/h * fe.N(b) * fe.detJxW();
										}
									}
								}
							}
						}
					}
				}
			}
		}


		///==============================================================================================================///
		///============================================ Weak BC (C+) ====================================================///
		///==============================================================================================================///
		if (input_data_->numWeakBCC1 != 0) {  /// loop for weak BC
			for (int num = 0; num < input_data_->numWeakBCC1; num++) {  /// loop over weakBC index
				if (input_data_->weakBCC1Idx[num] == sideInd) { /// if current node is on weakBC index
					//TODO set following variable as class members

					DENDRITE_REAL detJxW = fe.detJxW ();
					const DENDRITE_UINT nsd = DIM; /// Number of space dimensions
					const DENDRITE_UINT nbf = fe.nbf (); /// Number of basis function
					const DENDRITE_UINT ndof = NSPNPNodeData::PNP_DOF;

					const DENDRITE_UINT numSpecies = 2;  /// number of species
					std::vector<DENDRITE_REAL> coeff_C = input_data_->C_;  /// coefficient of diffusion term
					std::vector<DENDRITE_REAL> coeff_D = input_data_->D_;  /// coefficient of migration term
					//double coeff_E = input_data__->E_; /// non-dimensional Debye length (lambda / L)

					const ZEROPTV &normal = fe.surface ()->normal ();   /// normal to surface

					DENDRITE_REAL gamma = input_data_->gamma_NP;  /// coeff for adjoint term. either 1 or -1
					DENDRITE_REAL Cb = input_data_->Cb_NP;   /// coeff of penalty term

					DENDRITE_REAL h_b = BdElmSize (fe, normal);  /// boundary element size

					/// Define expression variable to link functions from 'config.txt' with actual value
					NSPNPExpr Expr;
					Expr.set_variable (t_, fe);
					//double phiBC = input_data_->weakBCPhi[num].value(); /// phi boundary condition
					DENDRITE_REAL C1BC = input_data_->weakBCC1[num].value ();

					DENDRITE_REAL c0[numSpecies]; /// previous concentration guess
					ZEROPTV gradC0[numSpecies]; /// grad{c0}
					DENDRITE_REAL gradCiDotN[numSpecies]; /// grad{c0}.vec{n}
					for (int species = 0; species < numSpecies; species++) {
						c0[species] = this->p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species);
						gradCiDotN[species] = 0.0;
						for (int i = 0; i < nsd; i++) {
							gradC0[species] (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species, i);
							gradCiDotN[species] += gradC0[species] (i) * normal (i);
						}
					}

					/// load fluid velocity
					const DENDRITE_REAL u = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
					const DENDRITE_REAL v = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
					DENDRITE_REAL w = 0.0;
					#if (DIM == 3)
					w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
					#endif
					ZEROPTV U (u, v, w);  /// non-dimensional velocity

					/// load gradient values from previous time step
					ZEROPTV dPhi_guess;
					for (int i = 0; i < nsd; i++) {
						dPhi_guess (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
					}
					double phi0 = this->p_data_->valueFEM (fe, NSPNPNodeData::POTENTIAL);

					/// SUPG calculation
					ZEROPTV conv = U - dPhi_guess;
					TezduyarUpwindFE adv;
					adv.calcSUPGWeightingFunction (fe, conv, 1.0); /// calculate coefficient tau

					for (int a = 0; a < fe.nbf (); a++) {

						double delta = adv.SUPG (a);
						double NaStar = fe.N (a) + delta;
						double gradNaDotN = .0;
						for (int k = 0; k < nsd; k++) {
							gradNaDotN += fe.dN (a, k) * normal (k);
						}
						double gradPhi0 = .0; /// grad {phi0}
						double gradPhi0DotN = 0.0;  /// grad{phi0} . vec{n}
						for (int i = 0; i < nsd; i++) {
							gradPhi0 += p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
							gradPhi0DotN += dPhi_guess (i) * normal (i);
						}

						for (int b = 0; b < fe.nbf (); b++) {
							double gradNbDotN = .0;
							for (int k = 0; k < nsd; k++) {
								gradNbDotN += fe.dN (b, k) * normal (k);
							}

							/// looping over sub-matrix j X k
							for (int j = 0; j < ndof; j++) {
								for (int k = 0; k < ndof; k++) {
									/// first row of sub-matrix is for Poisson equation
									if (j == 0) { /// for Poisson equation
										if (k == 0) { /// BC terms for Poisson equation
											/* /// consistency terms
											 Ae(a * ndof, b * ndof) += - coeff_E * NaStar * gradNbDotN * detJxW;
											 /// weak bc
											 Ae(a * ndof, b * ndof) += - gamma * coeff_E * gradNaDotN * fe.N(b) * detJxW;
											 /// penalty term
											 Ae(a * ndof, b * ndof) += fabs(coeff_E) * NaStar * Cb / h_b * fe.N(b) * detJxW;*/
										}
									} else {  /// for Nernst-Planck equation
										double zj = z[j - 1];
										double cj = p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + (j - 1));
										if (k == 0) { /// BC terms
											/// consistency term
											/// (w,grad{c}.n)
											Ae (a * ndof + j, b * ndof) += -coeff_D[j - 1] * NaStar * zj * cj * gradNbDotN * detJxW;
											/*Ae(a*ndof+j, b*ndof) += - gamma * coeff_D[j-1] * zj * (NaStar * gradCiDotN[j-1]
																																	 + cj * gradNaDotN) * fe.N(b) * fe.detJxW(); /// weak bc
											Ae(a*ndof+j, b*ndof) += - gamma *  coeff_D[j-1]* zj * cj * gradNaDotN * fe.N(b) * detJxW; /// weak bc
											Ae(a*ndof+j, b*ndof) += NaStar * Cb * fabs(zj*cj*coeff_D[j-1])/h_b * fe.N(b) * detJxW;*/  /// penalty term
										} else if (k == j) { /// diagonal terms, derivative with respect to each species
											/// consistency term
											/// C_ * (w,grad{del c}.n)
											Ae (a * ndof + j, b * ndof + k) += -coeff_C[j - 1] * NaStar * gradNbDotN * fe.detJxW ();
											/// D_ * (w,z_i*del(c_i)*grad{phi}.n)
											Ae (a * ndof + j, b * ndof + k) +=
													-coeff_D[j - 1] * NaStar * zj * gradPhi0DotN * fe.N (b) * fe.detJxW ();

											/// weak bc
											/// C_ * (grad{w}.n, del c)
											Ae (a * ndof + j, b * ndof + k) +=
													-coeff_C[j - 1] * gamma * gradNaDotN * fe.N (b) * fe.detJxW ();
											/// D_ * (z_i*phi*grad{w}.n, del(c_i))
											Ae (a * ndof + j, b * ndof + k) +=
													-coeff_D[j - 1] * gamma * zj * gradPhi0 * gradNaDotN * fe.N (b) * fe.detJxW ();

											/// penalty term
											/// C_ * (C/h*w, del(c))
											Ae (a * ndof + j, b * ndof + k) +=
													fabs (coeff_C[j - 1]) * NaStar * Cb / h_b * fe.N (b) * fe.detJxW ();
											//Ae(a*n_dof_+j, b*n_dof_+k) += NaStar * Cnp1 * fabs(zj*gradPhiGuess)/h * fe.N(b) * fe.detJxW();
										}
									}
								}
							}
						}
					}
				}
			}
		}

	}

	//todo
	void PNPIntegrands4sidebe(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be, const double dt, const double t) {
		DENDRITE_REAL dt_ = dt;
		DENDRITE_REAL t_ = t + dt_;
		///==============================================================================================================///
		///============================================ Weak BC (Phi) ===================================================///
		///==============================================================================================================///
		if (input_data_->numWeakBCPhi != 0) {  /// loop for weak BC
			for (int num = 0; num < input_data_->numWeakBCPhi; num++) {  /// loop over weakBC index
				if (input_data_->weakBCPhiIdx[num] == sideInd) { /// if current node is on weakBC index
					//TODO set following variable as class members

					double detJxW = fe.detJxW ();
					int nsd = DIM; /// Number of space dimensions
					int nbf = fe.nbf (); /// Number of basis function
					int ndof = NSPNPNodeData::PNP_DOF;

					const int numSpecies = 2;  /// number of species
					double coeff_E = input_data_->E_; /// non-dimensional Debye length (lambda / L)

					const ZEROPTV &normal = fe.surface ()->normal ();   /// normal to surface

					double gamma = input_data_->gamma_NP;  /// coeff for adjoint term. either 1 or -1
					double Cb = input_data_->Cb_Poisson;   /// coeff of penalty term

					double h_b = BdElmSize (fe, normal);  /// boundary element size

					/// Define expression variable to link functions from 'config.txt' with actual value
					NSPNPExpr Expr;
					Expr.set_variable (t_, fe);
					double phiBC = input_data_->weakBCPhi[num].value (); /// phi boundary condition


					double c0[numSpecies]; /// previous concentration guess
					ZEROPTV gradC0[numSpecies]; /// grad{c0}
					double gradCiDotN[numSpecies]; /// grad{c0}.vec{n}
					for (int species = 0; species < numSpecies; species++) {
						c0[species] = this->p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species);
						gradCiDotN[species] = 0.0;
						for (int i = 0; i < nsd; i++) {
							gradC0[species] (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species, i);
							gradCiDotN[species] += gradC0[species] (i) * normal (i);
						}
					}

					/// load fluid velocity
					const double u = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
					const double v = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
					double w = 0.0;
					#if (DIM == 3)
					w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
					#endif
					ZEROPTV U (u, v, w);  /// non-dimensional velocity


					/// load gradient values from previous time step
					ZEROPTV dPhi_guess;
					for (int i = 0; i < nsd; i++)
						dPhi_guess (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
					double phi0 = this->p_data_->valueFEM (fe, NSPNPNodeData::POTENTIAL);

					/// SUPG calculation
					ZEROPTV conv = U - dPhi_guess;
					TezduyarUpwindFE adv;
					adv.calcSUPGWeightingFunction (fe, conv, 1.0); /// calculate coefficient tau


					///=============================== Beginning Elemental assembly ===========================================///
					for (int a = 0; a < fe.nbf (); a++) {
						double delta = adv.SUPG (a);
						double NaStar = fe.N (a) + delta;
						double gradNaDotN = .0;
						for (int k = 0; k < fe.nsd (); k++)
							gradNaDotN += fe.dN (a, k) * normal (k);

						double gradPhi0 = .0; /// grad {phi0}
						double gradPhi0DotN = 0.0;  /// grad{phi0} . vec{n}
						for (int i = 0; i < nsd; i++) {
							gradPhi0 += p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
							gradPhi0DotN += dPhi_guess (i) * normal (i);
						}


						/// looping over sub-matrix j X k
						for (int j = 0; j < ndof; j++) {
							///NOTE: negative sign is already taken in SNES context
							if (j == 0) { /// BC terms in Residual of Poisson eq.
								be (a * ndof) += -coeff_E * NaStar * gradPhi0DotN * detJxW; /// consistency
								be (a * ndof) += -coeff_E * gamma * gradNaDotN * (phi0 - phiBC) * detJxW;  /// weak BCs
								be (a * ndof) += fabs (coeff_E) * NaStar * Cb / h_b * (phi0 - phiBC) * detJxW; /// penalty term

							} else {  /// BC terms in Residual of PNP eq.
								double zj = z[j - 1];
								double cj = p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + (j - 1));
								/// consistency term
								//be(a*n_dof_+j) += - fe.N(a) * gradCiDotN[j-1] * fe.detJxW();
								/*be(a*ndof+j) += - coeff_D[j-1] * fe.N(a) * zj*cj*gradPhi0DotN * fe.detJxW();
								/// weak bc
								//be(a*n_dof_+j) += - gamma * gradNaDotN * (cj - ciBC) * cj * fe.detJxW();
								//be(a*n_dof_+j) += - gamma * zj * fe.N(a) * gradPhi0DotN * (cj - ciBC) * fe.detJxW();
								*//*be(a*ndof+j) += - gamma * coeff_D[j-1] * zj * (NaStar * gradCiDotN[j-1]
                                                     + cj * gradNaDotN) * (phi0 - phiBC) * fe.detJxW();*//* /// weak bc
                  be(a*ndof+j) += - gamma *  coeff_D[j-1] * zj * cj * gradNaDotN * (phi0 - phiBC) * detJxW; /// weak bc
                  /// penalty term
                  //be(a*n_dof_+j) += NaStar * Cnp/h * (cj - ciBC) * fe.detJxW();
                  //be(a*n_dof_+j) += NaStar * Cnp1/h * fabs(zj*gradPhi0DotN)/h * (cj-ciBC) * fe.detJxW();
                  be(a*ndof+j) +=  NaStar * Cb * fabs(zj*cj*coeff_D[j-1])/h_b * (phi0 - phiBC) * detJxW;*/
							}
						}
					}
				}
			}
		}


		///==============================================================================================================///
		///============================================ Weak BC (C+) ====================================================///
		///==============================================================================================================///
		if (input_data_->numWeakBCC1 != 0) {  /// loop for weak BC
			for (int num = 0; num < input_data_->numWeakBCC1; num++) {  /// loop over weakBC index
				if (input_data_->weakBCC1Idx[num] == sideInd) { /// if current node is on weakBC index
					//TODO set following variable as class members

					DENDRITE_REAL detJxW = fe.detJxW ();
					const DENDRITE_UINT nsd = DIM; /// Number of space dimensions
					const DENDRITE_UINT nbf = fe.nbf (); /// Number of basis function
					const DENDRITE_UINT ndof = NSPNPNodeData::PNP_DOF;

					const DENDRITE_UINT numSpecies = 2;  /// number of species
					std::vector<DENDRITE_REAL> coeff_C = input_data_->C_;  /// coefficient of diffusion term
					std::vector<DENDRITE_REAL> coeff_D = input_data_->D_;  /// coefficient of migration term
					//double coeff_E = input_data__->E_; /// non-dimensional Debye length (lambda / L)

					const ZEROPTV &normal = fe.surface ()->normal ();   /// normal to surface

					DENDRITE_REAL gamma = input_data_->gamma_NP;  /// coeff for adjoint term. either 1 or -1
					DENDRITE_REAL Cb = input_data_->Cb_NP;   /// coeff of penalty term

					DENDRITE_REAL h_b = BdElmSize (fe, normal);  /// boundary element size

					/// Define expression variable to link functions from 'config.txt' with actual value
					NSPNPExpr Expr;
					Expr.set_variable (t_, fe);
					//double phiBC = input_data_->weakBCPhi[num].value(); /// phi boundary condition
					DENDRITE_REAL C1BC = input_data_->weakBCC1[num].value ();

					DENDRITE_REAL c0[numSpecies]; /// previous concentration guess
					ZEROPTV gradC0[numSpecies]; /// grad{c0}
					DENDRITE_REAL gradCiDotN[numSpecies]; /// grad{c0}.vec{n}
					for (int species = 0; species < numSpecies; species++) {
						c0[species] = this->p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species);
						gradCiDotN[species] = 0.0;
						for (int i = 0; i < nsd; i++) {
							gradC0[species] (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species, i);
							gradCiDotN[species] += gradC0[species] (i) * normal (i);
						}
					}

					/// load fluid velocity
					const DENDRITE_REAL u = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
					const DENDRITE_REAL v = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
					DENDRITE_REAL w = 0.0;
					#if (DIM == 3)
					w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
					#endif
					ZEROPTV U (u, v, w);  /// non-dimensional velocity

					/// load gradient values from previous time step
					ZEROPTV dPhi_guess;
					for (int i = 0; i < nsd; i++) {
						dPhi_guess (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
					}
					double phi0 = this->p_data_->valueFEM (fe, NSPNPNodeData::POTENTIAL);

					/// SUPG calculation
					ZEROPTV conv = U - dPhi_guess;
					TezduyarUpwindFE adv;
					adv.calcSUPGWeightingFunction (fe, conv, 1.0); /// calculate coefficient tau

					for (int a = 0; a < fe.nbf (); a++) {

						double delta = adv.SUPG (a);
						double NaStar = fe.N (a) + delta;
						double gradNaDotN = .0;
						for (int k = 0; k < nsd; k++) {
							gradNaDotN += fe.dN (a, k) * normal (k);
						}
						double gradPhi0 = .0; /// grad {phi0}
						double gradPhi0DotN = 0.0;  /// grad{phi0} . vec{n}
						for (int i = 0; i < nsd; i++) {
							gradPhi0 += p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
							gradPhi0DotN += dPhi_guess (i) * normal (i);
						}

						/// looping over sub-matrix j X k
						for (int j = 0; j < ndof; j++) {
							/**
							 * NOTE: negative sign is already taken in SNES context.
							 */
							if (j == 0) { /// BC terms in Residual of Poisson eq.
								/*be(a*ndof) += - coeff_E * NaStar * gradPhi0DotN * detJxW; /// consistency
								be(a*ndof) += - coeff_E * gamma * gradNaDotN * (phi0 - phiBC) * detJxW;  /// weak BCs
								be(a*ndof) +=  fabs(coeff_E) * NaStar * Cb/h_b * (phi0 - phiBC) * detJxW;*/ /// penalty term

							} else {  /// BC terms in Residual of PNP eq.
								double zj = z[j - 1];
								double cj = p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + (j - 1));
								/// consistency term
								/// C_ * (w,grad{c0}.n)
								be (a * ndof + j) += -coeff_C[j - 1] * NaStar * gradCiDotN[j - 1] * fe.detJxW ();
								/// D_ * (w,z_i*c0*grad{phi}.n)
								be (a * ndof + j) += -coeff_D[j - 1] * NaStar * zj * cj * gradPhi0DotN * fe.detJxW ();

								/// weak bc
								/// C_ * (grad{w}.n, c0 - g)
								be (a * ndof + j) += -coeff_C[j - 1] * gamma * gradNaDotN * (cj - C1BC) * fe.detJxW ();
								/// D_ * (z_i*phi*grad{w}.n, del(c_i))
								be (a * ndof + j) +=
										-coeff_D[j - 1] * gamma * zj * NaStar * gradPhi0DotN * (cj - C1BC) * fe.detJxW ();
								/*be(a*ndof+j) += - gamma * coeff_D[j-1] * zj * (NaStar * gradCiDotN[j-1]
																									 + cj * gradNaDotN) * (phi0 - phiBC) * fe.detJxW();*/ /// weak bc
								//be(a*ndof+j) += - gamma *  coeff_D[j-1] * zj * cj * gradNaDotN * (phi0 - phiBC) * detJxW; /// weak bc
								/// penalty term
								be (a * ndof + j) += fabs (coeff_C[j - 1]) * NaStar * Cb / h_b * (cj - C1BC) * fe.detJxW ();
								//be(a*n_dof_+j) += NaStar * Cnp1/h * fabs(zj*gradPhi0DotN)/h * (cj-ciBC) * fe.detJxW();
								//be(a*ndof+j) +=  NaStar * Cb * fabs(zj*cj*coeff_D[j-1])/h_b * (phi0 - phiBC) * detJxW;
							}
						}
					}

				}
			}
		}


		///==============================================================================================================///
		///=================================Charged Wall Boundary Condition=============================================////
		///==============================================================================================================///
		if (input_data_->numWCBC != 0) {  /// loop for Neuman BC
			for (int num = 0; num < input_data_->numWCBC; num++) {  /// loop over NBC index
				if (input_data_->WCBCIdx[num] == sideInd) { /// if current node is on NBC index


					DENDRITE_REAL detJxW = fe.detJxW ();
					const DENDRITE_UINT nsd = DIM; /// Number of space dimensions
					const DENDRITE_UINT nbf = fe.nbf (); /// Number of basis function
					const DENDRITE_UINT ndof = NSPNPNodeData::PNP_DOF;

					/// preparation for expression
					Expression::global_symbol_table ().set_variable ("t", t);

					double E = input_data_->E_;

					const int num_species_ = 2;   /// number of species
					double NBCValue = input_data_->WCBC[num].value ();                   /// [V/m]

					double Ci_guess[num_species_];
					for (int species = 0; species < num_species_; species++) {
						Ci_guess[species] = this->p_data_->valueFEM (fe, NSPNPNodeData::CONCENTRATION_1 + species);
					}
					const ZEROPTV &normal = fe.surface ()->normal ();   /// normal to surfacev

					/// load fluid velocity
					const DENDRITE_REAL u = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_X);  /// non-dimensional x-velcosity
					const DENDRITE_REAL v = this->p_data_->valueFEM (fe, NSPNPNodeData::VEL_Y);  /// non-dimensional y-velcosity
					DENDRITE_REAL w = 0.0;
					#if (DIM == 3)
					w = this->p_data_->valueFEM(fe,NSPNPNodeData::VEL_Z);  /// non-dimensional z-velcosity
					#endif
					ZEROPTV U (u, v, w);  /// non-dimensional velocity

					/// load gradient values from previous time step
					ZEROPTV dPhi_guess;
					for (int i = 0; i < fe.nsd (); i++) {
						dPhi_guess (i) = this->p_data_->valueDerivativeFEM (fe, NSPNPNodeData::POTENTIAL, i);
					}

					/// SUPG calculation
					ZEROPTV conv = U - dPhi_guess;
					TezduyarUpwindFE adv;
					adv.calcSUPGWeightingFunction (fe, conv, 1.0);   /// calculate coefficient tau

					/// loop over basis function
					for (int a = 0; a < fe.nbf (); a++) {
						double delta = adv.SUPG (a);
						double NaStar = fe.N (a) + delta;
						be (a * ndof) -= E * NaStar * NBCValue * fe.detJxW ();
					}
				}
			}

		}
	}

#pragma mark PNP helper functions
	double BdElmSize (const FEMElm &fe, const ZEROPTV &normal) {
		DENDRITE_UINT nsd = DIM;
		/// Get the cofactor matrix from FEELm class
		ZeroMatrix<double> ksiX;
		ksiX.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				ksiX(i, j) = fe.cof(j, i) / fe.volume_jacc();
			}
		}

		/// G_{ij} in equation 65 in Bazilevs et al. (2007)
		ZeroMatrix<double> Ge;
		Ge.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				Ge(i, j) = 0.0;
				for (int k = 0; k < nsd; k++)
					Ge(i, j) += ksiX(k, i) * ksiX(k, j);
			}
		}

		double h_b = 0.0;
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				h_b += normal(i) * Ge(i, j) * normal(j);
			}
		}

		h_b = 2. / sqrt(h_b);
		return h_b;
	}

	void load_param(){  /// load parameters
		A_ = this->input_data_->A_;
		B_ = this->input_data_->B_;
		C_ = this->input_data_->C_;
		D_ = this->input_data_->D_;
		E_ = this->input_data_->E_;
	}

	std::vector<double> calc_tauK (const FEMElm& fe, ZEROPTV U, ZEROPTV dPhi_g, double ddPhi_g){
		DENDRITE_UINT nsd = DIM;
		DENDRITE_REAL Ci_f_pnp = input_data_->Ci_f_pnp;
		DENDRITE_REAL dt_ = input_data_->timeStep;

		/// Get the cofactor matrix from FEELm class
		ZeroMatrix<double> ksiX;
		ksiX.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				ksiX(i, j) = fe.cof(j, i) / fe.jacc();
			}
		}
		// The formulation of the following three terms are the same as in Bazilevs et al. (2007) for NS

		/// the same G_{ij} as in equation 65 in Bazilevs et al. (2007)
		ZeroMatrix<double> Ge;
		Ge.redim(nsd, nsd);
		for (int i = 0; i < nsd; i++) {
			for (int j = 0; j < nsd; j++) {
				Ge(i, j) = 0.0;
				for (int k = 0; k < nsd; k++)
					Ge(i, j) += ksiX(k, i) * ksiX(k, j);
			}
		}

		/// eq.67 in Bazilevs et al. (2007), eq. 28 in Bauer et al.
		/*ZEROARRAY<double> u_Gu;
		u_Gu.redim(nSpcs);*/
		std::vector<double> u_Gu(nSpcs_);
		for(int k = 0; k < nSpcs_; k++){
			ZEROPTV Vk;
			Vk = U*B_[k] - dPhi_g * D_[k] * z[k];
			for (int i = 0; i < nsd; i++) {
				for (int j = 0; j < nsd; j++) {
					u_Gu[k] += Vk(i) * Ge(i, j) * Vk(j);
				}
			}
		}

		/// the same term as the last term in Equation 63 in Bazilevs et al. (2007)
		std::vector<double> G_G_u(nSpcs_);
		for(int k = 0; k < nSpcs_; k++){
			for (int i = 0; i < nsd; i++) {
				for (int j = 0; j < nsd; j++) {
					G_G_u[k] += Ci_f_pnp * C_[k] * C_[k] * Ge(i, j) * Ge(i, j);
				}
			}
		}

		std::vector<double> sig_mk(nSpcs_);
		for (int k = 0; k < nSpcs_; k++) {
			sig_mk[k] = D_[k] * z[k] * ddPhi_g;
		}

		//equaiton 28 in ref paper by Bauer et al.
		std::vector<double> tauM(nSpcs_);
		double tran; // transient term in tau
		double denom; // denominator in tau
		for(int k = 0; k < nSpcs_; k++){
			tran = 2.0*A_[k]/dt_ + fabs(sig_mk[k]);
			denom = tran*tran + u_Gu[k] + G_G_u[k];
			tauM[k] = 1./sqrt(denom);
		}
		return tauM;
	}


	#pragma mark Manufactured solutions helper PNP

	double Force_Poisson(double t, const ZEROPTV &pt, double E_) const {

		double x = pt.x();
		double y = pt.y();
		double force;

		/// forcing term for Potential in Poisson equation
		//force = - z[0]*cos(2.*pi*t)*cos(2.*pi*x)*sin(2.*pi*y) - z[1]*cos(2.*pi*t)*cos(2.*pi*y)*sin(2.*pi*x)
		//        - 8.*E_*pi*pi*cos(2.*pi*t)*cos(2.*pi*x)*sin(2.*pi*y);
		/*force = cos(2.*pi*t)*cos(2.*pi*y)*sin(2.*pi*x) - cos(2.*pi*t)*cos(2.*pi*x)*sin(2.*pi*y)
						- 8.*E_*pi*pi*cos(2.*pi*t)*cos(2.*pi*x)*sin(2.*pi*y);*/
		/*double t2 = sin(t);
		double t3 = pi * x * 2.0;
		double t4 = pi * y * 2.0;
		double t5 = cos (t3);
		double t6 = sin (t4);
		double t0 = t2 * cos (t4) * sin (t3) - t2 * t5 * t6 - E_ * (pi * pi) * t2 * t5 * t6 * 8.0;

		return t0;*/
		force = cos(2.*pi*y)*sin(2.*pi*x)*sin(t) - cos(2.*pi*x)*sin(2.*pi*y)*sin(t)
		        - 8.*E_*pi*pi*cos(2.*pi*x)*sin(2.*pi*y)*sin(t);
		return force;
	}

	double PNPCnForce(double t, const ZEROPTV &pt, double A2_, double B2_, double C2_, double D2_) const{
		/// Revert coefficients in NP equation to their original value regardless of 'ifMulDT'

		/*if (input_data_->ifMulDt)  {
			B2_ /= dt;
			C2_ /= dt;
			D2_ /= dt;
		} else {
			A2_ *= dt;
		}*/

		double x = pt.x();
		double y = pt.y();
		double cnForcing;
		double Re = input_data_->Re;

		/// forcing term for anion NP equation of full NSPNP equation
		if (!input_data_->ifPurePNP) {
//			cnForcing = B2_ * (2. * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * y)
//			                   * cos (2. * pi * y) * sin (2. * pi * x)
//			                   + 2. * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * sin (2. * pi * x)
//			                     * sin (2. * pi * y) * sin (2. * pi * y))
//			            + 8. * C2_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * y) * sin (2. * pi * x)
//			            - 2. * A2_ * pi * cos (2. * pi * y) * sin (2. * pi * t) * sin (2. * pi * x)
//			            + 16. * D2_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * y)
//			              * sin (2. * pi * x) * sin (2. * pi * y);
			double t2 = sin(t);
			double t3 = pi * pi;
			double t5 = pi * x * 2.0;
			double t6 = pi * y * 2.0;
			double t4 = t2 * t2;
			double t7 = cos (t5);
			double t8 = cos (t6);
			double t9 = sin (t5);
			double t10 = sin (t6);
			cnForcing = B2_ * (t3 * t4 * (t9 * t9) * t10 * pow (sin (pi * y), 2.0) * 2.0
			                   + t3 * t4 * t7 * t8 * t10 * pow (sin (pi * x), 2.0) * 2.0) + A2_ * t8 * t9 * cos (t)
			            + C2_ * t2 * t3 * t8 * t9 * 8.0 + D2_ * t3 * t4 * t7 * t8 * t9 * t10 * 1.6E+1;
		} else {  /// anion NP equation of full NSPNP equation
//			cnForcing = 8. * C2_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * y) * sin (2. * pi * x)
//			            - 2. * A2_ * pi * cos (2. * pi * y) * sin (2. * pi * t) * sin (2. * pi * x)
//			            - 16. * D2_ * z[1] * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x)
//			              * cos (2. * pi * y) * sin (2. * pi * x) * sin (2. * pi * y);
			/*double t2 = sin (t);
			double t3 = pi * pi;
			double t4 = pi * x * 2.0;
			double t5 = pi * y * 2.0;
			double t6 = cos (t5);
			double t7 = sin (t4);
			cnForcing = A2_ * t6 * t7 * cos (t) + C2_ * t2 * t3 * t6 * t7 * 8.0
			            + D2_ * (t2 * t2) * t3 * t6 * t7 * cos (t4) * sin (t5) * 1.6E+1;*/
			/*cnForcing = 8.*C2_*pi*pi*cos(2.*pi*t)*cos(2.*pi*y)*sin(2.*pi*x)
									- 2.*A2_*pi*cos(2.*pi*y)*sin(2.*pi*t)*sin(2.*pi*x)
									+ 16.*D2_*pi*pi*cos(2.*pi*t)*cos(2.*pi*t)*cos(2.*pi*x)*cos(2.*pi*y)*sin(2.*pi*x)*sin(2.*pi*y);*/
			cnForcing = A2_*cos(2.*pi*y)*sin(2.*pi*x)*cos(t) + 8.*C2_*pi*pi*cos(2.*pi*y)*sin(2.*pi*x)*sin(t)
			            + 16.*D2_*pi*pi*cos(2.*pi*x)*cos(2.*pi*y)*sin(2.*pi*x)*sin(2.*pi*y)*sin(t)*sin(t);
		}
		return cnForcing;
	}

	double PNPCpForce(double t, const ZEROPTV &pt, double A1_, double B1_, double C1_, double D1_) const{

		/// Revert coefficients in NP equation to their original value regardless of 'ifMulDT'
		/*if (input_data_->ifMulDt)  {
			B1_ /= dt;
			C1_ /= dt;
			D1_ /= dt;
		} else {
			A1_ *= dt;
		}*/

		double x = pt.x();
		double y = pt.y();
		double cpForcing;

		/// forcing term for cation NP equation of full NSPNP equation
		if (!input_data_->ifPurePNP) {
//			cpForcing = D1_ * (4. * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * x)
//			                   * cos (2. * pi * y) * cos (2. * pi * y)
//			                   + 4. * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * sin (2. * pi * x) * sin (2. * pi * x)
//			                     * sin (2. * pi * y) * sin (2. * pi * y))
//			            - B1_ * (2. * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * x)
//			                     * cos (2. * pi * y) * sin (2. * pi * y)
//			                     + 2. * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * y) * sin (2. * pi * x)
//			                       * sin (2. * pi * x) * sin (2. * pi * y))
//			            + 8. * C1_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * x) * sin (2. * pi * y)
//			            - 8. * D1_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * x)
//			              * sin (2. * pi * y) * sin (2. * pi * y)
//			            - 2. * A1_ * pi * cos (2. * pi * x) * sin (2. * pi * t) * sin (2. * pi * y);
			double t2 = sin(t);
			double t3 = pi * x;
			double t4 = pi * pi;
			double t7 = pi * y * 2.0;
			double t5 = t2 * t2;
			double t6 = t3 * 2.0;
			double t8 = sin (t3);
			double t10 = cos (t7);
			double t12 = sin (t7);
			double t9 = cos (t6);
			double t11 = sin (t6);
			double t13 = t8 * t8;
			double t15 = t12 * t12;
			double t14 = t9 * t9;
			cpForcing = D1_ * (t4 * t5 * (t10 * t10) * t14 * 4.0 + t4 * t5 * (t11 * t11) * t15 * 4.0)
			            - B1_ * (t4 * t5 * t11 * t13 * t15 * 2.0 - t4 * t5 * t9 * t10 * t12 * t13 * 2.0)
			            + A1_ * t9 * t12 * cos (t) + C1_ * t2 * t4 * t9 * t12 * 8.0 - D1_ * t4 * t5 * t14 * t15 * 8.0;
		} else {  /// if solving for PNP or P_NP equation
			/*cpForcing = D1_ * z[0]
			            * (4. * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x) * cos (2. * pi * x)
			               * cos (2. * pi * y) * cos (2. * pi * y)
			               + 4. * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * sin (2. * pi * x) * sin (2. * pi * x)
			                 * sin (2. * pi * y) * sin (2. * pi * y))
			            + 8. * C1_ * pi * pi * cos (2. * pi * t) * cos (2. * pi * x) * sin (2. * pi * y)
			            - 2. * A1_ * pi * cos (2. * pi * x) * sin (2. * pi * t) * sin (2. * pi * y)
			            - 8. * D1_ * z[0] * pi * pi * cos (2. * pi * t) * cos (2. * pi * t) * cos (2. * pi * x)
			              * cos (2. * pi * x) * sin (2. * pi * y) * sin (2. * pi * y);
			*cpForcing = D1_*(4.*pi*pi*cos(2.*pi*t)*cos(2.*pi*t)*cos(2.*pi*x)*cos(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*y)
											 + 4.*pi*pi*cos(2.*pi*t)*cos(2.*pi*t)*sin(2.*pi*x)*sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*y))
									+ 8.*C1_*pi*pi*cos(2.*pi*t)*cos(2.*pi*x)*sin(2.*pi*y)
									- 8.*D1_*pi*pi*cos(2.*pi*t)*cos(2.*pi*t)*cos(2.*pi*x)*cos(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*y)
									- 2.*A1_*pi*cos(2.*pi*x)*sin(2.*pi*t)*sin(2.*pi*y);*/
			/*double t2 = sin (t);
			double t3 = pi * pi;
			double t5 = pi * x * 2.0;
			double t6 = pi * y * 2.0;
			double t4 = t2 * t2;
			double t7 = cos (t5);
			double t8 = sin (t6);
			double t9 = t7 * t7;
			double t10 = t8 * t8;
			cpForcing = D1_ * (t3 * t4 * t10 * pow (sin (t5), 2.0) * 4.0 + t3 * t4 * t9 * pow (cos (t6), 2.0) * 4.0)
			            + A1_ * t7 * t8 * cos (t) + C1_ * t2 * t3 * t7 * t8 * 8.0 - D1_ * t3 * t4 * t9 * t10 * 8.0;*/
			cpForcing = D1_*(4.*pi*pi*cos(2.*pi*x)*cos(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*y)*sin(t)*sin(t)
			                 + 4.*pi*pi*sin(2.*pi*x)*sin(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*y)*sin(t)*sin(t))
			            + A1_*cos(2.*pi*x)*sin(2.*pi*y)*cos(t) + 8.*C1_*pi*pi*cos(2.*pi*x)*sin(2.*pi*y)*sin(t)
			            - 8.*D1_*pi*pi*cos(2.*pi*x)*cos(2.*pi*x)*sin(2.*pi*y)*sin(2.*pi*y)*sin(t)*sin(t);
		}
		return cpForcing;
	}



    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {
       getIntegrandsPNPAe(fe, Ae, dt_, t_);
    }

    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
        getIntegrandsPNPbe(fe, be, dt_, t_);
    }

    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT &boundaryType, const DENDRITE_UINT &
    boundaryID, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h) {
        //getIntegrands4sidePNPAe(fe, boundaryType, Ae, dt_, t_);
    }

    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT &boundaryType, const DENDRITE_UINT &
    boundaryID, TALYFEMLIB::ZEROARRAY<double> &be, const double *h)
    {
        //getIntegrands4sidePNPbe(fe, boundaryType, be, dt_, t_);
    }

};
