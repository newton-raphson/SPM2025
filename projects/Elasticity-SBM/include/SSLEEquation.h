#pragma once

#include <talyfem/fem/cequation.h>
#include <talyfem/stabilizer/tezduyar_upwind.h>
#include "LENodeData.h"
#include "LEInputData.h"
#include <talyfem/talyfem.h>
#include <link.h>
#include <cmath>

#ifdef DEEPTRACE
#include "SBMCalcDeepTrace.h"
#else
#include "SBMCalc.h"
#endif
#include "util.h"

class SSLEEquation : public TALYFEMLIB::CEquation<LENodeData>
{

  const bool IFSBM = true;
  bool ShortestDist = true;
  const IBM_METHOD method;
  const IMGA *imga_;
  my_kd_tree_t *kd_tree_;


public:
  explicit SSLEEquation(LEInputData *idata, const IBM_METHOD _method, const IMGA *imga)
      : TALYFEMLIB::CEquation<LENodeData>(false, TALYFEMLIB::kAssembleGaussPoints), method(_method), imga_(imga)
  {
    idata_ = idata;
  }

  void Solve(double dt, double t) override
  {
    assert(false);
  }

  void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                  TALYFEMLIB::ZEROARRAY<double> &be) override
  {
    assert(false);
  }

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae)
  {

    double thickness = 1;

    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = fe.nsd();

    // # of basis functions
    const int n_basis_functions = fe.nbf();
    // (determinant of J) cross W

    const double detJxW = fe.detJxW();

    // C matrix
      double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];
      CalcCmatrix(Cmatrix);

//
//
////        const double detJxW = fe.detJxW();
////    for (int a = 0; a < n_basis_functions; a++)
////    {
////      for (int b = 0; b < n_basis_functions; b++)
////      {
////        for (int i = 0; i < DIM; i++)
////        {
////          for (int j = 0; j < DIM; j++)
////          {
////            Ae(DIM * a + i, DIM * b + j) += fe.dN(b, j) * idata_->Cmatrix[i][j] * fe.dN(a, i) * detJxW + fe.dN(b, DIM - j - 1) * idata_->Cmatrix[DIM][DIM] * fe.dN(a, DIM - i - 1) * detJxW;
////          }
////        }
////      }
////    }
//
//for(int a = 0; a < n_basis_functions; a++)
//{ // loop over test basis functions (a)
//    for(int b = 0; b < n_basis_functions; b++)
//    { // loop over trial basis functions (b)
//        for(int i = 0; i < DIM; i++) { // loop over DOF test function i
//
//            for (int k = 0; k < DIM; k++) { // loop over DOF trial function k
//                for (int l = 0; l < DIM; l++) { // internal sum l
//                    for (int j = 0; j < DIM; j++) { // internal sum j
//                        Ae(DIM * a + i, DIM * b + k) +=
//                                fe.dN(b, l) * Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.dN(a, j) * detJxW * 0.5;
//                    } // end j
//
//                } // end l
//
//            } // end k
//
//            /////// symmetric part  //////
//            for(int l = 0; l < DIM; l++) { // loop over DOF trial function l
//                for (int k = 0; k < DIM; k++) { // internal sum k
//                    for (int j = 0; j < DIM; j++) { // internal sum j
//                        Ae(DIM * a + i, DIM * b + l) +=
//                                fe.dN(b, k) * Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.dN(a, j) * detJxW * 0.5;
//                    } // end j
//
//                } // end k
//            } // end l
//
//
//        } // end i
//    } // end b
//
//} // end a


#if (DIM == 3)

      //            std::vector<std::vector<double>> Ae_check(n_dimensions * n_basis_functions, std::vector<double>(n_dimensions * n_basis_functions));
      //            std::vector<std::vector<double>> Be(n_dimensions * n_basis_functions, std::vector<double>(6));
      //            CalcBe(fe,Be);
      //
      //            std::vector<std::vector<double>> BeCmatrix(n_dimensions * n_basis_functions, std::vector<double>(6));
      //            CalcBeCmatrix(fe, Be, idata_->Cmatrix, BeCmatrix);
      //
      ////             for final term -> previous implementation!
      //            for (int a = 0; a < n_dimensions * n_basis_functions; a++)
      //            {
      //                for (int b = 0; b < n_dimensions * n_basis_functions; b++)
      //                {
      ////                    Ae_check[a][b] = 0;
      //                    for (int k = 0; k < 6; k++)
      //                    {
      //                        Ae(a, b) += BeCmatrix[a][k] * Be[b][k] * detJxW;
      ////                        Ae_check[a][b] += BeCmatrix[a][k] * Be[b][k] * detJxW;
      //                    }
      //                }
      //            }

      /*
       * below is from CalcAe3D.m.
       * And, we try to know the pattern by ourselves
       */
      const double Emv = idata_->Cmatrix[0][0];
      const double Ev = idata_->Cmatrix[1][0];
      const double half = idata_->Cmatrix[3][3];
      //
      //    for (int i = 0; i < fe.nbf(); ++i) {
      //        for (int j = 0; j < fe.nbf(); ++j) {
      //            for (int k1 = 0; k1 < DIM; k1++) {
      //                for (int k2 = 0; k2 < DIM; k2++) {
      //                    // First type of element
      //                    if (k1 == k2) {
      //                        Ae(DIM * i + k1, DIM * j + k2) += Emv * fe.dN(i, k1) * fe.dN(j, k1) * detJxW;
      ////                        Ae_check[DIM * i + k1][ DIM * j + k2]-= Emv * fe.dN(i, k1) * fe.dN(j, k1) * detJxW;
      //                        if (k1 == 0) {
      //                            Ae(DIM * i + k1, DIM * j + k2) +=
      //                                    half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
      ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
      //                        } else if (k1 == 1) {
      //                            Ae(DIM * i + k1, DIM * j + k2) +=
      //                                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
      ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
      //                        } else if (k1 == 2) {
      //                            Ae(DIM * i + k1, DIM * j + k2) +=
      //                                    half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
      ////                            Ae_check[DIM * i + k1][ DIM * j + k2]-= half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
      //                        }
      //                    }
      //                        // Second type of element -> switch k1 and k2
      //                    else {
      //                        Ae(DIM * i + k1, DIM * j + k2) +=
      //                                (Ev * fe.dN(i, k1) * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
      ////                        Ae_check[DIM * i + k1][ DIM * j + k2]-= (Ev * fe.dN(i, k1) * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
      //                    }
      //                }
      //            }
      //        }
      //    }

      const int nbf = fe.nbf();
      for (int i = 0; i < nbf; ++i)
      {
          int DIM_i = DIM * i;

          for (int j = 0; j < nbf; ++j)
          {
              int DIM_j = DIM * j;

              for (int k1 = 0; k1 < DIM; k1++)
              {
                  double dN_ik1 = fe.dN(i, k1); // Cache this value for k1 loop

                  for (int k2 = 0; k2 < DIM; k2++)
                  {
                      if (k1 == k2)
                      {
                          Ae(DIM_i + k1, DIM_j + k2) += Emv * dN_ik1 * fe.dN(j, k1) * detJxW;

                          switch (k1)
                          {
                              case 0:
                                  Ae(DIM_i + k1, DIM_j + k2) +=
                                          half * (fe.dN(i, 1) * fe.dN(j, 1) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
                                  break;
                              case 1:
                                  Ae(DIM_i + k1, DIM_j + k2) +=
                                          half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 2) * fe.dN(j, 2)) * detJxW;
                                  break;
                              case 2:
                                  Ae(DIM_i + k1, DIM_j + k2) +=
                                          half * (fe.dN(i, 0) * fe.dN(j, 0) + fe.dN(i, 1) * fe.dN(j, 1)) * detJxW;
                                  break;
                          }
                      }
                      else
                      {
                          Ae(DIM_i + k1, DIM_j + k2) +=
                                  (Ev * dN_ik1 * fe.dN(j, k2) + half * fe.dN(i, k2) * fe.dN(j, k1)) * detJxW;
                      }
                  }
              }
          }
      }

#endif

  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be)
  {


    using namespace TALYFEMLIB;
    const ZEROPTV &p = fe.position();

    ZEROPTV BodyForce;
    bool ForceHaveSet = false;
    CalcForce(p, BodyForce, ForceHaveSet);
    if (!ForceHaveSet)
    {
      BodyForce = idata_->BodyForce;
    }

#if (DIM == 3)
    double body_z = idata_->BodyForce[2];
#endif
    double BR_V = idata_->radialbodyforce.br_v;
    int BR_POW = idata_->radialbodyforce.br_pow;

//      BodyForce.print();

    /*
     * please write something like this so that it can support both 2D and 3D
     */
    for (int a = 0; a < fe.nbf(); a++)
    {
      for (int dim = 0; dim < DIM; dim++)
      {
        be(DIM * a + dim) += fe.N(a) * BodyForce(dim) * fe.detJxW();
      }
    }

#if (DIM == 2)
    double x_min = idata_->mesh_def.physDomain.min[0];
    double y_min = idata_->mesh_def.physDomain.min[1];
    double x_max = idata_->mesh_def.physDomain.max[0];
    double y_max = idata_->mesh_def.physDomain.max[1];
    double x_mid = (x_min + x_max) / 2;
    double y_mid = (y_min + y_max) / 2;

    double x = p.x();
    double y = p.y();

    for (int a = 0; a < fe.nbf(); a++)
    {
      double r_value = 0;
      double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
      r_value = BR_V * pow(radius, BR_POW);
      // std::cout<<"r_value:"<<r_value<<"\n";
      double sin_x = (x - x_mid) / radius;
      double sin_y = (y - y_mid) / radius;
      be(2 * a) += fe.N(a) * r_value * sin_x * fe.detJxW();
      be(2 * a + 1) += fe.N(a) * r_value * sin_y * fe.detJxW();
    }
#endif

  }

  ///// ==================== ibm start====================================
#pragma mark Ae-be Surface IBM

  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {
    assert(method == IBM_METHOD::NITSCHE);
  } // end:ibm_Ae

  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {
    assert(method == IBM_METHOD::NITSCHE);
  } // end:ibm_be

  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                              const NodeAndValues<DENDRITE_REAL> &gpinfo,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h,
                              const std::vector<double> &surface_values)
  {
    ibm_Integrands4side_Ae(fe, Ae, gpinfo, position, h);
  }

  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &Ae,
                              const NodeAndValues<DENDRITE_REAL> &gpinfo,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h,
                              const std::vector<double> &surface_values)
  {
    ibm_Integrands4side_be(fe, Ae, gpinfo, position, h);
  }
  ///// ==================== ibm end====================================
#pragma mark Ae-be Surface Carved
  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZeroMatrix<double> &Ae) {



    assert(method == IBM_METHOD::SBM);

    double h = ElementSize(fe);
    double d[DIM];
    int geom_ID;
#ifdef DEEPTRACE
    //// this is of no use
    std::vector<my_kd_tree_t*> kd_tree_;
    //// add the kd_tree_ to the pointer

    SBMCalcDeepTrace sbmCalc(fe, idata_, imga_, kd_tree_, 0 /*CarvedOutGeomID*/);
#else

    SBMCalc sbmCalc(fe, idata_, imga_, kd_tree_);
#endif
      bool DirichletHaveSet = false;
      double BCValue[DIM];
      sbmCalc.Dist2Geo(d);
      SBMCalc::BCTypes bcType;
      sbmCalc.GetBC(d, BCValue, bcType);

    // penalty
    const double Cb_e = idata_->Cb_e.value_at(t_);
    const double TauN = idata_->TauN.value_at(t_);

    //  ====== parameters setting end ============

    const int nbf = fe.nbf();
    const int n_basis_functions = fe.nbf();
    const int nsd = DIM;
    const double detSideJxW = fe.detJxW();

    //////////////////////////////////////weak//////////////////////////////////////////

    double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];
    std::vector<std::vector<double>> Be(DIM * n_basis_functions);
    /// for middle term => B_T*C_T
    std::vector<std::vector<double>> BeCmatrix(DIM * n_basis_functions);
    double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
    /// for mid2 term => B_T*C_T*n
    std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);
    ////////////
    if (bcType == SBMCalc::BCTypes::DIRICHLET)
      {
          WriteCSV("dirichlet.csv", {fe.position().x(), fe.position().y(),fe.position().z(),BCValue[0],BCValue[1],BCValue[2]},"X,Y,Z,BCX,BCY,BCZ");
          CalcCmatrix(Cmatrix);
          CalcBe(fe, Be);
          CalcBeCmatrix(fe, Be, Cmatrix, BeCmatrix);
          CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);
          CalcStressDotNormal(fe, BeCmatrix, SurrogateNormalMatrix, StressDotSurrogateNormal); // (DIM*fe.nbf()) * DIM

          double Ne_[DIM][DIM * n_basis_functions];
          memset(Ne_, 0.0, sizeof Ne_);
          double Ne_con_[DIM][DIM * n_basis_functions];
          memset(Ne_con_, 0.0, sizeof Ne_con_);

          DENDRITE_REAL secondOrderTerm_a_(0);
          for (int j = 0; j < n_basis_functions; j++)
          {
              // secondOrderTerm_a = d[0] * (fe.d2N(j, 0, 0) * d[0] + fe.d2N(j, 0, 1) * d[1]) + d[1] * (fe.d2N(j, 1, 0) * d[0] + fe.d2N(j, 1, 1) * d[1]) / 2;
              double gradWdotd = 0.0;
              for (int dim = 0; dim < DIM; dim++)
              {
                  gradWdotd += fe.dN(j, dim) * d[dim];
              }

              for (int dim = 0; dim < DIM; dim++)
              {
                  for (int dim2 = 0; dim2 < DIM; dim2++)
                  {
                      if (dim2 == dim)
                      {
                          Ne_[dim][DIM * j + dim2] = fe.N(j) + gradWdotd + secondOrderTerm_a_
                                  /*TODO: fix small issue here when we use QBF*/;
                      }
                      else
                      {
                          Ne_[dim][DIM * j + dim2] = 0.0;
                      }
                  }
              }
          }
          for (int j = 0; j < n_basis_functions; j++)
          {
              for (int dim = 0; dim < DIM; dim++)
              {
                  for (int dim2 = 0; dim2 < DIM; dim2++)
                  {
                      if (dim2 == dim)
                      {
                          Ne_con_[dim][DIM * j + dim2] = fe.N(j);
                      }
                      else
                      {
                          Ne_con_[dim][DIM * j + dim2] = 0.0;
                      }
                  }
              }
          }

          // for Final term => B_T*C_T*n*N
          const double detJxW = fe.detJxW();
          for (int a = 0; a < DIM * n_basis_functions; a++)
          {
              for (int b = 0; b < DIM * n_basis_functions; b++)
              {
                  double N = 0;
                  double N_con = 0;
                  // StressDotSurrogateNormal -> (DIM*fe.nbf()) * DIM
                  for (int k = 0; k < DIM; k++)
                  {
                      N += StressDotSurrogateNormal[a][k] * Ne_[k][b] * detJxW;
                      N_con += StressDotSurrogateNormal[a][k] * Ne_con_[k][b] * detJxW;
                  }
                  // Ae is symmetric

                  Ae(a, b) += N; // adjoint consistency

                  Ae(b, a) -= N_con; // consistency
              }
          }

              DENDRITE_REAL  secondOrderTerm_b_(0);


          double weakBCpenaltyParameter_ = util_funcs::ReturnPenaltyParameters(idata_) * Cb_e / h;

              for (int a = 0; a < fe.nbf(); a++)
              {
#if (DIM == 2)
                  if (idata_->elemOrder == 2 && idata_->ifHessian)
          {
            secondOrderTerm_a_ = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                 d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                                2;
          }
          else
          {
            secondOrderTerm_a_ = 0;
          }
#endif

#if (DIM == 3)
                  if (idata_->elemOrder == 2 && idata_->ifHessian)
                  {

                      secondOrderTerm_a_ = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
                  }
                  else
                  {
                      secondOrderTerm_a_ = 0;
                  }
#endif

                  double gradWdotd = 0.0;
                  for (int k = 0; k < DIM; k++)
                  {
                      gradWdotd += fe.dN(a, k) * d[k];
                  }

                  for (int b = 0; b < fe.nbf(); b++)
                  {
#if (DIM == 2)
                      if (idata_->elemOrder == 2 && idata_->ifHessian)
            {
              secondOrderTerm_b_ = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1]) +
                                   d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1])) /
                                  2;
            }
            else
            {
              secondOrderTerm_b_ = 0;
            }
#endif

#if (DIM == 3)
                      if (idata_->elemOrder == 2 && idata_->ifHessian)
                      {

                          secondOrderTerm_b_ = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1] + fe.d2N(b, 0, 2) * d[2]) + d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1] + fe.d2N(b, 1, 2) * d[2]) + d[2] * (fe.d2N(b, 2, 0) * d[0] + fe.d2N(b, 2, 1) * d[1] + fe.d2N(b, 2, 2) * d[2])) / 2;
                      }
                      else
                      {
                          secondOrderTerm_b_ = 0;
                      }
#endif
                      double gradUdotd = 0.0;
                      for (int k = 0; k < DIM; k++)
                      {
                          gradUdotd += fe.dN(b, k) * d[k];
                      }

                      /*
                       *  make it j to match with what's inside Navier-Stokes
                       */
                      for (int j = 0; j < DIM; j++)
                      {
                          Ae(DIM * a + j, DIM * b + j) +=
                                  +weakBCpenaltyParameter_ * (fe.N(a) + gradWdotd + secondOrderTerm_a_)
                                  * (fe.N(b) + gradUdotd + secondOrderTerm_b_) * detSideJxW; // penalty
                      }
                  } // b loop`
              }   // a loop


          return;
      }
    if (bcType == SBMCalc::BCTypes::NEUMANN)
      {
        WriteCSV("neumann.csv", {fe.position().x(), fe.position().y(),fe.position().z(),BCValue[0],BCValue[1],BCValue[2]},"X,Y,Z,BCX,BCY,BCZ");

          const ZEROPTV &p = fe.position();
          const ZEROPTV &SurrogateNormal = fe.surface()->normal();

          ZEROPTV TrueNormal;

          sbmCalc.NormalofGeo(TrueNormal,d);

          double SurrogateDotTrueNormal = SurrogateNormal.innerProduct(TrueNormal);

          CalcCmatrix(Cmatrix);

          const double detJxW = fe.detJxW();


          //// In the loops below we try to perform the following operations basically:
          //// (w_i, n.\tilde{n} (C_ijkl \partial_q \partial_l u_k n_j + C_ijkl \partial_l u_k n_j) - C_ijkl \partial_l u_k \tilde{n_j}  ) dx
          for (int a = 0; a < fe.nbf(); a++) {

              for (int b = 0; b < fe.nbf(); b++) {

                  for (int i = 0; i < DIM; i++) {
                      for (int k = 0; k < DIM; k++) {
                          double Cijkl_d_Hessian_n = 0.0;
                          double Cijkl_gradu_n = 0.0;
                          double Cijkl_gradu_n_tilde = 0.0;

                          for (int j = 0; j < DIM; j++) {
                              double Cijkl_d_Hessian = 0.0;
                              double Cijkl_gradu = 0.0;

                              for (int l = 0; l < DIM; l++) {
                                  /// Compute the Hessian term correctly
                                  double hessian_dot_distance = 0.0;
                                  for (int q = 0; q < DIM; q++) {
                                      hessian_dot_distance += fe.d2N(b, l, q) * d[q];
                                  }


                                  Cijkl_d_Hessian += Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * hessian_dot_distance;
                                  Cijkl_gradu += Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(b, l);
                              }

                              // Correctly accumulate Hessian contribution
                              Cijkl_d_Hessian_n += Cijkl_d_Hessian * TrueNormal[j];
                              Cijkl_gradu_n += Cijkl_gradu * TrueNormal[j];
                              Cijkl_gradu_n_tilde += Cijkl_gradu * SurrogateNormal[j];
                          }
                          Ae(DIM * a + i, DIM * b + k) += fe.N(a) * (SurrogateDotTrueNormal *
                                                                     (Cijkl_d_Hessian_n + Cijkl_gradu_n) - // Hessian now contributes
                                                                     Cijkl_gradu_n_tilde) * detJxW;
                      }
                  }
              }
          }

          return;
      }



  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZEROARRAY<double> &be)
  {

    assert(method == IBM_METHOD::SBM);

    double h = ElementSize(fe);
    double d[DIM];
    int geom_ID;


#ifdef DEEPTRACE
      //// this is of no use
      std::vector<my_kd_tree_t*> kd_tree_;
      //// add the kd_tree_ to the pointer

      SBMCalcDeepTrace sbmCalc(fe, idata_, imga_, kd_tree_, 0 /*CarvedOutGeomID*/);
#else
    SBMCalc sbmCalc(fe, idata_, imga_, kd_tree_);
#endif
    bool DirichletHaveSet = false;
    double BCValue[DIM];
    sbmCalc.Dist2Geo(d);
    SBMCalc::BCTypes bcType;
    sbmCalc.GetBC(d, BCValue, bcType);
//    sbmCalc.GetDirichletBC(d, BCValue, DirichletHaveSet);

#ifndef  NDEBUG
ZEROPTV true_points = {fe.position().x()+d[0],fe.position().y()+d[1],fe.position().z()+d[2]};
ZEROPTV true_normal;
      sbmCalc.NormalofGeo(true_normal,d);
// write csv
WriteCSV("true_points.csv", {true_points.x(), true_points.y(),true_points.z(),true_normal.x(),true_normal.y(),true_normal.z()},"X,Y,Z,NX,NY,NZ");


#endif

    /// finish: calculate d vector ======================================================================


    // ====== parameters setting ============
    // penalty
    const double Cb_e = idata_->Cb_e.value_at(t_);
    const double TauN = idata_->TauN.value_at(t_);

    //  ====== parameters setting end ============

    const int nbf = fe.nbf();
    const int n_basis_functions = fe.nbf();
    const int nsd = DIM;
    const double detSideJxW = fe.detJxW();

    //////////////////////////////////////weak//////////////////////////////////////////

    ///
    double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];
    // for middle term => B_T*C_T
    std::vector<std::vector<double>> BeCmatrix(DIM * n_basis_functions);
    std::vector<std::vector<double>> Be(DIM * n_basis_functions);

      if (bcType == SBMCalc::BCTypes::DIRICHLET){

          WriteCSV("dirichlet_be.csv", {fe.position().x(), fe.position().y(),fe.position().z(),BCValue[0],BCValue[1],BCValue[2]},"X,Y,Z,BCX,BCY,BCZ");

          CalcCmatrix(Cmatrix);
          CalcBe(fe, Be);
          CalcBeCmatrix(fe, Be, Cmatrix, BeCmatrix);
          double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
          CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);

          // for mid2 term => B_T*C_T*n
          std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);
          CalcStressDotNormal(fe, BeCmatrix, SurrogateNormalMatrix, StressDotSurrogateNormal);

          double weakBCpenaltyParameter_ = util_funcs::ReturnPenaltyParameters(idata_) * Cb_e / h;



          DENDRITE_REAL secondOrderTerm_a_(0);
          for (int a = 0; a < fe.nbf(); a++)
          {
#if (DIM == 2)
              if (idata_->elemOrder == 2 && idata_->ifHessian)
        {
          secondOrderTerm_a_ = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                               d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                              2;
        }
        else
        {
          secondOrderTerm_a_ = 0;
        }
#endif

#if (DIM == 3)
              if (idata_->elemOrder == 2 && idata_->ifHessian)
              {

                  secondOrderTerm_a_ = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
              }
              else
              {
                  secondOrderTerm_a_ = 0;
              }
#endif
          }

          for (int a = 0; a < DIM * n_basis_functions; a++)
          {
              for (int dim = 0; dim < DIM; dim++)
              {
                  be(a) +=
                          StressDotSurrogateNormal[a][dim] * BCValue[dim] * detSideJxW; // adjoint
              }
          }

          for (int a = 0; a < n_basis_functions; a++)
          {
              double gradWdotd = 0.0;
              for (int k = 0; k < DIM; k++)
              {
                  gradWdotd += fe.dN(a, k) * d[k];
              }
              for (int i = 0; i < DIM; i++)
              {
                  be(DIM * a + i) +=
                          +weakBCpenaltyParameter_ * (fe.N(a) + gradWdotd + secondOrderTerm_a_) * BCValue[i] * detSideJxW; // penalty
              }
          }
          return;
      }

      if (bcType == SBMCalc::BCTypes::NEUMANN) {
            const ZEROPTV &p = fe.position();
            const ZEROPTV &SurrogateNormal = fe.surface()->normal();

            ZEROPTV TrueNormal;

            sbmCalc.NormalofGeo(TrueNormal,d);

            double SurrogateDotTrueNormal = SurrogateNormal.innerProduct(TrueNormal);


            const double detJxW = fe.detJxW();


            for (int a = 0; a < fe.nbf(); a++) {
                for (int i = 0; i < DIM; i++) {
                    be(DIM * a + i) += fe.N(a) * SurrogateDotTrueNormal * BCValue[i] *TrueNormal[i] * detJxW;
                }
            }

        }

      return;

  }

#pragma mark Ae-be to support more dendrite-kt version

  void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h)
  {
    Integrands_Ae(fe, Ae);
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const double *h)
  {
    Integrands_be(fe, be);
  }

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae,
                          const double *h)
  {
    Integrands4side_Ae(fe, side_idx, id, Ae);
  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be,
                          const double *h)
  {
    Integrands4side_be(fe, side_idx, id, be);
  }

    void setSbmCalc(my_kd_tree_t *kd_tree)
    {
        kd_tree_ = kd_tree;
    }

private:
  LEInputData *idata_;

  /**
   * This function is going to check the SBMGeo and find out the corresponding MMS force calculated from MMS solution
   * @param p IN
   * @param BodyForce OUT
   * @param ForceHaveSet OUT
   */
  void CalcForce(const TALYFEMLIB::ZEROPTV &p, ZEROPTV &BodyForce, bool &ForceHaveSet) const
  {
    switch (idata_->SbmGeo)
    {


    case LEInputData::SBMGeo::NONE:
    {
      ForceHaveSet = false;
      break;
    }
        case LEInputData::SBMGeo::SPHERE:
        {
            ///
            double pi = M_PI;
            double E = idata_->planeStress.young;
            double v = idata_->planeStress.poisson;

            ///
            double x = p.x();
            double y = p.y();
            double z = p.z();

            ///
            BodyForce.x() = (E * ((pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / 20 + (pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / 10) * (v - 0.5)) / ((2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / (20 * (2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z) * (v - 1)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z) * (v - 0.5)) / (5 * (2 * v - 1) * (v + 1));
            BodyForce.y() = (E * ((pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / 10 + (pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / 20) * (v - 0.5)) / ((2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z)) / (10 * (2 * v - 1) * (v + 1)) - (E * v * pow(pi, 2) * cos(pi * y) * sin(pi * x) * sin(pi * z)) / (20 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z) * (v - 1)) / (10 * (2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * x) * sin(pi * y) * sin(pi * z) * (v - 0.5)) / (5 * (2 * v - 1) * (v + 1));
            BodyForce.z() = (E * v * pow(pi, 2) * cos(pi * x) * cos(pi * y) * cos(pi * z)) / (5 * (2 * v - 1) * (v + 1)) - (2 * E * (v - 0.5) * ((pow(pi, 2) * cos(pi * x) * cos(pi * y) * cos(pi * z)) / 10 - (pow(pi, 2) * cos(pi * z) * sin(pi * x) * sin(pi * y)) / 20)) / ((2 * v - 1) * (v + 1)) + (E * pow(pi, 2) * cos(pi * z) * sin(pi * x) * sin(pi * y) * (v - 1)) / (20 * (2 * v - 1) * (v + 1));
            ForceHaveSet = true;
            break;
        }
        case LEInputData::SBMGeo::BUNNY:{
            BodyForce.x() = 0;
            BodyForce.y() = 0;
            BodyForce.z() = 0.0;
            ForceHaveSet = true;
            break;
        }

    default:
    {
      break;
    }
    }
  }

  DENDRITE_REAL normalDistance(const TALYFEMLIB::FEMElm &fe,
                               const ZEROPTV &normal,
                               const ZEROPTV &h)
  {
    double max_h = std::max(std::max(h.x(), h.y()), h.z());

    ZEROPTV root_node;
    fe.grid()->GetCoord(root_node, 0);
    int numNode = pow(idata_->elemOrder + 1, DIM);
    std::vector<ZEROPTV> elm_node(numNode);
#ifdef ENABLE_3D
    for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
    {
      for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
      {
        for (int p3 = 0; p3 < idata_->elemOrder + 1; p3++)
        {
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3] = root_node;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].x() += p3 * h.x() / idata_->elemOrder;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].y() += p2 * h.y() / idata_->elemOrder;
          elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].z() += p1 * h.y() / idata_->elemOrder;
        }
      }
    }
#else
    for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
    {
      for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
      {
        elm_node[p1 * (idata_->elemOrder + 1) + p2] = root_node;
        elm_node[p1 * (idata_->elemOrder + 1) + p2].x() += p2 * h.x() / idata_->elemOrder;
        elm_node[p1 * (idata_->elemOrder + 1) + p2].y() += p1 * h.y() / idata_->elemOrder;
      }
    }
#endif
    std::vector<double> hb(numNode);
    for (unsigned int i = 0; i < numNode; i++)
    {
      ZEROPTV vec1 = fe.position() - elm_node[i];
      hb[i] = vec1.innerProduct(normal) > 0 ? vec1.innerProduct(normal) : 0;
    }
    auto maxhb = std::max_element(hb.begin(), hb.end());
    if (*maxhb < max_h * 1e-2)
    {
      return max_h * 1e-2;
    }
    return *maxhb;
  }

#pragma mark normal matrix

  /*
   * this normal matrix calculation is based on t_i = stress_ij * n_j
   * the total formulation please check cheng-hau's LE presentation
   */
  void CalcSurrogateNormalMatrix(const TALYFEMLIB::FEMElm &fe, double (&SurrogateNormalMatrix)[DIM][3 * (DIM - 1)])
  {
#if (DIM == 2)
    SurrogateNormalMatrix[0][0] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[1][0] = 0;
    SurrogateNormalMatrix[0][1] = 0;
    SurrogateNormalMatrix[1][1] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[0][2] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[1][2] = fe.surface()->normal().data()[0];
#endif
#if (DIM == 3)
    SurrogateNormalMatrix[0][0] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[1][0] = 0;
    SurrogateNormalMatrix[2][0] = 0;
    SurrogateNormalMatrix[0][1] = 0;
    SurrogateNormalMatrix[1][1] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[2][1] = 0;
    SurrogateNormalMatrix[0][2] = 0;
    SurrogateNormalMatrix[1][2] = 0;
    SurrogateNormalMatrix[2][2] = fe.surface()->normal().data()[2];

    SurrogateNormalMatrix[0][3] = fe.surface()->normal().data()[1];
    SurrogateNormalMatrix[1][3] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[2][3] = 0;
    SurrogateNormalMatrix[0][4] = fe.surface()->normal().data()[2];
    SurrogateNormalMatrix[1][4] = 0;
    SurrogateNormalMatrix[2][4] = fe.surface()->normal().data()[0];
    SurrogateNormalMatrix[0][5] = 0;
    SurrogateNormalMatrix[1][5] = fe.surface()->normal().data()[2];
    SurrogateNormalMatrix[2][5] = fe.surface()->normal().data()[1];

#endif
  }

#pragma mark generating terms for integrations



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

  void CalcBe(const TALYFEMLIB::FEMElm &fe, std::vector<std::vector<double>> &Be)
  {
    assert(Be.size() == DIM * fe.nbf());
    for (int j = 0; j < fe.nbf(); j++)
    {
#if (DIM == 2)
      /// x-dir
      Be[2 * j].resize(3);
      /// y-dir
      Be[2 * j + 1].resize(3);

      Be[2 * j][0] = fe.dN(j, 0);
      Be[2 * j + 1][0] = 0;
      Be[2 * j][1] = 0;
      Be[2 * j + 1][1] = fe.dN(j, 1);
      Be[2 * j][2] = fe.dN(j, 1);
      Be[2 * j + 1][2] = fe.dN(j, 0);
#endif
#if (DIM == 3)
      /// x-dir
      Be[3 * j].resize(6);
      /// y-dir
      Be[3 * j + 1].resize(6);
      /// z-dir
      Be[3 * j + 2].resize(6);

      Be[3 * j][0] = fe.dN(j, 0);
      Be[3 * j + 1][0] = 0;
      Be[3 * j + 2][0] = 0;

      Be[3 * j][1] = 0;
      Be[3 * j + 1][1] = fe.dN(j, 1);
      Be[3 * j + 2][1] = 0;

      Be[3 * j][2] = 0;
      Be[3 * j + 1][2] = 0;
      Be[3 * j + 2][2] = fe.dN(j, 2);

      Be[3 * j][3] = fe.dN(j, 1);
      Be[3 * j + 1][3] = fe.dN(j, 0);
      Be[3 * j + 2][3] = 0;

      Be[3 * j][4] = fe.dN(j, 2);
      Be[3 * j + 1][4] = 0;
      Be[3 * j + 2][4] = fe.dN(j, 0);

      Be[3 * j][5] = 0;
      Be[3 * j + 1][5] = fe.dN(j, 2);
      Be[3 * j + 2][5] = fe.dN(j, 1);
#endif
    }
  }

  void CalcBeCmatrix(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &Be, const double (&Cmatrix)[3 * (DIM - 1)][3 * (DIM - 1)], std::vector<std::vector<double>> &BeCmatrix)
  {
    assert(BeCmatrix.size() == DIM * fe.nbf());
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      BeCmatrix[a].resize(3 * (DIM - 1));
      for (int b = 0; b < 3 * (DIM - 1); b++)
      {
        double N = 0;
        for (int k = 0; k < 3 * (DIM - 1); k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        BeCmatrix[a][b] = N;
      }
    }
  }

  void CalcStressDotNormal(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &BeCmatrix, const double (&NormalMatrix)[DIM][3 * (DIM - 1)], std::vector<std::vector<double>> &StressDotNormal)
  {
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      StressDotNormal[a].resize(DIM);
      for (int b = 0; b < DIM; b++)
      {
        for (int k = 0; k < 3 * (DIM - 1); k++)
        {
          StressDotNormal[a][b] += BeCmatrix[a][k] * NormalMatrix[b][k];
        }
      }
    }
  }

  void CalcBeCmatrix(const TALYFEMLIB::FEMElm &fe, const std::vector<std::vector<double>> &Be, const std::vector<std::vector<double>> Cmatrix, std::vector<std::vector<double>> &BeCmatrix)
  {
    assert(BeCmatrix.size() == DIM * fe.nbf());
    for (int a = 0; a < DIM * fe.nbf(); a++)
    {
      // BeCmatrix[a].resize(3*(DIM-1));
      for (int b = 0; b < 3 * (DIM - 1); b++)
      {
        double N = 0;
        for (int k = 0; k < 3 * (DIM - 1); k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        BeCmatrix[a][b] = N;
      }
    }
  }

  DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe)
  {
    return pow((pow(2, DIM) * fe.volume_jacc()), (double)1 / DIM);
  }


  int tensorToVoigt(int i, int j) {
        if (i == j) {
            // xx -> 0, yy -> 1
            return (i == 0) ? 0 : 1; // 0 for xx, 1 for yy
        } else {
            // xy or yx -> 2
            return 2;
        }
    }

};
