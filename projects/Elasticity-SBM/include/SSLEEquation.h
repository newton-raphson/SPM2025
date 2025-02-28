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

#ifndef NDEBUG
    std::ofstream fout("testdomaingp_Ae.txt", std::ios::app);
#if (DIM == 2)
    fout << fe.position().x() << "," << fe.position().y() << "\n";
#endif
#if (DIM == 3)
    fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << "\n";
#endif
    fout.close();
#endif
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


#if (DIM == 2)

//        const double detJxW = fe.detJxW();
//    for (int a = 0; a < n_basis_functions; a++)
//    {
//      for (int b = 0; b < n_basis_functions; b++)
//      {
//        for (int i = 0; i < DIM; i++)
//        {
//          for (int j = 0; j < DIM; j++)
//          {
//            Ae(DIM * a + i, DIM * b + j) += fe.dN(b, j) * idata_->Cmatrix[i][j] * fe.dN(a, i) * detJxW + fe.dN(b, DIM - j - 1) * idata_->Cmatrix[DIM][DIM] * fe.dN(a, DIM - i - 1) * detJxW;
//          }
//        }
//      }
//    }

for(int a = 0; a < n_basis_functions; a++)
{ // loop over test basis functions (a)
    for(int b = 0; b < n_basis_functions; b++)
    { // loop over trial basis functions (b)
        for(int i = 0; i < DIM; i++) { // loop over DOF test function i

            for (int k = 0; k < DIM; k++) { // loop over DOF trial function k
                for (int l = 0; l < DIM; l++) { // internal sum l
                    for (int j = 0; j < DIM; j++) { // internal sum j
                        Ae(DIM * a + i, DIM * b + k) +=
                                fe.dN(b, l) * Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.dN(a, j) * detJxW * 0.5;
                    } // end j

                } // end l

            } // end k

            /////// symmetric part  //////
            for(int l = 0; l < DIM; l++) { // loop over DOF trial function l
                for (int k = 0; k < DIM; k++) { // internal sum k
                    for (int j = 0; j < DIM; j++) { // internal sum j
                        Ae(DIM * a + i, DIM * b + l) +=
                                fe.dN(b, k) * Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.dN(a, j) * detJxW * 0.5;
                    } // end j

                } // end k
            } // end l


        } // end i
    } // end b

} // end a


#endif

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

#ifndef NDEBUG
      std::ofstream fout("testdomaingp_be.txt", std::ios::app);
#if (DIM == 2)
      fout << fe.position().x() << "," << fe.position().y() << "\n";
#endif
#if (DIM == 3)
      fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << "\n";
#endif
      fout.close();
#endif

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

      // TODO: some integration over the boundary except for the carved out region before this (if needed)
      if (idata_->ibm_geom_def.size() == 0) {
          return;
      }

      bool cantilever_Dirichlet_SBM = false;
      bool cantilever_Neumann_SBM = false;

      if (idata_->SbmGeo == LEInputData::cantilever and (fabs(idata_->DomainMax[0] - fe.position().x()) < 1e-12)) {
          cantilever_Dirichlet_SBM = true;
      }

      if (side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY) {
          // this is for 2 * DIM domain boundaries
          const auto &bc = idata_->boundary_def.at(side_idx);
          if (bc.disp_type == BoundaryDef::Disp_Type::SBM_NEUMANN_WALL) {
              cantilever_Neumann_SBM = true;
          }
      }

#ifndef NDEBUG
      if (cantilever_Dirichlet_SBM){
        std::cout << "Dirichlet normal  = "; fe.surface()->normal().print(); std::cout << "\n";
      }

      if (cantilever_Neumann_SBM){
          std::cout << "Neumann normal  = "; fe.surface()->normal().print(); std::cout << "\n";
      }
#endif


#ifndef NDEBUG
    std::ofstream fout("testsurgp_Ae.txt", std::ios::app);
    fout << fe.position().x() << "," << fe.position().y() << "\n";
    fout.close();
#endif

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

    sbmCalc.Dist2Geo(d);

#ifndef NDEBUG
    std::ofstream fout_dist("sur_dist_" + std::to_string(idata_->mesh_def.refine_lvl_base) + ".txt", std::ios::app);
#if (DIM == 2)
    fout_dist << sqrt(pow(d[0], 2) + pow(d[1], 2)) << "\n";
#endif
#if (DIM == 3)
    fout_dist << sqrt(pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2)) << "\n";
#endif
    fout_dist.close();
#endif

    /// finish: calculate d vector ======================================================================

    /// geometry
    auto &geo = idata_->ibm_geom_def.at(geom_ID);

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

    double Cmatrix[3 * (DIM - 1)][3 * (DIM - 1)];
    std::vector<std::vector<double>> Be(DIM * n_basis_functions);
    /// for middle term => B_T*C_T
    std::vector<std::vector<double>> BeCmatrix(DIM * n_basis_functions);
    double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
    /// for mid2 term => B_T*C_T*n
    std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);

    if (idata_->bccaseType == POSITION_DISPLACEMENT )
      {
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

                  // std::cout << " N_con = " << N_con << "\n";
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
                  } // b loop
              }   // a loop


          return;
      }

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL ||
        cantilever_Neumann_SBM ||
        cantilever_Dirichlet_SBM ||
        /*below only works for lambda = 1*/
            (idata_->RatioGPSBM == 1 and side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY and idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET))
    {


      CalcCmatrix(Cmatrix);
      CalcBe(fe, Be);
      CalcBeCmatrix(fe, Be, Cmatrix, BeCmatrix);
      CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);
      CalcStressDotNormal(fe, BeCmatrix, SurrogateNormalMatrix, StressDotSurrogateNormal); // (DIM*fe.nbf()) * DIM
    }

    //      for (int i = 0;i<DIM*fe.nbf();i++){
    //          for (int j=0;j<DIM;j++){
    //              std::cout << "StressDotSurrogateNormal[i][j]=" << StressDotSurrogateNormal[i][j] << "\n";
    //          }
    //      }

#ifndef NDEBUG
//      std::cout << "idata_->RatioGPSBM = " << idata_->RatioGPSBM << "\n ";
#endif

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL ||
        cantilever_Dirichlet_SBM ||
        /*below only works for lambda = 1*/
        (fabs(idata_->RatioGPSBM-1)< 1e-12 and side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY and idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET))
    {

      /// wall normal distance
      // std::cout<<"h:"<<h<<"\n";
      double weakBCpenaltyParameter = util_funcs::ReturnPenaltyParameters(idata_) * Cb_e / h;

      const ZEROPTV &SurrogateNormal = fe.surface()->normal();
////
////      //////// consisitency term (subtracted in LHS)//////
//        for(int a=0; a<fe.nbf(); a++)
//      { // loop over test basis functions (a)
//        for(int b=0; b<fe.nbf(); b++)
//        { // loop over trial basis functions (b)
//          for(int i=0; i<DIM; i++)
//          { // loop over DOF test function i
//            for(int k=0; k<DIM; k++)
//            { // loop over DOF trial function j
//                for(int j=0; j<DIM; j++)
//                { // loop over internal sum j
//                    for(int l=0; l<DIM; l++)
//                    { // loop over internal sum l
//                        Ae(DIM*a+i, DIM*b+k) -= Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.N(a) * fe.dN(b,l) * SurrogateNormal(j)* detSideJxW * 0.5;
//                    } // end l
//                } // end j
//            } // end k
//            for(int l=0; l<DIM; l++)
//            { // loop over DOF trial function j
//                for(int j=0; j<DIM; j++)
//                { // loop over internal sum j
//                    for(int k =0; k<DIM; k++)
//                    { // loop over internal sum k
//                        Ae(DIM*a+i, DIM*b+l) -= Cmatrix[tensorToVoigt(i,j)][tensorToVoigt(k,l)] * fe.N(a) * fe.dN(b,k) * SurrogateNormal(j)* detSideJxW * 0.5;
//                    } // end k
//
//
//                } // end j
//
//            } // end l
//
//
//          } // end i
//        } // end b
//      } // end a
//
////      //////// adjoint consistency term (added in LHS) //////

        //// it will include two parts to maintain the symmetry of the strain tensor
        //// we are trying to calculate the following term
        //// C_{ijkl}(\del_l wk + \del_k wl)) * \tilde{n}_j * (u_i + \del_i u_p * d_p)
//        for(int a=0; a<fe.nbf(); a++) /// test basis functions a loop begin
//        {
//            for (int b = 0; b < fe.nbf(); b++) /// trial basis functions b loop begin
//            {
//                double grad_u_dot_d = 0.0; //// u + \grad u * d
//                //// \sum_p N_b + \del_p N_b *d_p is same for all the test basis function
//                for (int p = 0; p < DIM; p++) { //// loop over the dimensions in p
//                    grad_u_dot_d += fe.dN(b, p) * d[p];
//                } //// loop over the dimensions in p ends
//
//                for (int k = 0; k < DIM; k++) { /// loop over the dimensions in k for test basis function dof
//                    for (int i = 0; i < DIM; i++) { /// loop over the dimensions in i    for trial basis function dof
//                        for (int j = 0; j < DIM; j++) { /// loop over the dimensions in j for summation
//                            for (int l = 0; l < DIM; l++) /// loop over the dimensions in l for summation
//                            {
//                                Ae(DIM * a + k, DIM * b + i) +=
//                                        Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, l) *
//                                        SurrogateNormal(j) *
//                                        (grad_u_dot_d+fe.N(b)) * detSideJxW *0.5;
//                            } /// l loop ends
//                        } /// j loop ends
//                    } /// i loop ends
//                } /// k loop ends
//                for (int l = 0; l < DIM; l++) { /// loop over the dimensions in l for test basis function
//                    for (int i = 0; i < DIM; i++) { /// loop over the dimensions in i for trial basis function
//                        for (int j = 0; j < DIM; j++) { /// loop over the dimensions in j for summation
//                            for (int k = 0; k < DIM; k++) /// loop over the dimensions in k for summation
//                            {
//                                Ae(DIM * a + l, DIM * b + i) +=
//                                        Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, k) *
//                                        SurrogateNormal(j) * (grad_u_dot_d+fe.N(b))  * detSideJxW * 0.5;
//                            } /// k loop ends
//                        } /// j loop ends
//                    } /// i loop ends
//                } /// l loop ends
//
//            } /// trial basis function b loop ends
//
//
//        } /// test basis function a loop ends

//////      //////// adjoint consistency term (added in LHS) //////
//
//      //// it will include two parts to maintain the symmetry of the strain tensor
//      //// we are trying to calculate the following term
//      //// C_{ijkl}(\del_l wk + \del_k wl)) * \tilde{n}_j * (u_i + \del_i u_p * d_p)
//    for(int a=0; a<fe.nbf(); a++) /// test basis functions a loop begin
//    {
//
//        for (int b = 0; b < fe.nbf(); b++) /// trial basis functions b loop begin
//        {
//            double grad_u_dot_d = 0.0; //// u + \grad u * d
//            //// \sum_p N_b + \del_p N_b *d_p is same for all the test basis function
//            for (int p = 0; p < DIM; p++) { //// loop over the dimensions in p
//                grad_u_dot_d += fe.dN(b, p) * d[p];
//            } //// loop over the dimensions in p ends
//
//                for (int i = 0; i < DIM; i++) { /// loop over the dimensions in i for test basis function dof
//                    for (int k = 0; k < DIM; k++) { /// loop over the dimensions in k for trial basis function dof
//                        for (int j = 0; j < DIM; j++) { /// loop over the dimensions in j for summation
//                            for (int l = 0; l < DIM; l++) /// loop over the dimensions in l for summation
//                            {
//                                Ae(DIM * a + i, DIM * b + k) +=
//                                        Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, j) *
//                                        SurrogateNormal(l) *
//                                                (grad_u_dot_d+fe.N(b)) * detSideJxW *0.5;
//                            } /// l loop ends
//                        } /// j loop ends
//                    } /// k loop ends
//                } /// i loop ends
//                for (int j = 0; j < DIM; j++) { /// loop over the dimensions in l for test basis function
//                    for (int k = 0; k < DIM; k++) { /// loop over the dimensions in i for trial basis function
//                        for (int i = 0; i < DIM; i++) { /// loop over the dimensions in j for summation
//                            for (int l = 0; l < DIM; l++) /// loop over the dimensions in k for summation
//                            {
//                                Ae(DIM * a + j, DIM * b + k) +=
//                                        Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, i) *
//                                        SurrogateNormal(l) * (grad_u_dot_d+fe.N(b))  * detSideJxW * 0.5;
//                            } /// k loop ends
//                        } /// j loop ends
//                    } /// i loop ends
//                } /// l loop ends
//
//        } /// trial basis function b loop ends
//
//
//    } /// test basis function a loop ends


//
//
//      /////// penalty term //////
//
//
//
//    for(int a =0; a<fe.nbf(); a++)
//    {
//        double grad_w_dot_d = 0.0; //// w + \grad w * d
//        for (int p = 0; p < DIM; p++)
//        { //// loop over the dimensions in p
//            grad_w_dot_d +=  fe.dN(a, p) * d[p];
//        } //// loop over the dimensions in p ends
//
//        for(int b=0; b<fe.nbf(); b++)
//        {
//            double grad_u_dot_d = 0.0; //// u + \grad u * d
//
//            for (int p = 0; p < DIM; p++)
//            { //// loop over the dimensions in p
//                grad_u_dot_d +=  fe.dN(b, p) * d[p];
//
//            } //// loop over the dimensions in p ends
//            for(int i=0; i<DIM; i++)
//            {
//
//                    Ae(DIM*a+i, DIM*b+i) += weakBCpenaltyParameter * (fe.N(a)+grad_w_dot_d)* (fe.N(b)+grad_u_dot_d)* detSideJxW;
//            }
//        }
//    }
    //////////////////// The below is the old code for the weak BCs
    //////////////////// based on matrix form enforcing direct symmetry

//
////      ////////////////////////////////
//
      double Ne[DIM][DIM * n_basis_functions];
      memset(Ne, 0.0, sizeof Ne);
      if (IFSBM)
      {
        DENDRITE_REAL secondOrderTerm_a(0);
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
                Ne[dim][DIM * j + dim2] = fe.N(j) + gradWdotd + secondOrderTerm_a;
              }
              else
              {
                Ne[dim][DIM * j + dim2] = 0.0;
              }
            }
          }

          // todo: print Ne to
          //  see whether it is correct

          //          Ne[0][2 * j] = fe.N(j) + gradWdotd + secondOrderTerm_a;
          //          Ne[0][2 * j + 1] = 0;
          //          Ne[1][2 * j] = 0;
          //          Ne[1][2 * j + 1] = fe.N(j) + gradWdotd + secondOrderTerm_a;
        }
      }


      else
      {
        for (int j = 0; j < n_basis_functions; j++)
        {
          for (int dim = 0; dim < DIM; dim++)
          {
            for (int dim2 = 0; dim2 < DIM; dim2++)
            {
              if (dim2 == dim)
              {
                Ne[dim][DIM * j + dim2] = fe.N(j);
              }
              else
              {
                Ne[dim][DIM * j + dim2] = 0.0;
              }
            }
          }
          //          Ne[0][2 * j] = fe.N(j);
          //          Ne[0][2 * j + 1] = 0;
          //          Ne[1][2 * j] = 0;
          //          Ne[1][2 * j + 1] = fe.N(j);
        }
      }
      double Ne_con[DIM][DIM * n_basis_functions];
      memset(Ne_con, 0.0, sizeof Ne_con);

      for (int j = 0; j < n_basis_functions; j++)
      {
        for (int dim = 0; dim < DIM; dim++)
        {
          for (int dim2 = 0; dim2 < DIM; dim2++)
          {
            if (dim2 == dim)
            {
              Ne_con[dim][DIM * j + dim2] = fe.N(j);
            }
            else
            {
              Ne_con[dim][DIM * j + dim2] = 0.0;
            }
          }
        }
        //        Ne_con[0][2 * j] = fe.N(j);
        //        Ne_con[0][2 * j + 1] = 0;
        //        Ne_con[1][2 * j] = 0;
        //        Ne_con[1][2 * j + 1] = fe.N(j);
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
            N += StressDotSurrogateNormal[a][k] * Ne[k][b] * detJxW;
            N_con += StressDotSurrogateNormal[a][k] * Ne_con[k][b] * detJxW;
          }
          // Ae is symmetric

          if (idata_->IfAdjointConsistency)
          {
            // std::cout << " N = " << N << "\n";
            Ae(a, b) += N; // adjoint consistency
          }

          // std::cout << " N_con = " << N_con << "\n";
          Ae(b, a) -= N_con; // consistency
        }
      }

      if (IFSBM)
      {
        DENDRITE_REAL secondOrderTerm_a(0), secondOrderTerm_b(0);

        for (int a = 0; a < fe.nbf(); a++)
        {
#if (DIM == 2)
          if (idata_->elemOrder == 2 && idata_->ifHessian)
          {
            secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                 d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                                2;
          }
          else
          {
            secondOrderTerm_a = 0;
          }
#endif

#if (DIM == 3)
          if (idata_->elemOrder == 2 && idata_->ifHessian)
          {

            secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
          }
          else
          {
            secondOrderTerm_a = 0;
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
              secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1]) +
                                   d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1])) /
                                  2;
            }
            else
            {
              secondOrderTerm_b = 0;
            }
#endif

#if (DIM == 3)
            if (idata_->elemOrder == 2 && idata_->ifHessian)
            {

              secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1] + fe.d2N(b, 0, 2) * d[2]) + d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1] + fe.d2N(b, 1, 2) * d[2]) + d[2] * (fe.d2N(b, 2, 0) * d[0] + fe.d2N(b, 2, 1) * d[1] + fe.d2N(b, 2, 2) * d[2])) / 2;
            }
            else
            {
              secondOrderTerm_b = 0;
            }
#endif
            double gradUdotd = 0.0;
            for (int k = 0; k < DIM; k++)
            {
              gradUdotd += fe.dN(b, k) * d[k];
            }

            for (int j = 0; j < DIM; j++)
            {
              Ae(DIM * a + j, DIM * b + j) +=
                  +weakBCpenaltyParameter * (fe.N(a) + gradWdotd) * (fe.N(b) + gradUdotd) * detSideJxW; // penalty
            }
          } // b loop
        }   // a loop
      }
      else
      {
        for (int a = 0; a < fe.nbf(); a++)
        {
          for (int b = 0; b < fe.nbf(); b++)
          {
            /*
             *  make it j to match with what's inside Navier-Stokes
             */
            for (int j = 0; j < DIM; j++)
            {
              Ae(DIM * a + j, DIM * b + j) +=
                  +weakBCpenaltyParameter * (fe.N(a)) * (fe.N(b)) * detSideJxW; // penalty
            }
          } // b loop
        }   // a loop
      }

        //////////////////// The below is the old code for the weak BCs
        //////////////////// based on matrix form enforcing direct symmetry

//
////      ////////////////////////////////
    }

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL || cantilever_Neumann_SBM)
    {

      const ZEROPTV &p = fe.position();
      const ZEROPTV &SurrogateNormal = fe.surface()->normal();

      ZEROPTV TrueNormal;

      if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL) {
          double x_mid = geo.InitialDisplacement.x();
          double y_mid = geo.InitialDisplacement.y();

          double x = p.x();
          double y = p.y();

          double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
          double sin_x = (x - x_mid) / radius;
          double sin_y = (y - y_mid) / radius;
          /// temporary using this -> need to implement SBMCalc inside!
          TrueNormal = ZEROPTV{sin_x, sin_y, 0.0};
      } else if (cantilever_Neumann_SBM){
          sbmCalc.NormalofGeo(TrueNormal,d);
      }

      double TrueNormalMatrix[DIM][3 * (DIM - 1)];
      CalcTrueNormalMatrix(TrueNormal, TrueNormalMatrix);

      // for mid2 term => B_T*C_T*n
      std::vector<std::vector<double>> StressDotTrueNormal(DIM * n_basis_functions);
      CalcStressDotNormal(fe, BeCmatrix, TrueNormalMatrix, StressDotTrueNormal);

      double SurrogateDotTrueNormal = SurrogateNormal.innerProduct(TrueNormal);


        std::vector<std::vector<double>> gradStrainDotDistance(DIM * n_basis_functions);
            CalcCmatrix(Cmatrix);
    //            Strain Dot Distance Matrix
            CalcGradStrainDotDistance(fe, gradStrainDotDistance, d);

    //            multiply Cmatrix with gradStrainDotDistance
            std::vector<std::vector<double>> CMatrixGradStrainDotDistance(DIM * n_basis_functions);

    //            just use THIS FUNCTION TO MULTIPLY THE TWO MAT IT'S ALREADY IMPLEMENTED
            CalcBeCmatrix(fe, gradStrainDotDistance, Cmatrix, CMatrixGradStrainDotDistance);
    //            multiply CMatrixGradStrainDotDistance with normal
            CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);
    //            multiply CMatrixGradStrainDotDistance with normal
            std::vector<std::vector<double>> Hessian_Dot_Normal(DIM * n_basis_functions);
    //            Hessian_Dot_Normal is DIM * (n_basis_functions*DIM)
            CalcStressDotNormal(fe, CMatrixGradStrainDotDistance, TrueNormalMatrix, Hessian_Dot_Normal);


    //            get the basis matrix , we are just doing this to make our life easier
            std::vector<std::vector<double>> basisMatrix(DIM * n_basis_functions);

      // for Final term => B_T*C_T*n*N

      const double detJxW = fe.detJxW();

//        for(int a = 0; a <n_basis_functions; a++){ // test function
//            for(int b=0; b<n_basis_functions; b++) { // trial function
//                for (int i = 0; i < DIM; i++) { // dof test function
//                    for (int k = 0; k < DIM; k++) { // dof trial function
//                        for (int j = 0; j < DIM; j++) {
//                            for (int l = 0; l < DIM; l++) {
//                                for (int q = 0; q < DIM; q++) { // inner loop to Hessian dot Distance q loop
//                                    Ae(a * DIM + i, b * DIM + k) +=
//                                            fe.N(a) * getCIJKLPlaneStress(i,j,k,l)*
//                                            /////  hessian(b,q,l) * distance(q) * n(j) * (n dot n tilde)
//                                            ( fe.d2N(b, l, q) *d[q] * TrueNormal(j) *SurrogateDotTrueNormal*0 +
//                                              ///// dN(b,l)* (n dot n tilde) * n(j)
//                                              fe.dN(b,l) * SurrogateDotTrueNormal * TrueNormal(j)
//                                              ///// dN(b,l)* tilde n(j)
//                                              -fe.dN(b,l) * SurrogateNormal[j] ) *detSideJxW;
//                                } // end q loop
//                            } // end l loop
//
//
//                        } // end J loop
//
//
//                    } // end trial function K loop in DOF
//                } // end test function I loop in DOF
//            } // end trial function b loop
//        } // end test function a loop



////     save Ae matrix in a file Ae.txt
//    std::ofstream fout("Ae.txt", std::ios::app);
//    for(int a = 0; a < DIM * n_basis_functions; a++){
//        for(int b=0; b<DIM * n_basis_functions; b++){
//            fout << Ae(a,b) << ",";
//        }
//        fout << "\n";
//    }
//    fout.close();
//// save stress dot surrogate normal
//
//for(int b=0; b<n_basis_functions; b++){
//    for(int k=0; k<DIM; k++){
//        double sigma_dot_n = 0;
//        double sigma_dot_n_tilde = 0;
//        for(int i=0; i<DIM; i++){
//            for(int j=0; j<DIM; j++){
//                for(int l=0; l<DIM; l++){
//                    sigma_dot_n+= getCIJKLPlaneStress(i,j,k,l)*fe.dN(b,l)*SurrogateNormal[i];
//                    sigma_dot_n_tilde+= getCIJKLPlaneStress(i,j,k,l)*fe.dN(b,l)*TrueNormal[i];
//                }
//            }
//        }
//
//        std::ofstream fout2("stress_dot_surrogate_normal.txt", std::ios::app);
//        fout2 << sigma_dot_n << ",\n";
//        fout2.close();
//
//        std::ofstream fout3("stress_dot_true_normal.txt", std::ios::app);
//        fout3 << sigma_dot_n_tilde << ",\n";
//        fout3.close();
//    }
//}
//
//
//    std::ofstream fout1("stress_dot_surrogate_normal_original.txt", std::ios::app);
//    for(int a = 0; a < DIM * n_basis_functions; a++){
//        for(int b=0; b<DIM; b++){
//            fout1 << StressDotSurrogateNormal[a][b] << ",";
//        }
//        fout1 << "\n";
//    }
//    fout1.close();
//
//    // save stress dot true normal
//    std::ofstream fout2("stress_dot_true_normal_original.txt", std::ios::app);
//    for(int a = 0; a < DIM * n_basis_functions; a++){
//        for(int b=0; b<DIM; b++){
//            fout2 << StressDotTrueNormal[a][b] << ",";
//        }
//        fout2 << "\n";
//    }
//    fout2.close();
//
//
//        exit(1);
//

      for (int a = 0; a < n_basis_functions; a++) // test function
      {
        for (int b = 0; b < n_basis_functions; b++) // trial function
        {
          for (int i = 0; i < DIM; i++)
          {
            for (int j = 0; j < DIM; j++)
            {

              Ae(a * DIM + i, b * DIM + j) -= fe.N(a) * StressDotSurrogateNormal[b * DIM + j][i] * detJxW;
              if(idata_->ifHessian)
              {
                  std::cout<<"Hessian is implemented"<<std::endl;
                  Ae(a * DIM + i, b * DIM + j) += fe.N(a) * Hessian_Dot_Normal[b * DIM + j][i] * SurrogateDotTrueNormal * detJxW; // Hessian Correction

              }
              Ae(a * DIM + i, b * DIM + j) += fe.N(a) * StressDotTrueNormal[b * DIM + j][i] * SurrogateDotTrueNormal * detJxW; // Area Correction
            }
          }
        }
      }



//      /// implement penalty
//      for (int a = 0; a < DIM * n_basis_functions; a++) // test function
//      {
//        for (int b = 0; b < DIM * n_basis_functions; b++) // trial function
//        {
//          for (int k = 0; k < DIM; k++)
//          {
//            Ae(a, b) += TauN * StressDotTrueNormal[a][k] * StressDotTrueNormal[b][k] * pow(SurrogateDotTrueNormal, 2) * detJxW; // penalty
//          }
//        }
//      }
      return;
    }

  }

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT side_idx, const DENDRITE_UINT id, TALYFEMLIB::ZEROARRAY<double> &be)
  {

//      std::cout << "be, position = "; fe.position().print(); std::cout << "\n";

#ifndef NDEBUG
      std::ofstream fout("testsurgp_be.txt", std::ios::app);
      fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << "\n";
      fout.close();
#endif

      /*
       * The reason that we implement below part is because that when we use side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY
       * we cannot get the GPs on the x_plus wall.
       */
      bool cantilever_Dirichlet_SBM = false;
      bool cantilever_Neumann_SBM = false;
      if (idata_->SbmGeo == LEInputData::cantilever and (fabs(idata_->DomainMax[0]-fe.position().x()) < 1e-12)){
          cantilever_Dirichlet_SBM = true;
      }


    /*
     * NOTE: this is essential
     * sides except for carved-out -> we need to implement neumann BC
     */
    if (side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY and idata_->bccaseType!= POSITION_DISPLACEMENT)
    {

#ifndef NDEBUG
        std::ofstream fout("testsurgp_be_wall.txt", std::ios::app);
        fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << "\n";
        fout.close();

        if (idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::NEUMANN) {
            std::cout << "idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::NEUMANN = "
                      << (idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::NEUMANN) << "\n";
        }

//        std::cout << "side_idx = "<< side_idx << "\n";

//        std::cout << "idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET = "
//                  << (idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET) << "\n";

#endif

        // this is for 2 * DIM domain boundaries
      const auto &bc = idata_->boundary_def.at(side_idx);
      if (bc.disp_type == BoundaryDef::Disp_Type::NEUMANN)
      {

          if (idata_->bccaseType== CSV_FORCE){

              ZEROPTV traction  = computeTraction(fe, idata_);
              for (int a = 0; a < fe.nbf(); a++) {
                  //  [bug fix] remove normal term
                  for (int dim = 0; dim < DIM; dim++) {
                      be(DIM * a + dim) += fe.N(a) * traction(dim) * fe.detJxW();
                  }
              }
          } else {
              for (int a = 0; a < fe.nbf(); a++) {
                  //  [bug fix] remove normal term
                  for (int dim = 0; dim < DIM; dim++) {
                      be(DIM * a + dim) += fe.N(a) * bc.Traction(dim) * fe.detJxW();
                  }
              }
          }

      }else if (bc.disp_type == BoundaryDef::Disp_Type::SBM_NEUMANN_WALL){
          cantilever_Neumann_SBM = true;
      }
    }

    if (idata_->ibm_geom_def.size() == 0)
    {
      return;
    }

    assert(method == IBM_METHOD::SBM);

    double h = ElementSize(fe);
    double d[DIM];
    int geom_ID;

    bool DirichletHaveSet = false;
    double DirichletBCValue[DIM];
#ifdef DEEPTRACE
      //// this is of no use
      std::vector<my_kd_tree_t*> kd_tree_;
      //// add the kd_tree_ to the pointer

      SBMCalcDeepTrace sbmCalc(fe, idata_, imga_, kd_tree_, 0 /*CarvedOutGeomID*/);
#else
    SBMCalc sbmCalc(fe, idata_, imga_, kd_tree_);
#endif

    sbmCalc.Dist2Geo(d);
    sbmCalc.GetDirichletBC(d, DirichletBCValue, DirichletHaveSet);

    /// finish: calculate d vector ======================================================================

    /// geometry
    auto &geo = idata_->ibm_geom_def.at(geom_ID);

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

      if (idata_->bccaseType == POSITION_DISPLACEMENT){

          CalcCmatrix(Cmatrix);
          CalcBe(fe, Be);
          CalcBeCmatrix(fe, Be, Cmatrix, BeCmatrix);
          double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
          CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);

          // for mid2 term => B_T*C_T*n
          std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);
          CalcStressDotNormal(fe, BeCmatrix, SurrogateNormalMatrix, StressDotSurrogateNormal);

          double weakBCpenaltyParameter_ = util_funcs::ReturnPenaltyParameters(idata_) * Cb_e / h;

          // Dirichlet
          ZEROPTV D_wall_;

          // x-dir
          for (int dim = 0; dim < DIM; dim++) {
                  if (DirichletHaveSet) {
                      D_wall_(dim) = DirichletBCValue[dim];
                  }
          }


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
                          StressDotSurrogateNormal[a][dim] * D_wall_(dim) * detSideJxW; // adjoint
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
                          +weakBCpenaltyParameter_ * (fe.N(a) + gradWdotd + secondOrderTerm_a_) * D_wall_(i) * detSideJxW; // penalty
              }
          }
          return;
      }

      if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL ||
        cantilever_Dirichlet_SBM ||
        cantilever_Neumann_SBM ||
        /*below only works for lambda = 1*/
        (fabs(idata_->RatioGPSBM-1)< 1e-12 and side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY and idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET))
    {
      CalcCmatrix(Cmatrix);
      CalcBe(fe, Be);
      CalcBeCmatrix(fe, Be, Cmatrix, BeCmatrix);
    }

    // x and y => below vector only has one value
    ZEROPTV bc_value;
#if (DIM == 2)
    bc_value = ZEROPTV{geo.getBC_D(0)[0], geo.getBC_D(1)[0], 0.0};
#endif
#if (DIM == 3)
    bc_value = ZEROPTV{geo.getBC_D(0)[0], geo.getBC_D(1)[0], geo.getBC_D(2)[0]};
#endif

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL ||
        cantilever_Dirichlet_SBM ||
        cantilever_Neumann_SBM ||
        /*below only works for lambda = 1*/
        (idata_->RatioGPSBM == 1 and side_idx < BoundaryTypes::MAX_WALL_TYPE_BOUNDARY and idata_->boundary_def.at(side_idx).disp_type == BoundaryDef::Disp_Type::SBM_ZERO_DIRICHLET))
    {
      double SurrogateNormalMatrix[DIM][3 * (DIM - 1)];
      CalcSurrogateNormalMatrix(fe, SurrogateNormalMatrix);

      // for mid2 term => B_T*C_T*n
      std::vector<std::vector<double>> StressDotSurrogateNormal(DIM * n_basis_functions);
      CalcStressDotNormal(fe, BeCmatrix, SurrogateNormalMatrix, StressDotSurrogateNormal);

      //      std::cout<<"dirichlet\n";

      // Dirichlet
      ZEROPTV D_wall;
      double DR;

      // SBM_ZERO_DIRICHLET just keep using D_wall.
        if (!cantilever_Dirichlet_SBM) {

            /// put into "be"
            // x-dir
            for (int dim = 0; dim < DIM; dim++) {
                if (geo.bc_type_D[dim] == IBMGeomDef::BCType::DIRICHLET) {
                    if (DirichletHaveSet) {
                        D_wall(dim) = DirichletBCValue[dim];
                        //                 std::cout << "DirichletBCValue[dim] = " << DirichletBCValue[dim] << "\n";
                    } else {
                        D_wall(dim) = bc_value(dim);
                    }
                }
            }

            if (geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL) {
                DR = bc_value(0);
            }
        }

      double weakBCpenaltyParameter = util_funcs::ReturnPenaltyParameters(idata_) * Cb_e / h;

      DENDRITE_REAL secondOrderTerm_a(0);
      for (int a = 0; a < fe.nbf(); a++)
      {
#if (DIM == 2)
        if (idata_->elemOrder == 2 && idata_->ifHessian)
        {
          secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                               d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                              2;
        }
        else
        {
          secondOrderTerm_a = 0;
        }
#endif

#if (DIM == 3)
        if (idata_->elemOrder == 2 && idata_->ifHessian)
        {

          secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
        }
        else
        {
          secondOrderTerm_a = 0;
        }
#endif
      }

      if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET
          or cantilever_Dirichlet_SBM)
      {
          const ZEROPTV &SurrogateNormal = fe.surface()->normal();

          ///// we imlement the adjoint consistency term
          //// \int_{\Gamma} -ud Cijkl (dw_i/dx_j+dw_j/dx_i) \tilde{n_l} dS

//          for(int a =0; a<n_basis_functions; a++){
//              for(int i=0; i<DIM; i++){
//                  for(int j=0; j<DIM; j++){
//                      for(int k=0; k<DIM; k++){
//                          for(int l=0; l<DIM; l++){
//                              be(a*DIM+i) -= D_wall[k]* Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, j)  * SurrogateNormal[l] * detSideJxW*0.5;
//                          }
//                      }
//                  }
//              }
//              for(int j=0; j<DIM; j++){
//                    for(int i=0; i<DIM; i++){
//                        for(int k=0; k<DIM; k++){
//                            for(int l=0; l<DIM; l++){
//                                be(a*DIM+j) -= D_wall[k]* Cmatrix[tensorToVoigt(j, k)][tensorToVoigt(i, l)] * fe.dN(a, i)  * SurrogateNormal[l] * detSideJxW*0.5;
//                            }
//                        }
//                    }
//              }
//
//          }
//
//          for (int a = 0; a < n_basis_functions; a++)
//          {
//            double gradWdotd = 0.0;
//            for (int k = 0; k < DIM; k++)
//            {
//              gradWdotd += fe.dN(a, k) * d[k];
//            }
//            for (int i = 0; i < DIM; i++)
//            {
//              be(DIM * a + i) +=
//                  +weakBCpenaltyParameter * (fe.N(a) + gradWdotd) * D_wall(i) * detSideJxW; // penalty
//            }
//          }


        if (idata_->IfAdjointConsistency)
        {

          for (int a = 0; a < DIM * n_basis_functions; a++)
          {
            for (int dim = 0; dim < DIM; dim++)
            {
              be(a) +=
                  StressDotSurrogateNormal[a][dim] * D_wall(dim) * detSideJxW; // adjoint
            }
          }
        }

        if (IFSBM)
        {
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
                  +weakBCpenaltyParameter * (fe.N(a) + gradWdotd + secondOrderTerm_a) * D_wall(i) * detSideJxW; // penalty
            }
          }
        }
        else
        {
          for (int a = 0; a < n_basis_functions; a++)
          {
            for (int i = 0; i < DIM; i++)
            {
              be(DIM * a + i) +=
                  +weakBCpenaltyParameter * (fe.N(a)) * D_wall(i) * detSideJxW; // penalty
            }
          }
        }
      }

      // r-dir
      if (geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
      {
        double x_mid = geo.InitialDisplacement.x();
        double y_mid = geo.InitialDisplacement.y();
        const ZEROPTV &p = fe.position();
        double x = p.x();
        double y = p.y();

        double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
        ZEROPTV sin_func;
        sin_func.x() = (x - x_mid) / radius;
        sin_func.y() = (y - y_mid) / radius;


          const ZEROPTV &SurrogateNormal = fe.surface()->normal();

          ///// we imlement the adjoint consistency term
          //// \int_{\Gamma} -ud Cijkl (dw_i/dx_j+dw_j/dx_i) \tilde{n_l} dS

//          for(int a =0; a<n_basis_functions; a++){
//              for(int i=0; i<DIM; i++){
//                  for(int j=0; j<DIM; j++){
//                      for(int k=0; k<DIM; k++){
//                          for(int l=0; l<DIM; l++){
//                              be(a*DIM+i) -= DR * sin_func(k)* Cmatrix[tensorToVoigt(i, j)][tensorToVoigt(k, l)] * fe.dN(a, j)  * SurrogateNormal[l] * detSideJxW*0.5;
//                          }
//                      }
//                  }
//              }
//              for(int j=0; j<DIM; j++){
//                  for(int i=0; i<DIM; i++){
//                      for(int k=0; k<DIM; k++){
//                          for(int l=0; l<DIM; l++){
//                              be(a*DIM+j) -= DR * sin_func(k)* Cmatrix[tensorToVoigt(j, k)][tensorToVoigt(i, l)] * fe.dN(a, i)  * SurrogateNormal[l] * detSideJxW*0.5;
//                          }
//                      }
//                  }
//              }
//
//          }
//
//          for (int a = 0; a < n_basis_functions; a++)
//          {
//              double gradWdotd = 0.0;
//              for (int k = 0; k < DIM; k++)
//              {
//                  gradWdotd += fe.dN(a, k) * d[k];
//              }
//              for (int i = 0; i < DIM; i++)
//              {
//                  be(DIM * a + i) +=
//                          +weakBCpenaltyParameter * (fe.N(a) + gradWdotd) * DR * sin_func(i)* detSideJxW; // penalty
//              }
//          }

        if (idata_->IfAdjointConsistency)
        {
          for (int a = 0; a < DIM * n_basis_functions; a++)
          {
            for (int dim = 0; dim < DIM; dim++)
            {
              be(a) += StressDotSurrogateNormal[a][dim] * DR * sin_func[dim] * detSideJxW; // adjoint
            }
          }
        }

        if (IFSBM)
        {
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
                  +weakBCpenaltyParameter * (fe.N(a) + gradWdotd + secondOrderTerm_a) * DR * sin_func(i) * detSideJxW; // penalty
            }
          }
        }

        else
        {
          for (int a = 0; a < n_basis_functions; a++)
          {
            for (int i = 0; i < DIM; i++)
            {
              be(DIM * a + i) +=
                  +weakBCpenaltyParameter * (fe.N(a)) * DR * sin_func(i) * detSideJxW; // penalty
            }
          }
        }
      }
    }

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::T_RADIAL)
    {
      for (int a = 0; a < fe.nbf(); a++)
      {

        const ZEROPTV &p = fe.position();
        const ZEROPTV &normal = fe.surface()->normal();

        double TR = bc_value(0);

        for (int dim = 0; dim < DIM; dim++)
        {
          be(DIM * a + dim) -= fe.N(a) * TR * normal(dim) * fe.detJxW();
        }
      }
    }

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL || cantilever_Neumann_SBM)
    {
      double x_mid = geo.InitialDisplacement.x();
      double y_mid = geo.InitialDisplacement.y();
      const ZEROPTV &p = fe.position();
      double x = p.x();
      double y = p.y();

      for (int a = 0; a < fe.nbf(); a++)
      {

        const ZEROPTV &SurrogateNormal = fe.surface()->normal();

          ZEROPTV TrueNormal;

          if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL) {
              double x_mid = geo.InitialDisplacement.x();
              double y_mid = geo.InitialDisplacement.y();

              double x = p.x();
              double y = p.y();

              double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
              double sin_x = (x - x_mid) / radius;
              double sin_y = (y - y_mid) / radius;
              /// temporary using this -> need to implement SBMCalc inside!
              TrueNormal = ZEROPTV{sin_x, sin_y, 0.0};
          } else if (cantilever_Neumann_SBM){
#ifndef NDEBUG
              std::cout << "calc true normal \n";
#endif
              sbmCalc.NormalofGeo(TrueNormal,d);
          }

        double TrueNormalMatrix[DIM][3 * (DIM - 1)];
        CalcTrueNormalMatrix(TrueNormal, TrueNormalMatrix);

        // for mid2 term => B_T*C_T*n
        std::vector<std::vector<double>> StressDotTrueNormal(DIM * n_basis_functions);
        CalcStressDotNormal(fe, BeCmatrix, TrueNormalMatrix, StressDotTrueNormal);

        double SurrogateDotTrueNormal = SurrogateNormal.innerProduct(TrueNormal);

        double TR = bc_value(0);

//        for(int dim=0; dim<DIM; dim++)
//        {
//                            be(DIM * a + dim)+=fe.N(a) * (TR * TrueNormal(dim)) * SurrogateDotTrueNormal * fe.detJxW(); // Area Correction
//        }

        for (int dim = 0; dim < DIM; dim++)
        {
          // be(DIM * a + dim) += fe.N(a) * (TR * SurrogateNormal(dim)) * SurrogateDotTrueNormal * fe.detJxW(); // Area Correction

          if (cantilever_Neumann_SBM){
#ifndef NDEBUG
              std::cout << "dim , idata_->boundary_def.at(side_idx).Traction(dim) = " << dim << ", " << idata_->boundary_def.at(side_idx).Traction(dim) << "\n ";
              std::cout << "SurrogateDotTrueNormal = " << SurrogateDotTrueNormal << "\n";
#endif

              be(DIM * a + dim) += fe.N(a) * ((-idata_->boundary_def.at(side_idx).Traction(dim)) * TrueNormal(dim)) * SurrogateDotTrueNormal * fe.detJxW(); // Area Correction

          } else if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL) {

              be(DIM * a + dim) +=
                      fe.N(a) * (TR * TrueNormal(dim)) * SurrogateDotTrueNormal * fe.detJxW(); // Area Correction
              // be(DIM * a + dim) += fe.N(a) * TR * SurrogateDotTrueNormal * fe.detJxW(); // Area Correction
              // be(DIM * a + dim) += fe.N(a) * TR * SurrogateNormal(dim) * fe.detJxW(); // Area Correction
          }

          for (int k = 0; k < DIM; k++)
          {
              if (cantilever_Neumann_SBM){
                  be(DIM * a + dim) += TauN * StressDotTrueNormal[DIM * a + dim][k] * ((-idata_->boundary_def.at(side_idx).Traction(dim)) * TrueNormal(k)) *
                                       pow(SurrogateDotTrueNormal, 2) * fe.detJxW();
              } else if (geo.bc_type_D[0] == IBMGeomDef::BCType::SBM_NEUMANN_RADIAL) {
                  be(DIM * a + dim) += TauN * StressDotTrueNormal[DIM * a + dim][k] * (TR * TrueNormal(k)) *
                                       pow(SurrogateDotTrueNormal, 2) * fe.detJxW();
              }
          }
        }
      }
    }

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::NEUMANN){
        if (idata_->bccaseType== CSV_FORCE){
#ifndef NDEBUG
            std::cout << "CSV_FORCE on carved-out places\n";
#endif
            ZEROPTV traction  = computeTraction(fe, idata_);
            for (int a = 0; a < fe.nbf(); a++) {
                //  [bug fix] remove normal term
                for (int dim = 0; dim < DIM; dim++) {
                    be(DIM * a + dim) += fe.N(a) * traction(dim) * fe.detJxW();
                }
            }
        }
    }

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

//    case LEInputData::SBMGeo::CIRCLE:
//    {
//      ///
//      double pi = M_PI;
//      double E = idata_->planeStress.young;
//      double poisson = idata_->planeStress.poisson;
//
//      ///
//      double x = p.x();
//      double y = p.y();
//
//      ///
//      BodyForce.x() = (E * pi * pi * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / 210 + (pi * pi * cos(pi * x) * sin(pi * y)) / 10)) / (2 * poisson + 2) + (E * poisson * pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / (210 * (pow(poisson, 2) - 1));
//      BodyForce.y() = (E * poisson * pi * pi * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1)) - (E * pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / (90 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos(pi * y) * sin(pi * x)) / 10 - (pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / 490)) / (2 * poisson + 2);
//      ForceHaveSet = true;
//      break;
//    }

    case LEInputData::SBMGeo::CIRCLE:
    case LEInputData::SBMGeo::STAR:
    {
      ///
      double pi = M_PI;
      double E = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;

      ///
      double x = p.x();
      double y = p.y();

      ///
      BodyForce.x() = (E * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (5 * (2 * poisson + 2)) - (E * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1)) - (E * poisson * pow(pi, 2) * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1));
      BodyForce.y() = (E * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (5 * (2 * poisson + 2)) - (E * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1)) - (E * poisson * pow(pi, 2) * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1));
      ForceHaveSet = true;
      break;
    }

    case LEInputData::SBMGeo::ROTATE:
    {
      ///
      double pi = M_PI;
      double E = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;

      ///
      double x_rotate = p.x() - 0.5 - 0.21;
      double y_rotate = p.y() - 0.5 - 0.21;

      /// rotate back to [0,0] to [1,1]
      double x = cos(M_PI / 4) * x_rotate - sin(M_PI / 4) * y_rotate + 0.5;
      double y = sin(M_PI / 4) * x_rotate + cos(M_PI / 4) * y_rotate + 0.5;

      // std::cout<< "x,y = " << x << "," << y << "\n";

      ///
      double BodyForceXp = (E * pi * pi * cos(pi * x) * sin(pi * y)) / (10 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / 210 + (pi * pi * cos(pi * x) * sin(pi * y)) / 10)) / (2 * poisson + 2) + (E * poisson * pi * pi * cos((pi * x) / 7) * cos((pi * y) / 3)) / (210 * (pow(poisson, 2) - 1));
      double BodyForceYp = (E * poisson * pi * pi * cos(pi * y) * sin(pi * x)) / (10 * (pow(poisson, 2) - 1)) - (E * pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / (90 * (pow(poisson, 2) - 1)) - (E * ((pi * pi * cos(pi * y) * sin(pi * x)) / 10 - (pi * pi * sin((pi * x) / 7) * sin((pi * y) / 3)) / 490)) / (2 * poisson + 2);

      ///
      BodyForce.x() = (BodyForceXp * cos(M_PI / 4) + BodyForceYp * cos(M_PI / 4));
      BodyForce.y() = (BodyForceYp * sin(M_PI / 4) - BodyForceXp * sin(M_PI / 4));
      ForceHaveSet = true;
      break;
    }

#if (DIM == 3)
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
#endif

    case LEInputData::SBMGeo::NONE:
    {
      ForceHaveSet = false;
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

  void CalcTrueNormalMatrix(const ZEROPTV TrueNormal, double (&TrueNormalMatrix)[DIM][3 * (DIM - 1)])
  {
#if (DIM == 2)
    TrueNormalMatrix[0][0] = TrueNormal.x();
    TrueNormalMatrix[1][0] = 0;
    TrueNormalMatrix[0][1] = 0;
    TrueNormalMatrix[1][1] = TrueNormal.y();
    TrueNormalMatrix[0][2] = TrueNormal.y();
    TrueNormalMatrix[1][2] = TrueNormal.x();
#endif
#if (DIM == 3)
    TrueNormalMatrix[0][0] = TrueNormal(0);
    TrueNormalMatrix[1][0] = 0;
    TrueNormalMatrix[2][0] = 0;
    TrueNormalMatrix[0][1] = 0;
    TrueNormalMatrix[1][1] = TrueNormal(1);
    TrueNormalMatrix[2][1] = 0;
    TrueNormalMatrix[0][2] = 0;
    TrueNormalMatrix[1][2] = 0;
    TrueNormalMatrix[2][2] = TrueNormal(2);

    TrueNormalMatrix[0][3] = TrueNormal(1);
    TrueNormalMatrix[1][3] = TrueNormal(0);
    TrueNormalMatrix[2][3] = 0;
    TrueNormalMatrix[0][4] = TrueNormal(2);
    TrueNormalMatrix[1][4] = 0;
    TrueNormalMatrix[2][4] = TrueNormal(0);
    TrueNormalMatrix[0][5] = 0;
    TrueNormalMatrix[1][5] = TrueNormal(2);
    TrueNormalMatrix[2][5] = TrueNormal(1);

#endif
  }

  void CalcGradStrainDotDistance(const TALYFEMLIB::FEMElm &fe, std::vector<std::vector<double>> &GradStressDotDistance, const double (&distance)[DIM])
    {
        // the shape of this matrix is similar as the BeCmatrix
//        this is the number of rows
        assert(GradStressDotDistance.size() == (DIM) * fe.nbf());
//        for each row have columns = 3*(DIM-1)

#if (DIM==2)
        for(int a=0; a<fe.nbf(); a++)
        {
            GradStressDotDistance[2*a].resize(3);  // 3 is the number of columns in BMatrix
            GradStressDotDistance[2*a+1].resize(3); // 3 is the number of columns in BMatrix
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y : for first column
            GradStressDotDistance[2*a][0] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1];
            GradStressDotDistance[2*a][1] = 0;
            GradStressDotDistance[2*a][2] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1];
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y : for second column
            GradStressDotDistance[2*a+1][0] = 0;
            GradStressDotDistance[2*a+1][1] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1];
            GradStressDotDistance[2*a+1][2] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1];

        }
#endif
#if (DIM==3)
        for(int a=0; a<fe.nbf(); a++)
    {
//        print the value of the second derivative of the shape function



        GradStressDotDistance[3*a].resize(6);  // 6 is the number of columns in BMatrix
        GradStressDotDistance[3*a+1].resize(6); // 6 is the number of columns in BMatrix
        GradStressDotDistance[3*a+2].resize(6); // 6 is the number of columns in BMatrix
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for first column
        GradStressDotDistance[3*a][0] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a][1] = 0;
        GradStressDotDistance[3*a][2] = 0;
        GradStressDotDistance[3*a][3] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1] + fe.d2N(a,1,2)*distance[2];
        GradStressDotDistance[3*a][4] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
        GradStressDotDistance[3*a][5] = 0;
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for second column
        GradStressDotDistance[3*a+1][0] = 0;
        GradStressDotDistance[3*a+1][1] = fe.d2N(a,1,0)*distance[0] + fe.d2N(a,1,1)*distance[1] + fe.d2N(a,1,2)*distance[2];
        GradStressDotDistance[3*a+1][2] = 0;
        GradStressDotDistance[3*a+1][3] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a+1][4] = 0;
        GradStressDotDistance[3*a+1][5] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
//      grad_wrt_x * distance_x + grad_wrt_y * distance_y + grad_wrt_z * distance_z : for third column
        GradStressDotDistance[3*a+2][0] = 0;
        GradStressDotDistance[3*a+2][1] = 0;
        GradStressDotDistance[3*a+2][2] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];
        GradStressDotDistance[3*a+2][3] = 0;
        GradStressDotDistance[3*a+2][4] = fe.d2N(a,0,0)*distance[0] + fe.d2N(a,0,1)*distance[1] + fe.d2N(a,0,2)*distance[2];
        GradStressDotDistance[3*a+2][5] = fe.d2N(a,2,0)*distance[0] + fe.d2N(a,2,1)*distance[1] + fe.d2N(a,2,2)*distance[2];

    }
#endif
// Check  if any of the value is nan

    }

#pragma mark generating terms for integrations


  double getCIJKLPlaneStress(int I, int J, int K, int L) {
        // Plane stress stiffness constants
        double E = idata_->planeStress.young;
        double nu = idata_->planeStress.poisson;
        double factor = E / (1 - nu * nu); // Common factor for plane stress
        double C_IJKL = 0.0;              // Initialize the tensor value

        // Handle the in-plane components (x=0, y=1)
        if (I == 0 && J == 0 && K == 0 && L == 0) { // C[xxxx]
            C_IJKL = factor;
        } else if (I == 1 && J == 1 && K == 1 && L == 1) { // C[yyyy]
            C_IJKL = factor;
        } else if ((I == 0 && J == 0 && K == 1 && L == 1) ||
                   (I == 1 && J == 1 && K == 0 && L == 0)) { // C[xx][yy] or C[yy][xx]
            C_IJKL = factor * nu;
        } else if ((I == 0 && J == 1 && K == 0 && L == 1) ||
                   (I == 1 && J == 0 && K == 1 && L == 0)) { // C[xy][xy] or C[yx][yx]
            C_IJKL = factor * (1 - nu) / 2.0;
        } else {
            // All other terms are zero in plane stress
            C_IJKL = 0.0;
        }

        return C_IJKL;
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

  ZEROPTV computeTraction(const TALYFEMLIB::FEMElm& fe, const LEInputData* idata_) {
        assert(idata_ != nullptr);  // Make sure input data is not null

        size_t num_results = 2;  // Searching for 2 nearest neighbors
        std::vector<uint32_t> ret_index(num_results);
        std::vector<double> out_dist_sqr(num_results);

#if (DIM == 2)
        const double query_pt[3] = {fe.position().x() - idata_->CsvForce.shift_pos[0], fe.position().y() - idata_->CsvForce.shift_pos[1], 0.0};
#endif

#if (DIM == 3)
        const double query_pt[3] = {fe.position().x() - idata_->CsvForce.shift_pos[0], fe.position().y() - idata_->CsvForce.shift_pos[1], fe.position().z() - idata_->CsvForce.shift_pos[2]};
#endif

        num_results = idata_->kdTree_->knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        auto closest_point_1 = idata_->traction_position_.pts[ret_index[0]];
        auto closest_point_2 = idata_->traction_position_.pts[ret_index[1]];

        ZEROPTV traction1 = idata_->traction_vector_[ret_index[0]];
        ZEROPTV traction2 = idata_->traction_vector_[ret_index[1]];

        ZEROPTV traction;

#if (DIM == 2)
        traction = util_funcs::linear_interpolation_extrapolation(ZEROPTV{closest_point_1.x, closest_point_1.y, 0.0}, ZEROPTV{closest_point_2.x, closest_point_2.y, 0.0}, traction1, traction2, ZEROPTV{query_pt[0], query_pt[1], 0.0});
#endif

#if (DIM == 3)
        traction = util_funcs::linear_interpolation_extrapolation(ZEROPTV{closest_point_1.x, closest_point_1.y, closest_point_1.z}, ZEROPTV{closest_point_2.x, closest_point_2.y, closest_point_2.z}, traction1, traction2, ZEROPTV{query_pt[0], query_pt[1], query_pt[2]});
#endif

#ifndef NDEBUG

#endif
              std::cout << "Query Point:" << std::endl;
              std::cout << "X: " << query_pt[0] << std::endl;
              std::cout << "Y: " << query_pt[1] << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << query_pt[2] << std::endl;
#endif
              std::cout << "Traction for Query Point:" << std::endl;
              std::cout << traction << std::endl;


              std::cout << "1st Closest Point:" << std::endl;
              std::cout << "X: " << closest_point_1.x << std::endl;
              std::cout << "Y: " << closest_point_1.y << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << closest_point_1.z << std::endl;
#endif
              std::cout << "Traction1 for 1st Closest Point:" << std::endl;
              std::cout << traction1 << std::endl;


              std::cout << "2nd Closest Point:" << std::endl;
              std::cout << "X: " << closest_point_2.x << std::endl;
              std::cout << "Y: " << closest_point_2.y << std::endl;
#if (DIM == 3)
              std::cout << "Z: " << closest_point_2.z << std::endl;
#endif
              std::cout << "Traction2 for 2nd Closest Point:" << std::endl;
              std::cout << traction2 << std::endl;

        return traction;
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
