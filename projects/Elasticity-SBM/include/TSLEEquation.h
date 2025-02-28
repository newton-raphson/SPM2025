#pragma once

#include <talyfem/fem/cequation.h>
#include <talyfem/stabilizer/tezduyar_upwind.h>
#include "LENodeData.h"
#include "LEInputData.h"
#include <talyfem/talyfem.h>
#include <link.h>
#include <cmath>

class TSLEEquation : public TALYFEMLIB::CEquation<LENodeData>
{

  double timespan = 0;

public:
  explicit TSLEEquation(LEInputData *idata /*, IBMSolver<NSHTNodeData> *ibmSolver*/)
      : TALYFEMLIB::CEquation<LENodeData>(false, TALYFEMLIB::kAssembleGaussPoints)
  {
    //    ibmSolver_ = ibmSolver;
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
    using namespace TALYFEMLIB;

    LEInputData inputData;

    double thickness = 1;

    using namespace TALYFEMLIB;
    // # of dimensions: 1D, 2D, or 3D
    const int n_dimensions = fe.nsd();

    // # of basis functions
    const int n_basis_functions = fe.nbf();

    double Cmatrix[3][3] = {{}};

    /// dynamic
    double dtime = idata_->dt[0];

    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      //C for plane stress
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
    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      //C for plane strain
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      //C for lame parameters
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
    }

    // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::vector<std::vector<double>> Be(2 * n_basis_functions);
    for (int j = 0; j < n_basis_functions; j++)
    {
      Be[2 * j].resize(3);
      Be[2 * j + 1].resize(3);

      Be[2 * j][0] = fe.dN(j, 0);
      Be[2 * j + 1][0] = 0;
      Be[2 * j][1] = 0;
      Be[2 * j + 1][1] = fe.dN(j, 1);
      Be[2 * j][2] = fe.dN(j, 1);
      Be[2 * j + 1][2] = fe.dN(j, 0);
    }

    // for middle term
    std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
    for (int a = 0; a < 2 * n_basis_functions; a++)
    {
      Matmid[a].resize(3);
      for (int b = 0; b < 3; b++)
      {
        double N = 0;
        for (int k = 0; k < 3; k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        Matmid[a][b] = N;
      }
    }

    // for final term

    double beta2 = 0.5;

    const double detJxW = fe.detJxW();
    for (int a = 0; a < 2 * n_basis_functions; a++)
    {
      for (int b = 0; b < 2 * n_basis_functions; b++)
      {
        double N = 0;
        for (int k = 0; k < 3; k++)
        {
          N += Matmid[a][k] * Be[b][k] * detJxW; // stiffness matrix
        }
        //Ae is symmetric
        Ae(a, b) += N * 0.5 * beta2 * dtime * dtime;
      }
    }

    std::vector<std::vector<double>> Nee(2 * n_basis_functions);

    // <change below to be simple>
    for (int j = 0; j < n_basis_functions; j++)
    {

      Nee[2 * j].resize(2);
      Nee[2 * j + 1].resize(2);

      Nee[2 * j][0] = fe.N(j);
      Nee[2 * j + 1][0] = 0;
      Nee[2 * j][1] = 0;
      Nee[2 * j + 1][1] = fe.N(j);
    }

    // adding "mass matrix"
    //double rho = 1;
    double rho = idata_->rho;
    for (int a = 0; a < 2 * n_basis_functions; a++)
    {
      for (int b = 0; b < 2 * n_basis_functions; b++)
      {
        double N = 0;
        for (int k = 0; k < 2; k++)
        {                                            // sum to achieve matrix multiply
          N += rho * Nee[a][k] * Nee[b][k] * detJxW; // mass matrix
        }
        Ae(a, b) += N;
      }
    }
  }

  void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be)
  {
    using namespace TALYFEMLIB;
    const int elmID = fe.elem()->elm_id();
    const ZEROPTV &p = fe.position();
    const int n_dimensions = fe.nsd();
    // # of basis functions
    const int n_basis_functions = fe.nbf();

    double body_x = idata_->BodyForce(0);
    double body_y = idata_->BodyForce(1);

    double BR_V = idata_->radialbodyforce.br_v;
    int BR_POW = idata_->radialbodyforce.br_pow;

    double x_min = idata_->mesh_def.physDomain.min[0];
    double y_min = idata_->mesh_def.physDomain.min[1];
    double x_max = idata_->mesh_def.physDomain.max[0];
    double y_max = idata_->mesh_def.physDomain.max[1];
    double x_mid = (x_min + x_max) / 2;
    double y_mid = (y_min + y_max) / 2;

    double x = p.x();
    double y = p.y();

    /// dynamic
    double beta2 = 0.5;
    double dtime = idata_->dt[0];

    std::vector<double> acc, vec, disp;
    const int nbf = fe.nbf();
    for (ElemNodeID a = 0; a < nbf; a++)
    {
      const int node_index = fe.elem()->ElemToLocalNodeID(a); // 0~3
      //std::cout<<"node_index:"<<node_index<<"\n";
      acc.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::AX));
      acc.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::AY));
      vec.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::VX));
      vec.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::VY));
      disp.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::UX));
      disp.push_back(this->p_data_->GetNodeData(node_index).value(LENodeData::UY));
    }

    std::vector<double> tot_term;
    for (int a = 0; a < acc.size(); a++)
    {
      tot_term.push_back(disp[a] + dtime * vec[a] + 0.5 * (1 - beta2) * dtime * dtime * acc[a]);
    }

    double thickness = 1;
    double Cmatrix[3][3] = {{}};

    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      //C for plane stress
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
    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      //C for plane strain
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      //C for lame parameters
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
    }

    // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::vector<std::vector<double>> Be(2 * n_basis_functions);
    for (int j = 0; j < n_basis_functions; j++)
    {
      Be[2 * j].resize(3);
      Be[2 * j + 1].resize(3);

      Be[2 * j][0] = fe.dN(j, 0);
      Be[2 * j + 1][0] = 0;
      Be[2 * j][1] = 0;
      Be[2 * j + 1][1] = fe.dN(j, 1);
      Be[2 * j][2] = fe.dN(j, 1);
      Be[2 * j + 1][2] = fe.dN(j, 0);
    }

    // for middle term
    std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
    for (int a = 0; a < 2 * n_basis_functions; a++)
    {
      Matmid[a].resize(3);
      for (int b = 0; b < 3; b++)
      {
        double N = 0;
        for (int k = 0; k < 3; k++)
        { // sum to achieve matrix multiply
          N += Be[a][k] * Cmatrix[k][b];
        }
        Matmid[a][b] = N;
      }
    }

    // for final term
    int matsize = 2 * n_basis_functions;
    const double detJxW = fe.detJxW();
    for (int a = 0; a < 2 * n_basis_functions; a++)
    {
      for (int b = 0; b < 2 * n_basis_functions; b++)
      {
        double N = 0;
        for (int k = 0; k < 3; k++)
        {
          N += Matmid[a][k] * Be[b][k] * detJxW; // stiffness matrix
        }
        be(a) -= N * tot_term[b];
      }
    }

    for (int a = 0; a < fe.nbf(); a++)
    {
      be(2 * a) += fe.N(a) * body_x * fe.detJxW();
      be(2 * a + 1) += fe.N(a) * body_y * fe.detJxW();
    }

    for (int a = 0; a < fe.nbf(); a++)
    {
      double r_value = 0;
      double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
      r_value=BR_V * pow(radius, BR_POW);
      double tan_x = (x - x_mid) / radius;
      double tan_y = (y - y_mid) / radius;
      be(2 * a) += fe.N(a) * r_value * tan_x * fe.detJxW();
      be(2 * a + 1) += fe.N(a) * r_value * tan_y * fe.detJxW();
    }
  }

  void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae)
  {
    const auto &bc = idata_->boundary_def.at(side_idx);
    Integrands4side_Ae_domain(fe, Ae, bc);
  }

  void Integrands4side_Ae_domain(const TALYFEMLIB::FEMElm &fe,
                                 TALYFEMLIB::ZeroMatrix<double> &Ae,
                                 const BoundaryDef &bc)
  {

    const double Cb_e = idata_->Cb_e.value_at(t_);
    const int nsd = DIM;
    const int n_basis_functions = fe.nbf();
    const double detSideJxW = fe.detJxW();

    double Cmatrix[3][3] = {{}};
    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      //C for plane stress
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
    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      //C for plane strain
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      //C for lame parameters
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
    }

    /*
    if (bc.disp_type == BoundaryDef::Disp_Type::WEAK ||
        bc.disp_type == BoundaryDef::Disp_Type::WEAK_DX ||
        bc.disp_type == BoundaryDef::Disp_Type::WEAK_DY)
    {
      TALYFEMLIB::ZeroMatrix<double> ksiX;
      ksiX.redim(nsd, nsd);

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          /// This is again correct, we need the volume jacabion on the surface of the box, not the same as in ibm
          ksiX(i, j) = fe.cof(j, i) / fe.volume_jacc();
        }
      }

      TALYFEMLIB::ZeroMatrix<double> Ge;
      Ge.redim(nsd, nsd);

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          Ge(i, j) = 0.0;
          for (int k = 0; k < nsd; k++)
            Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }

      double n_Gn = 0.0;

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          n_Gn += fe.surface()->normal()(j) * Ge(j, i) * fe.surface()->normal()(i);
        }
      }

      /// wall normal distance, check this to see if ksiX is right or not
      double hb = 2.0 / sqrt(n_Gn);
      double weakBCpenaltyParameter = Cb_e / pow(hb, idata_->orderOfhb);
      if (idata_->weakBCglobal)
      {
        weakBCpenaltyParameter = Cb_e;
      }
      if (weakBCpenaltyParameter < 4)
      {
        //        PrintWarning("WeakBC parameter less than 4");
      }

      // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      std::vector<std::vector<double>> Be(2 * n_basis_functions);
      for (int j = 0; j < n_basis_functions; j++)
      {
        Be[2 * j].resize(3);
        Be[2 * j + 1].resize(3);

        Be[2 * j][0] = fe.dN(j, 0);
        Be[2 * j + 1][0] = 0;
        Be[2 * j][1] = 0;
        Be[2 * j + 1][1] = fe.dN(j, 1);
        Be[2 * j][2] = fe.dN(j, 1);
        Be[2 * j + 1][2] = fe.dN(j, 0);
      }

      // for middle term
      std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
      for (int a = 0; a < 2 * n_basis_functions; a++)
      {
        Matmid[a].resize(3);
        for (int b = 0; b < 3; b++)
        {
          double N = 0;
          for (int k = 0; k < 3; k++)
          { // sum to achieve matrix multiply
            N += Be[a][k] * Cmatrix[k][b];
          }
          Matmid[a][b] = N;
        }
      }

      const ZEROPTV &normal = fe.surface()->normal();
      for (int a = 0; a < fe.nbf(); a++)
      {
        double supg_x = 0.0;
        double supg_y = 0.0;
        supg_x += +Matmid[2 * a][0] * fe.surface()->normal()(0)   //* n_x
                  + Matmid[2 * a][2] * fe.surface()->normal()(1); //* n_y

        supg_y += +Matmid[2 * a + 1][1] * fe.surface()->normal()(1)   //* n_y
                  + Matmid[2 * a + 1][2] * fe.surface()->normal()(0); //* n_x

        for (int b = 0; b < fe.nbf(); b++)
        {
          double supgb_x = 0.0;
          double supgb_y = 0.0;
          supgb_x += +Matmid[2 * b][0] * fe.surface()->normal()(0)   //* n_x
                     + Matmid[2 * b][2] * fe.surface()->normal()(1); //* n_y

          supgb_y += +Matmid[2 * b + 1][1] * fe.surface()->normal()(1)   //* n_y
                     + Matmid[2 * b + 1][2] * fe.surface()->normal()(0); //* n_x

          ////////////// fe.N(b) to wall_value in b //////////////////////
          if (bc.disp_type == BoundaryDef::Disp_Type::WEAK)
          {
            Ae(2 * a, 2 * b) +=
                -fe.N(a) * (supgb_x)*detSideJxW                            //consistency
                - supg_x * fe.N(b) * detSideJxW                            //adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
            Ae(2 * a + 1, 2 * b + 1) +=
                -fe.N(a) * (supgb_y)*detSideJxW                            //consistency
                - supg_y * fe.N(b) * detSideJxW                            //adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
          }
          if (bc.disp_type == BoundaryDef::Disp_Type::WEAK_DX)
          {
            Ae(2 * a, 2 * b) +=
                -fe.N(a) * (supgb_x)*detSideJxW                            //consistency
                - supg_x * fe.N(b) * detSideJxW                            //adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
          }
          if (bc.disp_type == BoundaryDef::Disp_Type::WEAK_DY)
          {
            Ae(2 * a + 1, 2 * b + 1) +=
                -fe.N(a) * (supgb_y)*detSideJxW                            //consistency
                - supg_y * fe.N(b) * detSideJxW                            //adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
          }
        } // b loop
      }   // a loop
    }     // end: weakBC
*/
  } // end: Integrands4side_Ae_domain

  void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be)
  {
    Integrands4side_be_domain(fe, be, side_idx);
  }

  void Integrands4side_be_domain(const TALYFEMLIB::FEMElm &fe,
                                 TALYFEMLIB::ZEROARRAY<double> &be,
                                 int side_idx)
  {
    using namespace TALYFEMLIB;
    const int nsd = DIM;
    const int n_basis_functions = fe.nbf();
    const double detSideJxW = fe.detJxW();
    const double Cb_e = idata_->Cb_e.value_at(t_);

    const auto boundary_def = idata_->boundary_def;

    double Cmatrix[3][3] = {{}};
    if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      //C for plane stress
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
    if (idata_->caseType == CaseType::PLANESTRAIN)
    {
      double young = idata_->planeStrain.young;
      double poisson = idata_->planeStrain.poisson;
      //C for plane strain
      Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
    }
    if (idata_->caseType == CaseType::LAME)
    {
      double lamda = idata_->lame.lamda;
      double mu = idata_->lame.mu;
      //C for lame parameters
      Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
    }

    ///////////////////////////////////neumann////////////////////////////////////
    if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::NEUMANN)
    {
      for (int a = 0; a < fe.nbf(); a++)
      {
        const ZEROPTV &normal = fe.surface()->normal();
        for (int dim = 0; dim < DIM; dim++)
        {
          be(2 * a + dim) += fe.N(a) * boundary_def[side_idx].Traction(dim) * fe.detJxW();
        }
      }
    }

    /*
    ///////////////////////////////////weak////////////////////////////////////

    if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK ||
        boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DX ||
        boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DY)
    {

      double DX_wall;
      double DY_wall;

      if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK)
      {
        DX_wall = boundary_def[side_idx].UX;
        DY_wall = boundary_def[side_idx].UY;
      }
      if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DX)
      {
        DX_wall = boundary_def[side_idx].UX;
      }
      if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DY)
      {
        DY_wall = boundary_def[side_idx].UY;
      }

      ZeroMatrix<double> ksiX;
      ksiX.redim(nsd, nsd);

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          /// This is again correct, we need the volume jacabion on the surface of the box, not the same as in ibm
          ksiX(i, j) = fe.cof(j, i) / fe.volume_jacc();
        }
      }

      ZeroMatrix<double> Ge;
      Ge.redim(nsd, nsd);

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          Ge(i, j) = 0.0;
          for (int k = 0; k < nsd; k++)
            Ge(i, j) += ksiX(k, i) * ksiX(k, j);
        }
      }

      double n_Gn = 0.0;

      for (int i = 0; i < nsd; i++)
      {
        for (int j = 0; j < nsd; j++)
        {
          n_Gn += fe.surface()->normal()(j) * Ge(j, i) * fe.surface()->normal()(i);
        }
      }

      /// wall normal distance, check this to see if ksiX is right or not
      double hb = 2.0 / sqrt(n_Gn);
      double weakBCpenaltyParameter = Cb_e / pow(hb, idata_->orderOfhb);
      if (idata_->weakBCglobal)
      {
        weakBCpenaltyParameter = Cb_e;
      }
      if (weakBCpenaltyParameter < 4)
      {
        //        PrintWarning("WeakBC parameter less than 4");
      }

      //std::cout <<"Penalty:" << weakBCpenaltyParameter << "\n";

      // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      std::vector<std::vector<double>> Be(2 * n_basis_functions);
      for (int j = 0; j < n_basis_functions; j++)
      {
        Be[2 * j].resize(3);
        Be[2 * j + 1].resize(3);

        Be[2 * j][0] = fe.dN(j, 0);
        Be[2 * j + 1][0] = 0;
        Be[2 * j][1] = 0;
        Be[2 * j + 1][1] = fe.dN(j, 1);
        Be[2 * j][2] = fe.dN(j, 1);
        Be[2 * j + 1][2] = fe.dN(j, 0);
      }

      // for middle term
      std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
      for (int a = 0; a < 2 * n_basis_functions; a++)
      {
        Matmid[a].resize(3);
        for (int b = 0; b < 3; b++)
        {
          double N = 0;
          for (int k = 0; k < 3; k++)
          { // sum to achieve matrix multiply
            N += Be[a][k] * Cmatrix[k][b];
          }
          Matmid[a][b] = N;
        }
      }

      for (int a = 0; a < fe.nbf(); a++)
      {
        double supg_x = 0.0;
        double supg_y = 0.0;
        supg_x += +Matmid[2 * a][0] * fe.surface()->normal()(0)   //* n_x
                  + Matmid[2 * a][2] * fe.surface()->normal()(1); //* n_y

        supg_y += +Matmid[2 * a + 1][1] * fe.surface()->normal()(1)   //* n_y
                  + Matmid[2 * a + 1][2] * fe.surface()->normal()(0); //* n_x

        if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK)
        {
          be(2 * a) +=
              -supg_x * DX_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DX_wall * detSideJxW; /// penalty
          be(2 * a + 1) +=
              -supg_y * DY_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DY_wall * detSideJxW; /// penalty
        }
        if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DX)
        {
          be(2 * a) +=
              -supg_x * DX_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DX_wall * detSideJxW; /// penalty
        }
        if (boundary_def[side_idx].disp_type == BoundaryDef::Disp_Type::WEAK_DY)
        {
          be(2 * a + 1) +=
              -supg_y * DY_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DY_wall * detSideJxW; /// penalty
        }
      }
    } //end:weak BC
*/
  } //end:Integrands4side_be_domain

  ///// ==================== ibm start====================================

  void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {

    auto &geo = idata_->ibm_geom_def.at(nodeAndValues.geomID);

    /// If the geometry is 2D plane, we should treat this differently.
    ZEROPTV normal = ZEROPTV{nodeAndValues.normal[0], nodeAndValues.normal[1], 0};
#if (DIM == 3)
    normal(2) = nodeAndValues.normal[2];
#endif

    // ====== parameters setting ============
    // penalty
    const double Cb_e = idata_->Cb_e.value_at(t_);

    //  ====== parameters setting end ============

    const int nbf = fe.nbf();
    const int n_basis_functions = fe.nbf();
    const int nsd = DIM;
    const double detSideJxW = fe.detJxW();

    //////////////////////////////////////weak//////////////////////////////////////////

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
    {

      /// wall normal distance
      //double hb = getNormalDistance(fe,nodeAndValues.location,nodeAndValues.normal,h,1E-2);
      double hb = normalDistance(fe, normal, h);
      double weakBCpenaltyParameter = Cb_e / hb;

      double Cmatrix[3][3] = {{}};
      if (idata_->caseType == CaseType::PLANESTRESS)
      {
        double young = idata_->planeStress.young;
        double poisson = idata_->planeStress.poisson;
        //C for plane stress
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
      if (idata_->caseType == CaseType::PLANESTRAIN)
      {
        double young = idata_->planeStrain.young;
        double poisson = idata_->planeStrain.poisson;
        //C for plane strain
        Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[0][2] = 0;
        Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[1][2] = 0;
        Cmatrix[2][0] = 0;
        Cmatrix[2][1] = 0;
        Cmatrix[2][2] = young / (1 + poisson) / 2;
      }
      if (idata_->caseType == CaseType::LAME)
      {
        double lamda = idata_->lame.lamda;
        double mu = idata_->lame.mu;
        //C for lame parameters
        Cmatrix[0][0] = lamda + 2 * mu;
        Cmatrix[0][1] = lamda;
        Cmatrix[0][2] = 0;
        Cmatrix[1][0] = lamda;
        Cmatrix[1][1] = lamda + 2 * mu;
        Cmatrix[1][2] = 0;
        Cmatrix[2][0] = 0;
        Cmatrix[2][1] = 0;
        Cmatrix[2][2] = mu;
      }

      // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      std::vector<std::vector<double>> Be(2 * n_basis_functions);
      for (int j = 0; j < n_basis_functions; j++)
      {
        Be[2 * j].resize(3);
        Be[2 * j + 1].resize(3);

        Be[2 * j][0] = fe.dN(j, 0);
        Be[2 * j + 1][0] = 0;
        Be[2 * j][1] = 0;
        Be[2 * j + 1][1] = fe.dN(j, 1);
        Be[2 * j][2] = fe.dN(j, 1);
        Be[2 * j + 1][2] = fe.dN(j, 0);
      }

      // for middle term
      std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
      for (int a = 0; a < 2 * n_basis_functions; a++)
      {
        Matmid[a].resize(3);
        for (int b = 0; b < 3; b++)
        {
          double N = 0;
          for (int k = 0; k < 3; k++)
          { // sum to achieve matrix multiply
            N += Be[a][k] * Cmatrix[k][b];
          }
          Matmid[a][b] = N;
        }
      }

      for (int a = 0; a < fe.nbf(); a++)
      {
        double supg_x = 0.0;
        double supg_y = 0.0;
        supg_x += +Matmid[2 * a][0] * nodeAndValues.normal[0]   //* n_x
                  + Matmid[2 * a][2] * nodeAndValues.normal[1]; //* n_y

        supg_y += +Matmid[2 * a + 1][1] * nodeAndValues.normal[1]   //* n_y
                  + Matmid[2 * a + 1][2] * nodeAndValues.normal[0]; //* n_x

        for (int b = 0; b < fe.nbf(); b++)
        {
          double supgb_x = 0.0;
          double supgb_y = 0.0;
          supgb_x += +Matmid[2 * b][0] * nodeAndValues.normal[0]   //* n_x
                     + Matmid[2 * b][2] * nodeAndValues.normal[1]; //* n_y

          supgb_y += +Matmid[2 * b + 1][1] * nodeAndValues.normal[1]   //* n_y
                     + Matmid[2 * b + 1][2] * nodeAndValues.normal[0]; //* n_x

          if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
              geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
          {
            Ae(2 * a, 2 * b) +=
                -fe.N(a) * (supgb_x)*detSideJxW                            ///consistency
                - supg_x * fe.N(b) * detSideJxW                            ///adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
          }
          if (geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
              geo.bc_type_D[1] == IBMGeomDef::BCType::W_RADIAL)
          {

            Ae(2 * a + 1, 2 * b + 1) +=
                -fe.N(a) * (supgb_y)*detSideJxW                            ///consistency
                - supg_y * fe.N(b) * detSideJxW                            ///adjoint
                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW; /// penalty, constant factor vs. calculated
          }

        } // b loop
      }   // a loop
    }     //end:weak BC
  }       // end:ibm_Ae

  void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const NodeAndValues<DENDRITE_REAL> &nodeAndValues,
                              const TALYFEMLIB::ZEROPTV &position,
                              const TALYFEMLIB::ZEROPTV &h)
  {

    auto &geo = idata_->ibm_geom_def.at(nodeAndValues.geomID);

    //nodeAndValues.normal does not have any value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /// If the geometry is 2D plane, we should treat this differently.
    ZEROPTV normal = ZEROPTV{nodeAndValues.normal[0], nodeAndValues.normal[1], 0};
#if (DIM == 3)
    normal(2) = nodeAndValues.normal[2];
#endif

    //x and y => below vector only has one value
    std::vector<double> bc_value_x = geo.getBC_D(0);
    std::vector<double> bc_value_y = geo.getBC_D(1);

    // ====== parameters setting ============
    // penalty
    const double Cb_e = idata_->Cb_e.value_at(t_);

    double x_mid = geo.InitialDisplacement.x();
    double y_mid = geo.InitialDisplacement.y();
    double x = fe.position().x();
    double y = fe.position().y();

    const int nbf = fe.nbf();
    const int n_basis_functions = fe.nbf();
    const int nsd = DIM;
    const double detSideJxW = fe.detJxW();

    ///////////////////////////////////neumann////////////////////////////////////

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::NEUMANN ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::NEUMANN ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::T_RADIAL)
    {

      for (int a = 0; a < fe.nbf(); a++)
      {

        // - => pressure/ + => traction
        // x-dir
        if (geo.bc_type_D[0] == IBMGeomDef::BCType::NEUMANN)
        {
          double TX = bc_value_x[0];
          be(2 * a) += fe.N(a) * TX * nodeAndValues.normal[0] * fe.detJxW();
        }
        // y-dir
        if (geo.bc_type_D[1] == IBMGeomDef::BCType::NEUMANN)
        {
          double TY = bc_value_y[0];
          be(2 * a + 1) += fe.N(a) * TY * nodeAndValues.normal[1] * fe.detJxW();
        }

        //r-dir  => this have problems!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (geo.bc_type_D[0] == IBMGeomDef::BCType::T_RADIAL)
        {
          double TR = bc_value_x[0];
          double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
          be(2 * a) += fe.N(a) * TR * nodeAndValues.normal[0] * fe.detJxW();
          be(2 * a + 1) += fe.N(a) * TR * nodeAndValues.normal[1] * fe.detJxW();
        }
      } // end: a loop
    }   // end: neumann
        //////////////////////////////////////weak//////////////////////////////////////////

    if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET ||
        geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
    {

      //std::cout<< "Dirichlet" <<"\n";

      // Dirichlet
      double DX_wall;
      double DY_wall;

      /// put into "be"
      //x-dir
      if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET)
      {
        DX_wall = bc_value_x[0];
      }
      //y-dir
      if (geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET)
      {
        DY_wall = bc_value_y[0];
      }

      double DR;
      if (geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
      {
        DR = bc_value_x[0];
      }

      //std::cout<< "x:" << x << ",y:" << y <<"\n";
      //std::cout<<"fe_x:" << fe.position().x() <<",fe_y:" << fe.position().y() << "\n";

      /// wall normal distance
      //double hb = getNormalDistance(fe,nodeAndValues.location,nodeAndValues.normal,h,1E-2);
      double hb = normalDistance(fe, normal, h);

      double weakBCpenaltyParameter = Cb_e / hb;

      double Cmatrix[3][3] = {{}};
      if (idata_->caseType == CaseType::PLANESTRESS)
      {
        double young = idata_->planeStress.young;
        double poisson = idata_->planeStress.poisson;
        //C for plane stress
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
      if (idata_->caseType == CaseType::PLANESTRAIN)
      {
        double young = idata_->planeStrain.young;
        double poisson = idata_->planeStrain.poisson;
        //C for plane strain
        Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[0][2] = 0;
        Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
        Cmatrix[1][2] = 0;
        Cmatrix[2][0] = 0;
        Cmatrix[2][1] = 0;
        Cmatrix[2][2] = young / (1 + poisson) / 2;
      }
      if (idata_->caseType == CaseType::LAME)
      {
        double lamda = idata_->lame.lamda;
        double mu = idata_->lame.mu;
        //C for lame parameters
        Cmatrix[0][0] = lamda + 2 * mu;
        Cmatrix[0][1] = lamda;
        Cmatrix[0][2] = 0;
        Cmatrix[1][0] = lamda;
        Cmatrix[1][1] = lamda + 2 * mu;
        Cmatrix[1][2] = 0;
        Cmatrix[2][0] = 0;
        Cmatrix[2][1] = 0;
        Cmatrix[2][2] = mu;
      }

      // "remember" start from zero!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      std::vector<std::vector<double>> Be(2 * n_basis_functions);
      for (int j = 0; j < n_basis_functions; j++)
      {
        Be[2 * j].resize(3);
        Be[2 * j + 1].resize(3);

        Be[2 * j][0] = fe.dN(j, 0);
        Be[2 * j + 1][0] = 0;
        Be[2 * j][1] = 0;
        Be[2 * j + 1][1] = fe.dN(j, 1);
        Be[2 * j][2] = fe.dN(j, 1);
        Be[2 * j + 1][2] = fe.dN(j, 0);
      }

      // for middle term
      std::vector<std::vector<double>> Matmid(2 * n_basis_functions);
      for (int a = 0; a < 2 * n_basis_functions; a++)
      {
        Matmid[a].resize(3);
        for (int b = 0; b < 3; b++)
        {
          double N = 0;
          for (int k = 0; k < 3; k++)
          { // sum to achieve matrix multiply
            N += Be[a][k] * Cmatrix[k][b];
          }
          Matmid[a][b] = N;
        }
      }

      for (int a = 0; a < fe.nbf(); a++)
      {
        double supg_x = 0.0;
        double supg_y = 0.0;
        supg_x += +Matmid[2 * a][0] * nodeAndValues.normal[0]   //* n_x
                  + Matmid[2 * a][2] * nodeAndValues.normal[1]; //* n_y

        supg_y += +Matmid[2 * a + 1][1] * nodeAndValues.normal[1]   //* n_y
                  + Matmid[2 * a + 1][2] * nodeAndValues.normal[0]; //* n_x

        // x-dir
        if (geo.bc_type_D[0] == IBMGeomDef::BCType::DIRICHLET)
        {
          be(2 * a) +=
              -supg_x * DX_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DX_wall * detSideJxW; /// penalty
        }

        // y-dir
        if (geo.bc_type_D[1] == IBMGeomDef::BCType::DIRICHLET)
        {
          be(2 * a + 1) +=
              -supg_y * DY_wall * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DY_wall * detSideJxW; /// penalty
        }

        // r-dir
        if (geo.bc_type_D[0] == IBMGeomDef::BCType::W_RADIAL)
        {
          double radius = sqrt(pow(x - x_mid, 2) + pow(y - y_mid, 2));
          double sin_x = (x - x_mid) / radius;
          double sin_y = (y - y_mid) / radius;

          be(2 * a) +=
              -supg_x * DR * sin_x * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DR * sin_x * detSideJxW; /// penalty
          be(2 * a + 1) +=
              -supg_y * DR * sin_y * detSideJxW                             /// adjoint
              + weakBCpenaltyParameter * fe.N(a) * DR * sin_y * detSideJxW; /// penalty
        }

      } // a loop

    } //end:weak BC
  }   // end:ibm_be
      ///// ==================== ibm end====================================

private:
  LEInputData *idata_;

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

  double calc_body_x(const TALYFEMLIB::ZEROPTV &pt) const
  {
    double r = sqrt(pow(pt.x() - 1, 2) + pow(pt.y() - 1, 2));
    double tan_x = (pt.x() - 1) / r;
    return 1 * tan_x / r / log(2);
  }

  double calc_body_y(const TALYFEMLIB::ZEROPTV &pt) const
  {
    double r = sqrt(pow(pt.x() - 1, 2) + pow(pt.y() - 1, 2));
    double tan_y = (pt.y() - 1) / r;
    return 1 * tan_y / r / log(2);
  }

  double calc_tract_x(const TALYFEMLIB::ZEROPTV &pt) const
  {
    double r = sqrt(pow(pt.x() - 1, 2) + pow(pt.y() - 1, 2));
    double tan_x = (pt.x() - 1) / r;
    return (log(0.25) + 1) * tan_x / 2 / log(2);
  }

  double calc_tract_y(const TALYFEMLIB::ZEROPTV &pt) const
  {
    double r = sqrt(pow(pt.x() - 1, 2) + pow(pt.y() - 1, 2));
    double tan_y = (pt.y() - 1) / r;
    return (log(0.25) + 1) * tan_y / 2 / log(2);
  }
};
