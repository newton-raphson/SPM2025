#ifndef LE_KT_CALCHOMOGENIZEDRESPONSE_H
#define LE_KT_CALCHOMOGENIZEDRESPONSE_H

#include <Traversal/Traversal.h>
#include <LEInputData.h>
#include "SSLEEquation.h"

// Collects domain averages of strain and stress:
//   <eps> = (1/|Ω|) ∫ eps dΩ
//   <sig> = (1/|Ω|) ∫ sig dΩ
//
// Then main can build C_eff by running 3 load cases:
//
//   Case 1: epsbar=[eps0,0,0]  => column1 = <sig>/eps0
//   Case 2: epsbar=[0,eps0,0]  => column2 = <sig>/eps0
//   Case 3: epsbar=[0,0,gam0]  => column3 = <sig>/gam0
//
class CalcHomogenizedResponse : public Traversal
{
public:
  struct AvgVec2D {
    double exx = 0.0, eyy = 0.0, gxy = 0.0; // engineering shear gamma_xy
    double sxx = 0.0, syy = 0.0, txy = 0.0; // tau_xy
    double vol = 0.0; // |Ω|
  };

  CalcHomogenizedResponse(DA *octDA,
                          const std::vector<TREENODE> &treePart,
                          const VecInfo &v,
                          const DomainExtents &domain,
                          const SubDomain *subDomain,
                          LEInputData *idata,
                          const TimeInfo ti,
                          const SSLEEquation *equation = nullptr)
      : Traversal(octDA, treePart, v, domain),
        subdomain_(subDomain),
        idata_(idata),
        ti_(ti),
        equation_(equation)
  {
    // run traversal immediately
    this->traverse();

    // reduce across MPI
    reduceMPI_();
  }

  AvgVec2D getAverages() const { return avg_; }

  // Convenience: just averaged stress vector in Voigt order [sxx, syy, txy]
  void getAvgStress(double out[3]) const {
    out[0] = avg_.sxx; out[1] = avg_.syy; out[2] = avg_.txy;
  }

  // Convenience: just averaged strain vector in Voigt order [exx, eyy, gxy]
  void getAvgStrain(double out[3]) const {
    out[0] = avg_.exx; out[1] = avg_.eyy; out[2] = avg_.gxy;
  }

private:
  const SubDomain *subdomain_ = nullptr;
  LEInputData *idata_ = nullptr;
  const TimeInfo ti_;
  const SSLEEquation *equation_ = nullptr;

  AvgVec2D avg_;      // final (MPI-reduced) averages
  AvgVec2D local_;    // local accumulators before MPI reduction

  void traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) override
  {
#if (DIM != 2)
    static_assert(DIM == 2, "CalcHomogenizedResponse currently implemented for DIM==2");
#endif

    const DENDRITE_UINT ndof = this->getNdof();
    fe.refill(0, 0);

    while (fe.next_itg_pt())
    {
      // Compute displacement gradients du_i/dx_j at this Gauss point
      DENDRITE_REAL duidj[LENodeData::LE_DOF * DIM];
      calcValueDerivativeFEM(fe, ndof, values, duidj);

      // Engineering strain (small strain) at GP:
      // exx = du/dx, eyy = dv/dy, gxy = du/dy + dv/dx
      const double exx = duidj[0 * DIM + 0];
      const double eyy = duidj[1 * DIM + 1];
      const double gxy = duidj[0 * DIM + 1] + duidj[1 * DIM + 0];

      double localC[3 * (DIM - 1)][3 * (DIM - 1)]{};
      const bool useLocalMaterial = (equation_ != nullptr);
      if (useLocalMaterial) {
        equation_->buildLocalCmatrix(fe.position(), localC);
      }
      const auto coeff = [&](int row, int col) -> double {
        return useLocalMaterial ? localC[row][col] : idata_->Cmatrix[row][col];
      };

      double sxx = coeff(0, 0) * exx + coeff(0, 1) * eyy + coeff(0, 2) * gxy;
      double syy = coeff(1, 0) * exx + coeff(1, 1) * eyy + coeff(1, 2) * gxy;
      double txy = coeff(2, 0) * exx + coeff(2, 1) * eyy + coeff(2, 2) * gxy;

      // GP weight (physical): detJ * w_gp
      const double w = fe.detJxW();

      // Accumulate integrals
      local_.exx += exx * w;
      local_.eyy += eyy * w;
      local_.gxy += gxy * w;

      local_.sxx += sxx * w;
      local_.syy += syy * w;
      local_.txy += txy * w;

      local_.vol += w;
    }
  }

  void reduceMPI_()
  {
    // Sum across ranks
    double send[7] = {local_.exx, local_.eyy, local_.gxy,
                      local_.sxx, local_.syy, local_.txy,
                      local_.vol};
    double recv[7] = {0,0,0,0,0,0,0};

    MPI_Allreduce(send, recv, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Normalize to averages
    const double vol = recv[6];
    avg_.vol = vol;

    if (vol > 0.0) {
      avg_.exx = recv[0] / vol;
      avg_.eyy = recv[1] / vol;
      avg_.gxy = recv[2] / vol;

      avg_.sxx = recv[3] / vol;
      avg_.syy = recv[4] / vol;
      avg_.txy = recv[5] / vol;
    }
  }
};

#endif // LE_KT_CALCHOMOGENIZEDRESPONSE_H
