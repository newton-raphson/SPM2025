//
// Created by maksbh on 2/1/20.
//

#ifndef DENDRITEKT_MAT_H
#define DENDRITEKT_MAT_H

#include "Basis.h"
#include <refel.h>
#include <DataTypes.h>
namespace TENSOROP {
  template<int eleOrder>
  class TensorMat {
    static constexpr int npe_ = intPow(eleOrder+1,DIM);
    const double *Q1d;
    const double *QT1d;
    const double *Dg;
    const double *DgT;
    const double *W1d;
    static constexpr int eleOrder_ = eleOrder;
    DENDRITE_REAL scale_[DIM];
    DENDRITE_REAL oneByScale;
    DENDRITE_REAL imMat1[npe_*npe_];
    DENDRITE_REAL imMat2[npe_*npe_];
    DENDRITE_REAL imMat3[npe_*npe_];
    DENDRITE_REAL Kx[npe_*npe_];
    DENDRITE_REAL Ky[npe_*npe_];
    DENDRITE_REAL Kz[npe_*npe_];
    DENDRITE_REAL UgradW[npe_*npe_];

    DENDRITE_REAL Nx[npe_*npe_];
    const int numEnteriesPerDof_;


    // Linear Kernels
    DENDRITE_REAL  kernel_W_IP_V[npe_*npe_];          /// (w,v)
    DENDRITE_REAL  kernel_gradW_IP_gradV[npe_*npe_];  /// (grad w, grad v)
    DENDRITE_REAL  kernel_W_IP_gradV[npe_*npe_];      /// (w,grad v)
    DENDRITE_REAL  kernel_gradW_IP_V[npe_*npe_];      /// (grad w,v)
    DENDRITE_REAL  kernel_gradWx_IP_V[npe_*npe_];     /// (wx,v)
    DENDRITE_REAL  kernel_gradWy_IP_V[npe_*npe_];     /// (wy,v)
    DENDRITE_REAL  kernel_gradWz_IP_V[npe_*npe_];     /// (wz,v)

    DENDRITE_REAL  kernel_W_IP_gradVx[npe_*npe_];      /// (w,vx)
    DENDRITE_REAL  kernel_W_IP_gradVy[npe_*npe_];      /// (w,vy)
    DENDRITE_REAL  kernel_W_IP_gradVz[npe_*npe_];     /// (w,vz)
    DENDRITE_REAL  weightTensorProduct[npe_];

    void computeWeightGP(){
      int counter = 0;
      for (int i = 0; i < eleOrder_ + 1; i++) {
        for (int j = 0; j < eleOrder_ + 1; j++) {
          for (int k = 0; k < eleOrder_ + 1; k++) {
            weightTensorProduct[counter] = W1d[i] * W1d[j] * W1d[k];
            counter++;
          }
        }
      }
    }


    /**
     * matVec = \alpha * mat*vec + \beta * matVec
     * @param mat
     * @param vec
     * @param matVec
     * @param alpha
     */
    void
    DGEMV(DENDRITE_REAL *mat, DENDRITE_REAL *vec, DENDRITE_REAL *matVec, const DENDRITE_UINT dof, DENDRITE_REAL alpha,
          DENDRITE_REAL beta = 0) {

      char TRANSA = 'N';
      int M = npe_;
      int N = npe_;
      int LDA = npe_;
      int LDB = npe_;
      int LDC = npe_;
      int INCX = 1;
      double BETA = 0;
      dgemv_(&TRANSA, &M, &N, &alpha, mat, &LDA, &vec[dof * npe_], &INCX, &beta, matVec, &INCX);
    }


    /**
     * mat = \alpha * mat1'mat2 + \beta * mat
     * @param mat1
     * @param mat2
     * @param mat
     * @param id_a
     * @param id_b
     * @param alpha
     * @param beta
     */
    void DGEMM(DENDRITE_REAL *mat1, DENDRITE_REAL *mat2, DENDRITE_REAL *mat, const DENDRITE_UINT id_a,
               const DENDRITE_UINT id_b, const int SolverDof, DENDRITE_REAL alpha, DENDRITE_REAL beta) {

      int numMatEntriesPerDof = npe_*npe_;
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numMatEntriesPerDof;
      char TRANSA = 'T';
      char TRANSB = 'N';
      int M = npe_;
      int K = npe_;
      int LDA = npe_;
      int LDB = npe_;
      int LDC = npe_;

      dgemm_(&TRANSA, &TRANSB, &M, &M, &M, &alpha, mat1, &LDA, mat2, &LDA, &beta, &mat[matOffset], &LDC);


    }


    void AXPY(DENDRITE_REAL *outputVec, DENDRITE_REAL *inputVec, DENDRITE_REAL alpha) {
      int INC = 1;
      int N = npe_ * npe_;
      daxpy_(&N, &alpha, inputVec, &INC, outputVec, &INC);
    }

    inline void performQuadratureLoop(DENDRITE_REAL *mat, const DENDRITE_REAL scale) {
      int counter = 0;
      for (int npe = 0; npe < npe_; npe++) {
        int startid = npe_ * npe;
#pragma ivdep
        for (counter = 0; counter < npe_; counter++) {
          mat[startid + counter] *= scale * weightTensorProduct[counter];
        }
        /*   for (int k = 0; k < (eleOrder_ + 1); k++) {
             for (int j = 0; j < (eleOrder_ + 1); j++) {
               for (int i = 0; i < (eleOrder_ + 1); i++) {
                 mat[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
                   scale * (W1d[i] * W1d[j] * W1d[k]);
               }
             }
           }*/
      }
    }

    inline void performQuadratureLoop(DENDRITE_REAL *mat, const DENDRITE_REAL *quadValues, const DENDRITE_REAL scale) {

      for (int npe = 0; npe < npe_; npe++) {
        int startid = npe_ * npe;
#pragma ivdep
        for (int counter = 0; counter < npe_; counter++) {
          mat[startid + counter] *= scale * weightTensorProduct[counter]*quadValues[counter];
        }
//        for (int k = 0; k < (eleOrder_ + 1); k++) {
//          for (int j = 0; j < (eleOrder_ + 1); j++) {
//            for (int i = 0; i < (eleOrder_ + 1); i++) {
//              mat[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
//                scale * (W1d[i] * W1d[j] * W1d[k]) * quadValues[counter++];
//            }
//          }
//        }
      }
    }

    inline void
    performQuadratureLoopValOnly(DENDRITE_REAL *mat, const DENDRITE_REAL *quadValues, const DENDRITE_REAL scale) {

      for (unsigned int npe = 0; npe < npe_; npe++) {
        int startid = npe_ * npe;
#pragma ivdep
        for (int counter = 0; counter < npe_; counter++) {
          mat[startid + counter] *= scale * quadValues[counter];
        }
      }
    }


    inline void preCompute_w_IP_v() {
      std::memset(kernel_W_IP_V,0.0, sizeof(double)*numEnteriesPerDof_);
      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1.0);
      DGEMM(Nx, imMat1, kernel_W_IP_V, 0, 0, 1, 1,0);
    }

    inline void preCompute_gradW_IP_gradV() {

      std::memset(kernel_gradW_IP_gradV,0.0, sizeof(double)*numEnteriesPerDof_);

      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1.0);
      DGEMM(Kx, imMat1, kernel_gradW_IP_gradV, 0, 0, 1,1, 0);

      std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat2, 1);
      DGEMM(Ky, imMat2, kernel_gradW_IP_gradV, 0, 0, 1, 1,1);

      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1);
      DGEMM(Kz, imMat1, kernel_gradW_IP_gradV, 0, 0, 1, 1,1);
    }

    inline void preCompute_gradW_IP_V() {
      std::memset(kernel_gradW_IP_V,0.0, sizeof(double)*numEnteriesPerDof_);

      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1);
      DGEMM(Nx, imMat1, kernel_gradWx_IP_V, 0, 0, 1,1, 0);

      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1);
      DGEMM(Nx, imMat1, kernel_gradWy_IP_V, 0, 0, 1, 1,0);

      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1);
      DGEMM(Nx, imMat1, kernel_gradWz_IP_V, 0, 0, 1, 1,0);

      memset(kernel_gradW_IP_V, 0, sizeof(DENDRITE_REAL) * npe_ * npe_);
      AXPY(kernel_gradW_IP_V, kernel_gradWx_IP_V, 1);
      AXPY(kernel_gradW_IP_V, kernel_gradWy_IP_V, 1);
      AXPY(kernel_gradW_IP_V, kernel_gradWz_IP_V, 1);


      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, 1);
      DGEMM(Kx, imMat1, kernel_W_IP_gradVx, 0, 0, 1,1, 0);
      DGEMM(Ky, imMat1, kernel_W_IP_gradVy, 0, 0, 1,1, 0);
      DGEMM(Kz, imMat1, kernel_W_IP_gradVz, 0, 0, 1,1, 0);


      memset(kernel_W_IP_gradV, 0, sizeof(DENDRITE_REAL) * npe_ * npe_);

      AXPY(kernel_W_IP_gradV, kernel_W_IP_gradVx, 1);
      AXPY(kernel_W_IP_gradV, kernel_W_IP_gradVy, 1);
      AXPY(kernel_W_IP_gradV, kernel_W_IP_gradVz, 1);

    }

    void preComputeLinearKernel() {
      preCompute_w_IP_v();
      preCompute_gradW_IP_gradV();
      preCompute_gradW_IP_V();
    }
  public:
    TensorMat(const RefElement *refElement)
      :numEnteriesPerDof_(npe_*npe_) {
      Q1d = refElement->getQ1d();
      QT1d = refElement->getQT1d();
      Dg = refElement->getDg1d();
      DgT = refElement->getDgT1d();
      W1d = refElement->getWgq();

      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder, Q1d, Q1d, Dg, Kx);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder, Q1d, Dg, Q1d, Ky);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder, Dg, Q1d, Q1d, Kz);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder, Q1d, Q1d, Q1d, Nx);

      computeWeightGP();
      preComputeLinearKernel();
    }

    inline void setScale(const DENDRITE_REAL *h) {
      for(int i = 0; i < DIM; i++){
        scale_[i] = h[i] * 0.5;
      }
      oneByScale = 1./scale_[0];
    }
    ~TensorMat() {
      /*delete[] imMat1;
      delete[] imMat2;
      delete[] Kx;
      delete[] Ky;
      delete[] Kz;
      delete[] Nx;


      delete[] kernel_W_IP_V;
      delete[] kernel_gradW_IP_gradV;
      delete[] kernel_W_IP_gradV;
      delete[] kernel_gradW_IP_V;
      delete[] kernel_gradWx_IP_V;
      delete[] kernel_gradWy_IP_V;
      delete[] kernel_gradWz_IP_V;

      delete[] kernel_W_IP_gradVx;
      delete[] kernel_W_IP_gradVy;
      delete[] kernel_W_IP_gradVz;*/
    }

    inline DENDRITE_UINT getEleOrder() const {
      return eleOrder_;
    }
    inline DENDRITE_REAL getScale() const {
      return scale_[0] / 0.5;
    }

#pragma mark "ValueFEM"
    inline void
    calcValueFEM(DENDRITE_REAL *val, DENDRITE_REAL *quad, const DENDRITE_UINT dof, const DENDRITE_REAL prefactor = 1) {
      DGEMV(Nx, val, quad, dof, prefactor, 0);
    }

    inline void calcValueDerivativeFEMX(DENDRITE_REAL *val, DENDRITE_REAL *quad, const DENDRITE_UINT dof,
                                        const DENDRITE_REAL prefactor = 1) {
      DENDRITE_REAL alpha = prefactor * oneByScale;
      DGEMV(Kx, val, quad, dof, alpha, 0);
    }

    inline void calcValueDerivativeFEMY(DENDRITE_REAL *val, DENDRITE_REAL *quad, const DENDRITE_UINT dof,
                                        const DENDRITE_REAL prefactor = 1) {
      DENDRITE_REAL alpha = prefactor * oneByScale;
      DGEMV(Ky, val, quad, dof, alpha, 0);
    }

    inline void calcValueDerivativeFEMZ(DENDRITE_REAL *val, DENDRITE_REAL *quad, const DENDRITE_UINT dof,
                                        const DENDRITE_REAL prefactor = 1) {
      DENDRITE_REAL alpha = prefactor * oneByScale;
      DGEMV(Kz, val, quad, dof, alpha, 0);
    }

    inline void computeUGradW(const DENDRITE_REAL *velx, const DENDRITE_REAL *vely, const DENDRITE_REAL *velz) {
      std::memcpy(UgradW, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoopValOnly(UgradW, velx, 1);
      std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoopValOnly(imMat2, vely, 1);
      AXPY(UgradW, imMat2, 1);
      std::memcpy(imMat3, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoopValOnly(imMat3, velz, 1);
      AXPY(UgradW, imMat3, 1);
    }
//
//    inline void computeUGradV(const DENDRITE_REAL *velx, const DENDRITE_REAL *vely, const DENDRITE_REAL *velz) {
//      std::memcpy(UgradV, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
//      performQuadratureLoopValOnly(UgradV, velx, 1);
//      std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
//      performQuadratureLoopValOnly(imMat2, vely, 1);
//      AXPY(UgradV, imMat2, 1);
//      std::memcpy(imMat3, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
//      performQuadratureLoopValOnly(imMat3, velz, 1);
//      AXPY(UgradV, imMat3, 1);
//    }


#pragma mark "Linear Kernel"
    inline void
    gradW_IP_V(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof, DENDRITE_REAL alpha = 1,
               DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_gradW_IP_V, scale);
    }

    inline void
    W_IP_gradVx(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof, DENDRITE_REAL alpha = 1,
                DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_W_IP_gradVx, scale);
    }

    inline void
    W_IP_gradVy(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof,
                DENDRITE_REAL alpha = 1,DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_W_IP_gradVy, scale);
    }

    inline void
    W_IP_gradVz(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b,const int SolverDof, DENDRITE_REAL alpha = 1,
                DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_W_IP_gradVz, scale);
    }


    inline void
    gradWx_IP_V(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof, DENDRITE_REAL alpha = 1,
                DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_gradWx_IP_V, scale);
    }

    inline void
    gradWy_IP_V(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof,DENDRITE_REAL alpha = 1,
                DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_gradWy_IP_V, scale);
    }

    inline void
    gradWz_IP_V(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b,const int SolverDof, DENDRITE_REAL alpha = 1,
                DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_gradWz_IP_V, scale);
    }

    /**
     * Linear version
     * @param [Out] mat The output matrix
     * @param alpha The prefactor
     * @param id_a
     * @param id_b
     * @param scale
     */

    inline void W_IP_V(DENDRITE_REAL *mat, const DENDRITE_UINT id_a, const DENDRITE_UINT id_b, const int SolverDof, DENDRITE_REAL alpha = 1,
                       DENDRITE_REAL beta = 1) {
      const DENDRITE_UINT matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      DENDRITE_REAL scale = scale_[0] * scale_[0] * scale_[0] * alpha;
      AXPY(&mat[matOffset], kernel_W_IP_V, scale);
    }

    /**
     * Linear version
     * @param mat
     * @param alpha
     * @param id_a
     * @param id_b
     * @param scale
     */
    inline void
    gradW_IP_gradV(DENDRITE_REAL *mat, const int id_a, const int id_b,  const int SolverDof, DENDRITE_REAL alpha = 1,
                   DENDRITE_REAL beta = 1) {
      DENDRITE_REAL scale = scale_[0] * alpha;
      int matOffset = (id_a * SolverDof + id_b) * numEnteriesPerDof_;
      AXPY(&mat[matOffset], kernel_gradW_IP_gradV, scale);
    }

#pragma mark "Non - Linear Kernel"

    inline void
    gradWx_IP_V(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b,const int SolverDof,
                DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Nx, imMat1, mat, id_a, id_b,SolverDof, alpha, beta);
    }

    inline void
    gradWy_IP_V(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b, const int SolverDof,
                DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Nx, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }


    inline void
    gradWz_IP_V(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b, const int SolverDof,
                DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Nx, imMat1, mat, id_a, id_b, SolverDof, alpha, beta);
    }

    /**
     * Non - linear version
     * @param [Out] mat The output matrix
     * @param alpha The prefactor
     * @param id_a
     * @param id_b
     * @param scale
     */

    inline void
    W_IP_V(DENDRITE_REAL *mat, DENDRITE_REAL *quadValues, const int id_a, const int id_b,const int SolverDof,
           DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {

      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      const DENDRITE_REAL scale = scale_[0] * scale_[0] * scale_[0];
      performQuadratureLoop(imMat1, quadValues, scale);
      DGEMM(Nx, imMat1, mat, id_a, id_b,SolverDof, alpha, beta);
    }




    /**
     * Non - linear version
     * @param mat
     * @param quadVal
     * @param alpha
     * @param id_a
     * @param id_b
     * @param scale
     */

    inline void gradW_IP_gradV(DENDRITE_REAL *mat, const DENDRITE_REAL *quadx, const DENDRITE_REAL *quady,
                               const DENDRITE_REAL *quadz, const int id_a, const int id_b, const int SolverDof,
                               DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];


      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadx, scale);
      DGEMM(Kx, imMat1, mat, id_a, id_b,SolverDof, alpha, beta);

      std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat2, quady, scale);
      DGEMM(Ky, imMat2, mat, id_a, id_b,SolverDof, alpha, beta);

      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadz, scale);
      DGEMM(Kz, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }

   /* inline void uDotGradW_IP_uDotGradV(DENDRITE_REAL *mat, const DENDRITE_REAL *velx, const DENDRITE_REAL *vely,
                                       const DENDRITE_REAL *velz, DENDRITE_REAL *quad, const int id_a,
                                       const int id_b, DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, velx, 1);
      std::memcpy(imMat2, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat2, vely, 1);
      AXPY(imMat1, imMat2, 1);
      std::memcpy(imMat2, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat2, velz, 1);
      AXPY(imMat1, imMat2, 1);

      std::memcpy(imMat3, imMat1, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, 1.0);
      DGEMM(imMat1, imMat3, imMat2, 0, 0, alpha, 0);
      const int matOffset = (id_a * ndof_ + id_b) * numEnteriesPerDof_;
      AXPY(&mat[matOffset], imMat2, scale);


    }*/

   /* inline void
    uDotGradW_IP_uDotGradV(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                           DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, UgradW, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, 1.0);
      DGEMM(imMat1, UgradW, imMat2, 0, 0, alpha, 0);
      const int matOffset = (id_a * ndof_ + id_b) * numEnteriesPerDof_;
      AXPY(&mat[matOffset], imMat2, scale);
    }*/

   /* inline void U1DotGradW_IP_U2DotGradV(DENDRITE_REAL *mat, const DENDRITE_REAL *u2gradW, const DENDRITE_REAL *quad,
                                         const int id_a, const int id_b, DENDRITE_REAL alpha = 1,
                                         DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, UgradW, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);
//    std::memcpy(imMat1,u2gradW, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quad,1.0);
//    DGEMM(imMat1,UgradW,imMat2,0,0,alpha,0);
//    const int matOffset = (id_a*ndof_ + id_b)*numEnteriesPerDof_;
//    AXPY(&mat[matOffset],imMat2,scale);
    }*/

    inline void
    uDotGradW_IP_V(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b, const int SolverDof,
                   DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(imMat1, UgradW, mat, id_a, id_b,SolverDof, alpha, beta);

    }
    inline void
    uDotGradW_IP_V(DENDRITE_REAL *mat, const int id_a, const int id_b, const int SolverDof, DENDRITE_REAL alpha = 1,
                   DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, scale);
      DGEMM(imMat1, UgradW, mat, id_a, id_b,SolverDof, alpha, beta);

    }
/*
    inline void
    uDotGradW_IP_gradWx(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(imMat1, UgradW, mat, id_a, id_b, alpha, beta);
//    const DENDRITE_REAL  scale = scale_[0];
//    std::memcpy(imMat1,UgradW, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quad,1.0);
//    DGEMM(Kx,imMat1,imMat2,0,0,alpha,0);
//    const int matOffset = (id_a*ndof_ + id_b)*numEnteriesPerDof_;
//    AXPY(&mat[matOffset],imMat2,scale);
    }

    inline void
    uDotGradW_IP_gradWy(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(imMat1, UgradW, mat, id_a, id_b, alpha, beta);
//    const DENDRITE_REAL  scale = scale_[0];
//    std::memcpy(imMat1,UgradW, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quad,1.0);
//    DGEMM(Ky,imMat1,imMat2,0,0,alpha,0);
//    const int matOffset = (id_a*ndof_ + id_b)*numEnteriesPerDof_;
//    AXPY(&mat[matOffset],imMat2,scale);
    }

    inline void
    uDotGradW_IP_gradWz(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(imMat1, UgradW, mat, id_a, id_b, alpha, beta);
//    const DENDRITE_REAL  scale = scale_[0];
//    std::memcpy(imMat1,UgradW, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quad,1.0);
//    DGEMM(Kz,imMat1,imMat2,0,0,alpha,0);
//    const int matOffset = (id_a*ndof_ + id_b)*numEnteriesPerDof_;
//    AXPY(&mat[matOffset],imMat2,scale);
    }

    inline void
    gradWx_IP_UDotgradV(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);

    }

    inline void
    gradWy_IP_UDotgradV(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);
    }

    inline void
    gradWz_IP_UDotgradV(DENDRITE_REAL *mat, DENDRITE_REAL *quad, const int id_a, const int id_b,
                        DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);
    }*/


    inline void
    gradWx_IP_gradVx(DENDRITE_REAL *mat, const DENDRITE_REAL *quadx, const int id_a, const int id_b,const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadx, scale);
      DGEMM(Kx, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }

    inline void
    gradWx_IP_gradVy(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b, const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Ky, imMat1, mat, id_a, id_b, SolverDof, alpha, beta);
    }


    inline void
    gradWx_IP_gradVz(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b, const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Kz, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }


    inline void
    gradWy_IP_gradVx(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b,const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(Kx, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }

    inline void
    gradWy_IP_gradVy(DENDRITE_REAL *mat, const DENDRITE_REAL *quady, const int id_a, const int id_b,const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quady, scale);
      DGEMM(Ky, imMat1, mat, id_a, id_b,SolverDof, alpha, beta);
    }

    inline void
    gradWy_IP_gradVz(DENDRITE_REAL *mat, const DENDRITE_REAL *quady, const int id_a, const int id_b,const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Ky, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quady, scale);
      DGEMM(Kz, imMat1, mat, id_a, id_b, SolverDof,alpha, beta);
    }

    inline void
    gradWz_IP_gradVx(DENDRITE_REAL *mat, const DENDRITE_REAL *quadz, const int id_a, const int id_b, const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadz, scale);
      DGEMM(Kx, imMat1, mat, id_a, id_b,SolverDof, alpha, beta);
    }

    inline void
    gradWz_IP_gradVy(DENDRITE_REAL *mat, const DENDRITE_REAL *quadz, const int id_a, const int id_b, const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadz, scale);
      DGEMM(Ky, imMat1, mat, id_a, id_b, SolverDof, alpha, beta);
    }


    inline void
    gradWz_IP_gradVz(DENDRITE_REAL *mat, const DENDRITE_REAL *quadz, const int id_a, const int id_b, const int SolverDof,
                     DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0];
      std::memcpy(imMat1, Kz, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quadz, scale);
      DGEMM(Kz, imMat1, mat, id_a, id_b, SolverDof, alpha, beta);
    }


/*    inline void
    W_IP_gradV(DENDRITE_REAL *mat, const DENDRITE_REAL *quadx, const DENDRITE_REAL *quady, const DENDRITE_REAL *quadz,
               const int id_a, const int id_b, DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];

      std::memcpy(imMat1, UgradV, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, scale);
//    performQuadratureLoop(imMat1,quadx,scale);
      DGEMM(imMat1, Nx, mat, id_a, id_b, alpha, beta);

//    std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quadx,scale);
//    DGEMM(imMat1,Nx,mat,id_a,id_b,alpha,beta);
//
//    std::memcpy(imMat2,Nx, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat2,quady,scale);
//    DGEMM(Ky,imMat2,mat,id_a,id_b,alpha,beta);
//
//    std::memcpy(imMat1,Nx, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quadz,scale);
//    DGEMM(Kz,imMat1,mat,id_a,id_b,alpha,beta);
    }

    inline void
    W_IP_UdotgradV(DENDRITE_REAL *mat, const int id_a, const int id_b, DENDRITE_REAL alpha = 1,
                   DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);
    }


    inline void
    W_IP_UdotgradV(DENDRITE_REAL *mat, const DENDRITE_REAL *quad, const int id_a, const int id_b,
                   DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];
      std::memcpy(imMat1, Nx, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, quad, scale);
      DGEMM(UgradV, imMat1, mat, id_a, id_b, alpha, beta);
    }


    inline void
    gradW_IP_V(DENDRITE_REAL *mat, const DENDRITE_REAL *quadx, const DENDRITE_REAL *quady, const DENDRITE_REAL *quadz,
               const int id_a, const int id_b, DENDRITE_REAL alpha = 1, DENDRITE_REAL beta = 1) {
      const DENDRITE_REAL scale = scale_[0] * scale_[0];

      this->computeUGradW(quadx, quady, quadz);
      std::memcpy(imMat1, UgradW, sizeof(DENDRITE_REAL) * npe_ * npe_);
      performQuadratureLoop(imMat1, scale);
      DGEMM(Nx, imMat1, mat, id_a, id_b, alpha, beta);
//    std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quadx,scale);
//    DGEMM(Nx,imMat1,mat,id_a,id_b,alpha,beta);
//
//    std::memcpy(imMat2,Ky, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat2,quady,scale);
//    DGEMM(Nx,imMat2,mat,id_a,id_b,alpha,beta);
//
//    std::memcpy(imMat1,Kz, sizeof(DENDRITE_REAL)*npe_*npe_);
//    performQuadratureLoop(imMat1,quadz,scale);
//    DGEMM(Nx,imMat1,mat,id_a,id_b,alpha,beta);
    }

    void copyUGradW(DENDRITE_REAL *ugradWVec) {
      std::memcpy(ugradWVec, UgradW, sizeof(DENDRITE_REAL) * npe_ * npe_);
    }

    void setUGradW(const DENDRITE_REAL *ugradWVec) {
      std::memcpy(UgradW, ugradWVec, sizeof(DENDRITE_REAL) * npe_ * npe_);
    }
*/


    void printMatToFile(const DENDRITE_REAL *mat, const int SolverDof) {
      std::ofstream fout("T.dat");
      for (int i = 0; i < npe_; i++) {
        for (int dofi = 0; dofi < SolverDof; dofi++) {
          for (int j = 0; j < npe_; j++) {

            for (int dofj = 0; dofj < SolverDof; dofj++) {
              const int tdof = dofi * SolverDof + dofj;
              if (fabs(mat[(npe_ * npe_ * tdof + i * npe_ + j)]) > 0) {
                fout << std::setprecision(7) << mat[(npe_ * npe_ * tdof + i * npe_ + j)] << " ";
              } else {
                fout << 0 << " ";
              }
            }
          }
          fout << "\n";
        }

      }
      fout.close();
//    std::ofstream fout("T.dat");
//    for(int i = 0; i < npe_*npe_; i++){
//      fout << mat[i] << "\n";
//      fout<< "\n";
//    }

    }

    void printMat(const DENDRITE_REAL *mat, const int SolverDof) {

      for (int i = 0; i < npe_; i++) {
        for (int dofi = 0; dofi < SolverDof; dofi++) {
          for (int j = 0; j < npe_; j++) {

            for (int dofj = 0; dofj < SolverDof; dofj++) {
              const int tdof = dofi * SolverDof + dofj;

              std::cout << std::setprecision(7) << mat[(npe_ * npe_ * tdof + i * npe_ + j)] << " ";
            }
          }
          std::cout << "\n";
//    std::ofstream fout("T.dat");
//    for(int i = 0; i < npe_*npe_; i++){
//      fout << mat[i] << "\n";
//      fout<< "\n";
//    }

        }
      }
    }

  };
}

#endif //DENDRITEKT_MAT_H
