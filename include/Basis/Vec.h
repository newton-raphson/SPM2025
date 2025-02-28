//
// Created by maksbh on 1/28/20.
//

#ifndef DENDRITEKT_VEC_H
#define DENDRITEKT_VEC_H
#include <DataTypes.h>
#include <refel.h>
#include "Basis.h"
namespace TENSOROP {
  template<int eleOrder>
  class TensorVec {
    static constexpr int npe_ = intPow(eleOrder+1,DIM);
    const double *Q1d;
    const double *QT1d;
    const double *Dg;
    const double *DgT;
    const double *W1d;
    static constexpr int eleOrder_ = eleOrder;
    DENDRITE_REAL oneByScale;
    DENDRITE_REAL imMat1[npe_*npe_];
    DENDRITE_REAL imMat2[npe_*npe_];
    DENDRITE_REAL imMat3[npe_*npe_];
    DENDRITE_REAL Kx[npe_*npe_];
    DENDRITE_REAL Ky[npe_*npe_];
    DENDRITE_REAL Kz[npe_*npe_];
    DENDRITE_REAL uGradW[npe_*npe_];
    DENDRITE_REAL  weightTensorProduct[npe_];

    DENDRITE_REAL Nx[npe_*npe_];
    static constexpr int numEnteriesPerDof_ = npe_*npe_;
    double scale_[DIM];

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
    void DGEMV(DENDRITE_REAL *mat, DENDRITE_REAL *vec, DENDRITE_REAL *matVec,const DENDRITE_UINT dof, DENDRITE_REAL alpha, DENDRITE_REAL beta = 0 ){

      char TRANSA = 'T';
      int M = npe_;
      int N = npe_;
      int LDA = npe_;
      int LDB = npe_;
      int LDC = npe_;
      int INCX = 1;
      double BETA = 0;
      dgemv_(&TRANSA,&M,&N,&alpha,mat,&LDA,vec,&INCX,&beta,&matVec[dof*npe_],&INCX);
    }

    void DGEMVN(DENDRITE_REAL *mat, DENDRITE_REAL *vec, DENDRITE_REAL *matVec,const DENDRITE_UINT dof, DENDRITE_REAL alpha, DENDRITE_REAL beta = 0 ){

      char TRANSA = 'N';
      int M = npe_;
      int N = npe_;
      int LDA = npe_;
      int LDB = npe_;
      int LDC = npe_;
      int INCX = 1;
      double BETA = 0;
      dgemv_(&TRANSA,&M,&N,&alpha,mat,&LDA,vec,&INCX,&beta,&matVec[dof*npe_],&INCX);
    }
    void AXPY(DENDRITE_REAL * outputVec,DENDRITE_REAL * inputVec, DENDRITE_REAL alpha){
      int INC = 1;
      int N = npe_*npe_;
      daxpy_(&N,&alpha,inputVec,&INC,outputVec,&INC);
    }
    void AXPYVec(DENDRITE_REAL * outputVec,DENDRITE_REAL * inputVec, DENDRITE_REAL alpha){
      int INC = 1;
      int N = npe_;
      daxpy_(&N,&alpha,inputVec,&INC,outputVec,&INC);
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
//    inline void performQuadratureLoopValOnly(DENDRITE_REAL * mat, const DENDRITE_REAL * quadValues, const DENDRITE_REAL scale){
//
//      for (unsigned int npe = 0; npe < npe_; npe++) {
//        DENDRITE_UINT  counter = 0;
//        int startid = npe_ * npe;
//        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
//          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
//            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
//              mat[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
//                scale*quadValues[counter++];
//            }
//          }
//        }
//      }
//    }
//    inline void performQuadratureLoop(DENDRITE_REAL * vec, const DENDRITE_REAL scale){
//      for (unsigned int npe = 0; npe < npe_; npe++) {
//        int startid = npe_ * npe;
//        for (unsigned int k = 0; k < (eleOrder_ + 1); k++) {
//          for (unsigned int j = 0; j < (eleOrder_ + 1); j++) {
//            for (unsigned int i = 0; i < (eleOrder_ + 1); i++) {
//              vec[startid + k * (eleOrder_ + 1) * (eleOrder_ + 1) + j * (eleOrder_ + 1) + i] *=
//                scale*(W1d[i] * W1d[j] * W1d[k]);
//            }
//          }
//        }
//      }
//    }
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
  public:
    TensorVec(const RefElement *refElement){
      Q1d = refElement->getQ1d();
      QT1d = refElement->getQT1d();
      Dg = refElement->getDg1d();
      DgT = refElement->getDgT1d();
      W1d = refElement->getWgq();
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Q1d,Dg,Kx);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Dg,Q1d,Ky);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Dg,Q1d,Q1d,Kz);
      TENSOROP::DENDRO_TENSOR_AAAX_APPLY_ELEM(eleOrder,Q1d,Q1d,Q1d,Nx);
      computeWeightGP();
    }
    inline void setScale(const double * h){
      for(int d = 0; d < DIM; d++){
        scale_[d] = 0.5*h[d];
      }
      oneByScale = 1/scale_[0];
    }


    inline void W_IP_Val(double * vec,  double * Quadvalues, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,Nx, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]*scale_[0]);
      DGEMV(imMat1,Quadvalues,vec,id,alpha,beta);
    }

    inline void gradW_IP_Val(double * vec,  double * QuadvaluesX, double * QuadvaluesY,double * QuadvaluesZ, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]);
      DGEMV(imMat1,QuadvaluesX,vec,id,alpha,beta);

      std::memcpy(imMat2,Ky, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat2,scale_[0]*scale_[0]);
      DGEMV(imMat2,QuadvaluesY,vec,id,alpha,beta);

      std::memcpy(imMat3,Kz, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat3,scale_[0]*scale_[0]);
      DGEMV(imMat3,QuadvaluesZ,vec,id,alpha,beta);

    }

    inline void gradWx_IP_Val(double * vec,  double * QuadvaluesX, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]);
      DGEMV(imMat1,QuadvaluesX,vec,id,alpha,beta);

    }

    inline void gradWy_IP_Val(double * vec,  double * QuadvaluesY, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,Ky, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]);
      DGEMV(imMat1,QuadvaluesY,vec,id,alpha,beta);

    }

    inline void gradWz_IP_Val(double * vec,  double * QuadvaluesZ, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,Kz, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]);
      DGEMV(imMat1,QuadvaluesZ,vec,id,alpha,beta);

    }

    inline void calcValueDerivativeFEMX(DENDRITE_REAL * val,DENDRITE_REAL * quad, const DENDRITE_UINT dof,const DENDRITE_REAL prefactor = 1){
      DENDRITE_REAL alpha = prefactor/scale_[0];
      DGEMVN(Kx, val, quad,dof, alpha,0);
    }

    inline void calcValueDerivativeFEMY(DENDRITE_REAL * val,DENDRITE_REAL * quad, const DENDRITE_UINT dof,const DENDRITE_REAL prefactor = 1){
      DENDRITE_REAL alpha = prefactor/scale_[0];
      DGEMVN(Ky, val, quad,dof,alpha,0);
    }

    inline void calcValueDerivativeFEMZ(DENDRITE_REAL * val,DENDRITE_REAL * quad,const DENDRITE_UINT dof, const DENDRITE_REAL prefactor = 1){
      DENDRITE_REAL alpha = prefactor/scale_[0];
      DGEMVN(Kz, val, quad,dof,alpha,0);
    }

    inline void calcValueFEM(DENDRITE_REAL * val,DENDRITE_REAL * quad, const DENDRITE_UINT dof, const DENDRITE_REAL prefactor = 1){
      DGEMVN(Nx, val, quad,dof, prefactor,0);
    }

    inline void computeUGradW( const double * velx,const double * vely,const double * velz){
      std::memcpy(uGradW,Kx, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoopValOnly(uGradW,velx,1);
      std::memcpy(imMat2,Ky, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoopValOnly(imMat2,vely,1);
      AXPY(uGradW,imMat2,1);
      std::memcpy(imMat3,Kz, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoopValOnly(imMat3,velz,1);
      AXPY(uGradW,imMat3,1);
    }
    inline void UGradW_IP_Val(double * vec,  double * Quadvalues, const DENDRITE_UINT id, DENDRITE_REAL alpha = 1, DENDRITE_REAL  beta = 1){
      std::memcpy(imMat1,uGradW, sizeof(DENDRITE_REAL)*npe_*npe_);
      performQuadratureLoop(imMat1,scale_[0]*scale_[0]);
      DGEMV(imMat1,Quadvalues,vec,id,alpha,beta);
    }

    /**
     * Convert from aaaaabbbb to abababab
     * @param in
     * @param val
     */
    inline void syncToDendro(const double * in, double * out, const int solverDof){
      for(int dof = 0; dof < solverDof; dof++){
        for(int i = 0; i < npe_; i++){
          out[i*solverDof + dof] = in[dof*npe_ + i];
        }
      }
    }

    inline void printVec(const double * out, const int solverDof){
      for(int i = 0; i < solverDof*npe_; i++){
        std::cout << out[i] << "  ";
      }
      std::cout << "\n";
    }

    bool detectNanORInf(const double * out, const int size){
      bool isNan = std::any_of(out, out + size, [](double x){
        return (std::isinf(x) or std::isnan(x));
      });
      return isNan;
    }

  };
}
#endif //DENDRITEKT_VEC_H
