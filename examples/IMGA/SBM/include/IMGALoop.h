//
// Created by maksbh on 4/1/21.
//

#ifndef DENDRITEKT_IMGALOOP_H
#define DENDRITEKT_IMGALOOP_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>



class IMGALoop: public IMGATraversal{
  double error = 0;
public:
  IMGALoop(DA * octDA, const IMGA * imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain);

  void imgaTraversalOperation(const TALYFEMLIB::FEMElm & fe,const NodeAndValues<DENDRITE_REAL> & gaussPoint, const TALYFEMLIB::ZEROPTV & h, const PetscScalar * values) override;

  void computeBoundaryError(DENDRITE_REAL * globalError);
};

IMGALoop::IMGALoop(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                   const DomainExtents &domain)
                   :IMGATraversal(octDA,imga,treePart,v,domain){
  this->imgaTraverse();
}

void IMGALoop::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV & h,
                                      const PetscScalar *values) {
  DENDRITE_REAL val_c,val_d[DIM];
  calcValueFEM(fe,1,values,&val_c);
  calcValueDerivativeFEM(fe,1,values,val_d);



  // TODO: Compute the error.
  // gaussPoint.normal = contains the position of point on surface.
  // gaussPoint.position = contains the position of point that we need to interpolate/

  // Use Taylor series: u(surface/unknown) = uo (val_c) + grad(u).d (valderivative)

  // d  can be formulated from the position.
}

void IMGALoop::computeBoundaryError(DENDRITE_REAL * globalError) {
    // TODO: add MPI Reduce for error
    // MPI_Reduce()
}
#endif //DENDRITEKT_IMGALOOP_H
