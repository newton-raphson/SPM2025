//
// Created by mehdi on 7/31/23.
//

#ifndef DENDRITEKT_NSPOST_H
#define DENDRITEKT_NSPOST_H

#ifndef LINNS_FORCECALC_H
#define LINNS_FORCECALC_H
#include <Traversal/SurfaceLoop.h>
#include <DendriteUtils.h>
#include "NSNodeData.h"

class ForceCalc: public SurfaceLoop{
    DENDRITE_REAL Force_[2*DIM]{};
    DENDRITE_REAL projectedArea_[DIM]{};
    const DENDRITE_REAL Re_;

public:
    ForceCalc(DA * octDA,const std::vector<TREENODE> & treePart,const VecInfo & nsSolution, const DomainExtents & domain,
              const DENDRITE_REAL Re, SubDomainBoundary *subDomainBoundary);
    void performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                            const BoundarySurface &boundarySurface, const PetscScalar *values) override;

    void calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe,const PetscScalar * values);

    void getTotalForce(DENDRITE_REAL * globalForce, DENDRITE_REAL * projectedArea);

};

ForceCalc::ForceCalc(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &nsSolution,
                     const DomainExtents &domain, const DENDRITE_REAL Re, SubDomainBoundary *subDomainBoundary)
        :SurfaceLoop(octDA,treePart,nsSolution,domain,subDomainBoundary),Re_(Re){
    std::memset(Force_,0,sizeof(DENDRITE_REAL)*DIM*2);
    std::memset(projectedArea_,0,sizeof(DENDRITE_REAL)*DIM);
    this->traverse();
}
void ForceCalc::performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                                   const BoundarySurface &boundarySurface, const PetscScalar *values) {
    if((boundarySurface.boundaryType == BoundaryTypes::CIRCLE) or (boundarySurface.boundaryType == BoundaryTypes::SPHERE)){
        while (surfaceFe.next_itg_pt()) {

            calcForceAtEachQuadraturePoint2D(surfaceFe,values);
        }
    }
}
void ForceCalc::calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe, const PetscScalar * values){
    using namespace TALYFEMLIB;
    const int nsd = DIM;
    double Coe_diff = 1.0 / Re_;
    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    DENDRITE_REAL vald[NSNodeData::NS_DOF*DIM];
    DENDRITE_REAL valc[NSNodeData::NS_DOF];
    calcValueDerivativeFEM(fe,NSNodeData::NS_DOF,values,vald);
    calcValueFEM(fe,NSNodeData::NS_DOF,values,valc);

    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            du(i, j) = vald[i*DIM+j] + vald[j*DIM + i];
        }
    }
    ZEROPTV diff;
    ZEROPTV surface_normal = fe.surface()->normal();

    for (int i = 0; i < nsd; i++) {
        diff(i) = 0;
        for (int j = 0; j < nsd; j++) {
            //diff(i) += Coe_diff*du(i,
            // j)*fe.normal(j);
            diff(i) += Coe_diff * du(i, j) * surface_normal(j);
        }
    }
    DENDRITE_REAL Force[DIM*2];
    const DENDRITE_REAL & p = valc[NSNodeData::PRESSURE];
//    std::cout << fe.surface()->normal() << " " << p << "\n";
    for (int j = 0; j < nsd; j++) {
        //Force(j) = (-p*fe.normal(j) + diff(j))*detJxW;
        Force[j*2 + 0] = (-p * surface_normal(j)) * fe.detJxW();
        Force[j*2 + 1] = (diff(j)) * fe.detJxW();
        projectedArea_[j]  += fabs(fe.detJxW()*surface_normal(j));

    }
    for (int d = 0; d < DIM; d++){
        Force_[2*d + 0] += Force[2*d+0];
        Force_[2*d + 1] += Force[2*d+1];
    }
}

void ForceCalc::getTotalForce(DENDRITE_REAL * globalForce, DENDRITE_REAL * projectedArea){
    MPI_Reduce(Force_,globalForce,DIM*2,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
    MPI_Reduce(projectedArea_,projectedArea,DIM,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());

}

#endif //LINNS_FORCECALC_H




#endif //DENDRITEKT_NSPOST_H
