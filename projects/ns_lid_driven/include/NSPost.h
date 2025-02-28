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



enum AnalyticType:short{
    L2ERROR=0,
    EVALFUNC=1
};

class Analytic: Traversal{
    DENDRITE_REAL * L2error_ = nullptr;
    DENDRITE_REAL *val_c_ = nullptr;
    AnalyticFunction analFunction_;
    AnalyticType analyticType;
    DENDRITE_REAL time_;
    DENDRITE_UINT sdof_;

    void calculateL2error(DENDRITE_REAL * globalError);
public:
    /**
     * @brief Constructor
     * @param octDA
     * @param v vector
     * @param f analytic function
     * @param time time
     */
    Analytic(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v,  const AnalyticFunction & f, const DomainExtents & domain, const DENDRITE_REAL time = 0);

    /**
     * @brief Overriding the travesal class operation
     * @param fe
     * @param values
     */
    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    /**
     * @brief Prints the  L2 error
     */
    void getL2error();

    /**
     * @brief returns the L2 error
     * @param error Retruns the L2 error. This must be allocated.
     */
    void getL2error(DENDRITE_REAL * error);

    ~Analytic();
};



Analytic::Analytic(DA *octDA, const std::vector<TREENODE> & treePart,const VecInfo & v, const AnalyticFunction &f,  const DomainExtents & domain, const DENDRITE_REAL time)
        :Traversal(octDA,treePart,v,domain){
    analyticType = AnalyticType::L2ERROR;
    analFunction_ = f;
    L2error_ = new DENDRITE_REAL[v.ndof];
    memset(L2error_,0, sizeof(DENDRITE_REAL)*v.ndof);
    val_c_ = new DENDRITE_REAL[v.ndof];
    time_ = time;
    sdof_ = v.nodeDataIndex;
    this->traverse();

}

void Analytic::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {

    const DENDRITE_UINT ndof = this->getNdof();
    while (fe.next_itg_pt()) {
        calcValueFEM(fe, ndof, values, val_c_);
        for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
            DENDRITE_REAL val_a = analFunction_(fe.position(), dof + sdof_, time_);
            L2error_[dof] += (val_c_[dof] - val_a) * (val_c_[dof] - val_a) * fe.detJxW();
        }
    }
}

void Analytic::getL2error() {
    DENDRITE_REAL * globalL2Eror = new DENDRITE_REAL[this->getNdof()];
    calculateL2error(globalL2Eror);
    if (not(TALYFEMLIB::GetMPIRank())) {
        std::cout << "L2Error : ";
        for (int dof = 0; dof < this->getNdof(); dof++) {
            std::cout << std::scientific << globalL2Eror[dof] << " ";
        }
        std::cout << std::defaultfloat << "\n";
    }

    delete [] globalL2Eror;

}

void Analytic::getL2error(DENDRITE_REAL *error) {
    calculateL2error(error);

}

Analytic::~Analytic() {
    if(analyticType == AnalyticType::L2ERROR){
        delete [] L2error_;
        delete [] val_c_;

    }
}

void Analytic::calculateL2error(DENDRITE_REAL *globalL2Error) {
    MPI_Reduce(L2error_,globalL2Error,this->getNdof(),MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
    if(TALYFEMLIB::GetMPIRank() == 0){
        for(int dof = 0; dof < this->getNdof(); dof++) {
            globalL2Error[dof] = sqrt(globalL2Error[dof]);
        }
    }
}






#endif //LINNS_FORCECALC_H




#endif //DENDRITEKT_NSPOST_H
