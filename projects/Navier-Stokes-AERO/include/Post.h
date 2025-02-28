//
// Created by mehdi on 10/23/23.
//

#ifndef DENDRITEKT_POST_H
#define DENDRITEKT_POST_H

#include "OctToPhysical.h"
#include <Traversal/Traversal.h>


class GPWriter2D : public Traversal {


    DENDRITE_UINT ndof_;

    std::string fPrefix = "trueGP.vtk";
    FILE *fp = fopen(fPrefix.c_str(), "w");

    DENDRITE_REAL * L2error_ = nullptr;
    DENDRITE_REAL *val_c_ = nullptr;
    AnalyticFunction analFunction_;
    DENDRITE_REAL time_;
    NSInputData

protected:





public:
    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values, const int elemID) override;

    GPWriter2D(DA *octDA, const std::vector<TREENODE> & treePart,const VecInfo & v, const AnalyticFunction &f,  const DomainExtents & domain, const DENDRITE_REAL time)
            :Traversal(octDA,treePart,v,domain),ndof_(v.ndof){
        if (TALYFEMLIB::GetMPIRank() == 0) { // if:master cpu, then:print
            fprintf(fp, "GP_x0,GP_y0\n");
        }

        analFunction_ = f;
        L2error_ = new DENDRITE_REAL[v.ndof];
        memset(L2error_,0, sizeof(DENDRITE_REAL)*v.ndof);
        val_c_ = new DENDRITE_REAL[v.ndof];
        time_ = time;


        this->traverse();
    }
};


void calc_force2D(NSInputData *idata, const int nodeIndicator, double t) {
    using namespace TALYFEMLIB;
    FEMElm fe(p_grid_, BASIS_ALL);
    std::vector<double> force(6, 0.0);
    for (int elmID = 0; elmID < p_grid_->n_elements(); elmID++) {
        if (!p_grid_->IsMyElement(elmID)) {
            continue;
        }
        fe.refill(elmID, 0);
        ELEM::SurfaceList_type::const_iterator it;
        for (it = fe.elem()->surface_indicator_.begin();
             it != fe.elem()->surface_indicator_.end(); it++) {
            fe.refill_surface(elmID, &*it, 0);
            fe.reset_element();
            ZEROPTV position;
            while (fe.next_itg_pt()) {
                if ((it->has_indicator(nodeIndicator))) {
                    calcForceAtEachQuadraturePoint2D(fe, nodeIndicator, force, idata);
                }
            }
        }
    }
    double globalForce[6];
    MPI_Reduce(force.data(), globalForce, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//        MPI_Allreduce(force.data(), globalForce, 6, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
//        PrintStatus("Force = " , 2*globalForce[0] , " " , 2*globalForce[1] );
    if (GetMPIRank() == 0) {
        std::ofstream forcefile("forces_time_series.txt", std::ios::app);
        forcefile << t << " " << 2 * globalForce[0] << " " << 2 * globalForce[1] << " " << 2 * globalForce[2]
                  << " " << 2 * globalForce[3] << " " << 2 * globalForce[4] << " " << 2 * globalForce[5] << "\n";
        forcefile.close();
    }
}

void calcForceAtEachQuadraturePoint2D(const FEMElm &fe, int sideInd, std::vector<double> &sumForce, NSInputData *idata) {
    using namespace TALYFEMLIB;
    const int nsd = idata->nsd;
    double Coe_diff = 1.0 / idata->ns_params.Re;
    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            du(i, j) = valueDerivativeFEM(fe, i, j) + valueDerivativeFEM(fe, j, i);
        }
    }
    double p = valueFEM(fe, NSNodeData::NS_PRESSURE_IDX);
    ZEROPTV diff;
    ZEROPTV surface_normal = fe.surface()->normal();
    double mag = fe.position().x() * fe.position().x() + fe.position().y() * fe.position().y();
    //surface_normal.x() = fe.position().x() / sqrt(mag);
    //surface_normal.y() = fe.position().y() / sqrt(mag);
    for (int i = 0; i < nsd; i++) {
        diff(i) = 0;
        for (int j = 0; j < nsd; j++) {
            //diff(i) += Coe_diff*du(i,j)*fe.normal(j);
            diff(i) += Coe_diff * du(i, j) * surface_normal(j);
        }
    }
    ZEROPTV pt = fe.position();
    ZEROPTV Force;
    ZEROPTV ForcePres;
    ZEROPTV ForceVisc;
    for (int j = 0; j < nsd; j++) {
        //Force(j) = (-p*fe.normal(j) + diff(j))*detJxW;
        Force(j) = (-p * surface_normal(j) + diff(j)) * fe.detJxW();
        ForcePres(j) = (-p * surface_normal(j)) * fe.detJxW();
        ForceVisc(j) = (diff(j)) * fe.detJxW();
    }
    sumForce[0] += (Force(0));
    sumForce[1] += (Force(1));
    sumForce[2] += (ForcePres(0));
    sumForce[3] += (ForcePres(1));
    sumForce[4] += (ForceVisc(0));
    sumForce[5] += (ForceVisc(1));
}











void AnalyticNew::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values, const int elemID) {
    fe.refill(0,0);
    const int ndof= 1;
    while (fe.next_itg_pt()) {
        calcValueFEM(fe, ndof, values, val_c_);
        for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
            DENDRITE_REAL val_a = analFunction_(fe.position(), dof , time_);
            L2error_[dof] += (val_c_[dof] - val_a) * (val_c_[dof] - val_a) * fe.detJxW();
        }
    }
}






void GPWriter2D::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values, const int elemID) {

    fe.refill(0,0);
    const int ndof= 1;
//    std::cout << "inside traverseOperation" << std::endl;
    while (fe.next_itg_pt()) {
        fprintf(fp, "%.10e,%.10e\n",
                fe.position().x(), fe.position().y());
    }
}







#endif //DENDRITEKT_POST_H
