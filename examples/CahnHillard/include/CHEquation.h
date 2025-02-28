//
// Created by maksbh on 1/21/19.
//

#pragma once
#include "CHNodeData.h"
#include <talyfem/fem/cequation.h>
#include "CHInputData.h"
class CHEquation : public TALYFEMLIB::CEquation<CHNodeData> {
    const CHInputData * inputData_;
public:
    void Solve(double delta_t, double current_time) override {
        assert(false);
    }
    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }
    CHEquation(const CHInputData * chInputData)
    :inputData_(chInputData){

    }
    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h);
    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const double *h);
    double calcd2f(double phi);
    double calcd1f(double phi);

};

void CHEquation::Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h) {
    using namespace TALYFEMLIB;

    const int spatial_dims = fe.nsd();
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();
    const double phi = this->p_data_->valueFEM(fe, 0);
    const double mobility = 1;
    const double epsilon = inputData_->epsilon;
    const int n_dof = 2;

    double N;
    for (int a = 0; a < nbf; a++) {
        for (int b = 0; b < nbf; b++) {
            double M = fe.N(a) * fe.N(b) * detJxW;
            N = 0;
            for (int k = 0; k < spatial_dims; k++) {
                N += mobility * fe.dN(a, k) * fe.dN(b, k) * detJxW;
            }

            Ae(a * n_dof + 0, b * n_dof + 0) += M / dt_;
            Ae(a * n_dof + 1, b * n_dof + 1) += M;
            Ae(a * n_dof + 0, b * n_dof + 1) += mobility * N;
            Ae(a * n_dof + 1, b * n_dof + 0) += -epsilon * N - M * calcd2f(phi);

        }
    }
}

void CHEquation::Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const double *h) {
    using namespace TALYFEMLIB;

    const int spatial_dims = fe.nsd();
    const int nbf = fe.nbf();
    const double detJxW = fe.detJxW();

    const double mobility = 1;
    const double epsilon = inputData_->epsilon;
    const double phi = this->p_data_->valueFEM(fe, 0);
    const double mu = this->p_data_->valueFEM(fe, 1);
    const double phi_prev = this->p_data_->valueFEM(fe, 2);
    const int n_dof = 2;

    for (int a = 0; a < nbf; a++) {
        be(a * n_dof + 0) += fe.N(a) * (phi - phi_prev) * detJxW / dt_;
        be(a * n_dof + 1) += fe.N(a) * mu * detJxW - fe.N(a) * calcd1f(phi) * detJxW;
        for (int k = 0; k < spatial_dims; k++) {
            be(a * n_dof + 0) += mobility * fe.dN(a, k) * this->p_data_->valueDerivativeFEM(fe, 1, k) * detJxW;
            be(a * n_dof + 1) += -epsilon * fe.dN(a, k) * this->p_data_->valueDerivativeFEM(fe, 0, k) * detJxW;
        }
    }

}

double CHEquation::calcd2f(double phi) {
    return ((3 * phi * phi - 1));
}

double CHEquation::calcd1f(double phi) {
    return ((phi * phi * phi - phi));
}
