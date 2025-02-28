//
// Created by makrand on 7/16/20.
//

#pragma once

#include <talyfem/talyfem.h>
#include <talyfem/stabilizer/tezduyar_upwind.h>
#include "NSPNPExpr.hpp"
#include "NSPNPNodeData.hpp"
#include "NSPNPInputData.h"
#include "NSEquation.h"


class PNPEquation1 : public NonlinearEquation<NSPNPNodeData> {
private:

    /// For generic integrands support
    NSIntegrandsGenForm *integrandsGenForm_;
    NSPNPInputData *input_data_;
public:
    PNPEquation1(NSPNPInputData *idata, NSIntegrandsGenForm *integrandsGenForm)
            : NonlinearEquation<NSPNPNodeData>(ENABLE_SURFACE_INTEGRATION, false, TALYFEMLIB::kAssembleGaussPoints),
              integrandsGenForm_(integrandsGenForm),
              input_data_(idata) {
    }


    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }

    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {
        /// Call the appropriate integrands from CHNSIntegrandsGenForm class;
        integrandsGenForm_->getIntegrandsPNPAe(fe, Ae, dt_, t_);
    }

    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
        /// Call the appropriate integrands from CHNSIntegrandsGenForm class;
        integrandsGenForm_->getIntegrandsPNPbe(fe, be, dt_, t_);
    }


//    virtual void
//    Integrands4side(const FEMElm &fe, int sideInd, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be, double *h) {
//        assert(false);
//    }



    virtual void
    Integrands4side(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT &boundaryType, const DENDRITE_UINT &
    boundaryID, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h) {
        assert(false);
    }


    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT &boundaryType, const DENDRITE_UINT &
    boundaryID, TALYFEMLIB::ZeroMatrix<double> &Ae, const double *h) {
        /// Call the appropriate integrands from CHNSIntegrandsGenForm class;
        //todo
        integrandsGenForm_->getIntegrands4sidePNPAe(fe, boundaryType, Ae, dt_, t_);

    }



//    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int sideInd, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {
//        /// Call the appropriate integrands from CHNSIntegrandsGenForm class;
//        //todo
//        integrandsGenForm_->getIntegrands4sidePNPAe(fe, sideInd, Ae, dt_, t_);
//    }





    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, const DENDRITE_UINT &boundaryType, const DENDRITE_UINT &
    boundaryID, TALYFEMLIB::ZEROARRAY<double> &be, const double *h)
    {
        /// Call the appropriate integrands from CHNSIntegrandsGenForm class;
        //todo
        integrandsGenForm_->getIntegrands4sidePNPbe(fe, boundaryType, be, dt_, t_);
    }



//	void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int sideInd, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
//		/// Call the appropriate integrands from CHNSIntegrandsGenForm class;
//		//todo
//		integrandsGenForm_->getIntegrands4sidePNPbe(fe, sideInd, be, dt_, t_);
//	}

};
