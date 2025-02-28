//
// Created by maksbh on 6/16/20.
//

#ifndef DENDRITEKT_PPEQUATION_H
#define DENDRITEKT_PPEQUATION_H
#include <talyfem/fem/cequation.h>
#include "NSPNPNodeData.hpp"
#include "NSPNPInputData.h"
#include "NSUtils.h"
#include <VMSparams.h>
#include <NSParams.h>
class PressurePoissonEquation: public TALYFEMLIB::CEquation<NSPNPNodeData> {
	TALYFEMLIB::ZEROPTV forcing_;
	TALYFEMLIB::ZEROPTV forcingPre_;

 private:
	const NSPNPInputData *idata_;
	VMSParams vmsParams;
	DENDRITE_REAL tauM_;
	const NSParams *nsParams_;
	NSIntegrandsGenForm *integrandsGenForm_;

 public:
	explicit PressurePoissonEquation(const NSPNPInputData *idata, const NSParams *nsParams, NSIntegrandsGenForm *integrandsGenForm)
			: TALYFEMLIB::CEquation<NSPNPNodeData> (false, TALYFEMLIB::kAssembleGaussPoints) {
		idata_ = idata;
		nsParams_ = nsParams;
		integrandsGenForm_ = integrandsGenForm;
	}


	void Solve(double dt, double t) override {
		assert(false);
	}

	void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
	                TALYFEMLIB::ZEROARRAY<double> &be) override {
		assert(false);
	}

	void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {
		integrandsGenForm_->getIntegrandsPPAe(fe, Ae, dt_, t_);
	}

	void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
		integrandsGenForm_->getIntegrandsPPbe(fe, be, dt_, t_);
	}
};
#endif //DENDRITEKT_PPEQUATION_H
