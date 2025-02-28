//
// Created by maksbh on 5/21/20.
//

#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSPNPNodeData.hpp"
#include "NSPNPInputData.h"
#include "NSUtils.h"
#include "NSParams.h"
#include "NSEquation.h"
#include <VMSparams.h>
class MometumEquation: public TALYFEMLIB::CEquation<NSPNPNodeData> {
	TALYFEMLIB::ZEROPTV forcing_;
	TALYFEMLIB::ZEROPTV forcingPre_;

 private:
	const NSPNPInputData *idata_;
	VMSParams vmsParms;
	const NSParams *nsparams_;
	DENDRITE_REAL tauM_;
	NSIntegrandsGenForm* integrandsGenForm_;

 public:
	explicit MometumEquation(const NSPNPInputData *idata, const NSParams *nsparams, NSIntegrandsGenForm
	*integrandsGenForm)
			: TALYFEMLIB::CEquation<NSPNPNodeData> (false, TALYFEMLIB::kAssembleGaussPoints) {
		idata_ = idata;
		nsparams_ = nsparams;
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
		/// Call the appropriate integrands from CHNSIntegrandsGenForm class;
		integrandsGenForm_->getIntegrandsNSAe(fe, Ae, dt_, t_);
	}

	void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
		/// Call the appropriate integrands from CHNSIntegrandsGenForm class;
		integrandsGenForm_->getIntegrandsNSbe(fe, be, dt_, t_);
	}

};

#endif //DENDRITEKT_NSEQUATION_H
