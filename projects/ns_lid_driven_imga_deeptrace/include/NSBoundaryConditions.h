//
// Created by maksbh on 6/15/20.
//

#ifndef DENDRITEKT_NSBCSETUP_H
#define DENDRITEKT_NSBCSETUP_H

#include <TimeInfo.h>
#include <PETSc/BoundaryConditions.h>
#include <DataTypes.h>
#include "NSInputData.h"
#include "NSNodeData.h"

class NSBoundaryConditions {
private:
    const NSInputData *inputData_;
    const TimeInfo *ti_;
    SubDomainBoundary *boundaries_;
    AnalyticFunction analyticFunction_ = nullptr;

    void returnmomentumMMSBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    void returnLDCBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

public:
    NSBoundaryConditions(const NSInputData *idata, const TimeInfo *ti, SubDomainBoundary *boundary);

    void setAnalyticalFunction(const AnalyticFunction &f);

    void getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
};

NSBoundaryConditions::NSBoundaryConditions(const NSInputData *idata, const TimeInfo *ti, SubDomainBoundary *boundary) {
    inputData_ = idata;
    ti_ = ti;
    boundaries_ = boundary;
}

void NSBoundaryConditions::setAnalyticalFunction(const AnalyticFunction &f) {

    analyticFunction_ = f;

}

void NSBoundaryConditions::getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    if (inputData_->ifMMS) {

        if ((analyticFunction_ == nullptr)) {
            throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
        }
        returnmomentumMMSBoundary(b, pos);
    } else {
        returnLDCBoundary(b, pos);
    }
}

void NSBoundaryConditions::returnLDCBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();
    Point<DIM> domainMin(inputData_->meshDef.min);
    Point<DIM> domainMax(inputData_->meshDef.max);

    DENDRITE_UINT objectID = -1;
    boundaries_->generateBoundaryFlags(pos, objectID);


        if (inputData_->InletBCType == NSInputData::LDC){
    #if(DIM == 3)
        double z = pos.z();

    //    if (FEQUALS(x, 0.0) and (FEQUALS(y, 0.0)) and (FEQUALS(z, 0.0))) {
    //        b.addDirichlet(NSNodeData::PRESSURE, 0.0);
    //    }

      bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
                     (fabs(x - domainMax.x(0)) < eps) ||
              (fabs(y - domainMin.x(1)) < eps) ||
              (fabs(z - domainMin.x(2)) < eps) ||
          (fabs(z - domainMax.x(2)) < eps) ||
          boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY);

        bool corner = (fabs(x - domainMin.x(0)) < eps) and
                       (fabs(y - domainMin.x(1)) < eps) and
                       (fabs(z - domainMin.x(2)) < eps);


        if (fabs(y - domainMax.x(1)) < eps) {
            b.addDirichlet(NSNodeData::VEL_X, 1);
            b.addDirichlet(NSNodeData::VEL_Y, 0);
            b.addDirichlet(NSNodeData::VEL_Z, 0);
        }
        else if (no_slip) {
            b.addDirichlet(NSNodeData::VEL_X, 0.0);
            b.addDirichlet(NSNodeData::VEL_Y, 0.0);
            b.addDirichlet(NSNodeData::VEL_Z, 0.0);
        }

        if (corner) {
            b.addDirichlet(NSNodeData::PRESSURE, 0.0);
        }

    //


    #endif
    #if(DIM == 2)
            if (FEQUALS(x, 0.0) and (FEQUALS(y, 0.0))) {
              b.addDirichlet(NSNodeData::PRESSURE, 0.0);
            }
            bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
                (fabs(y - domainMin.x(1)) < eps) ||
                (fabs(x - domainMax.x(0)) < eps);

            if (fabs(y - domainMax.x(1)) < eps) {
              b.addDirichlet(NSNodeData::VEL_X, 1);
              b.addDirichlet(NSNodeData::VEL_Y, 0);
            } else if (no_slip) {
              b.addDirichlet(NSNodeData::VEL_X, 0.0);
              b.addDirichlet(NSNodeData::VEL_Y, 0.0);
            }

    //        bool no_slip = boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY);

    #endif
        } else if (inputData_->InletBCType == NSInputData::FlowPassingObject) {

#if(DIM==3)
            throw TALYFEMLIB::TALYException() << "we only implemented for 2D, please fix it if you want to run 3D";
#endif
#if(DIM == 2)

//            if (fabs(x- inputData_->physDomain.min[0]) < eps and fabs(y - inputData_->physDomain.min[1]) < eps ) {
//
////                std::cout << "pos = " << x << ", " << y << "\n";
//                b.addDirichlet(NSNodeData::PRESSURE, 0.0);
//            }

//            std::cout << "x = " << domainMin.x(0) << ", " << domainMax.x(0) << "\n";
//            std::cout << "y = " << domainMin.x(1) << ", " << domainMax.x(1) << "\n";

            /// (Cheng-Hau) (2024/5/10) we should not use domainMin and domainMax here
            bool free_slip = (fabs(x - inputData_->physDomain.min[0]) < eps);

            if (free_slip) {
                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }

            bool no_penetraction = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                           (fabs(y - inputData_->physDomain.max[1]) < eps);

            if (no_penetraction){
//                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }

            if ((fabs(x - inputData_->physDomain.max[0]) < eps)){
                b.addDirichlet(NSNodeData::PRESSURE, 0);
            }

            // (TODO) haven't implemented the BACKFLOW STAB yet
#endif
        } else if (inputData_->InletBCType == NSInputData::PIPE_FLOW){
            bool free_slip = (fabs(x - inputData_->physDomain.min[0]) < eps);

            bool no_penetraction = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                                   (fabs(y - inputData_->physDomain.max[1]) < eps);

            if (free_slip) {
                b.addDirichlet(NSNodeData::VEL_X, (y - inputData_->physDomain.min[1])* (inputData_->physDomain.max[1] - y)
                /(0.25 * pow(inputData_->physDomain.max[1] - inputData_->physDomain.min[1],2)));
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }
            else if (no_penetraction){
                b.addDirichlet(NSNodeData::VEL_X, 0);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }


        }
    }
void NSBoundaryConditions::returnmomentumMMSBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
        double x = pos.x();
        double y = pos.y();
        if ((FEQUALS(x, 0.0)) or (FEQUALS(y, 0.0)) or (FEQUALS(x, 1.0)) or (FEQUALS(y, 1.0))) {
            b.addDirichlet(NSNodeData::VEL_X, analyticFunction_(pos, NSNodeData::VEL_X, ti_->getCurrentTime()));
            b.addDirichlet(NSNodeData::VEL_Y, analyticFunction_(pos, NSNodeData::VEL_Y, ti_->getCurrentTime()));
        }
        if ((FEQUALS(x, 0.0)) and (FEQUALS(y, 0.0))) {
            b.addDirichlet(NSNodeData::PRESSURE, analyticFunction_(pos, NSNodeData::PRESSURE, ti_->getCurrentTime()));
        }
    }

#endif //DENDRITEKT_NSBCSETUP_H
