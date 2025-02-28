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

    double InletParabolicVelocity(double y, double max_y, double max_velocity);
    double InletParabolicVelocity_3D(double y, double max_y, double z, double max_z,double max_velocity);
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
            if ((fabs(x - inputData_->physDomain.min[0]) < eps) and (fabs(y - inputData_->physDomain.min[1]) < eps)) {
              b.addDirichlet(NSNodeData::PRESSURE, 0.0);
            }
//        if (FEQUALS(x, -2.55) and (FEQUALS(y, -2.55))) {
//            b.addDirichlet(NSNodeData::PRESSURE, 0.0);
//        }
            bool no_slip = (fabs(x - inputData_->physDomain.min[0]) < eps) ||
                (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                (fabs(x - inputData_->physDomain.max[0]) < eps);

            if (fabs(y - inputData_->physDomain.max[1]) < eps) {
              b.addDirichlet(NSNodeData::VEL_X, 1);
              b.addDirichlet(NSNodeData::VEL_Y, 0);
            } else if (no_slip) {
              b.addDirichlet(NSNodeData::VEL_X, 0.0);
              b.addDirichlet(NSNodeData::VEL_Y, 0.0);
            }

    //        bool no_slip = boundaries_->checkBoundaryType(BoundaryTypes::VOXEL::GEOMETRY);

    #endif
        }
        else if (inputData_->InletBCType == NSInputData::FlowPassingObject) {

#if(DIM==3)
            double z = pos.z();
            if (inputData_->outletBcType == OutletBCType::PRESSUREOUTLET) {
//                THIS IS SAME AS X+
                if (fabs(x - inputData_->physDomain.max[0]) < eps) {
                    b.addDirichlet(NSNodeData::PRESSURE, 0);
                }
            } else {
                if (fabs(x- inputData_->physDomain.min[0]) < eps and fabs(y - inputData_->physDomain.min[1]) < eps and fabs(z- inputData_->physDomain.min[2])< eps) {

                    //                std::cout << "pos = " << x << ", " << y << "\n";
                    b.addDirichlet(NSNodeData::PRESSURE, 0.0);
                }
            }


            /// (Cheng-Hau) (2024/5/10) we should not use domainMin and domainMax here
            bool free_slip = (fabs(x - inputData_->physDomain.min[0]) < eps);

            if (free_slip) {
                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
                b.addDirichlet(NSNodeData::VEL_Z, 0);
            }

            bool no_penetraction_y = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                                   (fabs(y - inputData_->physDomain.max[1]) < eps);
            bool no_penetraction_z = (fabs(z - inputData_->physDomain.min[2]) < eps) ||
                                     (fabs(z - inputData_->physDomain.max[2]) < eps);


            if (no_penetraction_y){
                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }
            if (no_penetraction_z){
                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Z, 0);
            }
#endif
#if(DIM == 2)

            if (inputData_->outletBcType == OutletBCType::PRESSUREOUTLET) {
                if (fabs(x - inputData_->physDomain.max[0]) < eps) {
                    b.addDirichlet(NSNodeData::PRESSURE, 0);
                }
            } else {
                if (fabs(x- inputData_->physDomain.min[0]) < eps and fabs(y - inputData_->physDomain.min[1]) < eps ) {

    //                std::cout << "pos = " << x << ", " << y << "\n";
                    b.addDirichlet(NSNodeData::PRESSURE, 0.0);
                }
            }


//            if ((fabs(x - inputData_->physDomain.max[0]) < eps)){
//                b.addDirichlet(NSNodeData::PRESSURE, 0);
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
                b.addDirichlet(NSNodeData::VEL_X, 1);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }

            // (TODO) haven't implemented the BACKFLOW STAB yet
#endif
        }
        else if (inputData_->InletBCType == NSInputData::PIPE_FLOW){
#if (DIM==2)
            bool free_slip = (fabs(x - inputData_->physDomain.min[0]) < eps);

            bool no_penetraction = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                                   (fabs(y - inputData_->physDomain.max[1]) < eps);

            if (free_slip) {
                b.addDirichlet(NSNodeData::VEL_X,
                               InletParabolicVelocity(y,inputData_->physDomain.max[1] - inputData_->physDomain.min[1],1));
                b.addDirichlet(NSNodeData::VEL_Y, 0);

            }
            else if (no_penetraction){
                b.addDirichlet(NSNodeData::VEL_X, 0);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
            }

            if (fabs(x - inputData_->physDomain.max[0]) < eps) {
                b.addDirichlet(NSNodeData::PRESSURE, 0);
            }
#endif
#if (DIM==3)
            double z = pos.z();
            bool free_slip = (fabs(z - inputData_->physDomain.max[2]) < eps);

            bool no_penetraction = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
                                   (fabs(y - inputData_->physDomain.max[1]) < eps) ||
                                   (fabs(x - inputData_->physDomain.min[0]) < eps) ||
                                    (fabs(x- inputData_->physDomain.max[0]) < eps);

            if (free_slip) {
                b.addDirichlet(NSNodeData::VEL_X,
                               0);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
                b.addDirichlet(NSNodeData::VEL_Z, InletParabolicVelocity_3D(x,inputData_->physDomain.max[0],y,inputData_->physDomain.max[0],1));

            }
            else if (no_penetraction){
                b.addDirichlet(NSNodeData::VEL_X, 0);
                b.addDirichlet(NSNodeData::VEL_Y, 0);
                b.addDirichlet(NSNodeData::VEL_Z, 0);
            }

            if (fabs(z - inputData_->physDomain.min[2]) < eps) {
                b.addDirichlet(NSNodeData::PRESSURE, 0);
            }
#endif


    }
    else if (inputData_->InletBCType == NSInputData::GYROID_FLOW) {
        
// #if (DIM==2)
//             bool free_slip = (fabs(x - inputData_->physDomain.min[0]) < eps);

//             bool no_penetraction = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
//                                    (fabs(y - inputData_->physDomain.max[1]) < eps);

//             if (free_slip) {
//                 b.addDirichlet(NSNodeData::VEL_X,
//                                InletParabolicVelocity(y,inputData_->physDomain.max[1] - inputData_->physDomain.min[1],1));
//                 b.addDirichlet(NSNodeData::VEL_Y, 0);

//             }
//             else if (no_penetraction){
//                 b.addDirichlet(NSNodeData::VEL_X, 0);
//                 b.addDirichlet(NSNodeData::VEL_Y, 0);
//             }

//             if (fabs(x - inputData_->physDomain.max[0]) < eps) {
//                 b.addDirichlet(NSNodeData::PRESSURE, 0);
//             }
// #endif
#if (DIM == 3)

        double z = pos.z();
//        std::cout<<"free slip\n";
        bool free_slip = (fabs(z - inputData_->physDomain.min[2]) < eps);
//        std::cout<<"The min z is "<<inputData_->physDomain.min[2]<<"\n";
        bool z_max = (fabs(z - inputData_->physDomain.max[2]) < eps);
//        std::cout<<"The max z is "<<inputData_->physDomain.max[2]<<"\n";

//std::cout << "z = " << z << "\n";

        // bool no_slip = (fabs(y - inputData_->physDomain.min[1]) < eps) ||
        //        (fabs(y - inputData_->physDomain.max[1]) < eps) ||
        //        (fabs(x - inputData_->physDomain.min[0]) < eps) ||
        //        (fabs(x - inputData_->physDomain.max[0]) < eps);
        // compute the distance from the center of the domain so

        double distance_from_center = sqrt(pow(x, 2) + pow(y, 2));


        if(distance_from_center-0.6 > eps){
            b.addDirichlet(NSNodeData::VEL_X, 0);
            b.addDirichlet(NSNodeData::VEL_Y, 0);
            b.addDirichlet(NSNodeData::VEL_Z, 0);
        }
        else if (free_slip) {
//            std::cout<<"free slip, we are here chenghau\n";
            b.addDirichlet(NSNodeData::VEL_X, 0);
            b.addDirichlet(NSNodeData::VEL_Y, 0);
            b.addDirichlet(NSNodeData::VEL_Z, 1);
        }
        else if (z_max) {
            b.addDirichlet(NSNodeData::PRESSURE, 0);
        }



#endif
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

    double NSBoundaryConditions::InletParabolicVelocity(double y, double max_y, double max_velocity){
        return max_velocity * 4 * (y/max_y - pow(y/max_y,2));
    }
double NSBoundaryConditions::InletParabolicVelocity_3D(double y, double max_y, double z, double max_z, double max_velocity) {
    // Compute the parabolic velocity profile at a given point (y, z)
    double term_y = std::pow(y / max_y, 2)-1;  // Parabolic term for y
    double term_z = std::pow(z / max_z, 2)-1;  // Parabolic term for z

    // Calculate the profile value based on the parabolic profile formula
    double profile = term_y*term_z;

    return profile*max_velocity;
}


#endif //DENDRITEKT_NSBCSETUP_H
