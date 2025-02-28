//
// Created by maksbh on 6/15/20.
//

#ifndef DENDRITEKT_NSBCSETUP_H
#define DENDRITEKT_NSBCSETUP_H

#include <TimeInfo.h>
#include <PETSc/BoundaryConditions.h>
#include <DataTypes.h>
#include "NSPNPInputData.h"
#include "NSPNPNodeData.hpp"

/**
 * Class for NS Boundary conditions.
 */
class NSBoundaryConditions {
private:
    const NSPNPInputData *inputData_; /// Input Data
    const TimeInfo *ti_; /// Time Info
    AnalyticFunction analyticFunction_ = nullptr; /// Analytical function for BC setup for MMS
    /**
     * @brief returns the boundary condition for velocity for MMS case
     * @param [out] b Boundary
     * @param pos position
     */
    void returnMMSmomentumBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    /**
     * @brief returns the boundary condition for pressure for MMS case
     * @param [out] b Boundary
     * @param pos position
     */
    void returnMMSpressureBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    /// Not tested
    void returnLDCBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    /// PNP BC functions
    void returnMMSPNPBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    void returnPNPBoundaryEOF(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
    void returnNSBoundaryEOF(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);



public:
    /**
     * @brief constructor
     * @param idata input Data
     * @param ti Time Info
     */
    NSBoundaryConditions(const NSPNPInputData *idata, const TimeInfo *ti);

    /**
     * @brief sets the Analytical function for MMS case.
     * @param f Analytical function
     */
    void setAnalyticalFunction(const AnalyticFunction &f);

    /**
     * @brief  returns the boundary condition for velocity
     * @param [out] b  boundary condition for pressure
     * @param pos position
     */
    void getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    /**
     * @brief returns the boundary condition for all degrees of freedom for PNP
     * @param [out] b  boundary condition for all deg of freedom
     * @param pos position
     */
    void getPNPBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);

    /**
     * @brief  returns the boundary condition for pressure
     * @param [out] b  boundary condition for pressure
     * @param pos position
     */

    void getPressureBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos);
};

NSBoundaryConditions::NSBoundaryConditions(const NSPNPInputData *idata, const TimeInfo *ti) {
    inputData_ = idata;
    ti_ = ti;
}

void NSBoundaryConditions::setAnalyticalFunction(const AnalyticFunction &f) {

    analyticFunction_ = f;

}

void NSBoundaryConditions::getMomentumBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    if (inputData_->ifMMS) {

        if ((analyticFunction_ == nullptr)) {
            throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
        }
        returnMMSmomentumBoundary(b, pos);
    } else {
        returnNSBoundaryEOF(b, pos);
//        returnLDCBoundary(b, pos);
    }
}

void NSBoundaryConditions::getPNPBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    if (inputData_->ifMMS) {
        if ((analyticFunction_ == nullptr)) {
            throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
        }
        returnMMSPNPBoundary(b, pos);
    } else {
        returnPNPBoundaryEOF(b, pos);
    }
}

void NSBoundaryConditions::getPressureBoundaryCondition(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    if (inputData_->ifMMS) {

        if ((analyticFunction_ == nullptr)) {
            throw TALYFEMLIB::TALYException() << "Please set the analytical solution ";
        }
        returnMMSpressureBoundary(b, pos);
    } else {
        returnLDCBoundary(b, pos);
    }
}

void NSBoundaryConditions::returnLDCBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    //TODO: Fix this
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();
    Point<DIM> domainMin(inputData_->physDomain.min);
    Point<DIM> domainMax(inputData_->physDomain.max);
#if(DIM == 3)
    double z = pos.z();
      if(FEQUALS(x,0.0) and (FEQUALS(y,0.0)) and (FEQUALS(z,0.0))){
        b.addDirichlet(NSPNPNodeData::PRESSURE,0.0);
      }

      bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
          (fabs(y - domainMin.x(1)) < eps) ||
          (fabs(z - domainMin.x(2)) < eps) ||
          (fabs(x - domainMax.x(0)) < eps) ||
          (fabs(z - domainMax.x(2)) < eps);
      if(fabs(y - domainMax.x(1)) < eps){
        b.addDirichlet(NSPNPNodeData::VEL_X,-1);
        b.addDirichlet(NSPNPNodeData::VEL_Y,0);
        b.addDirichlet(NSPNPNodeData::VEL_Z,0);
      }
      else if(no_slip){
        b.addDirichlet(NSPNPNodeData::VEL_X,0.0);
        b.addDirichlet(NSPNPNodeData::VEL_Y,0.0);
        b.addDirichlet(NSPNPNodeData::VEL_Z,0);
      }
#endif
#if(DIM == 2)
    if (FEQUALS(x, 0.0) and (FEQUALS(y, 0.0))) {
        b.addDirichlet(NSPNPNodeData::PRESSURE, 0.0);
    }
    bool no_slip = (fabs(x - domainMin.x(0)) < eps) ||
                   (fabs(y - domainMin.x(1)) < eps) ||
                   (fabs(x - domainMax.x(0)) < eps);
    if (fabs(y - domainMax.x(1)) < eps) {

        b.addDirichlet(NSPNPNodeData::VEL_X, 1);
        b.addDirichlet(NSPNPNodeData::VEL_Y, 0);
    } else if (no_slip) {
        b.addDirichlet(NSPNPNodeData::VEL_X, 0.0);
        b.addDirichlet(NSPNPNodeData::VEL_Y, 0.0);
    }

#endif
}

void NSBoundaryConditions::returnMMSmomentumBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    double x = pos.x();
    double y = pos.y();
    if ((FEQUALS(x, 0.0)) or (FEQUALS(y, 0.0)) or (FEQUALS(x, 1.0)) or (FEQUALS(y, 1.0))) {
        b.addDirichlet(NSPNPNodeData::VEL_X, analyticFunction_(pos, NSPNPNodeData::VEL_X,
                                                               ti_->getCurrentTime() + ti_->getCurrentStep()));
        b.addDirichlet(NSPNPNodeData::VEL_Y, analyticFunction_(pos, NSPNPNodeData::VEL_Y,
                                                               ti_->getCurrentTime() + ti_->getCurrentStep()));
    }
    /// Apply corner condition for pressure for PSPG solver
    if ((FEQUALS(x, 0.0)) and (FEQUALS(y, 0.0))) {
        b.addDirichlet(NSPNPNodeData::PRESSURE, analyticFunction_(pos, NSPNPNodeData::PRESSURE,
                                                                  ti_->getCurrentTime() + ti_->getCurrentStep()));
    }
}

void NSBoundaryConditions::returnMMSpressureBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    double x = pos.x();
    double y = pos.y();
    if ((FEQUALS(x, 0.0)) and (FEQUALS(y, 0.0))) {
        b.addDirichlet(0, analyticFunction_(pos, NSPNPNodeData::PRESSURE, ti_->getCurrentTime()));
    }
}

void NSBoundaryConditions::returnMMSPNPBoundary(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    double x = pos.x();
    double y = pos.y();
    if ((FEQUALS(x, 0.0)) or (FEQUALS(y, 0.0)) or (FEQUALS(x, 1.0)) or (FEQUALS(y, 1.0))) {
        b.addDirichlet(0, analyticFunction_(pos, 0, ti_->getCurrentTime() + ti_->getCurrentStep()));
        b.addDirichlet(1, analyticFunction_(pos, 1, ti_->getCurrentTime() + ti_->getCurrentStep()));
        b.addDirichlet(2, analyticFunction_(pos, 2, ti_->getCurrentTime() + ti_->getCurrentStep()));
    }
}

void NSBoundaryConditions::returnPNPBoundaryEOF(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();
    Point<DIM> domainMin(inputData_->physDomain.min);
    Point<DIM> domainMax(inputData_->physDomain.max);


//    /// bottom
//    if (fabs(y - domainMin.x(1)) < eps){
//        /// Dirichlet boundary conditions for potential
//        b.addDirichlet(0, -2.32);
////        /// Dirichlet boundary conditions for C1
////        b.addDirichlet(1, 0.0);
////        /// Dirichlet boundary conditions for C2
////        b.addDirichlet(2, 0.0);
//    }
//
//    /// Top
//    else if (fabs(y - domainMax.x(1)) < eps){
//        /// Dirichlet boundary conditions for potential
//        b.addDirichlet(0, -2.32);
////                /// Dirichlet boundary conditions for C1
////        b.addDirichlet(1, 0.0);
////        /// Dirichlet boundary conditions for C2
////        b.addDirichlet(2, 0.0);
//    }
//
//    /// West
//    else if (fabs(x - domainMin.x(0)) < eps){
//        /// Dirichlet boundary conditions for potential
//        b.addDirichlet(0, 0.0);
//        /// Dirichlet boundary conditions for C1
//        b.addDirichlet(1, 1.0);
//        /// Dirichlet boundary conditions for C2
//        b.addDirichlet(2, 1.0);
//    }
//
//    /// East
//    else if (fabs(x - domainMax.x(0)) < eps){
//        /// Dirichlet boundary conditions for potential
//        b.addDirichlet(0, 0.39);
//        /// Dirichlet boundary conditions for C1
//        b.addDirichlet(1, 1.0);
//        /// Dirichlet boundary conditions for C2
//        b.addDirichlet(2, 1.0);
//    }


//    printf("Top %f\n", domainMax.x(1));



#if(DIM == 2)
    bool inlet = (fabs(x - domainMin.x(0)) < eps);
    bool outlet = (fabs(x - domainMax.x(0)) < eps);
    //bool walls = (fabs(y - domainMin.x(1)) < eps) || (fabs(y - domainMax.x(1)) < eps);
    bool walls = (fabs(y - domainMin.x(1)) < eps && !inlet && !outlet);
    bool wallSegmentInit = walls and ( x < 0.2);
    bool wallSegmentEnd = walls and ( x > (fabs(domainMax.x(0) - 0.2)) );
    bool bulk = (fabs(y - domainMax.x(1)) < eps);

    /// Distribution of wall in segments
    ///  --wallSegmentinit-------|---wallSegmentMiddle-----|----wallSegmentEnd----|
    bool wallSegmentMiddle = walls and ( x >= 0.3) and ( x <= (fabs(domainMax.x(0) - 0.3)) ) ;

    if (inlet) {
        if (ti_->getCurrentTime () < inputData_->restartTimeEof) {
        } else {
            b.addDirichlet (1, 1.0);
            b.addDirichlet (2, 1.0);
        }
    }
    if (outlet) {
        if (ti_->getCurrentTime () < inputData_->restartTimeEof) {
        } else {
            b.addDirichlet (1, 1.0); /// concentration1
            b.addDirichlet (2, 1.0); /// concentration2
        }
    }

    DENDRITE_REAL potentialOnWall = 1.0;
    if (ti_->getCurrentTime() < inputData_->restartTimeEof){
        potentialOnWall = inputData_->eofPNPonlybottomPotential;
        //potentialOnWall = -2.707;
    } else{
        //potentialOnWall = -0.03868*x-2.39756; /// Only for channel length of 8
        potentialOnWall = inputData_->eofbottomPotetialSlope * x + inputData_->eofbottomPotetialIntercept;
    }
    //potentialOnWall = -0.03868*x-2.39756;
    if(wallSegmentMiddle){
        b.addDirichlet(0, potentialOnWall);
    }
    if(walls){
        b.addDirichlet(0, potentialOnWall);
    }

    if(inlet){ /// Inlet with corners
        if (ti_->getCurrentTime() < inputData_->restartTimeEof){
        } else{
            b.addDirichlet (0, inputData_->eofInletPotential);
            //b.addDirichlet(0, 0.30944); /// This value chosen such a that we can compare to an analytical solution
        }
    }
    if(outlet){ /// Outlet with corners
        if (ti_->getCurrentTime() < inputData_->restartTimeEof){
        } else{
            b.addDirichlet(0, 0.0); /// potential
        }
    }

    if(bulk){
        //DENDRITE_REAL potentialBulk = -0.03868 * x + 0.30944; /// Only for channel length of 8
        DENDRITE_REAL potentialBulk = 0.0; /// Only for channel length of 8
        if (ti_->getCurrentTime() < inputData_->restartTimeEof){
            b.addDirichlet(0, potentialBulk);
            b.addDirichlet(1, 1.0);
            b.addDirichlet(2, 1.0);
        } else{
        }
    }
#endif





}

void NSBoundaryConditions::returnNSBoundaryEOF(PETSc::Boundary &b, const TALYFEMLIB::ZEROPTV &pos) {
    //TODO: Fix this
    static constexpr double eps = 1e-14;
    double x = pos.x();
    double y = pos.y();
    Point<DIM> domainMin(inputData_->physDomain.min);
    Point<DIM> domainMax(inputData_->physDomain.max);


#if(DIM == 2)

//    if (FEQUALS(x, 0.0) and (FEQUALS(y, 0.0))) {
//        b.addDirichlet(NSPNPNodeData::PRESSURE, 0.0);
//    }




    bool free_slip = (fabs(x - domainMin.x(0)) < eps) ||
                    (fabs(x - domainMax.x(0)) < eps);

    bool no_slip = (fabs(y - domainMin.x(1)) < eps) ||
                   (fabs(y - domainMax.x(1)) < eps);


     if (free_slip) {
        b.addDirichlet(NSPNPNodeData::VEL_Y, 0.0);
        b.addDirichlet(NSPNPNodeData::PRESSURE, 0.0);
    }

    else if (no_slip) {
        b.addDirichlet(NSPNPNodeData::VEL_X, 0.0);
        b.addDirichlet(NSPNPNodeData::VEL_Y, 0.0);
    }
#endif
}






#endif //DENDRITEKT_NSBCSETUP_H
