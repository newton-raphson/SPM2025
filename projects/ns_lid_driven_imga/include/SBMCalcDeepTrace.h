//
// Created by pc-14-2 on 3/22/24.
//

#ifndef NS_LD_CONS_IMGA_DEEPTRACE_SBMCALCDEEPTRACE_H
#define NS_LD_CONS_IMGA_DEEPTRACE_SBMCALCDEEPTRACE_H

#include "DeepTrace.h"
#include "SBMCalc.h"
#include <IMGA/IMGA.h>
#include "nanoflann.hpp"
#include "NSInputData.h"
#include "talyfem/grid/femelm.h"

// construct a kd-tree index:
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
        PointCloud<double>, 3 /* dim */>;

using my_kd_tree_ptr = std::unique_ptr<my_kd_tree_t>;

class SBMCalcDeepTrace: public SBMCalc
{
private:
    DeepTrace deepTrace_;
    TALYFEMLIB::FEMElm fe_ = TALYFEMLIB::FEMElm(nullptr);

public:
    SBMCalcDeepTrace( const TALYFEMLIB::FEMElm &fe, NSInputData *idata, const IMGA *imga, std::vector<my_kd_tree_t *> kd_trees, const int CarvedOutGeomID);

    void PtDist2Geo(const TALYFEMLIB::ZEROPTV &pt, double (&d)[DIM])
    {


    auto * coordinates = new double[3];

    // Assign x, y, and z values to the array
    coordinates[0] = pt.x();
    coordinates[1] = pt.y();
    #if (DIM==2)
        coordinates[2]   = 0.0;

    #endif
    #if (DIM==3)
        coordinates[2]   = pt.z();

    #endif


    deepTrace_.compute_distance_vector(coordinates,d);

    delete[] coordinates;

    }

    void Dist2Geo(double (&d)[DIM]);


};

SBMCalcDeepTrace::SBMCalcDeepTrace(const FEMElm &fe, NSInputData *idata, const IMGA *imga,
                                   std::vector<my_kd_tree_t *> kd_trees, const int CarvedOutGeomID)
        :SBMCalc(fe, idata, imga, kd_trees, CarvedOutGeomID), deepTrace_(idata->model_path,idata->scale), fe_(fe){}



//#if (DIM==2)
//void SBMCalcDeepTrace::Dist2Geo(double (&d)[2]) {
//    const TALYFEMLIB::ZEROPTV pt = fe_.position();
//    double d_temp[3]; // Array of doubles, not a reference to an array
//
//    PtDist2Geo(pt, d_temp); // Pass the array directly, not a reference
//
//    // Copy the first two values
//    d[0] = d_temp[0];
//    d[1] = d_temp[1];
//}
//
//
//#endif

//#if (DIM==3)
void SBMCalcDeepTrace::Dist2Geo(double (&d)[DIM]) {
    const TALYFEMLIB::ZEROPTV pt = fe_.position();

    PtDist2Geo(pt, d);

}

//#endif
#endif //NS_LD_CONS_IMGA_DEEPTRACE_SBMCALCDEEPTRACE_H
