#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"
#include "NSDataStructure.h"
#include "SBMCalc.h"
#include "SBMCalcDeepTrace.h"
#include <unordered_map>


using namespace TALYFEMLIB;

class NSEquation : public TALYFEMLIB::CEquation<NSNodeData> {
    TALYFEMLIB::ZEROPTV forcing_;
    TALYFEMLIB::ZEROPTV forcingPre_;
    DENDRITE_REAL theta_ = 1.0;
    NSInputData *idata_;

    const IMGA *imga_;
    std::vector<my_kd_tree_t*> kd_trees_;

    /// we hard code here, if want, we can change this inside the config.txt
    const double Cb_f;
    const double Ci_f = 36;
    const double orderOfhb = 1.0;
    const int BetaNS = 1.0;
    const double lambda = - 2.0 / 3.0;

    const double * rampingCoefficient_ = nullptr;

    DENDRITE_REAL normalDistance(const TALYFEMLIB::FEMElm &fe,
                                 const TALYFEMLIB::ZEROPTV &normal,
                                 const TALYFEMLIB::ZEROPTV &h)
    {
        double max_h = std::max(std::max(h.x(), h.y()), h.z());

        TALYFEMLIB::ZEROPTV root_node;
        fe.grid()->GetCoord(root_node, 0);
        int numNode = pow(idata_->elemOrder + 1, DIM);
        std::vector<TALYFEMLIB::ZEROPTV> elm_node(numNode);
#ifdef ENABLE_3D
        for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
        {
            for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
            {
                for (int p3 = 0; p3 < idata_->elemOrder + 1; p3++)
                {
                    elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3] = root_node;
                    elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].x() += p3 * h.x() / idata_->elemOrder;
                    elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].y() += p2 * h.y() / idata_->elemOrder;
                    elm_node[p1 * (idata_->elemOrder + 1) * (idata_->elemOrder + 1) + p2 * (idata_->elemOrder + 1) + p3].z() += p1 * h.y() / idata_->elemOrder;
                }
            }
        }
#else
        for (int p1 = 0; p1 < idata_->elemOrder + 1; p1++)
        {
            for (int p2 = 0; p2 < idata_->elemOrder + 1; p2++)
            {
                elm_node[p1 * (idata_->elemOrder + 1) + p2] = root_node;
                elm_node[p1 * (idata_->elemOrder + 1) + p2].x() += p2 * h.x() / idata_->elemOrder;
                elm_node[p1 * (idata_->elemOrder + 1) + p2].y() += p1 * h.y() / idata_->elemOrder;
            }
        }
#endif
        std::vector<double> hb(numNode);
        for (unsigned int i = 0; i < numNode; i++)
        {
            TALYFEMLIB::ZEROPTV vec1 = fe.position() - elm_node[i];
            hb[i] = vec1.innerProduct(normal) > 0 ? vec1.innerProduct(normal) : 0;
        }
        auto maxhb = std::max_element(hb.begin(), hb.end());
        if (*maxhb < max_h * 1e-2)
        {
            return max_h * 1e-2;
        }
        return *maxhb;
    }

    bool GpOnDomainWall(const TALYFEMLIB::FEMElm &fe)
    {
        /// All the possible cases
        static const double eps = 1e-7;

        bool x_minus_wall = fabs(fe.position().x() - idata_->physDomain.min[0]) < eps;
        bool y_minus_wall = fabs(fe.position().y() - idata_->physDomain.min[1]) < eps;
        bool x_max_wall = fabs(fe.position().x() - idata_->physDomain.max[0]) < eps;
        bool y_max_wall = fabs(fe.position().y() - idata_->physDomain.max[1]) < eps;
        bool z_minus_wall = false;
        bool z_max_wall = false;
#if(DIM == 3)
        z_minus_wall = fabs(fe.position().z() - idata_->physDomain.min[2]) < eps;
    z_max_wall = fabs(fe.position().z() - idata_->physDomain.max[2]) < eps;
#endif

        return x_minus_wall || y_minus_wall || x_max_wall || y_max_wall || z_minus_wall || z_max_wall;
    }

    void calcForcing(TALYFEMLIB::ZEROPTV &forcing, const TALYFEMLIB::ZEROPTV &location, const double time) {
        double Re(idata_->Re);
        forcing.x() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI)
                      + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI) * 2.0
                      + ((M_PI * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) *
                         2.0) / Re
                      + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
                        (cos(location.y() * M_PI)
                         * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0)
                      + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
                        (cos(location.x() * M_PI)
                         * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0);

        forcing.y() = M_PI * cos(time * M_PI * 2.0) * cos(location.x() * M_PI) * sin(location.y() * M_PI) * (-2.0)
                      + M_PI * cos(time * M_PI * 2.0) * cos(location.y() * M_PI) * sin(location.x() * M_PI)
                      - ((M_PI * M_PI) * cos(location.x() * M_PI) * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) *
                         2.0) / Re
                      + M_PI * cos(location.x() * M_PI) * cos(location.y() * M_PI) * sin(time * M_PI * 2.0) *
                        (cos(location.x() * M_PI)
                         * sin(time * M_PI * 2.0) * sin(location.y() * M_PI) - 2.0)
                      + M_PI * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) * sin(location.y() * M_PI) *
                        (cos(location.y() * M_PI)
                         * sin(time * M_PI * 2.0) * sin(location.x() * M_PI) + 2.0);
        forcing.z() = 0;

    }

#pragma mark Distance Calc

    void DistanceFuncChecking(const TALYFEMLIB::FEMElm &fe,double (&d)[DIM])
    {
        for (int dim = 0; dim < DIM; dim++)
        {
            if (fabs(d[dim]) > 2.0 * ElementSize(fe))
            {

//                        std::ofstream fout("GP_issue.txt", std::ios::app);
//#if (DIM == 2)
//        fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().x() + d[0] << "," << fe.position().y() + d[1]
//                << "," << d[0] << "," << d[1]
//        << "\n";
//#endif
//
//#if (DIM == 3)
//        fout << fe.position().x() << "," << fe.position().y() <<"," << fe.position().z() << "\n";
//#endif
//        fout.close();

                std::cout << "ElementSize(fe) = " << ElementSize(fe) << "\n";
                std::cout << "[Ae_Carved]there may be something wrong. Please check the config  and your change to the code again!\n";
#if (DIM == 3)
                std::cout << "d = " << d[0] << "," << d[1] << "," << d[2] << "\n";
#endif
#if (DIM == 2)
                std::cout << "d = " << d[0] << "," << d[1] << "\n";
#endif

                std::cout << "x_surrogate = " << fe.position().x() << ","
                          << "y_surrogate = " << fe.position().y() << ","
                          << "z_surrogate = " << fe.position().z() << "\n";
                std::cout << "x = " << fe.position().x() + d[0] << ","
                          << "y = " << fe.position().y() + d[1] << ","
                          << "z = " << fe.position().z() + d[2] << "\n";
            }
        }
    }

    void CalcDistFast(const int CarvedOutGeomID, const TALYFEMLIB::ZEROPTV &pt, const TALYFEMLIB::FEMElm &fe, double (&d)[DIM], bool CopyFromPrevious)
    {
        // std::cout << "CalcDistFast_NS\n";
        if (CopyFromPrevious)
        {

            std::string key = queryToStringKey(pt);


            auto it = queryToIndexMap_.find(key);

            if (it != queryToIndexMap_.end()) {
//                std::cout << "find the key\n";

                for (int dim = 0; dim <DIM;dim++) {
                    d[dim] = it->second[dim];
                }
            } else {
//                std::cout << "do not find the key\n";

//                SBMCalc sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
                SBMCalcDeepTrace sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
                // This is a new node. Compute things and store them.
                sbmCalc.Dist2Geo(d);
                std::vector<double> d_vector(DIM);
                for (int dim = 0; dim < DIM ; dim ++) {
                    d_vector[dim] = d[dim];
                }
                queryToIndexMap_[key] = d_vector;
            }

            // Only calculate distance and normal once
            // But, there is no node ID available, so need to find it.
//            unsigned int nodeID = 0;
//            bool found_it = false;
//
//            while (distanceSet[nodeID])
//            {
//                if (fabs(pt.x() - surfaceNodeCoords[nodeID * DIM]) < 1e-14)
//                {
//                    if (fabs(pt.y() - surfaceNodeCoords[nodeID * DIM + 1]) < 1e-14)
//                    {
//#if (DIM == 3)
//                        if (fabs(pt.z() - surfaceNodeCoords[nodeID * DIM + 2]) < 1e-14)
//#endif
//                        {
//                            // This is the node!
//                            found_it = true;
//                            break;
//                        }
//                    }
//                }
//                // Not the node.
//                nodeID++;
//            }
//
//            /*
//             * TODO: passing CarvedOutGeomID inside and the SBM geom is storing in imga_->geometry() to get different SBMGeo for different carved-out places
//             */
//            // Now we either found a set node, or need to compute a new one.
//            if (found_it)
//            {
//                for (int i = 0; i < DIM; i++)
//                {
//                    d[i] = distanceToBdry[nodeID * DIM + i];
//                }
//                // below is for neumann
//                // trueNormal = TALYFEMLIB::ZEROPTV(geoNormal[nodeID*DIM], geoNormal[nodeID*DIM+1], geoNormal[nodeID*DIM+2]);
//            }
//            else
//            {
//                SBMCalc sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
//                // This is a new node. Compute things and store them.
//                sbmCalc.Dist2Geo(d);
//                // below is for neumann
//                // sbmCalc.NormalofGeo(trueNormal, d);
//                // Save them
//                for (int i = 0; i < DIM; i++)
//                {
//                    distanceToBdry[nodeID * DIM + i] = d[i];
//                    // need to set for neumann
//                    // geoNormal[nodeID*DIM+i] = trueNormal.data()[i];
//                    surfaceNodeCoords[nodeID * DIM + i] = pt.data()[i];
//                }
//                distanceSet[nodeID] = true;
//            }
        }
        else
        {
            SBMCalcDeepTrace sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
            sbmCalc.Dist2Geo(d);
        }

//#ifndef NDEBUG
//        DistanceFuncChecking(fe,d);
//#endif
//        std::ofstream fout("GP_issue.txt", std::ios::app);
//#if (DIM == 2)
//        fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().x() + d[0] << "," << fe.position().y() + d[1]
//             << "," << d[0] << "," << d[1]
//             << "\n";
//#endif
//
//#if (DIM == 3)
//        fout << fe.position().x() << "," << fe.position().y() <<"," << fe.position().z() << "\n";
//#endif
//        fout.close();

    }


    int NearestGeo(const TALYFEMLIB::FEMElm &fe) {
        const TALYFEMLIB::ZEROPTV pt = fe.position();
        double Min_Dist = std::numeric_limits<double>::max();
        int PickGeomID = 0;

        for (int CarvedOutGeomID = 0;
             CarvedOutGeomID < imga_->getGeometries().size();
             ++CarvedOutGeomID) {

//            std::cout << "CarvedOutGeomID = " << CarvedOutGeomID << "\n";

            double d[DIM];
            CalcDistFast(CarvedOutGeomID, pt, fe, d, false);

            double dist_square = std::inner_product(d, d + DIM, d, 0.0);

//            std::cout << "CarvedOutGeomID = " << CarvedOutGeomID << "\n";
//            std::cout << "dist_square = " << dist_square << "\n";
//            std::cout << "pt = "; pt.print(); std::cout << "\n";
//            std::cout << "------------\n";

            if (dist_square < Min_Dist) {
                Min_Dist = dist_square;
                PickGeomID = CarvedOutGeomID;
            }
        }

        return PickGeomID;
    }

    void CalculateDistances(const TALYFEMLIB::FEMElm &fe,
                            const int CarvedOutGeomID,
                            double (&d)[DIM]) {
        const TALYFEMLIB::ZEROPTV pt = fe.position();

        bool UsePrevDistance = true;

        CalcDistFast(CarvedOutGeomID, pt, fe, d, UsePrevDistance);
    }


    DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe)
    {
        if (imga_->getIBMMethod() == SBM)
        {
            return pow((pow(2, DIM) * fe.volume_jacc()), (double)1 / DIM);
        }
        else
        {
            return pow((pow(2, DIM) * fe.jacc()), (double)1 / DIM);
        }
    }


#pragma mark SBM

    void SBM_Ae(TALYFEMLIB::ZeroMatrix<double> &Ae,
                const TALYFEMLIB::FEMElm &fe,
                const int nbf,
                const int ndof,
                const int nsd,
                const double detSideJxW,
                const double weakBCpenaltyParameter,
                const double Coe_diff,
                const double *d,
                const bool PenaltyOnly = false) {

//        std::cout <<"Coef_diff  (SBM_Ae)  = " << Coe_diff << "\n";
//        std::ofstream fout("DistanceFunc_Ae.txt",std::ios::app);
//        fout << fe.position().x() << " " << fe.position().y() << "\n";
//        fout.close();
//        fout << d[0] << " " << d[1] << " ";
//        fout << fe.position().x()+ d[0] << " " << fe.position().y()+d[1] << "\n";

        double RampingTerm_Coe = (rampingCoefficient_)?  *rampingCoefficient_ : idata_->RampingTerm_Coe.value_at(t_-dt_);

//                std::cout << "t_ = " << t_ << "\n";
//        std::cout << "RampingTerm_Coe = " << RampingTerm_Coe << "\n";

        TALYFEMLIB::ZEROPTV u;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = p_data_->valueFEM(fe, i);
            // u(i) = p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
        }

        double inflow_g = 0.0;
        for (int m = 0; m < nsd; m++)
        {
            inflow_g += u(m) * fe.surface()->normal()(m);
        }

        DENDRITE_REAL secondOrderTerm_a(0), secondOrderTerm_b(0);
        for (int a = 0; a < nbf; a++) {

#if (DIM == 2)
            if (idata_->elemOrder == 2) {
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                     d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                                    2;
            } else {
                secondOrderTerm_a = 0;
            }
#endif

#if (DIM == 3)
            if (idata_->elemOrder == 2)
                    {

                        secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
                    }
                    else
                    {
                        secondOrderTerm_a = 0;
                    }
#endif
            TALYFEMLIB::ZEROPTV gradWdotd;
            TALYFEMLIB::ZEROPTV gradWdotn;

            for (int dim_V = 0; dim_V < nsd; dim_V++) {
                gradWdotn(dim_V) = 0.0;
                gradWdotd(dim_V) = 0.0;
                for (int dim_Coord = 0; dim_Coord < nsd; dim_Coord++) {
                    gradWdotn(dim_V) += fe.dN(a, dim_Coord) * fe.surface()->normal()(dim_Coord);
                    gradWdotd(dim_V) += fe.dN(a, dim_Coord) * d[dim_Coord];
                }
            }

            for (int b = 0; b < nbf; b++) {

                TALYFEMLIB::ZEROPTV gradDELUdotd;
                TALYFEMLIB::ZEROPTV gradDELUdotn;

#if (DIM == 2)
                if (idata_->elemOrder == 2) {
                    secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1]) +
                                         d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1])) /
                                        2;
                } else {
                    secondOrderTerm_b = 0;
                }
#endif

#if (DIM == 3)
                if (idata_->elemOrder == 2)
                        {

                            secondOrderTerm_b = (d[0] * (fe.d2N(b, 0, 0) * d[0] + fe.d2N(b, 0, 1) * d[1] + fe.d2N(b, 0, 2) * d[2]) + d[1] * (fe.d2N(b, 1, 0) * d[0] + fe.d2N(b, 1, 1) * d[1] + fe.d2N(b, 1, 2) * d[2]) + d[2] * (fe.d2N(b, 2, 0) * d[0] + fe.d2N(b, 2, 1) * d[1] + fe.d2N(b, 2, 2) * d[2])) / 2;
                        }
                        else
                        {
                            secondOrderTerm_b = 0;
                        }
#endif

                for (int dim_V = 0; dim_V < nsd; dim_V++) {
                    gradDELUdotn(dim_V) = 0.0;
                    gradDELUdotd(dim_V) = 0.0;
                    for (int dim_Coord = 0; dim_Coord < nsd; dim_Coord++) {
                        gradDELUdotn(dim_V) += fe.dN(b, dim_Coord) * fe.surface()->normal()(dim_Coord);
                        gradDELUdotd(dim_V) += fe.dN(b, dim_Coord) * d[dim_Coord];
                    }
                }

                for (int j = 0; j < nsd; j++) {

                    // std::cout<<"gradDELUdotd="<<gradDELUdotd<<"\n";

                    // gradWdotd(j) = 0;
                    // gradDELUdotd(j) = 0;
                    //  gradUPREdotd(j) = 0;

                    Ae(ndof * a + j, ndof * b + j) +=
                            /// diffusion_1
                            -Coe_diff * fe.N(a) * gradDELUdotn(j) * detSideJxW
                            /// adjoint diffusion
                            - Coe_diff * gradWdotn(j) * (fe.N(b) + gradDELUdotd(j) + secondOrderTerm_b) *
                              detSideJxW
                            /// penalty
                            + weakBCpenaltyParameter * (fe.N(a) + gradWdotd(j) + secondOrderTerm_a) *
                              (fe.N(b) + gradDELUdotd(j) + secondOrderTerm_b) *
                              detSideJxW;

                    /// diffusion_2 (new term)
                    for (int k = 0; k < nsd; k++) {
                        Ae(ndof * a + j, ndof * b + k) +=
                                -Coe_diff * fe.N(a) * fe.dN(b, j) * fe.surface()->normal()(k) * detSideJxW;
                    }

                    /// pressure
                    Ae(ndof * a + j, ndof * b + nsd) +=
                            fe.N(a) * fe.N(b) * fe.surface()->normal()(j) * detSideJxW;

                    /// adjoint diffusion_2 --- grad w^T
                    for (int k = 0; k < nsd; k++) {
                        Ae(ndof * a + k, ndof * b + j) +=
                                -Coe_diff * fe.dN(a, j) * fe.surface()->normal()(k) *
                                (fe.N(b) + gradDELUdotd(j) + secondOrderTerm_b) * detSideJxW;
                    }

                    /// adjoint pressure -> we use symmetric form (same as Yuri Bazilevs)
                    Ae(ndof * a + nsd, ndof * b + j) +=
                            -fe.N(a) * (fe.N(b) + gradDELUdotd(j) + secondOrderTerm_b) *
                            fe.surface()->normal()(j) * detSideJxW;

                    if (inflow_g < 0) {
                        // NOTE: previously have bug here: missing "detSideJxW".
                        Ae(ndof * a + j, ndof * b + j) -=
                                fe.N(a) * (fe.N(b) + gradDELUdotd(j) + secondOrderTerm_b) * u(j) *
                                fe.surface()->normal()(j) *
                                detSideJxW;
                    }


                    if(idata_->BoundaryTerms) {
                        for (int i = 0; i < nsd; i++) {

                            /// advection term

                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                    RampingTerm_Coe * fe.N(a) * (u[i] * fe.N(b)) * fe.surface()->normal()(i) * detSideJxW;

                            Ae((nsd + 1) * a + j, (nsd + 1) * b + i) +=
                                    RampingTerm_Coe * fe.N(a) * (u[j] * fe.N(b)) * fe.surface()->normal()(i) * detSideJxW;

                            /// coarseDiffusion 1
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                    RampingTerm_Coe * Coe_diff * fe.N(a) * fe.dN(b, i)
                                    * fe.surface()->normal()(i) * detSideJxW;

                            /// coarseDiffusion 2
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + i) +=
                                    RampingTerm_Coe * Coe_diff * fe.N(a) * fe.dN(b, j)
                                    * fe.surface()->normal()(i) * detSideJxW;

                            if (idata_->SecondViscosity) {
                                /// coarseDiffusion 3
                                Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                        RampingTerm_Coe * Coe_diff * lambda * fe.N(a) * fe.dN(b, j)
                                        * fe.surface()->normal()(i) * detSideJxW;
                            }
                        }
                    }

                }
            }
        }


    }

    void SBM_be(TALYFEMLIB::ZEROARRAY<double> &be,
                const TALYFEMLIB::FEMElm &fe,
                const int nbf,
                const int ndof,
                const int nsd,
                const double detSideJxW,
                const double weakBCpenaltyParameter,
                const double Coe_diff,
                const double *d,
                const CarvedOutGeom &carved_out_def,
                const bool PenaltyOnly = false) {

//        std::ofstream fout("DistanceFunc_be.txt",std::ios::app);
//        fout << fe.position().x() << " " << fe.position().y() << "\n";
//        fout.close();

//        std::cout <<"Coef_diff  (SBM_be) = " << Coe_diff << "\n";

        double RampingTerm_Coe = (rampingCoefficient_)?  *rampingCoefficient_ : idata_->RampingTerm_Coe.value_at(t_-dt_);

//        std::ofstream fout("DistanceFunc_be.txt",std::ios::app);
//        fout << fe.position().x() << " " << fe.position().y() << " ";
//        fout << d[0] << " " << d[1] << " ";
//        fout << fe.position().x()+ d[0] << " " << fe.position().y()+d[1] << "\n";


        TALYFEMLIB::ZEROPTV u;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = p_data_->valueFEM(fe, i);
            // u(i) = p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
        }

        double p = p_data_->valueFEM(fe, NSNodeData::PRESSURE);
        // double p = p_data_->valueFEM(fe, NSNodeData::PRESSURE_PRE1);

        double inflow_g = 0.0;
        for (int m = 0; m < nsd; m++)
        {
            inflow_g += u(m) * fe.surface()->normal()(m);
        }

#if (DIM == 2)
        std::vector<double> bc_value;
        bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                    carved_out_def.getBC(NSNodeData::VEL_Y)[0]};

#elif (DIM == 3)
        std::vector<double> bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                                        carved_out_def.getBC(NSNodeData::VEL_Y)[0],
                                        carved_out_def.getBC(NSNodeData::VEL_Z)[0]};
#endif


        // Non-linear NS

        TALYFEMLIB::ZeroMatrix<double> du;
        du.redim(nsd, nsd);

        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
            }
        }


        double d2u[DIM][DIM][DIM];
        if (idata_->elemOrder == 2)
        {

            for (int dim = 0; dim < DIM; dim++)
            {
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        d2u[dim][i][j] = this->p_data_->value2DerivativeFEM(fe, dim, i, j);
                    }
                }
            }
        }

        TALYFEMLIB::ZEROPTV gradUPREdotn;
        TALYFEMLIB::ZEROPTV gradUPREdotnTransPose;
        TALYFEMLIB::ZEROPTV gradUPREdotd;
        TALYFEMLIB::ZEROPTV secondOrder_u0; // previous time step solution u0
        for (int i = 0; i < nsd; i++)
        {
            gradUPREdotn(i) = 0.0;
            gradUPREdotnTransPose(i) = 0.0;
            gradUPREdotd(i) = 0.0;
            if (idata_->elemOrder == 2)
            {
#if (DIM == 2)
                secondOrder_u0(i) = (d[0] * (d2u[i][0][0] * d[0] + d2u[i][0][1] * d[1]) +
                                     d[1] * (d2u[i][1][0] * d[0] + d2u[i][1][1] * d[1])) /
                                    2;
#endif
#if (DIM == 3)
                secondOrder_u0(i)= (d[0] * (d2u[i][0][0] * d[0] + d2u[i][0][1] * d[1] + d2u[i][0][2] * d[2])
                                            + d[1]* (d2u[i][1][0] * d[0] + d2u[i][1][1] * d[1] + d2u[i][1][2] * d[2])
                                            + d[2] * (d2u[i][2][0] * d[0] + d2u[i][2][1] * d[1] + d2u[i][2][2] * d[2])) / 2;

#endif
            }

            for (int j = 0; j < nsd; j++)
            {
                gradUPREdotn(i) += du(i, j) * fe.surface()->normal()(j);          //  (dui/dxj)nj
                gradUPREdotnTransPose(i) += du(j, i) * fe.surface()->normal()(j); // (duj/dxi)nj
                gradUPREdotd(i) += du(i, j) * d[j];
            }
        }

        DENDRITE_REAL secondOrderTerm_a(0);

        for (int a = 0; a < nbf; a++)
        {

            TALYFEMLIB::ZEROPTV gradWdotd;
            TALYFEMLIB::ZEROPTV gradWdotn;

#if (DIM == 2)
            if (idata_->elemOrder == 2)
            {
                secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1]) +
                                     d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1])) /
                                    2;
            }
            else
            {
                secondOrderTerm_a = 0;
            }
#endif

#if (DIM == 3)
            if (idata_->elemOrder == 2)
                    {
                        secondOrderTerm_a = (d[0] * (fe.d2N(a, 0, 0) * d[0] + fe.d2N(a, 0, 1) * d[1] + fe.d2N(a, 0, 2) * d[2]) + d[1] * (fe.d2N(a, 1, 0) * d[0] + fe.d2N(a, 1, 1) * d[1] + fe.d2N(a, 1, 2) * d[2]) + d[2] * (fe.d2N(a, 2, 0) * d[0] + fe.d2N(a, 2, 1) * d[1] + fe.d2N(a, 2, 2) * d[2])) / 2;
                    }
                    else
                    {
                        secondOrderTerm_a = 0;
                    }
#endif

            for (int dim_V = 0; dim_V < nsd; dim_V++)
            {
                gradWdotn(dim_V) = 0.0;
                gradWdotd(dim_V) = 0.0;
                for (int dim_Coord = 0; dim_Coord < nsd; dim_Coord++)
                {
                    gradWdotn(dim_V) += fe.dN(a, dim_Coord) * fe.surface()->normal()(dim_Coord);
                    gradWdotd(dim_V) += fe.dN(a, dim_Coord) * d[dim_Coord];
                }
            }

            for (int j = 0; j < nsd; j++)
            {
                // std::cout<<"gradUPREdotd(j)="<<gradUPREdotd(j)<<"\n";
                // gradUPREdotd(j) = 0;
                // gradWdotd(j) = 0;

                be(ndof * a + j) +=
                        /// pressure
                        fe.N(a) * p * fe.surface()->normal()(j) * detSideJxW
                        /// diffusion --- u
                        - Coe_diff * fe.N(a) * gradUPREdotn(j) * detSideJxW
                        /// diffusion --- u^T
                        - Coe_diff * fe.N(a) * gradUPREdotnTransPose(j) * detSideJxW
                        /// adjoint diffusion
                        - Coe_diff * (gradWdotn(j)) *
                          (u(j) + gradUPREdotd(j) + secondOrder_u0(j) - bc_value.at(j)) * detSideJxW
                        /// penalty
                        + weakBCpenaltyParameter * (fe.N(a) + gradWdotd(j) + secondOrderTerm_a) *
                          (u(j) + gradUPREdotd(j) + secondOrder_u0(j) - bc_value.at(j)) * detSideJxW;

                /// adjoint diffusion 2 -- grad w ^T
                for (int k = 0; k < nsd; k++)
                {
                    be(ndof * a + k) -=
                            +Coe_diff * fe.dN(a, j) * fe.surface()->normal()(k) * (u(j) + gradUPREdotd(j) + secondOrder_u0(j) - bc_value.at(j)) * detSideJxW;
                }

                /// adjoint pressure
                be(ndof * a + nsd) +=
                        -fe.N(a) * fe.surface()->normal()(j) * (u(j) + gradUPREdotd(j) + secondOrder_u0(j) - bc_value.at(j)) *
                        detSideJxW;

                if (inflow_g < 0)
                {
                    /// minus is because of the symmetric form
                    be(ndof * a + j) -= fe.N(a) * (u(j) + gradUPREdotd(j) + secondOrder_u0(j) - bc_value.at(j)) * u(j) * fe.surface()->normal()(j) * detSideJxW;
                }

                if(idata_->BoundaryTerms) {
                    for (int i = 0; i < nsd; i++) {
                        /// advention term
                        be((nsd + 1) * a + j) += RampingTerm_Coe* (fe.N(a) * u[j] * u[i])
                                                 * fe.surface()->normal()(i) * detSideJxW;
                        /// Diffusion term 1, 2
                        be((nsd + 1) * a + j) += RampingTerm_Coe* Coe_diff * fe.N(a) * (du(j, i) + du(i, j))
                                                 * fe.surface()->normal()(i) * detSideJxW;

                        if (idata_->SecondViscosity) {
                            /// Diffusion term 3
                            be((nsd + 1) * a + i) += RampingTerm_Coe*Coe_diff * lambda * fe.N(a) * du(j, j)
                                                     * fe.surface()->normal()(i) * detSideJxW;
                        }
                    }
                }


            }
        }
    }

    std::string queryToStringKey(const ZEROPTV &query_pt) {
        return std::to_string(query_pt[0]) + "," + std::to_string(query_pt[1]) + "," + std::to_string(query_pt[2]);
    }

protected:

    struct VMSParams {
        double tauM;
        double tauC;

        VMSParams() : tauM(1.0), tauC(1.0) {}
    };

    /// Eric add to do one time finding
    double *distanceToBdry;
    double *bdryNormal;
    double *surfaceNodeCoords;
    bool *distanceSet;
    unsigned int distanceArraySize;

    std::unordered_map<std::string, std::vector<double>> queryToIndexMap_;



public:
    explicit NSEquation(NSInputData *idata)
            : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints),Cb_f(idata->Cb_f) {
        idata_ = idata;
        if (idata->timeStepping == NSInputData::TIME_STEPPING::CRANK_NICHOLSON) {
            theta_ = 0.5;
        }
    }
    void Solve(double dt, double t) override {
        assert(false);
    }

    void Integrands(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae,
                    TALYFEMLIB::ZEROARRAY<double> &be) override {
        assert(false);
    }


    void calcAe_vms(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, const std::vector<double> &bdf) {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        double Ci_f = 36;
        const int nsd = DIM;
        const int nbf = fe.nbf();
        const double detJxW = fe.detJxW();
        const double dt = dt_;
        const double lambda = - 2.0 / 3.0;

        if (idata_->ifMMS) {
            calcForcing(forcing_, fe.position(), t_);
        }

        VMSParams vms_params = calc_tau(fe);


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

//    MPITimer timer("total Ae time");
//    timer.Start();
//
//
        const double temp_coeff = 1.0;
        const double adv_coeff = 1.0;
        const double press_coeff = 1.0;
        const double diff_coeff = 1.0 / Re;
        const double force_coeff = 1.0;


        /// time_scheme
        double force[nsd];

        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
        }


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
//    const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
//    const double tauC = vms_params.tauC;

        /// field variables
        /** Get u and u_pre from the NodeData
         * NOTE: ValueFEM is defined on NodeData object, therefore it follows
         * indices defined in NodeData subclass.
         * NOTE: Every iteration the solution is copied into NodeData, we use these
         * solved fields from the NodeData object to calculate the be (NS residual)
         * of current
         * iteration
         */


        ////////////////////////////// precalculate velocity terms to simplify the expressions


        double u[nsd];
        double u_pre[nsd];
        double u_pre_pre[nsd];

        double res_M[nsd];
        for (int i = 0; i < nsd; i++) {
            u[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
            u_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);

        }

        /// Define velocity gradient tensor
        double du[nsd][nsd];
        double dp[nsd];
        double ui_ujj[nsd];
        double uj_uij[nsd];
        for (int i = 0; i < nsd; i++) {
            dp[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
            for (int j = 0; j < nsd; j++) {
                du[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
            }
        }


        double d2u_ijj[nsd];
        double d2u_jij[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj[dof] = 0;
            d2u_jij[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dof, dir, dir);
                d2u_jij[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dir, dof, dir);
            }
        }


        double ukk;
        for (int i = 0; i < nsd; i++) {
            ui_ujj[i] = 0;
            uj_uij[i] = 0;
            ukk = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj[i] += u[i] * du[j][j];
                uj_uij[i] += u[j] * du[i][j];
                ukk += du[j][j];
            }
            res_M[i] = temp_coeff * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt + adv_coeff * (ui_ujj[i] + uj_uij[i]) +
                       press_coeff * dp[i] - diff_coeff * (d2u_ijj[i] + d2u_jij[i]) - force_coeff * force[i];

            /// since d2u_jij = d2u_kik we don't need to calculate it again
            if (idata_->SecondViscosity){
                res_M[i] -= diff_coeff * lambda * d2u_jij[i];
            }
        }


        ////////////////////////////////// end of precalculations




        ////////////////////////////////// calculate the derivative of the momentunm residual


        double dres_du_M_Diag[nsd][nbf];
        double dres_du_M_offDiag[nsd][nsd][nbf];

        for (int a = 0; a < nbf; a++) {

            for (int i = 0; i < nsd; i++) {

                dres_du_M_Diag[i][a] = temp_coeff * (fe.N(a) / dt);

                for (int j = 0; j < nsd; j++) {
                    dres_du_M_Diag[i][a] += adv_coeff * ((du[j][j] * fe.N(a)) + (u[j] * fe.dN(a, j)))  -
                                            (diff_coeff * fe.d2N(a, j, j));



                    dres_du_M_offDiag[i][j][a] =
                            adv_coeff * ((u[i] * fe.dN(a, j)) + (du[i][j] * fe.N(a))) -
                            (diff_coeff * fe.d2N(a, i, j));



                } // end of j loop

                if (idata_->SecondViscosity){
                    dres_du_M_Diag[i][a] -= diff_coeff * lambda * fe.d2N(a, i, i);
                }

            } // end of i loop
        } // end of a loop


        ////////////////////////////////// end of calculation of the derivative of the momentunm residual

//    printf("inside the integrands\n");


        for (int a = 0; a < nbf; a++) {
            for (int b = 0; b < nbf; b++) {
                for (int i = 0; i < nsd; i++) {


                    /// first part of momentum equation
                    /// Temporal term -- Term No. in the Document= 1
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += temp_coeff * fe.N(a) * (fe.N(b) / dt) * detJxW; ///ok!



                    for (int j = 0; j < nsd; j++) {


                        /// coarseAdvection 1 -- Term No. in the Document= 2
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * (u[j] * fe.N(b)) * detJxW; ///ok!

                        /// coarseAdvection 2 -- Term No. in the Document= 3
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * (u[i] * fe.N(b)) * detJxW; ///ok!

                        /// fineAdvection 1 -- Term No. in the Document= 4
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * (tauM * res_M[j] * fe.N(b)) * detJxW; ///ok!

                        /// fineAdvection 2 -- Term No. in the Document= 5
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * u[i] * tauM * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) * detJxW; ///ok!


                        /// fineAdvection 3a  -- Term No. in the Document= 6
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * tauM * dres_du_M_Diag[i][b] * u[j] * detJxW; ///ok!

                        /// fineAdvection 3b  -- Term No. in the Document= 7
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * tauM * dres_du_M_offDiag[i][j][b] * u[j] * detJxW; ///ok!


                        /// fineAdvection 4 -- Term No. in the Document= 8
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * (tauM * res_M[i] * fe.N(b)) * detJxW; ///ok!



                        /// fineAdvection 5a -- Term No. in the Document= 9
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * tauM * tauM * dres_du_M_Diag[i][b] * res_M[j] * detJxW; ///ok!


                        /// fineAdvection 5b -- Term No. in the Document= 10
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * tauM * tauM * dres_du_M_offDiag[i][j][b] * res_M[j] * detJxW; ///ok!


                        /// fineAdvection 6  -- Term No. in the Document= 11
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) -= adv_coeff * fe.dN(a,j) * tauM * tauM * res_M[i] * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) *  detJxW; ///ok!

                        /// finePressure 1 -- Term No. in the Document= 12
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += press_coeff * fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW; ///ok!

                        /// coarseDiffusion 1 -- Term No. in the Document= 13
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * fe.dN(a, j) * fe.dN(b, j) * detJxW; ///ok!

                        /// coarseDiffusion 2-- Term No. in the Document= 14
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.dN(a, j) * fe.dN(b, i) * detJxW; ///ok!


                        if (idata_->DiffFineTerms) {
                            /// fineDiffusion 1  -- Term No. in the Document= 15
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.d2N(a, j, i) * tauM * (dres_du_M_Diag[i][b] + dres_du_M_offDiag[j][j][b]) *  detJxW; ///ok!

                            /// second part of momentum equation
                            /// fineDiffusion 3-- Term No. in the Document= 23
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += diff_coeff * fe.d2N(a, i, j) * tauM * fe.dN(b, j) * detJxW;

                        }


                        /// second part of momentum equation
                        /// fineAdvection 7 -- Term No. in the Document= 18
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a, j) * u[i] * tauM * fe.dN(b, j) * detJxW; ///ok!


                        /// fineAdvection 8 -- Term No. in the Document= 19
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a, j) * u[j] * tauM * fe.dN(b, i) * detJxW; ///ok!

                        /// fineAdvection 9 -- Term No. in the Document= 20
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= adv_coeff * fe.dN(a, j) * tauM * tauM * fe.dN(b, i) *
                                                                      res_M[j] * detJxW; ///ok!

                        /// fineAdvection 10 -- Term No. in the Document= 21
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= adv_coeff * fe.dN(a, j) * tauM * tauM * fe.dN(b, j) *
                                                                      res_M[i] * detJxW; ///ok!


                        /// fineContinuityVelocity-a -- Term No. in the Document= 26
                        Ae((nsd + 1) * a + nsd, (nsd + 1) * b + j) += fe.dN(a, i) * tauM * dres_du_M_offDiag[i][j][b] * detJxW; ///ok!



                    } // end of j loop


                    /// coarsePressure 1 -- Term No. in the Document= 22
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= press_coeff * fe.dN(a, i) * fe.N(b) * detJxW; ///ok!

                    /// first part of continuity equation
                    /// coarseContinuityVelocity -- Term No. in the Document= 25
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.N(a) * fe.dN(b, i) * detJxW; ///ok!


                    /// fineContinuityVelocity-b -- Term No. in the Document= 27
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * dres_du_M_Diag[i][b] * detJxW; ///ok!


                    /// second part of continuity equation
                    /// fineContinuityPressure -- Term No. in the Document= 28
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) += fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW; ///ok!



                    if (idata_->SecondViscosity) {

                        /// coarseDiffusion 3-- Term No. in the Document= 16
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * lambda * fe.dN(a, i) * fe.dN(b, i) * detJxW; ///ok!

                        /// fineDiffusion 2-- Term No. in the Document= 17
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * lambda * fe.d2N(a, i, i) * tauM * dres_du_M_Diag[i][b] * detJxW; ///ok!

                        /// second part of continuity equation
                        /// fineDiffusion 4-- Term No. in the Document= 24
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += diff_coeff * lambda * fe.d2N(a, i, i) * tauM * fe.dN(b, i) * detJxW; ///ok!

                    }





                } // end of i loop
            } // end of b loop
        } // end of a loop

//
//    timer.Stop();
//    timer.PrintTotalTimeSeconds();





//    MPITimer timer("total Ae time");
//    timer.Start();



//    timer.Stop();
//    timer.PrintTotalTimeSeconds();


    }

    void calcbe_vms(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, const std::vector<double> &bdf) {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        double Ci_f = 36;
        const int nsd = DIM;
        const int nbf = fe.nbf();
        const double detJxW = fe.detJxW();
        const double dt = dt_;
        const double lambda = - 2.0 / 3.0;

        if (idata_->ifMMS) {
            calcForcing(forcing_, fe.position(), t_);
        }

        VMSParams vms_params = calc_tau(fe);

        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

//    MPITimer timer("total Ae time");
//    timer.Start();
//
//
        const double temp_coeff = 1.0;
        const double adv_coeff = 1.0;
        const double press_coeff = 1.0;
        const double diff_coeff = 1.0 / Re;
        const double force_coeff = 1.0;


        /// time_scheme
        double force[nsd];

        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
        }

        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
//    const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
//    const double tauC = vms_params.tauC;

        /// field variables
        /** Get u and u_pre from the NodeData
         * NOTE: ValueFEM is defined on NodeData object, therefore it follows
         * indices defined in NodeData subclass.
         * NOTE: Every iteration the solution is copied into NodeData, we use these
         * solved fields from the NodeData object to calculate the be (NS residual)
         * of current
         * iteration
         */


        ////////////////////////////// precalculate velocity terms to simplify the expressions


        double u[nsd];
        double u_pre[nsd];
        double u_pre_pre[nsd];
        double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
        for (int i = 0; i < nsd; i++) {
            u[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
            u_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre_pre[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);

        }

        /// Define velocity gradient tensor
        double du[nsd][nsd];
        double dp[nsd];
        double ui_ujj[nsd];
        double uj_uij[nsd];
        for (int i = 0; i < nsd; i++) {
            dp[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
            for (int j = 0; j < nsd; j++) {
                du[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
            }
        }


        double d2u_ijj[nsd];
        double d2u_jij[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj[dof] = 0;
            d2u_jij[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dof, dir, dir);
                d2u_jij[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + dir, dof, dir);

            }
        }






        double res_M[nsd];
        double res_C = 0;
        double ukk;
        for (int i = 0; i < nsd; i++) {
            ui_ujj[i] = 0;
            uj_uij[i] = 0;
            ukk = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj[i] += u[i] * du[j][j];
                uj_uij[i] += u[j] * du[i][j];
                ukk += du[j][j];
            }
            res_M[i] = temp_coeff * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt + adv_coeff * (ui_ujj[i] + uj_uij[i]) +
                       press_coeff * dp[i] - diff_coeff * (d2u_ijj[i] + d2u_jij[i]) - force_coeff * force[i];
            if (idata_->SecondViscosity){
                res_M[i] -= diff_coeff * lambda * d2u_jij[i];
            }
            res_C += du[i][i];
        }


        ////////////////////////////////// end of precalculations



        double R_NS[nsd];
        double R_C;
        for (int a = 0; a < nbf; a++) {
            R_C = 0;
            for (int i = 0; i < nsd; i++) {
                R_NS[i] = 0;
                for (int j = 0; j < nsd; j++) {
                    /// expression:  R_NS[i] += advection[i][j] + diffusion[i][j]
                    R_NS[i] += adv_coeff * ((-fe.dN(a, j) * u[i] * u[j]) + (fe.dN(a, j) * u[i] * tauM * res_M[j]) +
                                            (fe.dN(a, j) * tauM * res_M[i] * u[j]) - (fe.dN(a, j) * tauM * res_M[i] * tauM * res_M[j])) +
                               diff_coeff * (fe.dN(a, j) * (du[i][j] + du[j][i]));
                    if (idata_->DiffFineTerms){
                        R_NS[i] += diff_coeff *  (fe.d2N(a, j, i) * tauM * res_M[j]);
                    }

                    if (idata_->SecondViscosity) {
                        R_NS[i] += diff_coeff * lambda * (fe.dN(a, i) * du[j][j] + fe.d2N(a, i, j) * tauM * res_M[j]);
                    }

                } // end of j loop
                /// expression:  R_NS[i] += temporal[i] - pressure[i] - force[i]
                R_NS[i] += temp_coeff * (fe.N(a) * (bdf[0] * u[i] + bdf[1] * u_pre[i] + bdf[2] * u_pre_pre[i]) / dt) +
                           press_coeff * ((-fe.dN(a, i) * p) + (fe.dN(a, i) * tauC * res_C)) -
                           force_coeff * fe.N(a) * force[i];
                R_C += (fe.N(a) * du[i][i]) + (fe.dN(a, i) * tauM * res_M[i]);

                be((nsd + 1) * a + i) +=  R_NS[i] * detJxW;

            } // end of i loop
            be((nsd + 1) * a + nsd) += R_C * detJxW;
        } // end of a loop

    }



    VMSParams calc_tau(const TALYFEMLIB::FEMElm &fe) const {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;
        const double Ci_f = 36;
        const int nsd = DIM;

        VMSParams params;

        ZEROPTV u;
        for (int i = 0; i < nsd; i++) {
            u(i) = this->p_data_->valueFEM(fe, i);
        }

        ZeroMatrix<double> ksiX;
        ksiX.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                ksiX(i, j) = fe.cof(j, i) / fe.jacc();
            }
        }

        ZeroMatrix<double> Ge;
        Ge.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                Ge(i, j) = 0.0;
                for (int k = 0; k < nsd; k++)
                    Ge(i, j) += ksiX(k, i) * ksiX(k, j);
            }
        }

        double u_Gu = 0.0;
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                u_Gu += u(i) * Ge(i, j) * u(j);
            }
        }

        double G_G_u = 0.0;
        for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
                G_G_u += Ci_f * Coe_diff * Coe_diff * Ge(i, j) * Ge(i, j);
            }
        }

        params.tauM = 1.0 / sqrt(4.0 / (this->dt_ * this->dt_) + u_Gu + G_G_u);

        ZEROARRAY<double> ge;
        ge.redim(nsd);
        for (int i = 0; i < nsd; i++) {
            ge(i) = 0.0;
            for (int j = 0; j < nsd; j++)
                ge(i) += ksiX(j, i);
        }

        double g_g = 0.0;
        for (int i = 0; i < nsd; i++) {
            g_g += ge(i) * ge(i);
        }

        params.tauC = 1.0 / (params.tauM * g_g);

        return params;
    }

    void calcAe_vms_bdf2(const FEMElm &fe, ZeroMatrix<double> &Ae, std::vector<double> bdf2)
    {

        VMSParams vms_params = calc_tau(fe);

        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        double Re = idata_->Re;
        const double Coe_diff = 1/Re;
        const double Ci_f = 36;

        const int nsd = DIM;
        const int n_basis_functions = fe.nbf();
        const double detJxW = fe.detJxW();
        /// field variables
        /** Get u and u_pre from the NodeData
         * NOTE: ValueFEM is defined on NodeData object, therefore it follows
         * indices defined in NodeData subclass.
         * NOTE: Every iteration the solution is copied into NodeData, we use these
         * solved fields from the NodeData object to calculate the be (NS residual)
         * of current
         * iteration
         */
        ZEROPTV u, u_pre, u_pre_pre;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = this->p_data_->valueFEM(fe, i);
            u_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
        }

        /// Define velocity gradient tensor
        ZeroMatrix<double> du;
        du.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
            }
        }

        /** Calculate the laplacian of velocity
         * This is required for course scale residual of Navier-Stokes for
         * diffusion term
         */
        ZEROPTV d2u;
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++)
        {
            /// Summing over three directions of velocity as it is laplacian
            for (int dir = 0; dir < nsd; dir++)
            {
                d2u(dof) += this->p_data_->value2DerivativeFEM(fe, dof, dir, dir);
            }
        }

        /** Remember your defined on the ValueFEM is on NodeData object, so use
         * indices the way you defined in NodeData, we acquire pressure from the
         * NodeData.
         */
        double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
        // double theta = p_data_->valueFEM(fe, nsd + 1);

        /// Get gradient of pressure
        ZEROPTV dp;
        for (int i = 0; i < nsd; i++)
        {
            dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
        }

        // NS terms
        ZEROPTV NS_time_derv;
        for (int i = 0; i < nsd; i++)
        {
            NS_time_derv(i) = (bdf2[0] * u(i) + bdf2[1] * u_pre(i) + bdf2[2] * u_pre_pre(i)) / dt_;
        }
        /** We define the convection of Navier Stokes here from
         * using inner product of gradient tensor and fields we acquired above.
         **/
        ZEROPTV convec;
        ZEROPTV diffusion;
        for (int i = 0; i < nsd; i++)
        {
            convec(i) = 0.0;
            diffusion(i) = 0.0;
            for (int j = 0; j < nsd; j++)
            {
                convec(i) += du(i, j) * u(j);
            }
            diffusion(i) += Coe_diff * d2u(i);
        }

        /** Construct the Navier Stokes equation without diffusion
         * Here diffusion is not present as its contribution to stabilizers for
         * linear basis functions is zero
         * Diffusion term is added for higher order basis functions
         * Equation 61 of Bazilevs et al. (2007)
         */
        ZEROPTV NS, NS_pre;
        for (int i = 0; i < nsd; i++)
        {
            NS(i) = NS_time_derv(i) + convec(i) + dp(i) - diffusion(i);

        }

        /** Calculate continuity residual for PSPG stabilizer
         * This residual is a essential for calculating the PSPG stabilizer
         * Equation 62 of Bazilevs et al. (2007)
         */
        double cont = 0.0;
        for (int i = 0; i < nsd; i++)
        {
            cont += du(i, i);
        }

        /** Calculating the Elemental matrices requires loop over basis functions
         * We loop over basis functions and calculate the Jacobian and the Residual
         * The equation we are solving is \vect{J} \cdot \delta u = -E
         * Here \vect{J} is the Jacobian matrix operating on \delta u and E is the
         * residual.  We will be using stabilized forms of Jacobian matrix and
         * residual
         */

        for (int a = 0; a < n_basis_functions; a++)
        {
            ///---------------Calculating the Jacobian operator--------------------///
            // NS terms

            /** Course scale terms for the cross terms (terms 4 and 5 in equation 52
             * in Bazilevs et al. (2007)
             * PS: these are not calculated using the SUPG class, and this is an
             * approximation of (u, grad{u})
             */
            double crossTermVelocityPart = 0.0;
            for (int i = 0; i < nsd; i++)
            {
                crossTermVelocityPart += fe.dN(a, i) * u(i);
            }

            /// Actual fine scale contribution to cross terms from tauM(inverse
            /// estimate)
            double crossTermFineScalePart = 0.0;
            for (int i = 0; i < nsd; i++)
            {
                crossTermFineScalePart += fe.dN(a, i) * tauM * NS(i);
            }

            for (int b = 0; b < n_basis_functions; b++)
            {

                /// Convection term
                double conv = 0.0;
                for (int i = 0; i < nsd; i++)
                {
                    conv += fe.dN(b, i) * u(i);
                }

                /// Adding terms to the Jacobian matrix.
                for (int i = 0; i < nsd; i++)
                {

                    /// Transient terms and the convection term and part of the stress
                    /// tensor in the diagonals
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                            (fe.N(a) * (fe.N(b) * bdf2[0] / dt_ + conv)) * detJxW;

                    for (int j = 0; j < nsd; j++)
                    {
                        /// This term calculates (w, (delta_u . grad{u_n}))
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (fe.N(a) * du(i, j) * fe.N(b) * detJxW);
                        /** This term calculates (grad{w},grad{delta_u}), goes only in
                         * diagonals PS: In this case we are using the diffusion form not
                         * the stress tensor form of the momentun equation
                         */
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                                (Coe_diff * fe.dN(a, j) * fe.dN(b, j) * detJxW);

                        //            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += thetaTimeStepping *
                        //                (Coe_diff * fe.dN(a, j) * fe.dN(b, i) * detJxW);
                    }
                    /** This term calculates contribution of pressure to Jacobian matrix
                     * mathematically the term is:  -(div{w},delta_p). This term gets
                     * added to the last diagonal term of the matrix
                     */
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                            -fe.dN(a, i) * fe.N(b) * detJxW;
                }

                /** crossTerm 1: This term is written as, (w, grad{u u'}), which is
                 * weakened to (grad{w}, u\innerproduct u') which is then linearised.
                 * Here u' is fine scale velocity
                 * and u is resolved velocity. The fine scale velocity is approximated
                 * as -tau*Res{NS} (Equation 58 bazilevs et al. (2007)),
                 *
                 */
                for (int i = 0; i < nsd; i++)
                {

                    /// Contribution of laplacian of velocity(diffusion) to the diagonal
                    double diff_J = 0;
                    for (int j = 0; j < nsd; j++)
                    {
                        diff_J += Coe_diff * fe.d2N(b, j, j);
                    }
                    /** When you linearise (u.grad{w},u'), (B_2 from Equation 57 from
                     * Bazilevs et al. (2007)) you get two terms, one goes
                     * in diagonal and other in the non-diagonal parts of the matrix
                     * PS: crossTermVelocityPart is (u. grad{w} => u*dN(a))
                     */

                    /// Diagonal part
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                            (crossTermVelocityPart * tauM * (fe.N(b) * bdf2[0] / dt_ + conv - diff_J) * detJxW);
                    for (int j = 0; j < nsd; j++)
                    {
                        /// Off diagonal part
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (crossTermVelocityPart * tauM * du(i, j) * fe.N(b) * detJxW);
                        /** this term is essentially, (grad{w}*tauM*Residual(NS), delta u)
                         * Equation 57 Bazilevs et al. (2007)
                         * PS: u' = (tauM*residual(NS))
                         */
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (tauM * fe.dN(a, j) * NS(i) * fe.N(b) * detJxW);
                    }
                    /// Pressure part goes in the last element of the diagonal
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                            (crossTermVelocityPart * tauM * fe.dN(b, i) * detJxW);
                }

                /** Crossterm2:This term can be mathematically written as,
                 * (w, u'.grad{u}). In this term we do not further weaken it as done in
                 * cross term 1
                 */
                for (int i = 0; i < nsd; i++)
                {

                    /// Contribution of laplacian of velocity(diffusion) to the diagonal
                    double diff_J = 0;
                    for (int j = 0; j < nsd; j++)
                    {
                        diff_J += Coe_diff * fe.d2N(b, j, j);
                    }

                    /** This term represents the contribution of,
                     * (w, tau \delta{u} . grad{u}) term to the cross term, here \delta{} is
                     * the linearisation operator (basically delta{u'} represents the
                     * change in u').  This is the first term which arise due to
                     * linearisation of the nonliner term of the
                     * Navier-Stokes residual when u' is substitute for in
                     * (w,\delta{u'}.grad{u}). PS: u' = - tauM*res{NS}
                     */
                    for (int j = 0; j < nsd; j++)
                    {
                        /// k is the dummy index
                        for (int k = 0; k < nsd; k++)
                        {
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                    (-du(i, k) * fe.N(a) * tauM * fe.N(b) * du(k, j) * detJxW);
                        }
                    }

                    for (int j = 0; j < nsd; j++)
                    {
                        /** This term represents the contribution of,
                         * (w, \delta{u} . grad{u\delta{u'}})
                         */
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (-du(i, j) * fe.N(a) * tauM * (fe.N(b) * bdf2[0] / dt_ + conv - diff_J) * detJxW);
                        /** This term is the contribution of (w, u'.grad{delta{u}}) to the
                         * Jacobian operator.
                         */
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                                (-fe.N(a) * tauM * NS(j) * fe.dN(b, j) * detJxW);

                        /// Pressure part of the term which arise from (w,\delta{u'}
                        /// .grad{u}), always the last diagonal term
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                                (-du(i, j) * fe.N(a) * tauM * fe.dN(b, j) * detJxW);
                    }
                }

                /** The Reynolds Stress term: (w, u'.grad{u'}), we subsitute u' as
                 * -tau*Res(NS) and expand.
                 */
                for (int i = 0; i < nsd; i++)
                {
                    /// Contribution of laplacian of velocity(diffusion) to the diagonal
                    double diff_J = 0;
                    for (int j = 0; j < nsd; j++)
                    {
                        diff_J += Coe_diff * fe.d2N(b, j, j);
                    }

                    /** Three terms arising from (w,\delta(u').grad{u'})
                     * u' has to be expanded to -tau*res{NS}: when the linearisation acts
                     * on the res{NS} it gives three terms. First term which goes in the
                     * diagonals of the matrix is
                     * (-grad{w}* tauM*res(NS), tauM*(d_t{\delta{u}} + conv -diff_J))
                     * PS: grad{w}*tau*res(NS) is corssTermFineScalePart
                     */
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) +=
                            (-crossTermFineScalePart * tauM * (fe.N(b) * bdf2[0] / dt_ + conv - diff_J) * detJxW);

                    /** Second term from (w,\delta(u').grad{u'}) which goes to the off
                     * diagonals is: (-grad{w}* tauM*res(NS), tauM*(\delta{u}. grad{u}))
                     */
                    for (int j = 0; j < nsd; j++)
                    {
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (-crossTermFineScalePart * tauM * du(i, j) * fe.N(b) * detJxW);
                    }

                    /** Third term from (w,\delta(u').grad{u'}) which is the contribution
                     * of pressure and goes in the last diagonal term
                     */
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                            (-crossTermFineScalePart * tauM * fe.dN(b, i) * detJxW);

                    for (int j = 0; j < nsd; j++)
                    {
                        for (int k = 0; k < nsd; k++)
                        {
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                    (-tauM * NS(i) * fe.dN(a, k) * tauM * fe.N(b) * du(k, j) * detJxW);
                        }
                    }

                    /** Just as above terms which arise when(w,(u').grad{\delta{u'}}) is
                     * expanded is given below.
                     */
                    for (int j = 0; j < nsd; j++)
                    {
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (-tauM * NS(i) * fe.dN(a, j) * tauM * (fe.N(b) * bdf2[0] / dt_ + conv - diff_J) * detJxW);
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                                (-tauM * NS(i) * fe.dN(a, j) * tauM * fe.dN(b, j) * detJxW);
                    }
                }

                /** This term represents the fine scale pressure given by
                 * (grad{w}, \delta{tauC*res{Cont}}), where Cont is the continuity
                 * equation. See third term in equation (71) in Bazilev et al. (2007)
                 */
                for (int i = 0; i < nsd; i++)
                {
                    for (int j = 0; j < nsd; j++)
                    {
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                (fe.dN(a, i) * tauC * fe.dN(b, j) * detJxW);
                    }
                }
                /** pspg term from VMS theory: See second term in equation (71)
                 * from Bazilevs et al. (2007)
                 */
                for (int i = 0; i < nsd; i++)
                {
                    double diff_J = 0;
                    for (int j = 0; j < nsd; j++)
                    {
                        diff_J += Coe_diff * fe.d2N(b, j, j);
                    }
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
                            fe.N(a) * fe.dN(b, i) * detJxW +
                            fe.dN(a, i) * tauM * (fe.N(b) * bdf2[0] / dt_ + conv - diff_J) * detJxW;

                    for (int j = 0; j < nsd; j++)
                    {
                        Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) +=
                                fe.dN(a, j) * tauM * du(j, i) * fe.N(b) * detJxW;
                    }

                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) +=
                            fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW;
                }
            }
        }
    }

    void calcbe_vms_bdf2(const FEMElm &fe, ZEROARRAY<double> &be, std::vector<double> bdf2)
    {

        VMSParams vms_params = calc_tau(fe);

        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        double Re = idata_->Re;
        const double Coe_diff = 1/Re;
        const double Ci_f = 36;
        ;
        const int nsd = DIM;
        const int nbf = fe.nbf();
        const double detJxW = fe.detJxW();

        ZEROPTV u, u_pre, u_pre_pre;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X + i);
            u_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
        }


        ZeroMatrix<double> du;
        du.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                du(i, j) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X + i, j);
            }
        }

        ZEROPTV d2u(0, 0, 0);
        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                d2u(i) += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X + i, j, j);
            }
        }

        double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);
        ZEROPTV dp;
        for (int i = 0; i < nsd; i++)
        {
            dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
        }

        ZEROPTV convec;
        for (int i = 0; i < nsd; i++)
        {
            convec(i) = 0.0;
            for (int j = 0; j < nsd; j++)
            {
                convec(i) += du(i, j) * u(j);
            }
        }

        ZEROPTV diffusion;
        for (int i = 0; i < nsd; i++)
        {
            diffusion(i) = Coe_diff * d2u(i);
        }

        ZEROPTV NS_time_derv;
        for (int i = 0; i < nsd; i++)
        {
            NS_time_derv(i) = (bdf2[0] * u(i) + bdf2[1] * u_pre(i) + bdf2[2] * u_pre_pre(i)) / dt_;
        }

        ZEROPTV NS;
        for (int i = 0; i < nsd; i++)
        {
            NS(i) = NS_time_derv(i) + convec(i) + dp(i) - diffusion(i);
        }

        double cont = 0.0;
        for (int i = 0; i < nsd; i++)
        {
            cont += du(i, i);
        }


        for (int a = 0; a < nbf; a++)
        {
            double crossTermVelocityPart = 0.0;
            for (int i = 0; i < nsd; i++)
            {
                crossTermVelocityPart += fe.dN(a, i) * u(i);
            }

            /** Constructing the RHS vector of the residual of the equation
             * This involves just defining the terms involved in the equation without
             * any linearisation, as this is calculated from the existing data on the
             * NodeData structure
             */

            /// All the terms involved in the equation
            ZEROPTV time_der;
            ZEROPTV cross1, cross2, reystress, pressurep, diff, press, normal_NS, NScontrib;

            /// Diffusion term
            for (int i = 0; i < nsd; i++)
            {
                diff(i) = 0.0;
                for (int j = 0; j < nsd; j++)
                {
                    diff(i) += Coe_diff * fe.dN(a, j) * (du(i, j) /* + du(j, i)*/);
                }
            }

            /// Terms of Normal Navier Stokes without any additional VMS terms
            for (int i = 0; i < nsd; i++)
            {
                time_der(i) = fe.N(a) * NS_time_derv(i) * detJxW;
                press(i) = (-fe.dN(a, i) * p) * detJxW;
                normal_NS(i) = time_der(i) + press(i) + fe.N(a) * (convec(i)) * detJxW + (diff(i)) * detJxW;
                /// first cross term
                cross1(i) = (tauM * crossTermVelocityPart * NS(i)) * detJxW;

                cross2(i) = 0.0;
                reystress(i) = 0.0;

                for (int j = 0; j < nsd; j++)
                {
                    /// Second cross term
                    cross2(i) += fe.N(a) * (du(i, j) * tauM * NS(j)) * detJxW;

                    /// The Reynolds stress term
                    reystress(i) += fe.dN(a, j) * (tauM * NS(i) * tauM * NS(j)) * detJxW;
                }

                /// Contribution of pressure
                pressurep(i) = fe.dN(a, i) * (tauC * cont) * detJxW;

                /// Add all the terms to the vector
                NScontrib(i) = normal_NS(i) + cross1(i) - cross2(i) - reystress(i) + pressurep(i);
            }

            double pspg = 0.0;

            for (int i = 0; i < nsd; i++)
            {
                /// Fine scale contribution to pressure
                pspg += fe.dN(a, i) * tauM * NS(i);
            }
            /// Contribution of this term goes to last element of the vector
            double contContrib = fe.N(a) * (cont)*detJxW + pspg * detJxW;

            /// Adding all the calculated terms to the vector
            for (int i = 0; i < nsd; i++)
            {
                be((nsd + 1) * a + i) += NScontrib(i);
            }
            be((nsd + 1) * a + nsd) += contContrib;
        }
    }



    void ReRamping (double t) {
        idata_->Re = idata_->getReynolds(t);
    }


    void Integrands_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZeroMatrix<double> &Ae, double *h) {

        if (idata_->DoReRamping) {
            ReRamping(t_);
        }
        std::vector<double> bdf;
        if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
            t_ < 1.5 * dt_) {
            bdf = {1.0, -1.0, 0.0};
        } else if (idata_->timeStepping == NSInputData::BDF2) {
            bdf = {1.5, -2.0, 0.5};
        }

        if (idata_->OldNSFormulation){
            calcAe_vms_bdf2(fe, Ae, bdf);
        } else {
            calcAe_vms(fe, Ae, bdf);
        }
    }

    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
        if (idata_->DoReRamping) {

            ReRamping(t_);
        }

//        std::cout << "1/Re = " << 1/idata_->Re << "\n";
        std::vector<double> bdf;
        if (idata_->timeStepping == NSInputData::BACKWARD_EULER ||
            t_ < 1.5 * dt_) {
            bdf = {1.0, -1.0, 0.0};
        } else if (idata_->timeStepping == NSInputData::BDF2) {
            bdf = {1.5, -2.0, 0.5};
        }

        if (idata_->OldNSFormulation){
            calcbe_vms_bdf2(fe, be, bdf);
        } else {
            calcbe_vms(fe, be, bdf);
        }
    }



#ifdef IBM
    void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h)
    {
//        std::ofstream fout("IBM_Ae.txt",std::ios::app);
//        fout << fe.position().x() << " " << fe.position().y() << "\n";
//


        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        auto &geo = idata_->carved_out_geoms_def.at(gpinfo.geomID);
        TALYFEMLIB::ZEROPTV normal = TALYFEMLIB::ZEROPTV{gpinfo.normal[0], gpinfo.normal[1], 0};
#if (DIM == 3)
        normal(2) = gpinfo.normal[2];
#endif

        const int nsd = DIM;
        const int nbf = fe.nbf();

        const double detSideJxW = fe.detJxW();
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;

        /** Get u from the NodeData
         * NOTE: ValueFEM is defined on NodeData object, therefore it follows
         * indices defined in NodeData subclass.
         * NOTE: Every iteration the solution is copied into NodeData, we use these
         * solved fields from the NodeData object to calculate the be (NS residual)
         * of current iteration.
         */
        TALYFEMLIB::ZEROPTV u;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = this->p_data_->valueFEM(fe, i);
        }

        /// Define velocity gradient tensor ( grad{u} )
        TALYFEMLIB::ZeroMatrix<double> du;
        du.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
            }
        }

        /** Calculate the laplacian of velocity
         * This is required for course scale residual of Navier-Stokes for
         * diffusion term*/

        TALYFEMLIB::ZEROPTV d2u;
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++)
        {
            /// Summing over three directions of velocity as it is laplacian
            for (int dir = 0; dir < nsd; dir++)
            {
                d2u(dof) += this->p_data_->value2DerivativeFEM(fe, dof, dir, dir);
            }
        }

        TALYFEMLIB::ZEROPTV dp;
        for (int i = 0; i < nsd; i++)
        {
            dp(i) = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE, i);
        }

        /// wall normal distance
        double hb = normalDistance(fe, normal, h);

//        std::cout << "hb = " << hb << "\n";

        double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, orderOfhb);
//        double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);

        for (int a = 0; a < nbf; a++)
        {
            double supg = 0.0;
            for (int m = 0; m < nsd; m++)
            {
                supg += fe.N(a) * normal(m);
            }

            TALYFEMLIB::ZEROPTV supg1, supg2;
            for (int m = 0; m < nsd; m++)
            {
                supg1(m) = 0.0;
                supg2(m) = 0.0;
                for (int n = 0; n < nsd; n++)
                {
                    supg1(m) += fe.dN(a, n) * normal(n);
                    supg2(m) += fe.dN(a, m) * normal(n);
                }
            }

            for (int b = 0; b < nbf; b++)
            {
                double supgb = 0.0;
                for (int m = 0; m < nsd; m++)
                {
                    supgb += fe.dN(b, m) * normal(m);
                }
                for (int j = 0; j < nsd; j++)
                {
                    Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                            /// diffusion_1
                            -Coe_diff * fe.N(a) * (supgb)*detSideJxW
                            /// adjoint diffusion
                            - Coe_diff * (supg1(j) + supg2(j)) * fe.N(b) * detSideJxW
                            /// penalty
                            + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW;
                    //                    + Coe_weak * fe.N(a) * fe.N(b) * detSideJxW;

                    /// adjoint pressure
                    Ae((nsd + 1) * (a + 1) - 1, (nsd + 1) * b + j) +=
                            -fe.N(a) * normal(j) * fe.N(b) * detSideJxW;
                    /// pressure
                    Ae((nsd + 1) * a + j, (nsd + 1) * (b + 1) - 1) +=
                            fe.N(a) * fe.N(b) * normal(j) * detSideJxW;

                    if(idata_->BoundaryTerms) {
                        for (int i = 0; i < nsd; i++) {

                            /// advection term

                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                    fe.N(a) * (u[i] * fe.N(b)) * normal(i) * detSideJxW;

                            Ae((nsd + 1) * a + j, (nsd + 1) * b + i) +=
                                    fe.N(a) * (u[j] * fe.N(b)) * normal(i) * detSideJxW;

                            /// coarseDiffusion 1
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                    Coe_diff * fe.N(a) * fe.dN(b, i)
                                    * normal(i) * detSideJxW;

                            /// coarseDiffusion 2
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + i) +=
                                    Coe_diff * fe.N(a) * fe.dN(b, j)
                                    * normal(i) * detSideJxW;

                            if (idata_->SecondViscosity) {
                                /// coarseDiffusion 3
                                Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                        Coe_diff * lambda * fe.N(a) * fe.dN(b, j)
                                        * normal(i) * detSideJxW;
                            }
                        }
                    }

                }
            }
        }

    }

    void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h)
    {
//        std::ofstream fout("IBM_be.txt",std::ios::app);
//        fout << fe.position().x() << " " << fe.position().y() << "\n";

        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        auto &geo = idata_->carved_out_geoms_def.at(gpinfo.geomID);

        ZEROPTV normal = ZEROPTV{gpinfo.normal[0], gpinfo.normal[1], 0};
#if (DIM==3)
        normal(2) = gpinfo.normal[2];
#endif

//        std::cout << "normal = "; normal.print(); std::cout << "\n";
#if (DIM == 2)
        std::vector<double> bc_value = {geo.getBC(NSNodeData::VEL_X)[0],
                                        geo.getBC(NSNodeData::VEL_Y)[0]};

#elif (DIM == 3)
        std::vector<double> bc_value = {geo.getBC(NSNodeData::VEL_X)[0],
                                            geo.getBC(NSNodeData::VEL_Y)[0],
                                            geo.getBC(NSNodeData::VEL_Z)[0]};
#endif

        /// NOTE: previous bug is here!
        const int nsd = DIM;

        const int nbf = fe.nbf();
        const double detSideJxW = fe.detJxW();
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;


//        std::cout << "detSideJxW = " << detSideJxW << "\n";

        TALYFEMLIB::ZEROPTV u, u_pre;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = this->p_data_->valueFEM(fe, i);
            u_pre(i) = this->p_data_->valueFEM(fe, NSNodeData::NUM_VARS + i);
        }

//        std::cout << "u ="; u.print(); std::cout << "\n";
//        std::cout << "u_pre ="; u_pre.print(); std::cout << "\n";

        /// Define velocity gradient tensor ( grad{u} )
        TALYFEMLIB::ZeroMatrix<double> du;
        du.redim(nsd, nsd);
        for (int i = 0; i < nsd; i++)
        {
            for (int j = 0; j < nsd; j++)
            {
                du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
            }
        }

        double p = this->p_data_->valueFEM(fe, NSNodeData::PRESSURE);

        TALYFEMLIB::ZEROPTV diff1, diff2;
        for (int i = 0; i < nsd; i++)
        {
            diff1(i) = 0.0;
            diff2(i) = 0.0;
            for (int j = 0; j < nsd; j++)
            {
                diff1(i) += du(i, j) * normal(j);
                diff2(i) += du(j, i) * normal(j);
            }
        }

        /// wall normal distance (the maximum of 8 nodes)
        /// wall normal distance
        double hb = normalDistance(fe, normal, h);


        double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, orderOfhb);

//        double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);

        for (int a = 0; a < nbf; a++)
        {
            TALYFEMLIB::ZEROPTV supg1;
            for (int m = 0; m < nsd; m++)
            {
                supg1(m) = 0.0;
                for (int n = 0; n < nsd; n++)
                {
                    supg1(m) += fe.dN(a, n) * normal(n);
                }
            }

            for (int j = 0; j < nsd; j++)
            {

                be((nsd + 1) * a + j) +=
                        /// pressure
                        fe.N(a) * (p * normal(j)) * detSideJxW
                        /// diffusion
                        - fe.N(a) * (Coe_diff * (diff1(j) /*+diff2(j)*/)) * detSideJxW
                        /// adjoint diffusion
                        - Coe_diff * supg1(j) * (u(j) - bc_value.at(j)) * detSideJxW
                        /// penalty
                        + weakBCpenaltyParameter * fe.N(a) * (u(j) - bc_value.at(j)) * detSideJxW;
                /// adjoint pressure
                be((nsd + 1) * (a + 1) - 1) += -fe.N(a) * normal(j) * (u(j) - bc_value.at(j)) * detSideJxW;

                if(idata_->BoundaryTerms) {
                    for (int i = 0; i < nsd; i++) {
                        /// advention term
                        be((nsd + 1) * a + j) += (fe.N(a) * u[j] * u[i])
                                                 * normal(i) * detSideJxW;
                        /// Diffusion term 1, 2
                        be((nsd + 1) * a + j) += Coe_diff * fe.N(a) * (du(j, i) + du(i, j))
                                                 * normal(i) * detSideJxW;

                        if (idata_->SecondViscosity) {
                            /// Diffusion term 3
                            be((nsd + 1) * a + i) += Coe_diff * lambda * fe.N(a) * du(j, j)
                                                     * normal(i) * detSideJxW;
                        }
                    }
                }

//                std::cout << "j = " << j << std::endl;
//                std::cout <<"fe.N(a) = " << fe.N(a) << std::endl;
//                std::cout << "p * normal(j) = " << p * normal(j) << std::endl;
//                std::cout << "detSideJxW = " << detSideJxW << std::endl;
//                std::cout << fe.N(a) * (p * normal(j)) * detSideJxW << std::endl;
//                std::cout << fe.N(a) * (Coe_diff * (diff1(j) /*+diff2(j)*/)) * detSideJxW << std::endl;
//                std::cout << -Coe_diff * supg1(j) * (u(j) - bc_value.at(j)) * detSideJxW << std::endl;
            }
        }

    }
#endif

#pragma mark Ae-be Surface

    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae)
    {
        if (GpOnDomainWall(fe))
        {
            /// side_idx in the order of x-, x+, y-, y+, z-, z+
//            std::cout << "side_idx = " << side_idx << "\n";
            const auto &bc = idata_->boundary_def.at(side_idx);
//            std::cout << "bc.pressure_type = " << bc.pressure_type << "\n";
            Integrands4side_Ae_domain(fe, Ae,bc);
        }
        else
        {
            // this is for carved out boundaries
            id = NearestGeo(fe);

            const auto &carved_out_def = idata_->carved_out_geoms_def.at(id);
            Integrands4side_Ae_carved(fe, Ae,carved_out_def, id);
        }
    }

    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae,
                            const double *h)
    {
        Integrands4side_Ae(fe, side_idx, id, Ae);
    }

    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be)
    {
        if (GpOnDomainWall(fe))
        {
            /// side_idx in the order of x-, x+, y-, y+, z-, z+
            const auto &bc = idata_->boundary_def.at(side_idx);
            Integrands4side_be_domain(fe, be,bc);
        }
        else
        {
            // this is for carved out boundaries
            id = NearestGeo(fe);
            const auto &carved_out_def = idata_->carved_out_geoms_def.at(id);
            Integrands4side_be_carved(fe, be,carved_out_def, id);
        }
    }

    void Integrands4side_be(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZEROARRAY<double> &be,
                            const double *h)
    {
        Integrands4side_be(fe, side_idx, id, be);
    }

#pragma mark Ae-be Surface Domain

    void Integrands4side_Ae_domain(const TALYFEMLIB::FEMElm &fe,
                                   TALYFEMLIB::ZeroMatrix<double> &Ae,
                                   const BoundaryDef &bc)
    {
        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        if (bc.pressure_type == BoundaryDef::Pressure_Type::BACKFLOW_STAB)
        {
//            std::cout << "Ae backflow\n";
//            std::cout << "Ae, position = "; fe.position().print();std::cout << "\n";

            const int nbf = fe.nbf();
            const int nsd = DIM;
            const double detSideJxW = fe.detJxW();

            TALYFEMLIB::ZEROPTV u;
            for (int i = 0; i < nsd; i++)
            {
                u(i) = p_data_->valueFEM(fe, i);
            }

            double OutFlow = 0.0;
            for (int m = 0; m < nsd; m++)
            {
                OutFlow += u(m) * fe.surface()->normal()(m);
            }

            double Uneg = std::min(0.0, OutFlow);

            for (int a = 0; a < nbf; a++)
            {
                for (int b = 0; b < nbf; b++)
                {
                    for (int j = 0; j < nsd; j++)
                    {
                        Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                /// BackFlow Stabilization
                                -BetaNS * fe.N(a) * Uneg * fe.N(b) * detSideJxW;
                    }
                }
            }
        }

    }

    void Integrands4side_be_domain(const TALYFEMLIB::FEMElm &fe,
                                   TALYFEMLIB::ZEROARRAY<double> &be,
                                   const BoundaryDef &bc)
    {
        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        if (bc.pressure_type == BoundaryDef::Pressure_Type::BACKFLOW_STAB)
        {
//            std::cout << "be backflow\n";

//            std::cout << "be, position = "; fe.position().print();std::cout << "\n";
            const int nbf = fe.nbf();
            const int nsd = DIM;
            const double detSideJxW = fe.detJxW();

            TALYFEMLIB::ZEROPTV u;
            for (int i = 0; i < nsd; i++)
            {
                u(i) = p_data_->valueFEM(fe, i);
            }

            double OutFlow = 0.0;
            for (int m = 0; m < nsd; m++)
            {
                OutFlow += u(m) * fe.surface()->normal()(m);
            }

            double Uneg = std::min(0.0, OutFlow);

            for (int a = 0; a < nbf; a++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    be((nsd + 1) * a + j) +=
                            /// BackFlow Stabilization
                            -BetaNS * fe.N(a) * Uneg * u(j) * detSideJxW;
                }
            }
        }
    }

#pragma mark Ae-be Surface CarveOut

    void Integrands4side_Ae_carved(const TALYFEMLIB::FEMElm &fe,
                                   TALYFEMLIB::ZeroMatrix<double> &Ae,
                                   const CarvedOutGeom &carved_out_def,
                                   const int CarvedOutGeomID)
    {

        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        int ndof = NSNodeData::NS_DOF;

        double Re = idata_->Re;

//        std::cout << "Re = " << Re << "\n";

        const double Coe_diff = 1.0 / Re;

//        std::cout << "Coe_diff (Ae) = " << Coe_diff << "\n";
        const int nbf = fe.nbf();
        const int nsd = DIM;
        const double detSideJxW = fe.detJxW();




        const auto &bcType = carved_out_def.bc_type_V[NSNodeData::VEL_X];

        if (bcType == CarvedOutGeom::BCType::WEAK)
        {

            TALYFEMLIB::ZEROPTV u;
            for (int i = 0; i < nsd; i++)
            {
                u(i) = p_data_->valueFEM(fe, i);
            }

            double p = p_data_->valueFEM(fe, NSNodeData::PRESSURE);

            TALYFEMLIB::ZeroMatrix<double> ksiX;
            ksiX.redim(nsd, nsd);

            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    ksiX(i, j) = fe.cof(j, i) / fe.volume_jacc();
                }
            }

            TALYFEMLIB::ZeroMatrix<double> Ge;
            Ge.redim(nsd, nsd);
            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    Ge(i, j) = 0.0;
                    for (int k = 0; k < nsd; k++)
                    {
                        Ge(i, j) += ksiX(k, i) * ksiX(k, j);
                    }
                }
            }

            double n_Gn = 0.0;

            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    n_Gn += fe.surface()->normal()(j) * Ge(j, i) * fe.surface()->normal()(i);
                }
            }

            double hb = 2.0 / sqrt(n_Gn);

            double inflow_g = 0.0;
            for (int m = 0; m < nsd; m++)
            {
                inflow_g += u(m) * fe.surface()->normal()(m);
            }

            double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, orderOfhb);

            for (int a = 0; a < nbf; a++)
            {
                double normalComponentOfN = 0.0;
                for (int m = 0; m < nsd; m++)
                {
                    normalComponentOfN += fe.N(a) * fe.surface()->normal()(m);
                }

                TALYFEMLIB::ZEROPTV normalComponentOfdNa;
                for (int m = 0; m < nsd; m++)
                {
                    normalComponentOfdNa(m) = 0.0;
                    for (int n = 0; n < nsd; n++)
                    {
                        normalComponentOfdNa(m) += fe.dN(a, n) * fe.surface()->normal()(n);
                    }
                }

                for (int b = 0; b < nbf; b++)
                {
                    double normalComponentOfdNb = 0.0;
                    for (int m = 0; m < nsd; m++)
                    {
                        normalComponentOfdNb += fe.dN(b, m) * fe.surface()->normal()(m);
                    }
                    for (int j = 0; j < nsd; j++)
                    {

                        for (int i = 0; i < nsd;i++) {
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                    fe.N(a) * (u[i] * fe.N(b)) * detSideJxW;

                            Ae((nsd + 1) * a + j, (nsd + 1) * b + i) +=
                                    fe.N(a) * (u[j] * fe.N(b)) * detSideJxW;
                        }

                        Ae((nsd + 1) * a + j, (nsd + 1) * b + j) +=
                                /// diffusion_1
                                -Coe_diff * fe.N(a) * (normalComponentOfdNb)*detSideJxW
                                /// adjoint diffusion
                                - Coe_diff * normalComponentOfdNa(j) * fe.N(b) * detSideJxW
                                /// penalty
                                + weakBCpenaltyParameter * fe.N(a) * fe.N(b) * detSideJxW;
                        if (inflow_g < 0)
                        {
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + j) -=
                                    u(j) * fe.surface()->normal()(j) * fe.N(a) * fe.N(b);
                        }
                        /// adjoint pressure
                        Ae((nsd + 1) * (a + 1) - 1, (nsd + 1) * b + j) +=
                                -fe.N(a) * fe.N(b) * fe.surface()->normal()(j) * detSideJxW;
                        /// pressure
                        Ae((nsd + 1) * a + j, (nsd + 1) * (b + 1) - 1) +=
                                fe.N(a) * fe.N(b) * fe.surface()->normal()(j) * detSideJxW;
                        /// diffusion_2 (new term)
                        for (int k = 0; k < nsd; k++)
                        {
                            Ae((nsd + 1) * a + j, (nsd + 1) * b + k) +=
                                    -Coe_diff * fe.N(a) * fe.dN(b, j) * fe.surface()->normal()(k) * detSideJxW;
                        }
                    }
                }
            }

        }
#ifdef IBM
        else if (bcType == CarvedOutGeom::BCType::SBM)
        {
            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);

//            std::cout << "d (Ae) = " << d[0] << " " << d[1] << std::endl;
            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_Ae(Ae, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d);
        }

        else if (bcType == CarvedOutGeom::BCType::PENALTY_SBM)
        {

            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);
            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);

            SBM_Ae(Ae, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d, true /*PenaltyOnly*/);
        }

#endif
    }

    void Integrands4side_be_carved(const TALYFEMLIB::FEMElm &fe,
                                   TALYFEMLIB::ZEROARRAY<double> &be,
                                   const CarvedOutGeom &carved_out_def,
                                   const int CarvedOutGeomID)
    {

        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        int ndof = NSNodeData::NS_DOF;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;

        const int nbf = fe.nbf();
        // previous bug is here
//        const int nsd = idata_->nsd;
        const int nsd = DIM;
        const double detSideJxW = fe.detJxW();

        const auto &bcType = carved_out_def.bc_type_V[NSNodeData::VEL_X];

        if (bcType == CarvedOutGeom::BCType::WEAK)
        {
#if (DIM == 2)
            std::vector<double> bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                                            carved_out_def.getBC(NSNodeData::VEL_Y)[0]};

#elif (DIM == 3)
            std::vector<double> bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                                            carved_out_def.getBC(NSNodeData::VEL_Y)[0],
                                            carved_out_def.getBC(NSNodeData::VEL_Z)[0]};
#endif

            TALYFEMLIB::ZEROPTV u;
            for (int i = 0; i < nsd; i++)
            {
                u(i) = p_data_->valueFEM(fe, i);
            }
            //            double uDotn = 0.0;
            //            for (int dim = 0;dim<DIM;dim++){
            //                uDotn += u(dim) * fe.surface()->normal()(dim);
            //            }

            TALYFEMLIB::ZeroMatrix<double> du;
            du.redim(nsd, nsd);

            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
                }
            }

            double p = p_data_->valueFEM(fe, NSNodeData::PRESSURE);

            TALYFEMLIB::ZeroMatrix<double> ksiX;
            ksiX.redim(nsd, nsd);

            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    ksiX(i, j) = fe.cof(j, i) / fe.volume_jacc();
                }
            }

            TALYFEMLIB::ZeroMatrix<double> Ge;
            Ge.redim(nsd, nsd);
            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    Ge(i, j) = 0.0;
                    for (int k = 0; k < nsd; k++)
                    {
                        Ge(i, j) += ksiX(k, i) * ksiX(k, j);
                    }
                }
            }

            double n_Gn = 0.0;

            for (int i = 0; i < nsd; i++)
            {
                for (int j = 0; j < nsd; j++)
                {
                    n_Gn += fe.surface()->normal()(j) * Ge(j, i) * fe.surface()->normal()(i);
                }
            }

            double hb = 2.0 / sqrt(n_Gn);
            TALYFEMLIB::ZEROPTV normalComponentOfVelocity;
            for (int i = 0; i < nsd; i++)
            {
                normalComponentOfVelocity(i) = 0.0;
                for (int j = 0; j < nsd; j++)
                {
                    normalComponentOfVelocity(i) += du(i, j) * fe.surface()->normal()(j);
                }
            }

            double inflow_g = 0.0;
            for (int m = 0; m < nsd; m++)
            {
                // this is because the surface normal is pointing inward of the geometry
                inflow_g += u(m) * fe.surface()->normal()(m);
            }

            for (int a = 0; a < nbf; a++)
            {
                double normalComponentOfN = 0.0;
                for (int m = 0; m < nsd; m++)
                {
                    normalComponentOfN += fe.N(a) * fe.surface()->normal()(m);
                }

                TALYFEMLIB::ZEROPTV normalComponentOfdNa;
                for (int m = 0; m < nsd; m++)
                {
                    normalComponentOfdNa(m) = 0.0;
                    for (int n = 0; n < nsd; n++)
                    {
                        normalComponentOfdNa(m) += fe.dN(a, n) * fe.surface()->normal()(n);
                    }
                }
                double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, orderOfhb);

                for (int j = 0; j < nsd; j++)
                {

                    for (int i = 0; i<nsd;i++){
                        be((nsd + 1) * a + j) += (fe.N(a) * u[i] * u[j]) *detSideJxW;
                    }

                    be((nsd + 1) * a + j) +=
                            /// pressure
                            fe.N(a) * (p * fe.surface()->normal()(j)
                                       /// diffusion
                                       - Coe_diff * (normalComponentOfVelocity(j))) *
                            detSideJxW
                            /// adjoint diffusion
                            - Coe_diff * normalComponentOfdNa(j) * (u(j) - bc_value.at(j)) * detSideJxW
                            /// penalty
                            + weakBCpenaltyParameter * fe.N(a) * (u(j) - bc_value.at(j)) * detSideJxW;
                    if (inflow_g < 0)
                    {
                        be((nsd + 1) * a + j) -= u(j) * fe.surface()->normal()(j) * fe.N(a) *
                                                 (u(j) - bc_value.at(j)) * detSideJxW;
                    }

                    /// adjoint pressure
                    be((nsd + 1) * (a + 1) - 1) +=
                            -fe.N(a) * fe.surface()->normal()(j) * (u(j) - bc_value.at(j)) * detSideJxW;
                }
            }

        }
#ifdef IBM
        if (bcType == CarvedOutGeom::BCType::SBM)
        {
            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);

//            std::cout << "d (be) = " << d[0] << " " << d[1] << std::endl;

            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_be(be, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d,carved_out_def);
        }
        if (bcType == CarvedOutGeom::BCType::PENALTY_SBM)
        {
            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);
            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_be(be, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d,carved_out_def, true /*PenaltyOnly*/);
        }

#endif
    }

// change this function for DeepTrace Code
//    void setSbmCalc(const IMGA *imga, const std::vector<my_kd_tree_ptr> &kd_trees)
// we don't need kdtree for any computation
    void setSbmCalc(const IMGA *imga)
    {
        imga_ = imga;
//        std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
//                       [](const auto& uniquePtr) { return uniquePtr.get(); });
    }

    void setDistanceArrays(BoundaryDistanceData &boundaryDistanceData)
    {
        distanceToBdry = boundaryDistanceData.distanceToBdry;
        bdryNormal = boundaryDistanceData.bdryNormal;
        surfaceNodeCoords = boundaryDistanceData.surfaceNodeCoords;
        distanceSet = boundaryDistanceData.distanceSet;
        distanceArraySize = boundaryDistanceData.distanceArraySize;
    }

    void setQueryToIndexMap(std::unordered_map<std::string, std::vector<double>> &queryToIndexMap) {
        queryToIndexMap_ = queryToIndexMap;
    }

    void gainQueryToIndexMap(std::unordered_map<std::string, std::vector<double>> &queryToIndexMap) {
        queryToIndexMap = queryToIndexMap_;
    }



    void setRampingCoefficient(const double * rampingCoefficient){
        rampingCoefficient_= rampingCoefficient;
    }



};

#endif //DENDRITEKT_NSEQUATION_H
