#ifndef DENDRITEKT_NSEQUATION_H
#define DENDRITEKT_NSEQUATION_H

#include <talyfem/fem/cequation.h>
#include "NSNodeData.h"
#include "NSInputData.h"
#include "NSDataStructure.h"
#include "SBMCalc.h"

class NSEquation : public TALYFEMLIB::CEquation<NSNodeData> {
    TALYFEMLIB::ZEROPTV forcing_;
    TALYFEMLIB::ZEROPTV forcingPre_;
    DENDRITE_REAL theta_ = 1.0;
    NSInputData *idata_;

    const IMGA *imga_;
    std::vector<my_kd_tree_t*> kd_trees_;

    /// we hard code here, if want, we can change this inside the config.txt
    const double Cb_f = 200;
    const double Ci_f = 36;
    const double orderOfhb = 1.0;
    const int BetaNS = 0.5;

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
            // Only calculate distance and normal once
            // But, there is no node ID available, so need to find it.
            unsigned int nodeID = 0;
            bool found_it = false;

            while (distanceSet[nodeID])
            {
                if (fabs(pt.x() - surfaceNodeCoords[nodeID * DIM]) < 1e-14)
                {
                    if (fabs(pt.y() - surfaceNodeCoords[nodeID * DIM + 1]) < 1e-14)
                    {
#if (DIM == 3)
                        if (fabs(pt.z() - surfaceNodeCoords[nodeID * DIM + 2]) < 1e-14)
#endif
                        {
                            // This is the node!
                            found_it = true;
                            break;
                        }
                    }
                }
                // Not the node.
                nodeID++;
            }

            /*
             * TODO: passing CarvedOutGeomID inside and the SBM geom is storing in imga_->geometry() to get different SBMGeo for different carved-out places
             */
            SBMCalc sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
            // Now we either found a set node, or need to compute a new one.
            if (found_it)
            {
                for (int i = 0; i < DIM; i++)
                {
                    d[i] = distanceToBdry[nodeID * DIM + i];
                }
                // below is for neumann
                // trueNormal = TALYFEMLIB::ZEROPTV(geoNormal[nodeID*DIM], geoNormal[nodeID*DIM+1], geoNormal[nodeID*DIM+2]);
            }
            else
            {
                // This is a new node. Compute things and store them.
                sbmCalc.Dist2Geo(d);
                // below is for neumann
                // sbmCalc.NormalofGeo(trueNormal, d);
                // Save them
                for (int i = 0; i < DIM; i++)
                {
                    distanceToBdry[nodeID * DIM + i] = d[i];
                    // need to set for neumann
                    // geoNormal[nodeID*DIM+i] = trueNormal.data()[i];
                    surfaceNodeCoords[nodeID * DIM + i] = pt.data()[i];
                }
                distanceSet[nodeID] = true;
            }
        }
        else
        {
            SBMCalc sbmCalc(fe, idata_, imga_, kd_trees_, CarvedOutGeomID);
            sbmCalc.Dist2Geo(d);
        }

#ifndef NDEBUG
        DistanceFuncChecking(fe, d);
#endif

    }


    int NearestGeo(const TALYFEMLIB::FEMElm &fe) {
        const TALYFEMLIB::ZEROPTV pt = fe.position();
        double Min_Dist = std::numeric_limits<double>::max();
        int PickGeomID = 0;

        for (int CarvedOutGeomID = 0;
             CarvedOutGeomID < imga_->getGeometries().size();
             ++CarvedOutGeomID) {

            double d[DIM];
            CalcDistFast(CarvedOutGeomID, pt, fe, d, true);

            double dist_square = std::inner_product(d, d + DIM, d, 0.0);

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

#ifndef NDEBUG
        DistanceFuncChecking(fe, d);
#endif

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

        TALYFEMLIB::ZEROPTV u;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = p_data_->valueFEM(fe, i);
            // u(i) = p_data_->valueFEM(fe, NSHTNodeData::VEL_X_PRE1 + i);
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

        TALYFEMLIB::ZEROPTV u;
        for (int i = 0; i < nsd; i++)
        {
            u(i) = p_data_->valueFEM(fe, i);
            // u(i) = p_data_->valueFEM(fe, NSHTNodeData::VEL_X_PRE1 + i);
        }

        double p = p_data_->valueFEM(fe, NSNodeData::PRESSURE);
        // double p = p_data_->valueFEM(fe, NSHTNodeData::PRESSURE_PRE1);

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
            }
        }
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


public:
    explicit NSEquation(NSInputData *idata)
            : TALYFEMLIB::CEquation<NSNodeData>(false, TALYFEMLIB::kAssembleGaussPoints) {
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

        if (idata_->ifMMS) {
            calcForcing(forcing_, fe.position(), t_);
        }



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
        double force_pre[nsd];
        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
            force_pre[i] = forcing_(i);
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


        double u_pre1[nsd];
        double u_pre2[nsd];
        double u_pre3[nsd];


        for (int i = 0; i < nsd; i++) {
            u_pre1[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre2[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
            u_pre3[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE3 + i);

        }

        std::vector<double> u_extrapolated = calc_linear_advection(u_pre1, u_pre2, idata_->VelocityExtrapolationOrder);
        VMSParams vms_params = calc_tau(fe, u_extrapolated);


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        /// Define velocity gradient tensor
        double du_pre1[nsd][nsd];
        double du_pre2[nsd][nsd];
        double dp_pre1[nsd];
        double dp_pre2[nsd];
        for (int i = 0; i < nsd; i++) {
            dp_pre1[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
            dp_pre2[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
            for (int j = 0; j < nsd; j++) {
                du_pre1[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
                du_pre2[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + i, j);
            }
        }


        double d2u_ijj_pre1[nsd];
        double d2u_jij_pre1[nsd];
        double d2u_ijj_pre2[nsd];
        double d2u_jij_pre2[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj_pre1[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dof, dir, dir);
                d2u_jij_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dir, dof, dir);
                d2u_ijj_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dof, dir, dir);
                d2u_jij_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dir, dof, dir);

            }
        }

        double u_full_pre1[nsd];
        double u_full_pre2[nsd];
        double res_M_pre1[nsd], res_M_pre2[nsd];
        double ui_ujj, uj_uij, ui_ujj_pre, uj_uij_pre;
        for (int i = 0; i < nsd; i++) {
            ui_ujj = 0;
            uj_uij = 0;
            ui_ujj_pre = 0;
            uj_uij_pre = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj += u_pre1[i] * du_pre1[j][j];
                uj_uij += u_pre1[j] * du_pre1[i][j];
                ui_ujj_pre += u_pre2[i] * du_pre2[j][j];
                uj_uij_pre += u_pre2[j] * du_pre2[i][j];
            }
            res_M_pre1[i] = temp_coeff * (bdf[0] * u_pre1[i] + bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt + adv_coeff * (ui_ujj + 0.5 * uj_uij) +
                            press_coeff * dp_pre1[i] - diff_coeff * (d2u_ijj_pre1[i] + d2u_jij_pre1[i]) - force_coeff * force[i];

            /// for temporal term we consider backward euler scheme to avoid storing one more degree of freedom
            res_M_pre2[i] = temp_coeff * (u_pre2[i] - u_pre3[i]) / dt + adv_coeff * (ui_ujj_pre + 0.5 * uj_uij_pre) +
                            press_coeff * dp_pre2[i] - diff_coeff * (d2u_ijj_pre2[i] + d2u_jij_pre2[i]) - force_coeff * force_pre[i];

            u_full_pre1[i] = u_pre1[i] - tauM * res_M_pre1[i];
            u_full_pre2[i] = u_pre2[i] - tauM * res_M_pre2[i];

        }


        std::vector<double> u_star = calc_linear_advection(u_full_pre1, u_full_pre2, idata_->VelocityExtrapolationOrder);

//        if (u[0] != 0)
//        printf ("u_star = %f %f\n", u[0], Re);




        ////////////////////////////////// end of precalculations




        ////////////////////////////////// calculate the derivative of the momentunm residual


        double res_M_Diag[nsd][nbf];
        double res_M_offDiag[nsd][nsd][nbf];

        for (int a = 0; a < nbf; a++) {
            for (int i = 0; i < nsd; i++) {
                res_M_Diag[i][a] = temp_coeff * (bdf[0] * fe.N(a) / dt);

                for (int j = 0; j < nsd; j++) {
                    res_M_Diag[i][a] += adv_coeff * (u_star[j] *  fe.dN(a, j) + 0.5 * du_pre1[j][j] * fe.N(a)) - diff_coeff * fe.d2N(a, i, j);
                    res_M_offDiag[i][j][a] = - diff_coeff *  fe.d2N(a, j, i);
                } // end of j loop
            } // end of i loop
        } // end of a loop


        ////////////////////////////////// end of calculation of the momentunm residual


//    printf("inside the integrands\n");


        for (int a = 0; a < nbf; a++) {
            for (int b = 0; b < nbf; b++) {
                for (int i = 0; i < nsd; i++) {


                    /// first part of momentum equation
                    /// Temporal term -- Term No. in the Document= 1
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += temp_coeff * fe.N(a) * (fe.N(b) / dt) * detJxW; ///ok!



                    for (int j = 0; j < nsd; j++) {


                        /// coarseAdvection 1 -- Term No. in the Document= 2
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.dN(a,j) * (u_star[j] * fe.N(b)) * detJxW; ///ok!



                        /// fineAdvection 1a -- Term No. in the Document= 3
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * res_M_Diag[i][b] * detJxW; ///ok!


                        /// fineAdvection 1b -- Term No. in the Document= 4
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * res_M_offDiag[i][j][b] * detJxW; ///ok!


                        /// fineAdvection 1c -- Term No. in the Document= 5
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.dN(a,j) * u_star[j] * tauM * fe.dN(b, i) * detJxW; ///ok!


                        /// coarseAdvection 2 -- Term No. in the Document= 6
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) -= adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * fe.N(b) * detJxW; ///ok!


                        /// fineAdvection 2a -- Term No. in the Document= 7
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * res_M_Diag[i][b] * detJxW; /// ok!


                        /// fineAdvection 2b -- Term No. in the Document= 8
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * res_M_offDiag[i][j][b] * detJxW; /// ok!


                        /// fineAdvection 2c -- Term No. in the Document= 9
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) += adv_coeff * fe.N(a) * 0.5 * du_pre1[j][j] * tauM * fe.dN(b, i) * detJxW; /// ok!


                        /// coarseDiffusion 1 -- Term No. in the Document= 12
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += diff_coeff * fe.dN(a,j) * fe.dN(b, j) * detJxW; ///ok!

                        /// coarseDiffusion 2 -- Term No. in the Document= 13
                        Ae((nsd + 1) * a + i, (nsd + 1) * b + j) += diff_coeff * fe.dN(a,j) * fe.dN(b, i) * detJxW; ///ok!

                        if (idata_->DiffFineTerms) {
                            /// fineDiffusion 1 -- Term No. in the Document= 14
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + j) +=
                                    diff_coeff * fe.d2N(a, j, i) * tauM * (res_M_Diag[j][b] + res_M_offDiag[j][j][b]) *
                                    detJxW; ///ok!

                            /// fineDiffusion 1 -- Term No. in the Document= 15
                            Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) +=
                                    diff_coeff * fe.d2N(a, j, i) * tauM * fe.dN(b, i) * detJxW; ///ok!

                        }

                        /// fineContinuity 1b -- Term No. in the Document= 18
                        Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * res_M_offDiag[i][j][b] * detJxW; ///ok!



                    } // end of j loop



                    /// coarsePressure 1 -- Term No. in the Document= 10
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + nsd) -= press_coeff * fe.dN(a, i) * fe.N(b) * detJxW; ///ok!

                    /// finePressure 1 -- Term No. in the Document= 11
                    Ae((nsd + 1) * a + i, (nsd + 1) * b + i) += press_coeff * fe.dN(a, i) * tauC * fe.dN(b, i) * detJxW; ///ok!


                    /// coarseContinuity 1 -- Term No. in the Document= 16
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.N(a) * fe.dN(b, i) * detJxW; ///ok!


                    /// fineContinuity 1a -- Term No. in the Document= 17
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + i) += fe.dN(a, i) * tauM * res_M_Diag[i][b] * detJxW; ///ok!


                    /// fineContinuity 1c -- Term No. in the Document= 19
                    Ae((nsd + 1) * a + nsd, (nsd + 1) * b + nsd) += fe.dN(a, i) * tauM * fe.dN(b, i) * detJxW; ///ok!


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

        if (idata_->ifMMS) {
            calcForcing(forcing_, fe.position(), t_);
        }


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
        double force_pre[nsd];
        for (int i = 0; i < nsd; i++) {
            force[i] = forcing_(i);
            force_pre[i] = forcing_(i);
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


        double u_pre1[nsd];
        double u_pre2[nsd];
        double u_pre3[nsd];


        for (int i = 0; i < nsd; i++) {
            u_pre1[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE1 + i);
            u_pre2[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE2 + i);
            u_pre3[i] = this->p_data_->valueFEM(fe, NSNodeData::VEL_X_PRE3 + i);

        }

        std::vector<double> u_extrapolated = calc_linear_advection(u_pre1, u_pre2, idata_->VelocityExtrapolationOrder);
        VMSParams vms_params = calc_tau(fe, u_extrapolated);


        /// Calculate the tau_M based on the projections calculated, Equation 63 in
        /// Bazilevs et al. (2007), here we use scale as a parameter to tune tauM,
        /// in the paper it is set to 1.
        const double tauM = vms_params.tauM;
        /// Calculate continuity residual based on
        const double tauC = vms_params.tauC;

        /// Define velocity gradient tensor
        double du_pre1[nsd][nsd];
        double du_pre2[nsd][nsd];
        double dp_pre1[nsd];
        double dp_pre2[nsd];
        for (int i = 0; i < nsd; i++) {
            dp_pre1[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE1, i);
            dp_pre2[i] = this->p_data_->valueDerivativeFEM(fe, NSNodeData::PRESSURE_PRE2, i);
            for (int j = 0; j < nsd; j++) {
                du_pre1[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + i, j);
                du_pre2[i][j]= this->p_data_->valueDerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + i, j);
            }
        }


        double d2u_ijj_pre1[nsd];
        double d2u_jij_pre1[nsd];
        double d2u_ijj_pre2[nsd];
        double d2u_jij_pre2[nsd];
        /// loop for three directions of velocity (MAX NSD is no of velocity
        /// directions in the problem)
        /// PS: This assumes the first three degrees are always velocity, it is
        /// wise to keep it that way
        for (int dof = 0; dof < nsd; dof++) {
            /// Summing over three directions of velocity as it is laplacian
            d2u_ijj_pre1[dof] = 0;
            for (int dir = 0; dir < nsd; dir++) {
                d2u_ijj_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dof, dir, dir);
                d2u_jij_pre1[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE1 + dir, dof, dir);
                d2u_ijj_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dof, dir, dir);
                d2u_jij_pre2[dof] += this->p_data_->value2DerivativeFEM(fe, NSNodeData::VEL_X_PRE2 + dir, dof, dir);

            }
        }

        double u_full_pre1[nsd];
        double u_full_pre2[nsd];
        double res_M_pre1[nsd], res_M_pre2[nsd];
        double ui_ujj, uj_uij, ui_ujj_pre, uj_uij_pre;
        for (int i = 0; i < nsd; i++) {
            ui_ujj = 0;
            uj_uij = 0;
            ui_ujj_pre = 0;
            uj_uij_pre = 0;
            for (int j = 0; j < nsd; j++) {
                ui_ujj += u_pre1[i] * du_pre1[j][j];
                uj_uij += u_pre1[j] * du_pre1[i][j];
                ui_ujj_pre += u_pre2[i] * du_pre2[j][j];
                uj_uij_pre += u_pre2[j] * du_pre2[i][j];
            }
            res_M_pre1[i] = temp_coeff * (bdf[0] * u_pre1[i] + bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt + adv_coeff * (ui_ujj + 0.5 * uj_uij) +
                            press_coeff * dp_pre1[i] - diff_coeff * (d2u_ijj_pre1[i] + d2u_jij_pre1[i]) - force_coeff * force[i];

            /// for temporal term we consider backward euler scheme to avoid storing one more degree of freedom
            res_M_pre2[i] = temp_coeff * (u_pre2[i] - u_pre3[i]) / dt + adv_coeff * (ui_ujj_pre + 0.5 * uj_uij_pre) +
                            press_coeff * dp_pre2[i] - diff_coeff * (d2u_ijj_pre2[i] + d2u_jij_pre2[i]) - force_coeff * force_pre[i];

            u_full_pre1[i] = u_pre1[i] - tauM * res_M_pre1[i];
            u_full_pre2[i] = u_pre2[i] - tauM * res_M_pre2[i];

        }


        std::vector<double> u_star = calc_linear_advection(u_full_pre1, u_full_pre2, idata_->VelocityExtrapolationOrder);



        ////////////////////////////////// end of precalculations








        double res_M_RHS[nsd];
        for (int i = 0; i < nsd; i++) {
            res_M_RHS[i] = temp_coeff * (bdf[1] * u_pre2[i] + bdf[2] * u_pre3[i]) / dt - force_coeff * force[i];
        }


        for (int a = 0; a < nbf; a++) {
            for (int i = 0; i < nsd; i++) {
                for (int j = 0; j < nsd; j++) {

                    /// fineAdvection RHS 1 -- Term No. in the Document= 1
                    be((nsd + 1) * a + i) -= fe.dN(a, j) * adv_coeff * u_star[i] * tauM * res_M_RHS[i] * detJxW; ///ok!

                    /// fineAdvection RHS 2 -- Term No. in the Document= 2
                    be((nsd + 1) * a + i) -= fe.N(a) * adv_coeff * tauM * 0.5 * du_pre1[j][j] * res_M_RHS[i] * detJxW; ///ok!

                    if (idata_->DiffFineTerms) {
                        /// fineDiffusion RHS 1 -- Term No. in the Document= 3
                        be((nsd + 1) * a + i) -= diff_coeff * fe.d2N(a, j, i) * tauM * res_M_RHS[i] * detJxW; ///ok!
                    }

                    /// force -- Term No. in the Document= 4
                    be((nsd + 1) * a + i) += force_coeff * fe.N(a) * force[i] * detJxW; ///ok!


                } // end of j loop

                /// fineContinuity RHS 1 -- Term No. in the Document= 5
                be((nsd + 1) * a + nsd) -= fe.dN(a, i) * tauM * res_M_RHS[i] * detJxW; ///ok!

            } // end of i loop

        } // end of a loop





    }


    VMSParams calc_tau(const TALYFEMLIB::FEMElm &fe, std::vector<double> u) const {
        using namespace TALYFEMLIB;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;
        const double Ci_f = 36;
        const int nsd = DIM;

        VMSParams params;

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
                u_Gu += u[i] * Ge(i, j) * u[j];
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


    std::vector<double> calc_linear_advection(double u[DIM], double u_prev[DIM], int option) {

        std::vector<double> u_star(DIM);

        switch (option) {
            case 1:
                for (int i = 0; i < DIM; i++) {
                    u_star[i] = u[i];
                }
                break;
            case 2:
                for (int i = 0; i < DIM; i++) {
                    u_star[i] = 2.0 * u[i] - u_prev[i];
                }
                break;
        }
        return u_star;
    }



    void ReRamping (double t) {
        if (t <= idata_->RampingTime) {
            idata_->Re = ((idata_->RampingRe - 10) / idata_->RampingTime) * t + 10;
        } else {
            idata_->Re = idata_->Re_tmp;
        }
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
        calcAe_vms(fe, Ae, bdf);
    }

    void Integrands_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be, double *h) {
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
        calcbe_vms(fe, be, bdf);
    }


#ifdef IBM
    void ibm_Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROMATRIX<double> &Ae,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h)
    {

        throw std::runtime_error("You should set true normal before using this function.");

        if (idata_->DoReRamping) {
            ReRamping(t_);
        }

        TALYFEMLIB::ZEROPTV normal;

        const int nsd = fe.nsd();
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

        double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, orderOfhb);

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
                }
            }
        }

    }

    void ibm_Integrands4side_be(const TALYFEMLIB::FEMElm &fe, TALYFEMLIB::ZEROARRAY<double> &be,
                                const NodeAndValues<DENDRITE_REAL> &gpinfo,
                                const TALYFEMLIB::ZEROPTV &position,
                                const TALYFEMLIB::ZEROPTV &h)
    {
//        auto &geo = idata_->carved_out_geoms_def.at(gpinfo.geomID);
//
//        /// If the geometry is 2D plane, we should treat this differently.
//        ZEROPTV normal = ZEROPTV{gpinfo.normal[0], gpinfo.normal[1], 0};
//#if (DIM == 3)
//        normal(2) = gpinfo.normal[2];
//#endif
//
//        const ZEROPTV &vel = geo.ic_vel; // todo this would be different for "rotating geometry" -> IMGA moving code is different
//        const int nsd = fe.nsd();
//        const int nbf = fe.nbf();
//        const double Cb_f = idata_->Cb_f.value_at(t_);
//        const double Coe_weak = Cb_f;
//        const double detSideJxW = fe.detJxW();
//        double Coe_diff = idata_->getDiffusionNS(t_);
//        double Coe_bdyforce = idata_->getBodyForce(t_);
//        /** Get u from the NodeData
//         * NOTE: ValueFEM is defined on NodeData object, therefore it follows
//         * indices defined in NodeData subclass.
//         * NOTE: Every iteration the solution is copied into NodeData, we use these
//         * solved fields from the NodeData object to calculate the be (NS residual)
//         * of current iteration.*/
//
//        ZEROPTV u, u_pre;
//        for (int i = 0; i < nsd; i++)
//        {
//            u(i) = this->p_data_->valueFEM(fe, i);
//            u_pre(i) = this->p_data_->valueFEM(fe, NSHTNodeData::NUM_VARS + i);
//        }
//
//        /// Define velocity gradient tensor ( grad{u} )
//        ZeroMatrix<double> du;
//        du.redim(nsd, nsd);
//        for (int i = 0; i < nsd; i++)
//        {
//            for (int j = 0; j < nsd; j++)
//            {
//                du(i, j) = this->p_data_->valueDerivativeFEM(fe, i, j);
//            }
//        }
//
//        double p = this->p_data_->valueFEM(fe, NSHTNodeData::PRESSURE);
//
//        ZEROPTV diff1, diff2;
//        for (int i = 0; i < nsd; i++)
//        {
//            diff1(i) = 0.0;
//            diff2(i) = 0.0;
//            for (int j = 0; j < nsd; j++)
//            {
//                diff1(i) += du(i, j) * normal(j);
//                diff2(i) += du(j, i) * normal(j);
//            }
//        }
//
//        /// wall normal distance (the maximum of 8 nodes)
//        /// wall normal distance
//        double hb = normalDistance(fe, normal, h);
////    double hb = h.x();
//#ifndef NDEBUG
//        if (hb < 0)
//        {
//            //      PrintStatus("Elm node:\t", ZEROPTV{gpinfo.node., gpinfo.node.y(), gpinfo.node.z()});
//            PrintStatus("gp_pos:\t", fe.position());
//            PrintStatus("normal:\t", normal);
//            //      for (int i = 0; i < 3; ++i) {
//            //        PrintStatus("triangle:\t", gpinfo.tri_nodes[i]);
//            //      }
//        }
//#endif
//
//        double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(hb, idata_->orderOfhb);
//        if (idata_->weakBCglobal)
//        {
//            weakBCpenaltyParameter = Cb_f;
//        }
//        //    if (weakBCpenaltyParameter < 4) {
//        //      PrintWarning("WeakBC parameter less than 4");
//        //    }
//
//        /// Based on what kind of solver chosen, choose the corresponding solver (SUPG or VMS)
//        if (idata_->solverTypeNS == NSHTInputData::STABILIZED_NS)
//        {
//            for (int a = 0; a < nbf; a++)
//            {
//                double supg = 0.0;
//                for (int m = 0; m < nsd; m++)
//                {
//                    supg += fe.N(a) * normal(m);
//                }
//
//                ZEROPTV supg1, supg2;
//                for (int m = 0; m < nsd; m++)
//                {
//                    supg1(m) = 0.0;
//                    supg2(m) = 0.0;
//                    for (int n = 0; n < nsd; n++)
//                    {
//                        supg1(m) += fe.dN(a, n) * normal(n);
//                        supg2(m) += fe.dN(a, m) * normal(n);
//                    }
//                }
//
//                for (int j = 0; j < nsd; j++)
//                {
//                    be((nsd + 1) * a + j) +=
//                            /// pressure
//                            fe.N(a) * (p * normal(j)) * detSideJxW
//                            /// diffusion
//                            - fe.N(a) * (Coe_diff * (diff1(j) + diff2(j))) * detSideJxW
//                            /// adjoint diffusion
//                            - Coe_diff * (supg1(j) + supg2(j)) * (u(j) - vel(j)) * detSideJxW
//                            /// penalty
//                            + weakBCpenaltyParameter * fe.N(a) * (u(j) - vel(j)) * detSideJxW;
//                    //                  + Coe_weak * fe.N(a) * (u(j) - vel(j)) * detSideJxW;
//
//                    /// adjoint pressure
//                    be((nsd + 1) * (a + 1) - 1) += -fe.N(a) * normal(j) * (u(j) - vel(j)) * detSideJxW;
//                }
//            }
//        }
//        else if (idata_->solverTypeNS == NSHTInputData::RBVMS_NS)
//        {
//            for (int a = 0; a < nbf; a++)
//            {
//                ZEROPTV supg1;
//                for (int m = 0; m < nsd; m++)
//                {
//                    supg1(m) = 0.0;
//                    for (int n = 0; n < nsd; n++)
//                    {
//                        supg1(m) += fe.dN(a, n) * normal(n);
//                    }
//                }
//
//                for (int j = 0; j < nsd; j++)
//                {
//                    be((nsd + 1) * a + j) +=
//                            /// pressure
//                            fe.N(a) * (p * normal(j)) * detSideJxW
//                            /// diffusion
//                            - fe.N(a) * (Coe_diff * (diff1(j) /*+diff2(j)*/)) * detSideJxW
//                            /// adjoint diffusion
//                            - Coe_diff * supg1(j) * (u(j) - vel(j)) * detSideJxW
//                            /// penalty
//                            + weakBCpenaltyParameter * fe.N(a) * (u(j) - vel(j)) * detSideJxW;
//                    //                  + Coe_weak * fe.N(a) * (u(j) - vel(j)) * detSideJxW;
//
//                    /// adjoint pressure
//                    be((nsd + 1) * (a + 1) - 1) += -fe.N(a) * normal(j) * (u(j) - vel(j)) * detSideJxW;
//                }
//            }
//        }
//        else if (idata_->solverTypeNS == NSHTInputData::LINEAR_NS)
//        {
//            for (int a = 0; a < nbf; a++)
//            {
//                ZEROPTV supg1, supg2;
//                for (int m = 0; m < nsd; m++)
//                {
//                    supg1(m) = 0.0;
//                    supg2(m) = 0.0;
//                    for (int n = 0; n < nsd; n++)
//                    {
//                        supg1(m) += fe.dN(a, n) * normal(n);
//                        supg2(m) += fe.dN(a, m) * normal(n);
//                    }
//                }
//
//                for (int j = 0; j < nsd; j++)
//                {
//                    be((nsd + 1) * a + j) +=
//                            /// adjoint diffusion
//                            +Coe_diff * (supg1(j) + supg2(j)) * (-vel(j)) * detSideJxW
//                            /// penalty
//                            - weakBCpenaltyParameter * fe.N(a) * (-vel(j)) * detSideJxW;
//                    ;
//
//                    /// adjoint pressure (inside formulation, it is with adjoint diffusion)
//                    be((nsd + 1) * (a + 1) - 1) += fe.N(a) * normal(j) * (-vel(j)) * detSideJxW;
//                }
//            }
//        }
//#ifndef NDEBUG
//        if (idata_->Debug_Integrand)
//        {
//            int interchangedmat[8]{0, 1, 3, 2, 4, 5, 7, 6};
//            FILE *fp = fopen("be_ibm_dendro5.dat", "w");
//            for (int i = 0; i < 8; ++i)
//            {
//                for (int dof1 = 0; dof1 < 4; ++dof1)
//                {
//                    fprintf(fp, "%.10f\t", be(interchangedmat[i] * 4 + dof1));
//                }
//            }
//            fprintf(fp, "\n");
//            fclose(fp);
//        }
//#endif
    }
#endif

#pragma mark Ae-be Surface

    void Integrands4side_Ae(const TALYFEMLIB::FEMElm &fe, int side_idx, int id, TALYFEMLIB::ZeroMatrix<double> &Ae)
    {
        if (GpOnDomainWall(fe))
        {
            /// side_idx in the order of x-, x+, y-, y+, z-, z+
            const auto &bc = idata_->boundary_def.at(side_idx);
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
        if (bc.pressure_type == BoundaryDef::Pressure_Type::BACKFLOW_STAB)
        {

            const int nbf = fe.nbf();
            const int nsd = idata_->nsd;
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
        if (bc.pressure_type == BoundaryDef::Pressure_Type::BACKFLOW_STAB)
        {
            const int nbf = fe.nbf();
            const int nsd = idata_->nsd;
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
                            -bc.BetaNS * fe.N(a) * Uneg * u(j) * detSideJxW;
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

        int ndof = NSNodeData::NS_DOF;

        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;
        const int nbf = fe.nbf();
        const int nsd = idata_->nsd;
        const double detSideJxW = fe.detJxW();

        const auto &bcType = carved_out_def.bc_type_V[NSNodeData::VEL_X];
        if (bcType == CarvedOutGeom::BCType::WEAK)
        {
//            std::cout << "weak Ae \n";

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
//            std::cout << "SBM Ae \n";

            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);

            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_Ae(Ae, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d);
        }

        else if (bcType == CarvedOutGeom::BCType::PENALTY_SBM)
        {
//            std::cout << "P SBM Ae \n";

            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);
            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_Ae(Ae, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d, true /*PenaltyOnly*/);
        }

#endif
    }

    void Integrands4side_be_carved(const TALYFEMLIB::FEMElm &fe,
                                   TALYFEMLIB::ZEROARRAY<double> &be,
                                   const CarvedOutGeom &carved_out_def,
                                   const int CarvedOutGeomID)
    {
        const auto &bcType = carved_out_def.bc_type_V[NSNodeData::VEL_X];
        int ndof = NSNodeData::NS_DOF;
        double Re = idata_->Re;
        const double Coe_diff = 1.0 / Re;

        const int nbf = fe.nbf();
        const int nsd = idata_->nsd;
        const double detSideJxW = fe.detJxW();

        if (bcType == CarvedOutGeom::BCType::WEAK)
        {
//            std::cout << "weak be \n";

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
//            std::cout << "be SBM \n";

            double d[DIM];
            CalculateDistances(fe,CarvedOutGeomID,d);
            /*
             * minor change 1: instead of using hb, we use h because of SBM
             */
            double weakBCpenaltyParameter = Cb_f * Coe_diff / ElementSize(fe);
            // double weakBCpenaltyParameter = Cb_f * Coe_diff / pow(CalcHb(fe), idata_->orderOfhb);

            SBM_be(be, fe, nbf, ndof, nsd, detSideJxW, weakBCpenaltyParameter, Coe_diff, d,carved_out_def);
        }
        if (bcType == CarvedOutGeom::BCType::PENALTY_SBM)
        {
//            std::cout << "P SBM be \n";

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


    void setSbmCalc(const IMGA *imga, const std::vector<my_kd_tree_ptr> &kd_trees)
    {
        imga_ = imga;
        std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
                       [](const auto& uniquePtr) { return uniquePtr.get(); });
    }

    void setDistanceArrays(BoundaryDistanceData &boundaryDistanceData)
    {
        distanceToBdry = boundaryDistanceData.distanceToBdry;
        bdryNormal = boundaryDistanceData.bdryNormal;
        surfaceNodeCoords = boundaryDistanceData.surfaceNodeCoords;
        distanceSet = boundaryDistanceData.distanceSet;
        distanceArraySize = boundaryDistanceData.distanceArraySize;
    }



};

#endif //DENDRITEKT_NSEQUATION_H
