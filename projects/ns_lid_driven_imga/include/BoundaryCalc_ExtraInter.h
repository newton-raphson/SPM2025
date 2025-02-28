//
// Created by chenghau on 7/30/23.
//

#ifndef NSIBM_KT_BOUNDARYCALC_EXTRAINTER_H
#define NSIBM_KT_BOUNDARYCALC_EXTRAINTER_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include <NSInputData.h>
#include <unordered_map>
#include "NSNodeData.h"
#include "TimeInfo.h"
#include <unordered_map>
#ifdef DEEPTRACE
#include "SBMCalcDeepTrace.h"
#else
#include "SBMCalc.h"
#endif

using namespace TALYFEMLIB;


// construct a kd-tree index:
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
        PointCloud<double>, 3 /* dim */>;


struct GPValuesWithDist{
    ZEROPTV position;
    DENDRITE_REAL value;
    DENDRITE_REAL distance;
};

struct GpZEROPTV{
    ZEROPTV position;
    ZEROPTV ZeroPtV_values;
    DENDRITE_REAL value;
};

class BoundaryCalc_ExtraInter : public IMGATraversal
{

    TimeInfo ti_;

    Point<DIM> shift_;
    double valueWJ = 0.0;
    double Length = 0.0;
    double totalGP = 0.0;
    double totalGPGlobal = 0.0;
    double ErrorSum = 0.0;
    double ErrorSumGlobal = 0.0;
    double MaxError = 0.0;
    double MaxErrorGlobal = 0.0;

    double x_true, y_true;

    std::vector<ZEROPTV> DomainPtLocal;
    std::vector<ZEROPTV> DomainPtGlobal;
    std::vector<ZEROPTV> GaussPtLocal;
    std::vector<ZEROPTV> GaussPtGlobal;
    std::vector<ZEROPTV> NormalLocal;
    std::vector<ZEROPTV> NormalGlobal;
    std::vector<ForceInfo> localforce_infos;

    const IMGA *imga_;
    NSInputData &idata_;
    const my_kd_tree_t *kd_tree_trueGP_;
    const my_kd_tree_t *kd_tree_surrogateGP_;
    std::vector<ZEROPTV> TrueGPnormalGlobal_;
    std::vector<double> TrueGPAreaGlobal_;

    std::vector<std::vector<GPValuesWithDist>> GPValuesWithDist_Geo;
    std::vector<std::vector<GPValuesWithDist>> GlobalGPValuesWithDist_Geo;
    std::vector<std::vector<GpZEROPTV>> GPBoundaryErrorWithDist_Geo;
    std::vector<std::vector<GpZEROPTV>> GlobalGPBoundaryErrorWithDist_Geo;

    const VecInfo v_;

    std::vector<double> velocityWJ;

    std::vector<my_kd_tree_t *> kd_trees_;

    // Map to store the mapping from query to ret_index
    std::unordered_map<std::string, uint32_t> queryToIndexMap_;
    // Utility function to convert a query point to a string key
    std::string queryToStringKey(const double query_pt[3]);

    // Function to handle the retrieval or computation of ret_index for a given query point
    uint32_t getRetIndexForQuery(const double query_pt[3]);

    double ElementSize(const FEMElm &fe);
    double linear_extrapolation(ZEROPTV p1, ZEROPTV p2, double value1, double value2, ZEROPTV p3);
    void calculateDiff(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values, const double Coe_diff,
                       const ZEROPTV &normal, ZEROPTV &diff);

    void calculatePressure(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values,
                           double &p);

    void calculateVelocity(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values,
                           ZEROPTV &vel);

    void NSBoundaryError(std::vector<double> &boundaryError);

    int NearestGeo(const FEMElm &fe);

public:
    /**
     * @brief constructor: loop over the surface-based Gauss Points to calculate the L2 norm
     * @param octDA
     * @param imga imga context. we use this to access Gauss Points on a surface
     * @param treePart treePartition
     * @param v vecInfo for syncing vector
     * @param domain Domain information
     * @param sbmgeo the specific SBM geometry
     */
    BoundaryCalc_ExtraInter(DA *octDA, const IMGA *imga, const IMGA *imga_TruePt_InterpPT, const std::vector<TREENODE> &treePart, const VecInfo &v,
                            const DomainExtents &domain, NSInputData &idata,
                            const TimeInfo &ti, /*const std::unordered_map<PointData, PointData, PointHash>& in_point_map */
                            const my_kd_tree_t *kdTree_TrueGP,
                            const my_kd_tree_t *kdTree_SurrogateGP,
                            const std::vector<my_kd_tree_ptr> & kd_trees,
                            std::vector<ZEROPTV> &TrueGPnormalGlobal,
                            std::vector<double> &TrueGPAreaGlobal, std::unordered_map<std::string, uint32_t> &queryToIndexMap);

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;
    void computeBoundaryError(DENDRITE_REAL *boundaryValue);
    void computeForce(std::vector<ForceInfo> &globalForces);
    void forceCalcIBM();


    void WriteTrueGPToFile();
    void WriteGPValue();

    void WriteGPBoundaryError();

    void WriteNSBoundaryError();

    // Function to gain queryTOIndexMap_
    void gainQueryToIndexMap(std::unordered_map<std::string, uint32_t> &queryToIndexMap);

};



BoundaryCalc_ExtraInter::BoundaryCalc_ExtraInter(DA *octDA, const IMGA *imga, const IMGA *imga_TruePt_InterpPT, const std::vector<TREENODE> &treePart, const VecInfo &v,
                                                 const DomainExtents &domain, NSInputData &idata,
                                                 const TimeInfo &ti, /*const std::unordered_map<PointData, PointData, PointHash>& in_point_map */
                                                 const my_kd_tree_t *kdTree_TrueGP,
                                                 const my_kd_tree_t *kdTree_SurrogateGP,
                                                 const std::vector<my_kd_tree_ptr> & kd_trees,
                                                 std::vector<ZEROPTV> &TrueGPnormalGlobal,
                                                 std::vector<double> &TrueGPAreaGlobal,
                                                 std::unordered_map<std::string, uint32_t> &queryToIndexMap)
        : IMGATraversal(octDA, imga_TruePt_InterpPT, treePart, v, domain), shift_(Point<DIM>(0,0,0)), imga_(imga), idata_(idata),
          kd_tree_trueGP_(kdTree_TrueGP), kd_tree_surrogateGP_(kdTree_SurrogateGP),
          TrueGPnormalGlobal_(TrueGPnormalGlobal), TrueGPAreaGlobal_(TrueGPAreaGlobal), v_(v), ti_(ti),
          queryToIndexMap_(queryToIndexMap)
{
#ifndef DEEPTRACE
    std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
                   [](const auto& uniquePtr) { return uniquePtr.get(); });

#endif

    velocityWJ.resize(DIM);
    GPValuesWithDist_Geo.resize(imga->getGeometries().size());
    GlobalGPValuesWithDist_Geo.resize(imga->getGeometries().size());
    GPBoundaryErrorWithDist_Geo.resize(imga->getGeometries().size());
    GlobalGPBoundaryErrorWithDist_Geo.resize(imga->getGeometries().size());
    localforce_infos.resize(imga->getGeometries().size());
    this->imgaTraverse();
}

void BoundaryCalc_ExtraInter::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                                     const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values)
{

    // gaussPoint.normal = contains the position of point on surface. => important
    // gaussPoint.location = contains the position of point that we need to interpolate/

    const double query_pt[3] = {gaussPoint.normal[0], gaussPoint.normal[1], gaussPoint.normal[2]}; // shift_
    uint32_t savedIndex = getRetIndexForQuery(query_pt);

//    std::cout << "savedIndex: " << savedIndex << std::endl;
    ZEROPTV normal = TrueGPnormalGlobal_[savedIndex];

    static constexpr DENDRITE_REAL weight = (DIM==2) ? 2.0 : 3.0;
    double jacc_x_w=TrueGPAreaGlobal_[savedIndex]/weight;


    double Re = idata_.Re;
    const double Coe_diff = 1.0 / Re;

    DomainPtLocal.push_back(fe.position());
    TALYFEMLIB::ZEROPTV TrueGPloc;
    std::memcpy(TrueGPloc.data(), gaussPoint.normal, sizeof(double)*DIM);
    GaussPtLocal.push_back(TrueGPloc);

    /*
     * fe: "move" inside the surrogate domain
     * fe2: "move more" inside the surrogate domain
     * We use these two points to extrapolate to get the values on the true boundary
     */
    TALYFEMLIB::ZEROPTV InsideTrueloc;
    std::memcpy(InsideTrueloc.data(), gaussPoint.location, sizeof(double)*DIM);
    TALYFEMLIB::ZEROPTV ptvg, ptvl;
    ZEROPTV Vector = InsideTrueloc - TrueGPloc;
    Vector.SafeNormalize();
    double scale = 0.05;
    ptvg=  InsideTrueloc + Vector * ElementSize(fe) *scale;
    TALYFEMLIB::FEMElm fe2 = fe;
    GetLocalPtv(fe2, ptvg, ptvl);
    fe2.calc_at(ptvl);

    double p1, p2, p;
    ZEROPTV diff1,diff2,diff,vel, vel1, vel2,val_a /* analytical solution*/;
    if (fe.position().distanceTo(TrueGPloc) < 1e-12/*Evalulated at True Boundary*/){
        calculatePressure(fe, values, p);
        calculateDiff(fe,values,Coe_diff,normal,diff);
        calculateVelocity(fe, values, vel);
    } else {
        calculatePressure(fe, values, p1);
        calculatePressure(fe2, values, p2);
        calculateVelocity(fe, values, vel1);
        calculateVelocity(fe2, values, vel2);

        p = linear_extrapolation(fe2.position(), fe.position(), p2, p1, TrueGPloc);
        calculateDiff(fe,values,Coe_diff,normal,diff1);
        calculateDiff(fe2,values,Coe_diff,normal,diff2);
        for (int i=0;i<DIM;i++){
            diff(i) = linear_extrapolation(fe2.position(),fe.position(),diff2(i),diff1(i),TrueGPloc);
            vel(i) = linear_extrapolation(fe2.position(),fe.position(),vel2(i),vel1(i),TrueGPloc);
        }
    }

#ifndef  NDEBUG
    if (std::isnan(p)) {
        std::cout << "The value is NaN" << std::endl;
    }
#endif



    // TODO: support multiple geometries
    GPValuesWithDist GPValuesWithDist_local;
    for (int dim = 0; dim < DIM; dim++) {
        GPValuesWithDist_local.position = fe.position();
    }

    GPValuesWithDist_local.value = p;
    GPValuesWithDist_local.distance = 0.0; /*TODO: we skip the distance here for faster calculations*/

    GPValuesWithDist_Geo[0].push_back(GPValuesWithDist_local);



    ZEROPTV Force_pressure;
    ZEROPTV Force_viscous;
    for (int dim = 0; dim < DIM; dim++) {
        Force_pressure(dim) = (p * normal(dim)) * jacc_x_w;
        Force_viscous(dim) = -(diff(dim)) * jacc_x_w;
    }

    // TODO: support multiple geometries
    localforce_infos[0].Force_pressure += Force_pressure;
    localforce_infos[0].Force_viscous += Force_viscous;

    GpZEROPTV gpZeroptv;
    gpZeroptv.position = fe.position();

    /*
     * [NOTE] we do not set any value to val_a because the boundary error case we test is flow passing cylinder.
     * User should set other val_a if it is not no-slip at the boundary.
     */
    /*
     * We test integrating only detSideJxW and the results make sense
     */

    for (int dim = 0; dim < DIM; dim++) {
        velocityWJ[dim] += (vel[dim] - val_a[dim]) * (vel[dim] - val_a[dim]) * jacc_x_w;
        /// hack distance be be error
        gpZeroptv.ZeroPtV_values(dim) = fabs(vel[dim] - val_a[dim]);
    }

    gpZeroptv.value = p;


    GPBoundaryErrorWithDist_Geo[0].push_back(gpZeroptv);

}

void BoundaryCalc_ExtraInter::computeBoundaryError(DENDRITE_REAL *boundaryValue)
{
    MPI_Reduce(&valueWJ, &boundaryValue[0], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
    boundaryValue[0] = sqrt(boundaryValue[0]);
    MPI_Reduce(&totalGP, &totalGPGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
    boundaryValue[1] = totalGPGlobal;
    MPI_Reduce(&ErrorSum, &ErrorSumGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
    MPI_Reduce(&Length, &boundaryValue[2], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
    boundaryValue[3] = ErrorSumGlobal/totalGPGlobal;
    MPI_Reduce(&MaxError, &boundaryValue[4], 1, MPI_DOUBLE, MPI_MAX, 0, m_octDA->getCommActive());

}

void BoundaryCalc_ExtraInter::computeForce(std::vector<ForceInfo> &globalForces)
{
    assert(localforce_infos.size() == globalForces.size());
    for (int i = 0; i < localforce_infos.size(); i++)
    {
        MPI_Allreduce(&localforce_infos[i].Force_pressure, &globalForces[i].Force_pressure, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
        MPI_Allreduce(&localforce_infos[i].Force_viscous, &globalForces[i].Force_viscous, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
        MPI_Allreduce(&localforce_infos[i].Force_penalty, &globalForces[i].Force_penalty, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
        MPI_Allreduce(&localforce_infos[i].mass_conservation, &globalForces[i].mass_conservation, 1, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
        MPI_Allreduce(&localforce_infos[i].Torque, &globalForces[i].Torque, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
    }
}

void BoundaryCalc_ExtraInter::forceCalcIBM()
{
    std::vector<ForceInfo> globalForces;
    globalForces.resize(imga_->getGeometries().size());
    computeForce(globalForces);

    for (int i = 0; i < imga_->getGeometries().size(); i++)
    {
        // print to screen
        PrintStatus(idata_.carved_out_geoms_def[i].name, " force: ", globalForces[i].Force_all(true));

        // write to file
        if (TALYFEMLIB::GetMPIRank() == 0)
        {
            std::string fname = "Force_" + std::to_string(i) + "_" + idata_.carved_out_geoms_def[i].name + ".dat";
            FILE *fp = fopen(fname.c_str(), "a");
            if (!fp)
            {
                throw TALYException() << "Cannot create file: " << fname;
            }
#if (DIM == 3)
            fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Fx_pre = % .10e, Fy_pre = % .10e, Fz_pre = % .10e\n"
                "Fx_vis = % .10e, Fy_vis = % .10e, Fz_vis = % .10e\n"
                "Fx_pen = % .10e, Fy_pen = % .10e, Fz_pen = % .10e\n"
                "mass_conservation = % .10e\n"
                "Tx = % .10e, Ty = % .10e, Tz = % .10e\n",
                ti_.getTimeStepNumber(),
                ti_.getCurrentTime(),
                globalForces[i].Force_pressure.x(), globalForces[i].Force_pressure.y(), globalForces[i].Force_pressure.z(),
                globalForces[i].Force_viscous.x(), globalForces[i].Force_viscous.y(), globalForces[i].Force_viscous.z(),
                globalForces[i].Force_penalty.x(), globalForces[i].Force_penalty.y(), globalForces[i].Force_penalty.z(),
                globalForces[i].mass_conservation,
                globalForces[i].Torque.x(), globalForces[i].Torque.y(), globalForces[i].Torque.z());
#endif
#if (DIM == 2)
            fprintf(fp,
                    "Timestep: %1d, Time: %.5f\n"
                    "Fx_pre = % .10e, Fy_pre = % .10e\n"
                    "Fx_vis = % .10e, Fy_vis = % .10e\n"
                    "Fx_pen = % .10e, Fy_pen = % .10e\n"
                    "mass_conservation = % .10e\n"
                    "Tz = % .10e\n",
                    ti_.getTimeStepNumber(),
                    ti_.getCurrentTime(),
                    globalForces[i].Force_pressure.x(), globalForces[i].Force_pressure.y(),
                    globalForces[i].Force_viscous.x(), globalForces[i].Force_viscous.y(),
                    globalForces[i].Force_penalty.x(), globalForces[i].Force_penalty.y(),
                    globalForces[i].mass_conservation,
                    globalForces[i].Torque.z());
#endif
            fclose(fp);
        }
    }
}

void BoundaryCalc_ExtraInter::WriteTrueGPToFile() {
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = DomainPtLocal.size();
    std::vector<int> eachProcData(nProc);
    MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> disp(nProc, 0);
    for (int i = 1; i < disp.size(); i++)
    {
        disp[i] = disp[i - 1] + eachProcData[i - 1];
    }

    int totalProcData = 0;
    for (int i = 0; i < nProc; i++)
    {
        totalProcData += eachProcData[i];
    }

    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        DomainPtGlobal.resize(totalProcData);
        GaussPtGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(DomainPtLocal.data(), DomainPtLocal.size(), ZEROPTVtype, DomainPtGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(GaussPtLocal.data(), GaussPtLocal.size(), ZEROPTVtype, GaussPtGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    std::string fPrefix = "trueGP_Step_" + std::to_string(ti_.getTimeStepNumber())+".vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "Do_x0,Do_y0,GP_x0,GP_y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "Do_x0,Do_y0,Do_z0,GP_x0,GP_y0,GP_z0\n");
#endif

        for (int i = 0; i < DomainPtGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e\n",
                    DomainPtGlobal[i](0), DomainPtGlobal[i](1),GaussPtGlobal[i](0), GaussPtGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                    DomainPtGlobal[i](0), DomainPtGlobal[i](1), DomainPtGlobal[i](2)
                    ,GaussPtGlobal[i](0), GaussPtGlobal[i](1) , GaussPtGlobal[i](2));
#endif
        }

        fclose(fp);
    }

}

void BoundaryCalc_ExtraInter::WriteGPValue() {

    for (int geomID = 0;geomID < imga_->getGeometries().size();geomID++) {
        int nProc = TALYFEMLIB::GetMPISize(); // be careful
        int numNodes = GPValuesWithDist_Geo[geomID].size();
        std::vector<int> eachProcData(nProc);
        MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> disp(nProc, 0);
        for (int i = 1; i < disp.size(); i++) {
            disp[i] = disp[i - 1] + eachProcData[i - 1];
        }

        int totalProcData = 0;
        for (int i = 0; i < nProc; i++) {
            totalProcData += eachProcData[i];
        }

        if (TALYFEMLIB::GetMPIRank() == 0) {
            GlobalGPValuesWithDist_Geo[geomID].resize(totalProcData);
        }


        MPI_Datatype ZEROPTVplus2type;
        MPI_Type_contiguous(5, MPI_DOUBLE, &ZEROPTVplus2type);
        MPI_Type_commit(&ZEROPTVplus2type);
        MPI_Gatherv(GPValuesWithDist_Geo[geomID].data(), GPValuesWithDist_Geo[geomID].size(), ZEROPTVplus2type, GlobalGPValuesWithDist_Geo[geomID].data(), eachProcData.data(),
                    disp.data(), ZEROPTVplus2type, 0, MPI_COMM_WORLD);

        if (TALYFEMLIB::GetMPIRank() == 0) {
            std::string fname = (v_.ndof == NSNodeData::NS_DOF) ?
                                "Local_Pressure_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat" :
                                (false) ? // TODO: separate temperature and FLux
                                "Local_Flux_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat" :
                                "Local_Temperature_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat";
            FILE *fp = fopen(fname.c_str(), "a");
            if (!fp) {
                throw TALYException() << "Cannot create file: " << fname;
            }

            for (int i = 0; i < GlobalGPValuesWithDist_Geo[geomID].size(); i++) {
#if (DIM == 3)
                fprintf(fp,
                    "x = % .10e, y = % .10e, z = %.10e, value = %.10e\n",
                    GlobalGPValuesWithDist_Geo[geomID][i].position[0],GlobalGPValuesWithDist_Geo[geomID][i].position[1]
                    ,GlobalGPValuesWithDist_Geo[geomID][i].position[2],GlobalGPValuesWithDist_Geo[geomID][i].value);
#endif
#if (DIM == 2)
                fprintf(fp,
                        "x = % .10e, y = % .10e, value = %.10e\n",
                        GlobalGPValuesWithDist_Geo[geomID][i].position[0], GlobalGPValuesWithDist_Geo[geomID][i].position[1], GlobalGPValuesWithDist_Geo[geomID][i].value);
#endif
            }
            fclose(fp);


            std::string fname2 = "Dist2SurrogateGP_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat";
            FILE *fp2 = fopen(fname2.c_str(), "a");
            if (!fp2) {
                throw TALYException() << "Cannot create file: " << fname2;
            }

            for (int i = 0; i < GlobalGPValuesWithDist_Geo[geomID].size(); i++) {
#if (DIM == 3)
                fprintf(fp2,
                    "x = % .10e, y = % .10e, z = %.10e, value = %.10e\n",
                    GlobalGPValuesWithDist_Geo[geomID][i].position[0],GlobalGPValuesWithDist_Geo[geomID][i].position[1]
                    ,GlobalGPValuesWithDist_Geo[geomID][i].position[2],GlobalGPValuesWithDist_Geo[geomID][i].distance);
#endif
#if (DIM == 2)
                fprintf(fp2,
                        "x = % .10e, y = % .10e, value = %.10e\n",
                        GlobalGPValuesWithDist_Geo[geomID][i].position[0], GlobalGPValuesWithDist_Geo[geomID][i].position[1], GlobalGPValuesWithDist_Geo[geomID][i].distance);
#endif
            }
            fclose(fp2);
        }
    }

}

void BoundaryCalc_ExtraInter::WriteGPBoundaryError() {

    for (int geomID = 0;geomID < imga_->getGeometries().size();geomID++) {
        int nProc = TALYFEMLIB::GetMPISize(); // be careful
        int numNodes = GPBoundaryErrorWithDist_Geo[geomID].size();
        std::vector<int> eachProcData(nProc);
        MPI_Allgather(&numNodes, 1, MPI_INT, eachProcData.data(), 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> disp(nProc, 0);
        for (int i = 1; i < disp.size(); i++) {
            disp[i] = disp[i - 1] + eachProcData[i - 1];
        }

        int totalProcData = 0;
        for (int i = 0; i < nProc; i++) {
            totalProcData += eachProcData[i];
        }

        if (TALYFEMLIB::GetMPIRank() == 0) {
            GlobalGPBoundaryErrorWithDist_Geo[geomID].resize(totalProcData);
        }


        MPI_Datatype DouobleZEROPTVtype;
        MPI_Type_contiguous(6 + 1, MPI_DOUBLE, &DouobleZEROPTVtype);
        MPI_Type_commit(&DouobleZEROPTVtype);
        MPI_Gatherv(GPBoundaryErrorWithDist_Geo[geomID].data(), GPBoundaryErrorWithDist_Geo[geomID].size(), DouobleZEROPTVtype, GlobalGPBoundaryErrorWithDist_Geo[geomID].data(), eachProcData.data(),
                    disp.data(), DouobleZEROPTVtype, 0, MPI_COMM_WORLD);

        if (TALYFEMLIB::GetMPIRank() == 0) {
            std::string fname = "Local_NSValues_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat";
            FILE *fp = fopen(fname.c_str(), "a");
            if (!fp) {
                throw TALYException() << "Cannot create file: " << fname;
            }

            for (int i = 0; i < GlobalGPBoundaryErrorWithDist_Geo[geomID].size(); i++) {
#if (DIM == 3)
                fprintf(fp,
                    "x = % .10e, y = % .10e, z = %.10e, u = %.10e, v = %.10e, w = %.10e\n",
                        GlobalGPBoundaryErrorWithDist_Geo[geomID][i].position[0],GlobalGPBoundaryErrorWithDist_Geo[geomID][i].position[1]
                    ,GlobalGPBoundaryErrorWithDist_Geo[geomID][i].position[2]
                    , GlobalGPBoundaryErrorWithDist_Geo[geomID][i].ZeroPtV_values[0], GlobalGPBoundaryErrorWithDist_Geo[geomID][i].ZeroPtV_values[1],
                    GlobalGPBoundaryErrorWithDist_Geo[geomID][i].ZeroPtV_values[2]);
#endif
#if (DIM == 2)
                fprintf(fp,
                        "x = % .10e, y = % .10e, u = %.10e, v = %.10e, p = %.10e\n",
                        GlobalGPBoundaryErrorWithDist_Geo[geomID][i].position[0], GlobalGPBoundaryErrorWithDist_Geo[geomID][i].position[1], GlobalGPBoundaryErrorWithDist_Geo[geomID][i].ZeroPtV_values[0], GlobalGPBoundaryErrorWithDist_Geo[geomID][i].ZeroPtV_values[1],GlobalGPBoundaryErrorWithDist_Geo[geomID][i].value);
#endif
            }

            fclose(fp);


            std::string fname2 = "Dist2SurrogateGP_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat";
            FILE *fp2 = fopen(fname2.c_str(), "a");
            if (!fp2) {
                throw TALYException() << "Cannot create file: " << fname2;
            }

            for (int i = 0; i < GlobalGPValuesWithDist_Geo[geomID].size(); i++) {
#if (DIM == 3)
                fprintf(fp2,
                    "x = % .10e, y = % .10e, z = %.10e, value = %.10e\n",
                    GlobalGPValuesWithDist_Geo[geomID][i].position[0],GlobalGPValuesWithDist_Geo[geomID][i].position[1]
                    ,GlobalGPValuesWithDist_Geo[geomID][i].position[2],GlobalGPValuesWithDist_Geo[geomID][i].distance);
#endif
#if (DIM == 2)
                fprintf(fp2,
                        "x = % .10e, y = % .10e, value = %.10e\n",
                        GlobalGPValuesWithDist_Geo[geomID][i].position[0], GlobalGPValuesWithDist_Geo[geomID][i].position[1], GlobalGPValuesWithDist_Geo[geomID][i].distance);
#endif
            }
            fclose(fp2);
        }
    }

}


void BoundaryCalc_ExtraInter::calculateDiff(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values, const double Coe_diff,
                                            const ZEROPTV &normal, ZEROPTV &diff) {

    ZeroMatrix<double> du;
    du.redim(DIM, DIM);
    std::vector<double> vald(NSNodeData::NS_DOF * DIM);
    calcValueDerivativeFEM(fe, NSNodeData::NS_DOF, values, vald.data());

    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            du(i, j) = vald[i * DIM + j] + vald[j * DIM + i];
        }
    }

    for (int i = 0; i < DIM; i++) {
        diff(i) = 0;
        for (int j = 0; j < DIM; j++) {
            diff(i) += Coe_diff * du(i, j) * normal[j];
        }
    }
}


void BoundaryCalc_ExtraInter::calculatePressure(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values, double &p) {

    std::vector<double> u(NSNodeData::NS_DOF);
    calcValueFEM(fe, NSNodeData::NS_DOF, values, u.data());

    p = u[NSNodeData::PRESSURE];
}

void BoundaryCalc_ExtraInter::calculateVelocity(const FEMElm &fe, const PetscScalar *values, ZEROPTV &vel) {
    std::vector<double> u(NSNodeData::NS_DOF);
    calcValueFEM(fe, NSNodeData::NS_DOF, values, u.data());

#if (DIM==3)
    vel = {u[NSNodeData::VEL_X],u[NSNodeData::VEL_Y],u[NSNodeData::VEL_Z]};
#endif
#if (DIM==2)
    vel = {u[NSNodeData::VEL_X],u[NSNodeData::VEL_Y],0.0};
#endif

}


DENDRITE_REAL BoundaryCalc_ExtraInter::ElementSize(const FEMElm &fe) {
    return pow((pow(2, DIM) * fe.jacc()), (double)1 / DIM);
}

double BoundaryCalc_ExtraInter::linear_extrapolation(ZEROPTV p1, ZEROPTV p2, double value1, double value2, ZEROPTV p3) {
    double slope = (value2 - value1) / p1.distanceTo(p2);
    double value3 = slope * p1.distanceTo(p3) + value1;
    return value3;
}

void BoundaryCalc_ExtraInter::NSBoundaryError(std::vector<double> &boundaryError)
{
    boundaryError.resize(DIM);
    for (int dim = 0; dim < DIM; dim++)
    {
        MPI_Reduce(&velocityWJ[dim], &boundaryError[dim], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
        boundaryError[dim] = sqrt(boundaryError[dim]);
    }
}

void BoundaryCalc_ExtraInter::WriteNSBoundaryError() {

    std::vector<double> boundaryError;

    NSBoundaryError(boundaryError);

    // print to screen
    for (int dim = 0;dim<DIM;dim++) {
        PrintStatus("NS boundary error:", dim," direction = ", boundaryError[dim]);
    }

    // write to file
    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        std::string fname = "NSBoundaryError.dat";
        FILE *fp = fopen(fname.c_str(), "a");
        if (!fp)
        {
            throw TALYException() << "Cannot create file: " << fname;
        }
#if (DIM == 3)
        fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Ex = % .10e, Ey = % .10e, Ez = % .10e\n",
                ti_.getTimeStepNumber(),
                ti_.getCurrentTime(),
                boundaryError[0], boundaryError[1], boundaryError[2]);
#endif
#if (DIM == 2)
        fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Ex = % .10e, Ey = % .10e\n",
                ti_.getTimeStepNumber(),
                ti_.getCurrentTime(),
                boundaryError[0], boundaryError[1]);
#endif
        fclose(fp);
    }

}

std::string BoundaryCalc_ExtraInter::queryToStringKey(const double query_pt[3]) {
    return std::to_string(query_pt[0]) + "," + std::to_string(query_pt[1]) + "," + std::to_string(query_pt[2]);
}

uint32_t BoundaryCalc_ExtraInter::getRetIndexForQuery(const double query_pt[3]) {
    std::string key = queryToStringKey(query_pt);

    // Check if the query is already in the map
    auto it = queryToIndexMap_.find(key);
    if (it == queryToIndexMap_.end()) {
//        std::cout << "do not find the pt\n";
        // If not found, perform the knnSearch and store the result
        size_t num_results = 1;
        std::vector<uint32_t> ret_index(num_results);
        std::vector<double> out_dist_sqr(num_results);

        num_results = kd_tree_trueGP_->knnSearch(
                &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        uint32_t resultIndex = ret_index[0];

        queryToIndexMap_[key] = resultIndex;
//        std::cout << "Result for query " << key << " computed and stored: " << resultIndex << std::endl;
        return resultIndex;
    } else {
//        std::cout << "find the pt\n";

        // If found, retrieve the result without recomputing
//        std::cout << "Retrieved saved result for query " << key << ": " << it->second << std::endl;
        return it->second;
    }
}

void BoundaryCalc_ExtraInter::gainQueryToIndexMap(std::unordered_map<std::string, uint32_t> &queryToIndexMap) {
    queryToIndexMap = queryToIndexMap_;
}

int BoundaryCalc_ExtraInter::NearestGeo(const TALYFEMLIB::FEMElm &fe) {
    const TALYFEMLIB::ZEROPTV pt = fe.position();
    double Min_Dist = std::numeric_limits<double>::max();
    int PickGeomID = 0;

    for (int CarvedOutGeomID = 0;
         CarvedOutGeomID < imga_->getGeometries().size();
         ++CarvedOutGeomID) {

//            std::cout << "CarvedOutGeomID = " << CarvedOutGeomID << "\n";

        double d[DIM];
#ifdef DEEPTRACE
        SBMCalcDeepTrace sbmCalc(fe, &idata_, imga_, kd_trees_, 0);
#else
        SBMCalc sbmCalc(fe,&idata_,imga_,kd_trees_,0);
#endif

//        SBMCalc sbmCalc(fe, &idata_, imga_, kd_trees_, CarvedOutGeomID);
//        sbmCalc.Dist2Geo(d);

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
//        std::cout << "PickGeomID = " << PickGeomID << "\n";

    return PickGeomID;
}


#endif //NSIBM_KT_BOUNDARYCALC_EXTRAINTER_H
