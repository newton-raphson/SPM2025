//
// Created by chenghau on 3/2/24.
//

#ifndef DENDRITEKT_EXTRAPOLATION_H
#define DENDRITEKT_EXTRAPOLATION_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include "NSNodeData.h"
#include "NSInputData.h"

class Extrapolation : public IMGATraversal
{
    Point<DIM> shift_;
    std::vector<double> velocityWJ;

    double x_true, y_true;

    std::vector<ZEROPTV> DomainPtLocal;
    std::vector<ZEROPTV> DomainPtGlobal;
    std::vector<ZEROPTV> GaussPtLocal;
    std::vector<ZEROPTV> GaussPtGlobal;
    std::vector<ZEROPTV> NormalLocal;
    std::vector<ZEROPTV> NormalGlobal;

    const IMGA *imga_;
    std::vector<std::vector<GpZEROPTV>> GPBoundaryErrorWithDist_Geo;
    std::vector<std::vector<GpZEROPTV>> GlobalGPBoundaryErrorWithDist_Geo;

    const TimeInfo ti_;

    NSInputData *idata_;

    std::vector<ForceInfo> localforce_infos;

    std::vector<ZEROPTV> TrueGPnormalGlobal_;

    const my_kd_tree_t *kd_tree_trueGP_;


    // Map to store the mapping from query to ret_index
    std::unordered_map<std::string, uint32_t> queryToIndexMap_;
    // Utility function to convert a query point to a string key
    std::string queryToStringKey(const double query_pt[3]);

    // Function to handle the retrieval or computation of ret_index for a given query point
    uint32_t getRetIndexForQuery(const double query_pt[3]);

    void computeForce(std::vector<ForceInfo> &globalForces);


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
    Extrapolation(DA *octDA, const IMGA *imga, const IMGA *imga1, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                  const TimeInfo &ti, NSInputData *idata,
                  std::vector<ZEROPTV> &TrueGPnormalGlobal,
                  const my_kd_tree_t *kdTree_TrueGP,
                  std::unordered_map<std::string, uint32_t> &queryToIndexMap);

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;

    void WriteTrueGPToFile();

    void WriteNSBoundaryError();

    void NSBoundaryError(std::vector<double> &boundaryError);

    void WriteGPBoundaryError();

    // Function to gain queryTOIndexMap_
    void gainQueryToIndexMap(std::unordered_map<std::string, uint32_t> &queryToIndexMap);

    void forceCalc();
};

Extrapolation::Extrapolation(DA *octDA, const IMGA *imga, const IMGA *imga1, const std::vector<TREENODE> &treePart, const VecInfo &v,
                             const DomainExtents &domain,const TimeInfo &ti, NSInputData *idata,
                             std::vector<ZEROPTV> &TrueGPnormalGlobal,
                             const my_kd_tree_t *kdTree_TrueGP,
                             std::unordered_map<std::string, uint32_t> &queryToIndexMap)
        : IMGATraversal(octDA, imga1, treePart, v, domain), shift_(imga->getGeometries()[0]->getTranslations()[0]),
        imga_(imga),ti_(ti),idata_(idata),
        queryToIndexMap_(queryToIndexMap), TrueGPnormalGlobal_(TrueGPnormalGlobal),kd_tree_trueGP_(kdTree_TrueGP)
{
    velocityWJ.resize(DIM);
    GPBoundaryErrorWithDist_Geo.resize(imga->getGeometries().size());
    GlobalGPBoundaryErrorWithDist_Geo.resize(imga->getGeometries().size());
    localforce_infos.resize(imga->getGeometries().size());
    this->imgaTraverse();
}

void Extrapolation::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                           const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values)
{
    static constexpr DENDRITE_REAL weight = (DIM==2) ? 2.0 : 3.0;
    double jacc_x_w=imga_->getSurfaceGaussPoints()[0].elemArea/weight;

    DENDRITE_REAL val_c_OnInterpolate[NSNodeData::NS_DOF], val_c[NSNodeData::NS_DOF], val_a[NSNodeData::NS_DOF],vald[NSNodeData::NS_DOF*DIM];

    calcValueFEM(fe, NSNodeData::NS_DOF, values, val_c_OnInterpolate);
    calcValueDerivativeFEM(fe,NSNodeData::NS_DOF,values,vald);


    DomainPtLocal.push_back(fe.position());
    TALYFEMLIB::ZEROPTV ptvn;
    std::memcpy(ptvn.data(), gaussPoint.normal, sizeof(double)*DIM);
    GaussPtLocal.push_back(ptvn);

    //std::cout<< "val_c = " << val_c <<"\n";
    //std::cout <<"gaussPoint.normal = " << gaussPoint.normal[0] << "," << gaussPoint.normal[1] << "\n";
    
    x_true = gaussPoint.normal[0];
    y_true = gaussPoint.normal[1];

    // TODO: Compute the error.
    // gaussPoint.normal = contains the position of point on surface. => important
    // gaussPoint.location = contains the position of point that we need to interpolate/

    // Use Taylor series: u(surface/unknown) = uo (val_c) + grad(u).d (valderivative)

    double GradUdotD[NSNodeData::NS_DOF];
    for (int i = 0; i < NSNodeData::NS_DOF; i++) {
        for (int j = 0; j < DIM; j++) {
            GradUdotD[i] = vald[i*DIM+j]*(gaussPoint.normal[j]-fe.position()[j]);
//            if (fabs(gaussPoint.normal[j]-gaussPoint.location[j])>0) {
//                std::cout << "gaussPoint.normal[j]-fe.position()[j] = " << gaussPoint.normal[j] - gaussPoint.location[j]
//                          << "\n";
//            }
        }
    }

    for (int i = 0; i < NSNodeData::NS_DOF; i++) {
        val_c[i] = val_c_OnInterpolate[i] + GradUdotD[i];
        val_a[i] = 0.0;
    }

    GpZEROPTV gpZeroptv;
    gpZeroptv.position = ZEROPTV{gaussPoint.normal[0], gaussPoint.normal[1], gaussPoint.normal[2]}; // position on surface (see above)

    for (int dim = 0; dim < DIM; dim++) {
        velocityWJ[dim] += (val_c[dim] - val_a[dim]) * (val_c[dim] - val_a[dim]) * jacc_x_w;
        gpZeroptv.ZeroPtV_values(dim) = val_c[dim] - val_a[dim];
    }

    gpZeroptv.value = val_c[NSNodeData::PRESSURE];

    const double p = val_c[NSNodeData::PRESSURE];

    GPBoundaryErrorWithDist_Geo[0].push_back(gpZeroptv);

    const double query_pt[3] = {gaussPoint.normal[0], gaussPoint.normal[1], gaussPoint.normal[2]}; // position on surface (see above)
    uint32_t savedIndex = getRetIndexForQuery(query_pt);

//    std::cout << "savedIndex = " << savedIndex << "\n";
    ZEROPTV normal = TrueGPnormalGlobal_[savedIndex];


    unsigned int geomID = (imga_->getGeometries().empty())? 0: gaussPoint.geomID;
    const int nsd = DIM;
    const double Coe_diff = 1.0/idata_->Re;

    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            du(i, j) = vald[i * DIM + j] + vald[j * DIM + i];
        }
    }

    ZEROPTV diff;
    for (int i = 0; i < nsd; i++) {
        diff(i) = 0;
        for (int j = 0; j < nsd; j++) {
            diff(i) += Coe_diff * du(i, j) * normal(j);
        }
    }

    // std::cout << "hb = " << hb << "\n";

    /// Force of object on the fluids, because the normals are pointed INSIDE (in the loading process)
    /// This needs extra case when you are calculating the kinematic info of the object
    ZEROPTV Force_pressure;
    ZEROPTV Force_viscous;
    ZEROPTV Force_penalty;

//    std::cout << "normal = "; normal.print(); std::cout << "\n";
//    std::cout << "p = " << p << "\n";
    for (int dim = 0; dim < nsd; dim++) {
        Force_pressure(dim) = (p * normal(dim)) * jacc_x_w;
        Force_viscous(dim) = -(diff(dim)) * jacc_x_w;

        // TODO: fix in future
        //    if (true) {
        //      if (idata_.ibm_geom_def[geomID].penalty_method == 0) {
        //        Force_penalty(dim) = Cb_f/hb*Coe_diff*(u[dim] - vel(dim))*detSideJxW * Coe_diff;
        //      } else if (idata_.ibm_geom_def[geomID].penalty_method == 1) {
        //        Force_penalty(dim) = Cb_f*Coe_diff*(u[dim] - vel(dim))*detSideJxW * Coe_diff;
        //      } else if (idata_.ibm_geom_def[geomID].penalty_method == 2) {
        //        Force_penalty(dim) = Cb_f/*/hb*Coe_diff*/*(u[dim] - vel(dim))*detSideJxW * Coe_diff;
        //      }
        //    }

        localforce_infos[geomID].mass_conservation += du(0, 0) + du(1, 1);
#if (DIM == 3)
        localforce_infos[geomID].mass_conservation += du(2, 2);
#endif
    }

    localforce_infos[geomID].Force_pressure += Force_pressure;
    localforce_infos[geomID].Force_viscous += Force_viscous;
    localforce_infos[geomID].Force_penalty += Force_penalty;


}

void Extrapolation::NSBoundaryError(std::vector<double> &boundaryError)
{
    boundaryError.resize(DIM);
    for (int dim = 0; dim < DIM; dim++)
    {
        MPI_Reduce(&velocityWJ[dim], &boundaryError[dim], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
        boundaryError[dim] = sqrt(boundaryError[dim]);
    }
}

void Extrapolation::WriteNSBoundaryError() {

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

void Extrapolation::WriteTrueGPToFile() {
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

    std::string fPrefix = "trueGP.vtk";
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

void Extrapolation::WriteGPBoundaryError() {

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
            std::string fname = "Local_NSValues_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + "_r" +std::to_string(idata_->RelativeOrderIntercepted) + ".dat";
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

        }

    }

}


std::string Extrapolation::queryToStringKey(const double query_pt[3]) {
    return std::to_string(query_pt[0]) + "," + std::to_string(query_pt[1]) + "," + std::to_string(query_pt[2]);
}

uint32_t Extrapolation::getRetIndexForQuery(const double query_pt[3]) {
    std::string key = queryToStringKey(query_pt);

    // Check if the query is already in the map
    auto it = queryToIndexMap_.find(key);
    if (it == queryToIndexMap_.end()) {
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
        // If found, retrieve the result without recomputing
//        std::cout << "Retrieved saved result for query " << key << ": " << it->second << std::endl;
        return it->second;
    }
}

void Extrapolation::gainQueryToIndexMap(std::unordered_map<std::string, uint32_t> &queryToIndexMap) {
    queryToIndexMap = queryToIndexMap_;
}


void Extrapolation::computeForce(std::vector<ForceInfo> &globalForces)
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

void Extrapolation::forceCalc()
{
    std::vector<ForceInfo> globalForces;
    globalForces.resize(imga_->getGeometries().size());
    computeForce(globalForces);

    for (int i = 0; i < imga_->getGeometries().size(); i++)
    {
        // print to screen
        PrintStatus(idata_->carved_out_geoms_def[i].name, " force: ", globalForces[i].Force_all(true));

        // write to file
        if (TALYFEMLIB::GetMPIRank() == 0)
        {
            std::string fname = "Force_" + std::to_string(i) + "_" + idata_->carved_out_geoms_def[i].name + ".dat";
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


#endif //DENDRITEKT_EXTRAPOLATION_H
