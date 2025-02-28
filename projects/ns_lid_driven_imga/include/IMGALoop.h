//
// Created by chenghau on 2/7/24.
//

#ifndef DENDRITEKT_IMGALOOP_H
#define DENDRITEKT_IMGALOOP_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include "NSInputData.h"
#include "NSNodeData.h"
#include "SBMCalc.h"
#include "TimeInfo.h"

using namespace TALYFEMLIB;


struct GPValues
{
    ZEROPTV position;
    DENDRITE_REAL values[DIM+1];
};


/// Force info
struct ForceInfo
{
    ZEROPTV Force_pressure;
    ZEROPTV Force_viscous;
    ZEROPTV Force_penalty;
    DENDRITE_REAL mass_conservation = 0.0;
    ZEROPTV Torque;
    double GradUDotN_PlusP = 0.0;
    const ZEROPTV Force_all(bool includePenalty = true)
    {
        ZEROPTV force_all;
        if (includePenalty)
        {
            force_all = Force_pressure + Force_viscous + Force_penalty;
        }
        else
        {
            force_all = Force_pressure + Force_viscous;
        }
        return force_all;
    }

    /**
     * Adds rhs to *this.
     * @param rhs other vector to add to this one
     * @returns *this
     */
    inline ForceInfo &operator+=(const ForceInfo &rhs)
    {
        this->Force_pressure += rhs.Force_pressure;
        this->Force_viscous += rhs.Force_viscous;
        this->Force_penalty += rhs.Force_penalty;
        this->mass_conservation += rhs.mass_conservation;
        this->Torque += rhs.Torque;
        return *this;
    }
};

class IMGALoop : public IMGATraversal {

    NSInputData &idata_;
    const TimeInfo &ti_;
    const IMGA *imga_;
    std::vector<double> valueWJ;

    const VecInfo v_;

    std::vector<std::vector<GPValues>> gpValues_Geo;
    std::vector<std::vector<GPValues>> GlobalgpValues_Geo;


    std::vector<ForceInfo> localforce_infos;

public:


    IMGALoop(DA *octDA,
             const IMGA *imga,
             const std::vector<TREENODE> &treePart,
             const VecInfo &v,
             const DomainExtents &domain,
             NSInputData &idata_,
             const TimeInfo &ti_);

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;

    void NSBoundaryError(std::vector<double> &boundaryError);

    void WriteGPValue();

    void WriteTrueGPToFile();

    void computeForce(std::vector<ForceInfo> &globalForces);
};

IMGALoop::IMGALoop(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                   const DomainExtents &domain, NSInputData &idata, const TimeInfo &ti)
        : IMGATraversal(octDA, imga, treePart, v, domain), idata_(idata), ti_(ti), imga_(imga),v_(v)
{
    valueWJ.resize(DIM);
    gpValues_Geo.resize(imga->getGeometries().size());
    GlobalgpValues_Geo.resize(imga->getGeometries().size());
    localforce_infos.resize(imga->getGeometries().size());
    this->imgaTraverse();
}

void IMGALoop::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                      const ZEROPTV &h, const PetscScalar *values) {

    DENDRITE_REAL u[NSNodeData::NS_DOF], val_a[NSNodeData::NS_DOF];
    calcValueFEM(fe, NSNodeData::NS_DOF, values, u);

    const double detSideJxW = fe.detJxW();


    unsigned int geomID = gaussPoint.geomID;

    GPValues gpValues_local;
    for (int dim = 0; dim < DIM; dim++) {
        gpValues_local.position = fe.position();
    }

    for (int dof = 0; dof < DIM +1 ; dof++) {
        gpValues_local.values[dof] = u[dof];
    }

    gpValues_Geo[geomID].push_back(gpValues_local);



    /// TODO: multiple geometries
    auto carved_out_def = idata_.carved_out_geoms_def.at(0 /*We hard code here, it should be GeomID*/);
#if (DIM == 2)
    std::vector<double> bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                                    carved_out_def.getBC(NSNodeData::VEL_Y)[0]};

#elif (DIM == 3)
    std::vector<double> bc_value = {carved_out_def.getBC(NSNodeData::VEL_X)[0],
                                            carved_out_def.getBC(NSNodeData::VEL_Y)[0],
                                            carved_out_def.getBC(NSNodeData::VEL_Z)[0]};
#endif

    for (int dim = 0; dim < DIM; dim++) {
        val_a[dim] = 0.0; // initialization
        valueWJ[dim] += (u[dim] - val_a[dim]) * (u[dim] - val_a[dim]) * detSideJxW;
    }

    const int nsd = DIM;
    const double Coe_diff = 1.0/idata_.Re;
    ZEROPTV normal{gaussPoint.normal[0], gaussPoint.normal[1], 0};
#if (DIM == 3)
    normal.z() = gaussPoint.normal[2];
#endif

    double p = u[NSNodeData::PRESSURE];


    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    DENDRITE_REAL vald[NSNodeData::NS_DOF * DIM];
    calcValueDerivativeFEM(fe, NSNodeData::NS_DOF, values, vald);
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


    /// Force of object on the fluids, because the normals are pointed INSIDE (in the loading process)
    /// This needs extra case when you are calculating the kinematic info of the object
    ZEROPTV Force_pressure;
    ZEROPTV Force_viscous;
    ZEROPTV Force_penalty;
    for (int dim = 0; dim < nsd; dim++) {
        Force_pressure(dim) = (p * normal(dim)) * detSideJxW;
        Force_viscous(dim) = -(diff(dim)) * detSideJxW;

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

    // this is the same for 2D and 3D
//            for (int dim = 0; dim < 3; dim++) {
//              localforce_infos[geomID].Torque(dim) +=
//                  (Force_pressure + Force_viscous + Force_penalty)((dim + 2)%3)
//                      *(fe.position()((dim + 1)%3) - kI.center((dim + 1)%3))
//                      - (Force_pressure + Force_viscous + Force_penalty)((dim + 1)%3)
//                          *(fe.position()((dim + 2)%3) - kI.center((dim + 2)%3));
//            }
#if (DIM == 2)
    localforce_infos[geomID].Torque.x() = 0.0;
    localforce_infos[geomID].Torque.y() = 0.0;
#endif
#if (DIM == 3)
    localforce_infos[geomID].Torque.x() = 0.0;
    localforce_infos[geomID].Torque.y() = 0.0;
    localforce_infos[geomID].Torque.z() = 0.0;
#endif

}

void IMGALoop::NSBoundaryError(std::vector<double> &boundaryError)
{
    boundaryError.resize(DIM);
    for (int dim = 0; dim < DIM; dim++)
    {
        MPI_Reduce(&valueWJ[dim], &boundaryError[dim], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
        boundaryError[dim] = sqrt(boundaryError[dim]);
    }
}

void IMGALoop::computeForce(std::vector<ForceInfo> &globalForces)
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


void IMGALoop::WriteGPValue() {

    for (int geomID = 0;geomID < imga_->getGeometries().size();geomID++) {
        int nProc = TALYFEMLIB::GetMPISize(); // be careful
        int numNodes = gpValues_Geo[geomID].size();
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
            GlobalgpValues_Geo[geomID].resize(totalProcData);
        }


        MPI_Datatype ZEROPTVplusNSvalues;
        MPI_Type_contiguous(4+DIM, MPI_DOUBLE, &ZEROPTVplusNSvalues);
        MPI_Type_commit(&ZEROPTVplusNSvalues);
        MPI_Gatherv(gpValues_Geo[geomID].data(), gpValues_Geo[geomID].size(), ZEROPTVplusNSvalues, GlobalgpValues_Geo[geomID].data(), eachProcData.data(),
                    disp.data(), ZEROPTVplusNSvalues, 0, MPI_COMM_WORLD);

        if (TALYFEMLIB::GetMPIRank() == 0) {
            std::string fname = "Local_NSValues_" + std::to_string(geomID) + "_Step_" + std::to_string(ti_.getTimeStepNumber()) + ".dat";
             ;
            FILE *fp = fopen(fname.c_str(), "a");
            if (!fp) {
                throw TALYException() << "Cannot create file: " << fname;
            }

            for (int i = 0; i < GlobalgpValues_Geo[geomID].size(); i++) {
#if (DIM == 3)
                fprintf(fp,
                    "x = % .10e, y = % .10e, z = %.10e, u = %.10e, v = %.10e, w = %.10e, p = %.10e\n",
                    GlobalgpValues_Geo[geomID][i].position[0],GlobalgpValues_Geo[geomID][i].position[1]
                    ,GlobalgpValues_Geo[geomID][i].position[2]
                    , GlobalgpValues_Geo[geomID][i].values[0], GlobalgpValues_Geo[geomID][i].values[1], GlobalgpValues_Geo[geomID][i].values[2],
                    GlobalgpValues_Geo[geomID][i].values[3]);
#endif
#if (DIM == 2)
                fprintf(fp,
                        "x = % .10e, y = % .10e, u = %.10e, v = %.10e, p = %.10e\n",
                        GlobalgpValues_Geo[geomID][i].position[0], GlobalgpValues_Geo[geomID][i].position[1], GlobalgpValues_Geo[geomID][i].values[0], GlobalgpValues_Geo[geomID][i].values[1], GlobalgpValues_Geo[geomID][i].values[2]);
#endif
            }
            fclose(fp);

        }
    }

}


#endif //DENDRITEKT_IMGALOOP_H


