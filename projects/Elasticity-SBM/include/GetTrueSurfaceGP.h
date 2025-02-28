//
// Created by chenghau on 10/31/22.
//

#ifndef LE_KT_GETTRUESURFACEGP_H
#define LE_KT_GETTRUESURFACEGP_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include <LEInputData.h>

class GetTrueSurfaceGP : public IMGATraversal
{
    LEInputData::SBMGeo sbmgeo_;
    LEInputData *idata_;
    Point<DIM> shift_;
    std::vector<ZEROPTV> surfaceGPposLocal;
    std::vector<ZEROPTV> surfaceGPposGlobal;
    std::vector<ZEROPTV> surfaceGPvalueLocal;
    std::vector<ZEROPTV> surfaceGPvalueGlobal;
    std::vector<ZEROPTV> surfaceGPErrorLocal;
    std::vector<ZEROPTV> surfaceGPErrorGlobal;

    int iter;

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
    GetTrueSurfaceGP(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                     LEInputData *idata);

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;

    /**
     * @brief print True GPs to file
     */
    void WriteTrueGPToFile();
};

GetTrueSurfaceGP::GetTrueSurfaceGP(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                                   const DomainExtents &domain, LEInputData *idata)
        : IMGATraversal(octDA, imga, treePart, v, domain), sbmgeo_(idata->SbmGeo), shift_(imga->getGeometries()[0]->getTranslations()[0]),idata_(idata)
{
    MPI_Barrier(MPI_COMM_WORLD);
    this->imgaTraverse();
}

void GetTrueSurfaceGP::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                              const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values)
{
    DENDRITE_REAL *val_c;
    DENDRITE_REAL val_a;
    ZEROPTV PTV_calculated, PTV_error;
    const DENDRITE_UINT ndof = this->getNdof();

    calcValueFEM(fe, ndof, values, val_c);

    switch (sbmgeo_) {
        case LEInputData::SBMGeo::CIRCLE:
        {
            for (DENDRITE_UINT dof = 0; dof < ndof; dof++)
            {
                double x_true = fe.position().x();
                double y_true = fe.position().y();
                val_a = (dof == 0)? -cos(M_PI*x_true)*sin(M_PI*y_true)/10.0: sin(M_PI*x_true/7)*sin(M_PI*y_true/3)/10.0;
                PTV_calculated(dof) = val_c[dof];
                PTV_error(dof) = fabs(val_c[dof]-val_a);
            }
            break;
        }

        default:
        {
            break;
        }
    }

    surfaceGPvalueLocal.push_back(PTV_calculated);
    surfaceGPErrorLocal.push_back(PTV_error);
    surfaceGPposLocal.push_back(fe.position());
}

void GetTrueSurfaceGP::WriteTrueGPToFile() {
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = surfaceGPposLocal.size();
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
        surfaceGPposGlobal.resize(totalProcData);
        surfaceGPvalueGlobal.resize(totalProcData);
        surfaceGPErrorGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(surfaceGPposLocal.data(), surfaceGPposLocal.size(), ZEROPTVtype, surfaceGPposGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(surfaceGPvalueLocal.data(), surfaceGPvalueLocal.size(), ZEROPTVtype, surfaceGPvalueGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(surfaceGPErrorLocal.data(), surfaceGPErrorLocal.size(), ZEROPTVtype, surfaceGPErrorGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    std::string fPrefix = "trueGP.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "x0,y0,Ux,Uy,Error_x,Error_y\n");
#endif

#if (DIM == 3)
        fprintf(fp, "x0,y0,z0,Ux,Uy,Uz,Error_x,Error_y,Error_z\n");
#endif

        for (int i = 0; i < surfaceGPposGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                    surfaceGPposGlobal[i](0), surfaceGPposGlobal[i](1), surfaceGPvalueGlobal[i](0),surfaceGPvalueGlobal[i](1),surfaceGPErrorGlobal[i](0),surfaceGPErrorGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                    surfaceGPposGlobal[i](0), surfaceGPposGlobal[i](1), surfaceGPposGlobal[i](2), surfaceGPvalueGlobal[i](0),surfaceGPvalueGlobal[i](1),surfaceGPvalueGlobal[i](2)
                    ,surfaceGPErrorGlobal[i](0),surfaceGPErrorGlobal[i](1),surfaceGPErrorGlobal[i](2));
#endif
        }

        fclose(fp);
    }

}


#endif //LE_KT_GETTRUESURFACEGP_H
