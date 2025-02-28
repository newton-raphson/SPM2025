//
// Created by chenghau on 9/14/23.
//

#ifndef LE_KT_CALCSTRESS_H
#define LE_KT_CALCSTRESS_H

#include <Traversal/Traversal.h>
#include <LEInputData.h>

class CalcStress : public Traversal
{

    const SubDomain *subdomain_;
    LEInputData *idata_;
    std::vector<DENDRITE_REAL> stress_;
    std::vector<DENDRITE_REAL>::iterator it;

    std::vector<GEOMETRY::MSH *> mshs_;
    std::vector<GEOMETRY::STL *> stls_;
    ZEROPTV shift_;
    const IMGA *imga_;
    std::vector<ZEROPTV> GPpos;
    std::vector<ZEROPTV> GPposGlobal;

    /// for setting diff the False Intercepted Element
    int localElemID = 0;

    const TimeInfo ti_;

public:
    CalcStress(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
               const DomainExtents &domain, const SubDomain *subDomain, LEInputData *idata, const TimeInfo ti);

    /**
     * @brief Overriding the traversal class operation
     * @param fe the element we use to access gauss points
     * @param values nodal value on the cell
     */
    void traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) override;

    /**
     * @brief return L2 residual
     * @param residual L2 residual
     */
    // TODO: implementation
    //    void getResidual(ZEROPTV *residual);

    /**
     * @brief return the stress on element cells
     * @return residual vector on element cells
     */
    void getElementalstress(std::vector<DENDRITE_REAL> &stress)
    {
        stress = stress_;
    }

    void WritePostProcessVolumeGPToFile();
};

CalcStress::CalcStress(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                       const SubDomain *subDomain, LEInputData *idata, const TimeInfo ti)
    : Traversal(octDA, treePart, v, domain), subdomain_(subDomain), idata_(idata), ti_(ti)
{
    stress_.resize((6 * (DIM - 1) + 1) * treePart.size(), 0.0);
    it = stress_.begin();

    this->traverse();
}

void CalcStress::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values)
{
    const DENDRITE_UINT ndof = this->getNdof();
    double VonMisesStress = 0.0;

    fe.refill(0, 0);

    std::vector<double> StrainVector(3 * (DIM - 1));
#ifndef NDEBUG
//    std::vector<std::vector<double>> StrainVector_eachGP(fe.nbf());
#endif

#ifndef NDEBUG
//    int GPID = 0;
#endif

    while (fe.next_itg_pt())
    {
#ifndef NDEBUG
//        StrainVector_eachGP[GPID].resize(3 * (DIM - 1));
#endif

        DENDRITE_REAL duidj[LENodeData::LE_DOF * DIM];
        calcValueDerivativeFEM(fe, ndof, values, duidj);

        for (int i = 0; i < DIM; i++)
        {
            StrainVector[i] += duidj[i * DIM + i] / fe.nbf();
#ifndef NDEBUG
//            StrainVector_eachGP[GPID][i] = duidj[i * DIM + i] / fe.nbf();
#endif
        }

        StrainVector[DIM] += (duidj[0 * DIM + 1] + duidj[1 * DIM + 0]) / fe.nbf();
#ifndef NDEBUG
//        StrainVector_eachGP[GPID][DIM] = duidj[0 * DIM + 1] + duidj[1 * DIM + 0];
#endif

#if (DIM == 3)
        StrainVector[4] += (duidj[0 * DIM + 2] + duidj[2 * DIM + 0]) / fe.nbf();
        StrainVector[5] += (duidj[1 * DIM + 2] + duidj[2 * DIM + 1]) / fe.nbf();
#ifndef NDEBUG
//        StrainVector_eachGP[GPID][4] = duidj[0 * DIM + 2] + duidj[2 * DIM + 0];
//        StrainVector_eachGP[GPID][5] = duidj[1 * DIM + 2] + duidj[2 * DIM + 1];
#endif
#endif

#ifndef NDEBUG
//        std::cout << "GPID = " << GPID << ", strain = ";
//        for (int i = 0; i < 3 * (DIM - 1); i++)
//        {
//            std::cout << StrainVector_eachGP[GPID][i] << " ";
//        }
//        std::cout << "\n";
//        GPID++;
#endif
    }
#ifndef NDEBUG
//    std::cout << " ------------------------- \n";
#endif

    std::vector<double> StressVector(3 * (DIM - 1));

    for (int row = 0; row < StrainVector.size(); row++)
    {
        for (int col = 0; col < StrainVector.size(); col++)
        {
            StressVector[row] += idata_->Cmatrix[row][col] * StrainVector[col];
        }
    }

    // von-mises stress
#if (DIM == 3)
    VonMisesStress = sqrt((pow(StressVector[0] - StressVector[1], 2) + pow(StressVector[1] - StressVector[2], 2) + pow(StressVector[2] - StressVector[0], 2) + 6 * (pow(StressVector[3], 2) + pow(StressVector[4], 2) + pow(StressVector[5], 2))) / 2);
#endif

#if (DIM == 2)
    // TODO: please make sure it is correct!
    VonMisesStress = sqrt(pow(StressVector[0], 2) + pow(StressVector[1], 2) - StressVector[0] * StressVector[1] + 3 * pow(StressVector[2], 2));
#endif

    for (int i = 0; i < 3 * (DIM - 1); i++)
    {
        *it = StrainVector[i];
        it = std::next(it);
    }

    for (int i = 0; i < 3 * (DIM - 1); i++)
    {
        *it = StressVector[i];
        it = std::next(it);
    }

    *it = VonMisesStress;
    it = std::next(it); // do forget this (previous mistake I made)
    //    std::cout << "*it = " << *it << "\n";
}

void CalcStress::WritePostProcessVolumeGPToFile()
{
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = GPpos.size();
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
        GPposGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(GPpos.data(), GPpos.size(), ZEROPTVtype, GPposGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    std::string fPrefix = "PostProcessvolumeGP.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "x0,y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "x0,y0,z0\n");
#endif

        for (int i = 0; i < GPposGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e\n",
                    GPposGlobal[i](0), GPposGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e\n",
                    GPposGlobal[i](0), GPposGlobal[i](1), GPposGlobal[i](2));
#endif
        }

        fclose(fp);
    }
}

#endif // LE_KT_CALCSTRESS_H
