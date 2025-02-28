//
// Created by chenghau on 9/12/22.
//

#ifndef NSHTIBM_KT_GETSURFACEGP_H
#define NSHTIBM_KT_GETSURFACEGP_H
#include <Traversal/SurfaceLoop.h>
#include "util.h"
#include "SBMCalc.h"

class GetSurfaceGP : public SurfaceLoop
{

    LEInputData &inputData_;
    std::vector<ZEROPTV> surfaceGPlocal;
    std::vector<ZEROPTV> surfaceGPglobal;
    std::vector<double> Adjointconsistencylocal;
    std::vector<double> Consistencylocal;
    std::vector<double> Penaltylocal;
    std::vector<double> AdjointconsistencyGlocal;
    std::vector<double> ConsistencyGlocal;
    std::vector<double> PenaltyGlocal;
    std::vector<double> gradUdotnlocal;
    std::vector<double> gradUdotnGlocal;
    const IMGA *imga_;
    double LocalDistSquare = 0.0;
    double LocalDistForth = 0.0;
    double MaxDist = 0.0;
    int NumberOfGp = 0;

    ZEROPTV MaxDomain;
    ZEROPTV MinDomain;
    ZEROPTV GlobalMaxPos;
    ZEROPTV GlobalMinPos;

public:
    /**
     * @brief constructor: loop over the surface to obtain GPs on surrogate boundary
     * @param octDA
     * @param treePart tree partition
     * @param v vector information
     * @param domain domain information
     * @param subDomainBoundary boundary information
     * @param inputData input data
     */
    GetSurfaceGP(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                 SubDomainBoundary *subDomainBoundary, LEInputData &inputData, const IMGA *imga);
    /**
     * @brief override performSurfaceOperation
     * @param fe element
     * @param surfaceCoords
     * @param boundarySurface
     * @param values
     */
    void performSurfaceOperation(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                                 const BoundarySurface &boundarySurface, const PetscScalar *values) override;

    void getEuclideanDistance(double &EuclideanDistance);

    void getMaxDistance(double &MaxDistance);

    void getRMSDistance(double &RMSDistance);

    void getSpecialDistance(double &SpecialDistance);

    void GetRegion(const MPI_Comm &comm, TALYFEMLIB::ZEROPTV &MaxPos_);

    void GetMinMaxRegion(const MPI_Comm &comm, TALYFEMLIB::ZEROPTV &MinPos_, TALYFEMLIB::ZEROPTV &MaxPos_);

};

GetSurfaceGP::GetSurfaceGP(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                           SubDomainBoundary *subDomainBoundary, LEInputData &inputData,const IMGA *imga)
        : SurfaceLoop(octDA, treePart, v, domain, subDomainBoundary), inputData_(inputData),imga_(imga)
{

    for (int dim = 0; dim < DIM; dim++)
    {
        MaxDomain(dim) = -10000;
        MinDomain(dim) = 10000;
    }
    this->traverse();
}

void GetSurfaceGP::performSurfaceOperation(TALYFEMLIB::FEMElm &fe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
                                           const BoundarySurface &boundarySurface, const PetscScalar *values)
{

    while (fe.next_itg_pt())
    {

//        std::ofstream fout("testGP.txt", std::ios::app);
//#if (DIM==2)
//        fout << fe.position().x() << "," << fe.position().y() << "\n";
//#endif
//#if (DIM==3)
//        fout << fe.position().x() << "," << fe.position().y() << "," << fe.position().z() << "\n";
//#endif
//        fout.close();

        for (int dim = 0 ; dim <DIM ;dim ++) {
            if (MaxDomain(dim) < fe.position()(dim)) {
                MaxDomain(dim) = fe.position()(dim);
            }
            if (MinDomain(dim) > fe.position()(dim)){
                MinDomain(dim) = fe.position()(dim);
            }
        }

        const ZEROPTV pt = fe.position();

        double d[DIM];
        int geomID;

        SBMCalc sbmCalc(fe,&inputData_,imga_);
        sbmCalc.Dist2Geo(d,geomID);

        double DistMag = 0.0;
        for (int dim = 0; dim < DIM; dim++)
        {
            DistMag += pow(d[dim],2);
            LocalDistSquare += pow(d[dim],2);
            LocalDistForth += pow(d[dim],16);
        }

        NumberOfGp++;

        DistMag = sqrt(DistMag);

        if (DistMag > MaxDist)
        {
            MaxDist = DistMag;
        }
    }

}

void GetSurfaceGP::getEuclideanDistance(double &EuclideanDistance)
{
    MPI_Reduce(&LocalDistSquare,&EuclideanDistance,1,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
}

void GetSurfaceGP::getMaxDistance(double &MaxDistance)
{
    MPI_Reduce(&MaxDist,&MaxDistance,1,MPI_DOUBLE,MPI_MAX,0,this->m_octDA->getCommActive());
}

void GetSurfaceGP::getRMSDistance(double &RMSDistance)
{
    double EuclideanDistance;
    MPI_Reduce(&LocalDistSquare,&EuclideanDistance,1,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
    int NumberOfGpSum;
    MPI_Reduce(&NumberOfGp,&NumberOfGpSum,1,MPI_INT,MPI_SUM,0,this->m_octDA->getCommActive());

    RMSDistance = sqrt(EuclideanDistance/NumberOfGpSum);
}


void GetSurfaceGP::getSpecialDistance(double &SpecialDistance) {
    double SpDist;
    MPI_Reduce(&LocalDistForth,&SpDist,1,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
    int NumberOfGpSum;
    MPI_Reduce(&NumberOfGp,&NumberOfGpSum,1,MPI_INT,MPI_SUM,0,this->m_octDA->getCommActive());

    SpecialDistance = sqrt(sqrt(sqrt(sqrt(SpDist/NumberOfGpSum))));
}

void GetSurfaceGP::GetRegion(const MPI_Comm &comm, TALYFEMLIB::ZEROPTV &MaxPos_)
{
    for (int dim = 0; dim < DIM; dim++)
    {
        MPI_Allreduce(&MaxDomain(dim), &GlobalMaxPos(dim), 1, MPI_DOUBLE, MPI_MAX, comm);
    }

    MaxPos_ = GlobalMaxPos;
}

void GetSurfaceGP::GetMinMaxRegion(const MPI_Comm &comm, TALYFEMLIB::ZEROPTV &MinPos_, TALYFEMLIB::ZEROPTV &MaxPos_)
{
    for (int dim = 0; dim < DIM; dim++)
    {
        MPI_Allreduce(&MinDomain(dim), &GlobalMinPos(dim), 1, MPI_DOUBLE, MPI_MIN, comm);
        MPI_Allreduce(&MaxDomain(dim), &GlobalMaxPos(dim), 1, MPI_DOUBLE, MPI_MAX, comm);
    }

    MinPos_ = GlobalMinPos;
    MaxPos_ = GlobalMaxPos;
}


#endif // NSHTIBM_KT_GETSURFACEGP_H
