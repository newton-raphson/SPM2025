//
// Created by chenghau.
//

#ifndef DENDRITEKT_IMGALOOP_H
#define DENDRITEKT_IMGALOOP_H

#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include <LEInputData.h>

class IMGALoop : public IMGATraversal
{
    LEInputData *idata_;
    Point<DIM> shift_;
    double valueWJ = 0.0;
    double totalGP = 0.0;
    double x_true, y_true;
    std::vector<NodeAndValues<DENDRITE_REAL>> gpInfo_;
    std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoIt_;
    std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoEnd_;

    DENDRITE_UINT lclElemID = 0;

    int iter;

public:
    IMGALoop(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
             LEInputData *idata);

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;
    void computeBoundaryError(DENDRITE_REAL *boundaryValue);
};

IMGALoop::IMGALoop(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                   const DomainExtents &domain, LEInputData *idata)
        : IMGATraversal(octDA, imga, treePart, v, domain), idata_(idata), shift_(imga->getGeometries()[0]->getTranslations()[0])
{
    gpInfo_ = imga->getSurfaceGaussPoints();
    gpInfoIt_ = gpInfo_.cbegin();
    gpInfoEnd_ = gpInfo_.cend();
    MPI_Barrier(MPI_COMM_WORLD);
    this->imgaTraverse();
}

void IMGALoop::imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                      const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values)
{

    DENDRITE_REAL val_c, val_a;

    x_true = fe.position().x();
    y_true = fe.position().y();

    switch (idata_->SbmGeo)
    {

        case LEInputData::SBMGeo::RING:
        {
            auto &geo_tmp = idata_->ibm_geom_def.at(0);
            double r = sqrt(pow(x_true - geo_tmp.InitialDisplacement.x(), 2) + pow(y_true - geo_tmp.InitialDisplacement.y(), 2));

            double cosx = (x_true - geo_tmp.InitialDisplacement.x()) / r;
            double sinx = (y_true - geo_tmp.InitialDisplacement.y()) / r;

            val_a = (idata_->CalcUxError)? -r * log(r) / 2 / log(2) * cosx: -r * log(r) / 2 / log(2) * sinx;

            break;
        }

        case LEInputData::SBMGeo::STAR:
        {
            val_a = (idata_->CalcUxError)? sin(M_PI*x_true)*cos(M_PI*y_true)/10.0: cos(M_PI*x_true)*sin(M_PI*y_true)/10.0;

            break;

        }
        case LEInputData::SBMGeo::CIRCLE:
        {

            val_a = (idata_->CalcUxError)? -cos(M_PI*x_true)*sin(M_PI*y_true)/10.0: sin(M_PI*x_true/7)*sin(M_PI*y_true/3)/10.0;

            break;
        }


        default:
        {
            break;
        }
    }
    calcValueFEM(fe, 1, values, &val_c);

    valueWJ += (val_c - val_a) * (val_c - val_a) * fe.detJxW();
}

void IMGALoop::computeBoundaryError(DENDRITE_REAL *boundaryValue)
{
    MPI_Reduce(&valueWJ, &boundaryValue[0], 1, MPI_DOUBLE, MPI_SUM, 0, m_octDA->getCommActive());
    boundaryValue[0] = sqrt(boundaryValue[0]);
}
#endif // DENDRITEKT_IMGALOOP_H
