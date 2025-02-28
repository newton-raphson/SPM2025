//
// Created by Cheng-Hau.
//

#ifndef DENDRITEKT_CALCERROR_H
#define DENDRITEKT_CALCERROR_H
#include <Traversal/Traversal.h>
#include <math.h>

class CalcError : public Traversal
{

    DENDRITE_REAL L2error_;
    DENDRITE_REAL LInferror_ = -INFINITY;
    const SubDomain *subdomain_;
    std::vector<DENDRITE_REAL> error_;
    std::vector<DENDRITE_REAL>::iterator it;
    ZEROPTV shift_;
    const std::vector<ElementMarker> ElementMarker_;
    int localElemID = 0;


public:
    CalcError(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
              const DomainExtents &domain, const SubDomain *subDomain, LEInputData *idata,const std::vector<ElementMarker> &ElementMarker);
    void traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) override;

    void getL2error(double *error);

    const std::vector<DENDRITE_REAL> &getElementalError() const
    {
        return error_;
    }

private:
    const LEInputData *idata_;
};

CalcError::CalcError(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &v,
                     const DomainExtents &domain, const SubDomain *subDomain, LEInputData *idata,const std::vector<ElementMarker> &ElementMarker)
        : Traversal(octDA, treePart, v, domain), subdomain_(subDomain), idata_(idata),
          shift_((idata->ibm_geom_def.size() == 0) ? ZEROPTV(0, 0, 0)
                                                   : idata->ibm_geom_def[0].InitialDisplacement),ElementMarker_(ElementMarker)
{
    L2error_ = 0;
    error_.resize(treePart.size(), 0.0);
    it = error_.begin();
    this->traverse();
}

void CalcError::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values)
{
    const DENDRITE_UINT ndof = this->getNdof();
    double localError = 0;
    auto &geo_tmp = idata_->ibm_geom_def.at(0);
    DENDRITE_REAL val_c_,val_a;

    if (ElementMarker_[localElemID] == ElementMarker::OUT_ELEMENT) {
        fe.refill(0, 0);
    } else if (ElementMarker_[localElemID] == ElementMarker::INTERCEPTED_ELEMENT) {
        fe.refill(0, idata_->RelOrderInterceptedElement);
    }

    while (fe.next_itg_pt())
    {

        bool inside = false;
        int InsideNumber = 0;

        DENDRITE_REAL coords[DIM];
        for (DENDRITE_UINT dim = 0; dim < DIM; dim++)
        {
            coords[dim] = fe.position()[dim] - shift_[dim];
        }

        if (subdomain_->getGeometry().size() != 0)
        {
#if (DIM == 2)

            for (int i = 0; i < subdomain_->getGeometry().size(); i++)
            {
                if (subdomain_->getGeometry()[i]->getRetainSide() == RetainSide::IN)
                {
                    if (subdomain_->getGeometry()[i]->getMSH()->ifInside(coords))
                    {
                        InsideNumber++;
                    }
                } else {
                    if (!subdomain_->getGeometry()[i]->getMSH()->ifInside(coords))
                    {
                        InsideNumber++;
                    }
                }

            }
#endif

#if (DIM == 3)
            for (int i = 0; i < subdomain_->getGeometry().size(); i++)
            {
              if (subdomain_->getGeometry()[i]->getRetainSide() == RetainSide::IN)
              {
                  if (subdomain_->getGeometry()[i]->getSTL()->ifInside(coords))
                  {
                      InsideNumber++;
                  }
              } else {
                  if (!subdomain_->getGeometry()[i]->getSTL()->ifInside(coords))
                  {
                      InsideNumber++;
                  }
              }
            }
#endif
        }
        if (InsideNumber == subdomain_->getGeometry().size())
        {
            inside = true;
        }
        else
        {
            inside = false;
        }


        if (inside)
        {
            calcValueFEM(fe, ndof, values, &val_c_);

            switch (idata_->SbmGeo) {

                case LEInputData::SBMGeo::RING:
                {
                    for (DENDRITE_UINT dof = 0; dof < ndof; dof++)
                    {
                        const double x_mid = geo_tmp.InitialDisplacement.x();
                        const double y_mid = geo_tmp.InitialDisplacement.y();
                        const DENDRITE_REAL r = sqrt((fe.position().x() - x_mid) * (fe.position().x() - x_mid) + (fe.position().y() - x_mid) * (fe.position().y() - x_mid));
                        double cosx = (fe.position().x() - x_mid) / r;
                        double sinx = (fe.position().y() - y_mid) / r;

                        val_a = (idata_->CalcUxError)? -r * log(r) / 2 / log(2) * cosx: -r * log(r) / 2 / log(2) * sinx;

                        localError += (val_c_ - val_a) * (val_c_ - val_a) * fe.detJxW();

                        if (LInferror_ < fabs(val_c_ - val_a))
                        {
                            LInferror_ = fabs(val_c_ - val_a);
                        }
                    }
                    break;
                }

                case LEInputData::SBMGeo::CIRCLE:
                case LEInputData::SBMGeo::STAR:
                {
                    double x = fe.position().x();
                    double y = fe.position().y();

                    for (DENDRITE_UINT dof = 0; dof < ndof; dof++)
                    {
                        val_a = (idata_->CalcUxError)? sin(M_PI*x)*cos(M_PI*y)/10.0:  cos(M_PI*x)*sin(M_PI*y)/10.0;
                        localError += (val_c_ - val_a) * (val_c_ - val_a) * fe.detJxW();

                        if (LInferror_ < fabs(val_c_ - val_a))
                        {
                            LInferror_ = fabs(val_c_ - val_a);
                        }
                    }
                    break;
                }

//                case LEInputData::SBMGeo::CIRCLE:
//                {
//                    double x = fe.position().x();
//                    double y = fe.position().y();
//
//                    for (DENDRITE_UINT dof = 0; dof < ndof; dof++)
//                    {
//                        val_a = (idata_->CalcUxError)? -cos(M_PI*x)*sin(M_PI*y)/10.0:  sin(M_PI*x/7)*sin(M_PI*y/3)/10.0;
//                        localError += (val_c_ - val_a) * (val_c_ - val_a) * fe.detJxW();
//
//                        if (LInferror_ < fabs(val_c_ - val_a))
//                        {
//                            LInferror_ = fabs(val_c_ - val_a);
//                        }
//                    }
//                    break;
//                }

                case LEInputData::SBMGeo::ROTATE:
                {
                    double x = fe.position().x() - 0.5 - shift_[0];
                    double y = fe.position().y() - 0.5 - shift_[1];

                    double x_true = cos(M_PI / 4) * x - sin(M_PI / 4) * y + 0.5;
                    double y_true = sin(M_PI / 4) * x + cos(M_PI / 4) * y + 0.5;


                    for (DENDRITE_UINT dof = 0; dof < ndof; dof++)
                    {
                        val_a = (idata_->CalcUxError)? -cos(M_PI*x_true)*sin(M_PI*y_true)/10.0:  sin(M_PI*x_true/7)*sin(M_PI*y_true/3)/10.0;
                        localError += (val_c_ - val_a) * (val_c_ - val_a) * fe.detJxW();

                        if (LInferror_ < fabs(val_c_ - val_a))
                        {
                            LInferror_ = fabs(val_c_ - val_a);
                        }
                    }
                    break;
                }

                default:
                {
                    break;
                }
            }



        }
    }

    localElemID++;
    *it = sqrt(localError);
    it = std::next(it);
    L2error_ += localError;
}

void CalcError::getL2error(double *error)
{
    DENDRITE_REAL globalL2Error, globalLInfError;
    MPI_Reduce(&L2error_, &globalL2Error, this->getNdof(), MPI_DOUBLE, MPI_SUM, 0, this->m_octDA->getCommActive());
    MPI_Reduce(&LInferror_, &globalLInfError, this->getNdof(), MPI_DOUBLE, MPI_MAX, 0, this->m_octDA->getCommActive());
    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        for (int dof = 0; dof < this->getNdof(); dof++)
        {
            globalL2Error = sqrt(globalL2Error);
        }
    }
    error[0] = globalL2Error;
    error[1] = globalLInfError;
}

#endif // DENDRITEKT_CALCERROR_H
