//
// Created by maksbh on 11/6/22.
//

#ifndef DENDRITEKT_ELEMENTMARKER_H
#define DENDRITEKT_ELEMENTMARKER_H
#include <Traversal/Traversal.h>
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
#include "NSInputData.h"


class SBMMarker : public Traversal {

    const double RatioGPSBM_;
    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    SubDomainBoundary * boundary_;
    const int RelativeOrder_;
    const NSInputData *idata_;

    std::vector<TALYFEMLIB::ZEROPTV> PosLocal;
    std::vector<TALYFEMLIB::ZEROPTV> PosGlobal;
    std::vector<int> InOutLocal;
    std::vector<int> InOutGlobal;

public:
    SBMMarker(DA * octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domain,
              std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
              const SubDomain * subDomain, const int &RelativeOrder, const NSInputData *idata);
    void traverseOperation(TALYFEMLIB::FEMElm & fe) override;

    void WriteDomainGPToFile();

};

SBMMarker::SBMMarker(DA * octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domain,
                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
                     const SubDomain * subDomain, const int &RelativeOrder, const NSInputData *idata)
        : Traversal(octDA, treePart, domain),eleMarkers_(eleMarkers),subdomain_(subDomain),RelativeOrder_(RelativeOrder),idata_(idata),RatioGPSBM_(idata->RatioGPSBM){

    std::bitset<ElementMarker::MAX_ELMENT_TYPE> init;
    init.reset();
    eleMarkers_.resize(octDA->getLocalElementSz(),init);
    this->traverse();

}

void SBMMarker::traverseOperation(TALYFEMLIB::FEMElm & fe) {
    /// out means ACTIVE here
    int outPts = 0;

    int npe = m_octDA->getNumNodesPerElement();
    const auto coords = this->m_coords;
    for(int i = 0; i < npe; i++){
        if (subdomain_->functionToRetainPhysNodes(&coords[DIM*i]) == ibm::Partition::OUT) {
            outPts++;
        }
    }

        if (idata_->ImmersedMethodType == NSInputData::ShiftedBoundary){

            if(outPts == npe){
                eleMarkers_[localElemID].set(ElementMarker::OUT_ELEMENT,true);
            }
            else if(outPts == 0){
                eleMarkers_[localElemID].set(ElementMarker::IN_ELEMENT, true);
            }
            else {

                eleMarkers_[localElemID].set(ElementMarker::INTERCEPTED_ELEMENT, true);

                /// SBM
                int ActiveGP = 0;
                fe.refill(0, RelativeOrder_);
                while (fe.next_itg_pt()) {
                    if (subdomain_->functionToRetainPhysNodes(fe.position().data()) == ibm::Partition::OUT) {
                        ActiveGP++; // active GPs
                        //InOutLocal.push_back(1);
                    } else {
                        //InOutLocal.push_back(0);
                    }
                    //PosLocal.push_back(fe.position());
                }

                if ((1 - (double) ActiveGP / fe.n_itg_pts()) > RatioGPSBM_) {
                    eleMarkers_[localElemID].set(ElementMarker::SBM_FALSE_INTERCEPTED, true); // Cut cell
                }
            }
        } else if (idata_->ImmersedMethodType == NSInputData::ImmersedBoundary) {

            if(outPts == npe){
                eleMarkers_[localElemID].set(ElementMarker::OUT_GP,true);
            }
            else if(outPts == 0){
                eleMarkers_[localElemID].set(ElementMarker::IN_GP, true);
            }
            else {
                eleMarkers_[localElemID].set(ElementMarker::INTERCEPTED_GP, true);
                /// Pure IMGA
                int ActiveGP = 0;
                fe.refill(0, idata_->RelativeOrderIntercepted);
                while (fe.next_itg_pt()) {
                    if (subdomain_->functionToRetainPhysNodes(fe.position().data()) == ibm::Partition::OUT) {
                        ActiveGP++; // active GPs
                    }
                }

                if (ActiveGP == 0) {
                    eleMarkers_[localElemID].set(ElementMarker::SBM_FALSE_INTERCEPTED, true); // Cut cell
                }
            }
        } else{
            std::cout << "Error: ImmersedMethodType not defined" << std::endl;
        }


    localElemID++;
}

void SBMMarker::WriteDomainGPToFile()
{
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = PosLocal.size();
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
        PosGlobal.resize(totalProcData);
        if (PosLocal.size() == InOutLocal.size()) {
            InOutGlobal.resize(totalProcData);
        }
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(PosLocal.data(), PosLocal.size(), ZEROPTVtype, PosGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    if (PosLocal.size() == InOutLocal.size()) {
        MPI_Gatherv(InOutLocal.data(), InOutLocal.size(), MPI_INT, InOutGlobal.data(), eachProcData.data(), disp.data(),
                    MPI_INT, 0, MPI_COMM_WORLD);
    }

    std::string fPrefix = "InterCeptedDomainGP.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");
        if (PosLocal.size() == InOutLocal.size()) {

#if (DIM == 2)
            fprintf(fp, "x0,y0,InOrOut\n");
#endif

#if (DIM == 3)
            fprintf(fp, "x0,y0,z0,InOrOut\n");
#endif

            for (int i = 0; i < PosGlobal.size(); i++) {
#if (DIM == 2)
                fprintf(fp, "%.10e,%.10e,%d\n",
                        PosGlobal[i](0), PosGlobal[i](1), InOutGlobal[i]);
#endif
#if (DIM == 3)
                fprintf(fp, "%.10e,%.10e,%.10e,%d\n",
                        PosGlobal[i](0), PosGlobal[i](1), PosGlobal[i](2), InOutGlobal[i]);
#endif
            }
        } else {
#if (DIM == 2)
            fprintf(fp, "x0,y0\n");
#endif

#if (DIM == 3)
            fprintf(fp, "x0,y0,z0\n");
#endif

            for (int i = 0; i < PosGlobal.size(); i++) {
#if (DIM == 2)
                fprintf(fp, "%.10e,%.10e\n",
                        PosGlobal[i](0), PosGlobal[i](1));
#endif
#if (DIM == 3)
                fprintf(fp, "%.10e,%.10e,%.10e\n",
                        PosGlobal[i](0), PosGlobal[i](1), PosGlobal[i](2));
#endif
            }
        }
        fclose(fp);
    }
}

#endif //DENDRITEKT_ELEMENTMARKER_H
