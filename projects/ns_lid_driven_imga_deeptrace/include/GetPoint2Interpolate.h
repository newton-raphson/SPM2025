//
// Created by chenghau on 12/1/22.
//

#ifndef NSSBM_GETPOINT2INTERPOLATE_H
#define NSSBM_GETPOINT2INTERPOLATE_H
#include <IMGA/IMGATraversal.h>
#include <DendriteUtils.h>
#include <NSInputData.h>
#include <Boundary/SubDomainBoundary.h>

using namespace TALYFEMLIB;

using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
        PointCloud<double>, 3 /* dim */>;


class GetPoint2Interpolate : public IMGATraversal
{
    ///
    static constexpr  int numPoints = 1u << DIM;
    static constexpr  int numFaces = 2*DIM;
    std::vector<std::bitset<numFaces>> faceMarker_;

    bool hasLocalGaussPoint = true;
    std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoIt_;
    std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoEnd_;
    const SubDomain * subdomain_;
    const NSInputData *idata_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    SubDomainBoundary * subDomainBoundary_;

    DENDRITE_UINT lclElemID = 0;

    //Point<DIM> shift_;
    double valueWJ = 0.0;
    double totalGP = 0.0;
    double x_true, y_true;

    std::vector<ZEROPTV> Pt2Interpolate;
    std::vector<ZEROPTV> GP;
    std::vector<ZEROPTV> Pt2InterpolateGlobal;
    std::vector<ZEROPTV> GPGlobal;

    std::vector<ZEROPTV> SurrogateGPGlobal_;

    DENDRITE_REAL ElementSize(const TALYFEMLIB::FEMElm &fe);

    const my_kd_tree_t *kd_tree_;
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
    GetPoint2Interpolate(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v, const DomainExtents &domain,
                         const SubDomain * subDomain, const NSInputData *idata,
                         std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers
            ,SubDomainBoundary *subDomainBoundary,std::vector<ZEROPTV> &GPposVector,
            const my_kd_tree_t *kd_tree);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    void imgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values) override;

    void ImgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint, const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values, bool TrueIntercepted
            ,ZEROPTV &GPpos);

    void WriteInterGPToFile();

    void GetpositionOfpointToInterpolateOnAndpositionOnSurface(std::vector<TALYFEMLIB::ZEROPTV> &positionOfpointToInterpolateOn,std::vector<TALYFEMLIB::ZEROPTV> &positionOnSurface);
};

GetPoint2Interpolate::GetPoint2Interpolate(DA *octDA, const IMGA *imga, const std::vector<TREENODE> &treePart, const VecInfo &v,
                                           const DomainExtents &domain,
                                           const SubDomain * subDomain, const NSInputData *idata,
                                           std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,SubDomainBoundary *subDomainBoundary
        ,std::vector<ZEROPTV> &GPposVector,const my_kd_tree_t *kd_tree)
        : IMGATraversal(octDA, imga, treePart, v, domain)
        ,idata_(idata),subdomain_(subDomain),eleMarkers_(eleMarkers), subDomainBoundary_(subDomainBoundary),SurrogateGPGlobal_(GPposVector),kd_tree_(kd_tree)
{
    std::bitset<numFaces> mark;
    mark.reset();

    faceMarker_.resize(octDA->getLocalElementSz());

    const auto & gpInfo_ = imga->getSurfaceGaussPoints();
    if(gpInfo_.empty()){
        hasLocalGaussPoint = false;
    }
    else {
        gpInfoIt_ = gpInfo_.cbegin();
        gpInfoEnd_ = gpInfo_.cend();
    }
    MPI_Barrier(octDA->getCommActive());
    this->imgaTraverse();
}

void GetPoint2Interpolate::ImgaTraversalOperation(const TALYFEMLIB::FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                                  const TALYFEMLIB::ZEROPTV &h, const PetscScalar *values,bool TrueIntercepted
        ,ZEROPTV &GPpos)
{
    GP.push_back(fe.position());


    /// If it is on TrueIntercepted element, return True GPs.
    /// If it is on FalseIntercepted element, return Surrogate GPs with moving a little bit toward ACTIVE region.
    if (TrueIntercepted){
        Pt2Interpolate.push_back(fe.position());
    } else{
        Pt2Interpolate.push_back(GPpos);
    }

}

void GetPoint2Interpolate::traverseOperation(FEMElm &fe, const PetscScalar *values) {

    bool TrueIntercepted;

    if (eleMarkers_[lclElemID].test(ElementMarker::SBM_FALSE_INTERCEPTED))
    {
        TrueIntercepted = false;

    } else {
        TrueIntercepted = true;
    }

    ZEROPTV SurrogateGPposWithSmallMove;
    //std::cout <<" nodePTVs.size() = "<< nodePTVs.size() <<"\n";
    if(hasLocalGaussPoint) {
        while ((gpInfoIt_->localElemID==lclElemID) and (gpInfoIt_ < gpInfoEnd_)) {
            const auto gaussPoint = *gpInfoIt_;

            /// TODO Fix it. Seems stupid;
            TALYFEMLIB::ZEROPTV ptvg, ptvl;
            std::memcpy(ptvg.data(), gaussPoint.location, sizeof(double)*DIM);
            GetLocalPtv(fe, ptvg, ptvl);
            fe.calc_at(ptvl);
            const auto &coords = this->m_coords;
            TALYFEMLIB::ZEROPTV h;
            const DENDRITE_UINT nPe = this->m_octDA->getNumNodesPerElement();
            for (int d = 0; d < DIM; d++) {
                h.data()[d] = coords[(nPe - 1)*DIM + d] - coords[d];
            }
            static constexpr DENDRITE_REAL weight = (DIM==2) ? 2.0 : 3.0;
            fe.set_jacc_x_w(gaussPoint.elemArea/weight);

            /// Obtaining the closest SurrogateGPpos
            size_t num_results = 1;
            std::vector<uint32_t> ret_index(num_results);
            std::vector<double> out_dist_sqr(num_results);

            double query_pt[3];
            std::copy_n(fe.position().data(), 3, query_pt);

            num_results = kd_tree_->knnSearch(
                    &query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

            int CloseGPIndex = ret_index[0];

            ZEROPTV TrueGP2SurrogateGP = SurrogateGPGlobal_[CloseGPIndex] - fe.position();
            TrueGP2SurrogateGP.SafeNormalize();
            SurrogateGPposWithSmallMove= SurrogateGPGlobal_[CloseGPIndex] + TrueGP2SurrogateGP * ElementSize(fe) * 0.01;


            ImgaTraversalOperation(fe, *gpInfoIt_, h, values, TrueIntercepted, SurrogateGPposWithSmallMove);
            gpInfoIt_ = std::next(gpInfoIt_);
        }
    }
    lclElemID++;
}

void GetPoint2Interpolate::WriteInterGPToFile() {
    int nProc = TALYFEMLIB::GetMPISize(); // be careful

    int numNodes = GP.size();
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
        GPGlobal.resize(totalProcData);
        Pt2InterpolateGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(GP.data(), GP.size(), ZEROPTVtype, GPGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(Pt2Interpolate.data(), Pt2Interpolate.size(), ZEROPTVtype, Pt2InterpolateGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    std::string fPrefix = "InterceptedPtWithGP.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "Do_x0,Do_y0,GP_x0,GP_y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "Do_x0,Do_y0,Do_z0,GP_x0,GP_y0,GP_z0\n");
#endif

        for (int i = 0; i < GPGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e\n",
                    Pt2InterpolateGlobal[i](0), Pt2InterpolateGlobal[i](1),GPGlobal[i](0), GPGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                    Pt2InterpolateGlobal[i](0), Pt2InterpolateGlobal[i](1), Pt2InterpolateGlobal[i](2)
                    ,GPGlobal[i](0), GPGlobal[i](1) , GPGlobal[i](2));
#endif
        }

        fclose(fp);
    }

}

void GetPoint2Interpolate::imgaTraversalOperation(const FEMElm &fe, const NodeAndValues<DENDRITE_REAL> &gaussPoint,
                                                  const ZEROPTV &h, const PetscScalar *values) {

}

void GetPoint2Interpolate::GetpositionOfpointToInterpolateOnAndpositionOnSurface(
        std::vector<TALYFEMLIB::ZEROPTV> &positionOfpointToInterpolateOn,
        std::vector<TALYFEMLIB::ZEROPTV> &positionOnSurface) {

    positionOfpointToInterpolateOn= Pt2InterpolateGlobal;
    positionOnSurface = GPGlobal;
}

DENDRITE_REAL GetPoint2Interpolate::ElementSize(const FEMElm &fe) {
    return pow((pow(2, DIM) * fe.jacc()), (double)1 / DIM);
}

#endif //NSSBM_GETPOINT2INTERPOLATE_H
