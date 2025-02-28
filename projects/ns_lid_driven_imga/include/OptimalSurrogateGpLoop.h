//
// Created by chenghau on 12/15/22.
//

#ifndef NSSBM_OPTIMALSURROGATEGPLOOP_H
#define NSSBM_OPTIMALSURROGATEGPLOOP_H
#include "Traversal/Traversal.h"
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"
#ifdef DEEPTRACE
#include "SBMCalcDeepTrace.h"
#else
#include "SBMCalc.h"
#endif
// Global counter variable


class OptimalSurrogateGpLoop:public Traversal{



    static constexpr  int numPoints = 1u << DIM;
    static constexpr  int numFaces = 2*DIM;
    const SubDomain * subdomain_;
    std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
    int localElemID = 0;
    SubDomainBoundary * subDomainBoundary_;

    std::vector<TALYFEMLIB::ZEROPTV> PosLocal;
    std::vector<TALYFEMLIB::ZEROPTV> PosGlobal;
    std::vector<double> PosDist;

    std::vector<TALYFEMLIB::ZEROPTV> NormalLocal;
    std::vector<TALYFEMLIB::ZEROPTV> NormalGlobal;
#ifdef DEEPTRACE
    std::vector<TALYFEMLIB::ZEROPTV> DISTANCEVEC;
    bool postprocess;
#endif
    std::vector<double> PosDistGlobal;
    std::vector<TALYFEMLIB::ZEROPTV> DISTANCEVECGLOBAL;

    double SumDistSquare = 0;
    int SumOfGP = 0;
    double SumDistSquareGlobal = 0;
    int SumOfGPGlobal = 0;

    /// For SBMCalc
    NSInputData *inputData_;
    const IMGA *imga_;
    std::vector <my_kd_tree_t *> kd_trees_;

    int NearestGeo(const FEMElm &fe);

public:
    OptimalSurrogateGpLoop(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
                           std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
                           const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary
            , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees);


    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    virtual void optimalSurrogatetraversaloperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values);

    void WriteOptimalGPToFile();
    void WriteOptimalGPDistToFile();

    void GetPos(std::vector<ZEROPTV> &GPpos);
    void GetNormal(std::vector<ZEROPTV> &GPnormal);

    void GetRMSDandSumOfGP(double &RMSD,int &SumOfGP_);
    void GetDistanceVector(std::vector<ZEROPTV> &DISTVEC);


    void write_distance_vector();

    void set_postprocess(void);

};

OptimalSurrogateGpLoop::OptimalSurrogateGpLoop(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
                                               std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary
        , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees)
        :eleMarkers_(eleMarkers), Traversal(octDA,treePart,v,domain),subdomain_(subDomain),subDomainBoundary_(subDomainBoundary)
        ,imga_(imga),inputData_(inputData){

    std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
                   [](const auto& uniquePtr) { return uniquePtr; });

    std::bitset<numFaces> mark;
    mark.reset();

    this->traverse();
}
//OptimalSurrogateGpLoop::OptimalSurrogateGpLoop(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
//                                               std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary
//        , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees,VecInfo v)
//        :eleMarkers_(eleMarkers), Traversal(octDA,treePart,v,domain),subdomain_(subDomain),subDomainBoundary_(subDomainBoundary)
//        ,imga_(imga),inputData_(inputData){
//
//    std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
//                   [](const auto& uniquePtr) { return uniquePtr; });
//
//    postprocess = false;
//
//    std::bitset<numFaces> mark;
//    mark.reset();
//
//    this->traverse();
//}

void OptimalSurrogateGpLoop::optimalSurrogatetraversaloperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {

//    THROW_NOT_IMPLEMENTED TALYFEMLIB::TALYException() << "You need to overload this" << __func__ << "\n";
    TALYFEMLIB::TALYException() << "You need to overload this" << __func__ << "\n";

}
void OptimalSurrogateGpLoop::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values)
{


#if (DIM == 3)
    static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0,1,2,3,4,5,6,7},{0,2,6,8,18,20,24,26}};
    static constexpr DENDRITE_UINT numNodesPerFace = 4;
    static constexpr DENDRITE_UINT faceID[2*DIM][4]
            {
                    {0,2,4,6}, // Left
                    {1,3,5,7}, // Right
                    {0,1,4,5}, // Bottom
                    {2,3,6,7}, // Top
                    {0,1,2,3}, // Back
                    {4,5,6,7},  // Front
            };
#elif(DIM == 2)
    static constexpr DENDRITE_UINT cornerMap[2][numPoints]{{0,1,2,3},{0,2,6,8}};
    static constexpr DENDRITE_UINT numNodesPerFace = 2;
    static constexpr DENDRITE_UINT faceID[2*DIM][2]
            {
                    {2,0}, // Left
                    {1,3}, // Right
                    {0,1}, // Bottom
                    {3,2}, // Top
            };

#else
    throw std::runtime_error("Not implemented\n");
#endif
    const int eleOrder = m_octDA->getElementOrder();
    if(eleMarkers_[localElemID].test(ElementMarker::INTERCEPTED_ELEMENT) or eleMarkers_[localElemID].test(ElementMarker::OUT_ELEMENT))
    {
        const int *surf_arr = m_elem->GetSurfaceCheckArray();
        const int surf_row_len = m_elem->GetSurfaceCheckArrayRowLength();

        std::vector<TALYFEMLIB::ZEROPTV> coords(m_octDA->getNumNodesPerElement());
        std::vector<bool> isBoundary(m_octDA->getNumNodesPerElement(),false);
        //std::array<std::bitset<MAX_BOUNDARY_TYPES>, numPoints> boundaryBits;

        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            std::memcpy(coords[i].data(), &this->m_coords[i*DIM],sizeof(double)*DIM);
        }
        for(int i = 0; i < m_octDA->getNumNodesPerElement(); i++) {
            unsigned int id;
            subDomainBoundary_->generateBoundaryFlags(coords[i],id);
            //boundaryBits[i] = subDomainBoundary_->getBoundary();
        }


        for(int i = 0; i < numFaces; i++) {
            bool isBoundaryFaceinterceptedElement = true;
            int surf_id = surf_arr[i * surf_row_len];  // id determines which nodes on the surface

            for(int j = 0; j < numNodesPerFace; j++) {
                unsigned int cornerNodesOnFace = faceID[i][j];
                unsigned int nodeID = cornerMap[eleOrder - 1][cornerNodesOnFace];
                isBoundaryFaceinterceptedElement = isBoundaryFaceinterceptedElement and isBoundary[nodeID];
            }

            bool SBMCheck = true;


//                std::cout<<"Pre_counter: "<<functionCallCounter_pre<<std::endl;

                for (int j = 0; j < numNodesPerFace; j++) {
                    SBMCheck = SBMCheck and ((subdomain_->functionToRetainPhysNodes(
                            &this->m_coords[cornerMap[eleOrder - 1][faceID[i][j]] * DIM]) ==
                                              ibm::IN)// This node belong to inactive region
                                             or ((values[cornerMap[eleOrder - 1][faceID[i][j]]] ==
                                                  NodeMarker::SBM_FALSE_INTERCEPTED_NODES)));

                }




            if(SBMCheck) {

                TALYFEMLIB::SurfaceIndicator surf(surf_id);
                TALYFEMLIB::ZEROPTV normal = m_elem->CalculateNormal(&m_grid, surf_id);
                surf.set_normal(normal);
                fe.refill_surface(m_elem, &surf, 0);

                    while (fe.next_itg_pt()) {


                        bool OnWall = false;
                        static const double eps = 1e-10;
                        auto position = fe.position();
                        bool x_minus_wall = fabs(position.x() - inputData_->physDomain.min[0]) < eps;
                        bool y_minus_wall = fabs(position.y() - inputData_->physDomain.min[1]) < eps;
                        bool x_max_wall = fabs(position.x() - inputData_->physDomain.max[0]) < eps;
                        bool y_max_wall = fabs(position.y() - inputData_->physDomain.max[1]) < eps;
                        bool z_minus_wall = false;
                        bool z_max_wall = false;
#if(DIM == 3)
                        z_minus_wall = fabs(position.z() - inputData_->physDomain.min[2]) < eps;
        z_max_wall = fabs(position.z() - inputData_->physDomain.max[2]) < eps;
#endif
#if (DIM == 2)
                        if (x_minus_wall || y_minus_wall || x_max_wall || y_max_wall)
#endif
#if (DIM == 3)
                            if (x_minus_wall || y_minus_wall || x_max_wall || y_max_wall || z_minus_wall || z_max_wall)
#endif
                        {
                            //TODO: fix hard code
                            OnWall = (inputData_->carved_out_geoms_def[0].outer_boundary == true) ? false : true;
                        }
//                    std::cout << "OnWall = " << OnWall << "\n";

                        if (!OnWall) {
                            PosLocal.push_back(fe.position());
                            NormalLocal.push_back(fe.surface()->normal());

                            const double query_pt[3] = {fe.position().x(), fe.position().y(),
                                                        fe.position().z()}; // shift_

//                        std::cout << "Surrogate GP pos = ";fe.position().print(); std::cout << "\n";
#ifdef DEEPTRACE
                            SBMCalcDeepTrace sbmCalc(fe, inputData_, imga_, kd_trees_, 0);
#else
                            SBMCalc sbmCalc(fe,inputData_,imga_,kd_trees_,0);
#endif
                            double d[DIM];
                            sbmCalc.Dist2Geo(d);
                            TALYFEMLIB::ZEROPTV distance_vector;

                            double distMag = 0.0;
                            for (int dim = 0; dim < DIM; dim++) {
                                distMag += pow(d[dim], 2);

                            }
                            distance_vector.x() = d[0];
                            distance_vector.y() = d[1];
#if DIM == 3
                            distance_vector.z()=d[2];
#endif
#ifdef DEEPTRACE
                            DISTANCEVEC.push_back(distance_vector);
#endif

                            PosDist.push_back(std::sqrt(distMag));

//                        if (distMag > 1e-2) {
//                            std::cout << "std::sqrt(distMag) = " << std::sqrt(distMag) << "\n";
//                        }



                            // TODO: fix below bug in future!
//                    double d[DIM];
//
//                    // TODO: multiple geometries
//                    SBMCalc sbmCalc(fe,inputData_,imgaTrueGP_,0);
//                    sbmCalc.Dist2Geo(d);
//#if (DIM ==2)
//                    SumDistSquare += pow(d[0],2) + pow(d[1],2);
//#endif
//#if (DIM ==3)
//                    SumDistSquare += pow(d[0],2) + pow(d[1],2) + pow(d[2],2);
//#endif
//                    SumOfGP++;
                        }  //end if !ONWALL
                    } //end while
            } //end if SBMCheck
        } //end for

    } //end if INTERCEPTED_ELEMENT or OUT_ELEMENT
    localElemID++;
} //end traverseOperation

void OptimalSurrogateGpLoop::WriteOptimalGPToFile() {
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
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(PosLocal.data(), PosLocal.size(), ZEROPTVtype, PosGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    std::string fPrefix = "optimalGP.vtk";


//    std::string fPrefix = "optimalGP_ppo.vtk";
    if (TALYFEMLIB::GetMPIRank() == 0)
    { // if:master cpu, then:print
        FILE *fp = fopen(fPrefix.c_str(), "w");

#if (DIM == 2)
        fprintf(fp, "x0,y0\n");
#endif

#if (DIM == 3)
        fprintf(fp, "x0,y0,z0\n");
#endif

        for (int i = 0; i < PosGlobal.size(); i++)
        {
#if (DIM == 2)
            fprintf(fp, "%.10e,%.10e\n",
                    PosGlobal[i](0), PosGlobal[i](1));
#endif
#if (DIM == 3)
            fprintf(fp, "%.10e,%.10e,%.10e\n",
                    PosGlobal[i](0), PosGlobal[i](1), PosGlobal[i](2));
#endif
        }

        fclose(fp);
    }

}

void OptimalSurrogateGpLoop::WriteOptimalGPDistToFile() {
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
        PosDistGlobal.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(PosLocal.data(), PosLocal.size(), ZEROPTVtype, PosGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(PosDist.data(), PosDist.size(), MPI_DOUBLE, PosDistGlobal.data(), eachProcData.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (TALYFEMLIB::GetMPIRank() == 0) {

            std::string fname2 = "Dist2TrueBoundary.dat";
            FILE *fp2 = fopen(fname2.c_str(), "a");
            if (!fp2) {
                throw TALYException() << "Cannot create file: " << fname2;
            }

            for (int i = 0; i < PosDistGlobal.size(); i++) {

//                std::cout << "PosDistGlobal = " << PosDistGlobal[i] << "\n";
#if (DIM == 3)
                fprintf(fp2,
                    "x = % .10e, y = % .10e, z = %.10e, value = %.10e\n",
                    PosGlobal[i].x(),PosGlobal[i].y()
                    ,PosGlobal[i].z(),PosDistGlobal[i]);
#endif
#if (DIM == 2)
                fprintf(fp2,
                        "x = % .10e, y = % .10e, value = %.10e\n",
                        PosGlobal[i].x(),PosGlobal[i].y()
                        ,PosDistGlobal[i]);
#endif
            }
            fclose(fp2);
        }


}
//write the distance vector as well as the position of the optimal GP
void OptimalSurrogateGpLoop::write_distance_vector() {
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = PosLocal.size();
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
        PosGlobal.resize(totalProcData);
        PosDistGlobal.resize(totalProcData);
        DISTANCEVECGLOBAL.resize(totalProcData);
    }

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Gatherv(PosLocal.data(), PosLocal.size(), ZEROPTVtype, PosGlobal.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);
    MPI_Gatherv(PosDist.data(), PosDist.size(), MPI_DOUBLE, PosDistGlobal.data(), eachProcData.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(DISTANCEVEC.data(), DISTANCEVEC.size(), ZEROPTVtype, DISTANCEVECGLOBAL.data(), eachProcData.data(), disp.data(), ZEROPTVtype, 0, MPI_COMM_WORLD);

    if (TALYFEMLIB::GetMPIRank() == 0) {
        static bool headerWritten = false;
        std::string fname2 = "Dist2TrueBoundaryVector.csv";
        FILE *fp2 = fopen(fname2.c_str(), "a");
        if (!fp2) {
            throw TALYException() << "Cannot create file: " << fname2;
        }

        if (!headerWritten) {
#if (DIM == 3)
            fprintf(fp2, "x, y, z, value, dx, dy, dz\n");
#endif
#if (DIM == 2)
            fprintf(fp2, "x, y, value, dx, dy\n");
#endif
            headerWritten = true;
        }

        for (int i = 0; i < PosDistGlobal.size(); i++) {
#if (DIM == 3)
            fprintf(fp2,
                    "x = % .10e, y = % .10e, z = %.10e, value = %.10e, dx = %.10e, dy = %.10e, dz = %.10e\n",
                    PosGlobal[i].x(), PosGlobal[i].y(), PosGlobal[i].z(), PosDistGlobal[i],
                    DISTANCEVECGLOBAL[i].x(), DISTANCEVECGLOBAL[i].y(), DISTANCEVECGLOBAL[i].z());
#endif
#if (DIM == 2)
            fprintf(fp2,
                    "% .10e,% .10e,%.10e,%.10e,%.10e\n",
                    PosGlobal[i].x(), PosGlobal[i].y(), PosDistGlobal[i],
                    DISTANCEVECGLOBAL[i].x(), DISTANCEVECGLOBAL[i].y());
#endif
        }
        fclose(fp2);
    }
    if(TALYFEMLIB::GetMPIRank() == 0) {
        TALYFEMLIB::PrintInfo("Distance vector written to file\n");
    }
}
void OptimalSurrogateGpLoop::GetPos(std::vector<ZEROPTV> &GPpos) {
    GPpos = PosLocal;
}
#ifdef DEEPTRACE
void OptimalSurrogateGpLoop::GetDistanceVector(std::vector<ZEROPTV> &DISTVEC) {
    DISTVEC = DISTANCEVEC;
}
#endif


void OptimalSurrogateGpLoop::GetRMSDandSumOfGP(double &RMSD,int &SumOfGP_) {

    MPI_Reduce(&SumOfGP, &SumOfGPGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&SumDistSquare, &SumDistSquareGlobal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (TALYFEMLIB::GetMPIRank() == 0) {
        RMSD = sqrt(SumDistSquareGlobal/SumOfGPGlobal);
        SumOfGP_ = SumOfGPGlobal;
    }
}

void OptimalSurrogateGpLoop::GetNormal(std::vector<ZEROPTV> &GPnormal) {
    GPnormal = NormalLocal;
}

int OptimalSurrogateGpLoop::NearestGeo(const TALYFEMLIB::FEMElm &fe) {
    const TALYFEMLIB::ZEROPTV pt = fe.position();
    double Min_Dist = std::numeric_limits<double>::max();
    int PickGeomID = 0;

    for (int CarvedOutGeomID = 0;
         CarvedOutGeomID < imga_->getGeometries().size();
         ++CarvedOutGeomID) {

//            std::cout << "CarvedOutGeomID = " << CarvedOutGeomID << "\n";

        double d[DIM];
        SBMCalc sbmCalc(fe, inputData_, imga_, kd_trees_, CarvedOutGeomID);
        sbmCalc.Dist2Geo(d);

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


#endif //NSSBM_OPTIMALSURROGATEGPLOOP_H
