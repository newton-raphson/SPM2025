//
// Created by mehdi on 7/31/23.
//

#ifndef DENDRITEKT_NSPOST_H
#define DENDRITEKT_NSPOST_H

#ifndef LINNS_FORCECALC_H
#define LINNS_FORCECALC_H

#include <DendriteUtils.h>
#include "NSNodeData.h"
#include "Traversal/Traversal.h"
#include "SubDA/SubDomain.h"
#include "Boundary/SubDomainBoundary.h"
#include "IMGA/IMGA.h"

#include <fstream>
#include <iostream>

#include "SBMCalcDeepTrace.h"


void calcSecondDerivativeFEM(const TALYFEMLIB::FEMElm &fe, DENDRITE_UINT ndof,
                             const DENDRITE_REAL *value, DENDRITE_REAL *val_dd) {

    for (int dof = 0; dof < ndof; dof++) {
        for (int dir1 = 0; dir1 < DIM; dir1++) {
            for (int dir2 = 0; dir2 < DIM; dir2++) {
                val_dd[DIM * DIM * dof + DIM * dir1 + dir2] = 0;
                for (ElemNodeID a = 0; a < fe.nbf(); a++) {
                    val_dd[DIM * DIM * dof + DIM * dir1 + dir2] += fe.d2N(a, dir1, dir2) * value[a * ndof + dof];
                }
            }
        }
    }
}
 class ForceCalc: public Traversal{

     static constexpr  int numPoints = 1u << DIM;
     static constexpr  int numFaces = 2*DIM;
     const SubDomain * subdomain_;
     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers_;
     int localElemID = 0;
     SubDomainBoundary * subDomainBoundary_;

     const IMGA *imga_;


     ForceInfo localforce_info;
     VecInfo ns_solution_;
     std::vector<TALYFEMLIB::ZEROPTV> PosLocal;
     std::vector<TALYFEMLIB::ZEROPTV> PosGlobal;
     NSInputData *inputData_;
     std::vector<my_kd_tree_t *> kd_trees_;
     const TimeInfo ti_;
public:
    ForceCalc(DA * octDA, const std::vector<TREENODE> &treePart,  const VecInfo & nsSolution, const DomainExtents &domain,
              std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, SubDomain* subDomain,SubDomainBoundary *subDomainBoundary
            , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees,const TimeInfo &ti);
    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
    void WriteOptimalGPToFile();
    void calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe,const PetscScalar * values,TALYFEMLIB::ZEROPTV distance_vector);
    void getTotalForce(ForceInfo &globalForce);
    void writeTotalForce();

 };

ForceCalc::ForceCalc(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo & nsSolution,
                     const DomainExtents &domain, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> &eleMarkers,
                     SubDomain* subDomain, SubDomainBoundary *subDomainBoundary, NSInputData *inputData,
                     const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees,const TimeInfo &ti)
                     :eleMarkers_(eleMarkers),imga_(imga),kd_trees_(kd_trees), Traversal(octDA,treePart,nsSolution,domain),subdomain_(subDomain),subDomainBoundary_(subDomainBoundary),inputData_(inputData),ti_(ti)
{

    std::bitset<numFaces> mark;
    mark.reset();

    this->traverse();

}
void ForceCalc::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values)
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


                for (int j = 0; j < numNodesPerFace; j++) {
                    SBMCheck = SBMCheck and ((subdomain_->functionToRetainPhysNodes(
                            &this->m_coords[cornerMap[eleOrder - 1][faceID[i][j]] * DIM]) ==
                                              ibm::IN)// This node belong to inactive region
                                             or
                                             ((values[cornerMap[eleOrder - 1][faceID[i][j]] * (DIM + 2) + (DIM + 1)] ==
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

                            //#ifdef DEEPTRACE
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

                            calcForceAtEachQuadraturePoint2D(fe, values,distance_vector);

//                            NormalLocal.push_back(fe.surface()->normal());
//
//                            const double query_pt[3] = {fe.position().x(), fe.position().y(),
//                                                        fe.position().z()}; // shift_
//
////                        std::cout << "Surrogate GP pos = ";fe.position().print(); std::cout << "\n";
//#ifdef DEEPTRACE
//                            SBMCalcDeepTrace sbmCalc(fe, inputData_, imga_, kd_trees_, 0);
//#else
//                            SBMCalc sbmCalc(fe,inputData_,imga_,kd_trees_,0);
//#endif
//                            double d[DIM];
//                            sbmCalc.Dist2Geo(d);
//                            TALYFEMLIB::ZEROPTV distance_vector;
//
//                            double distMag = 0.0;
//                            for (int dim = 0; dim < DIM; dim++) {
//                                distMag += pow(d[dim], 2);
//
//                            }
//                            distance_vector.x() = d[0];
//                            distance_vector.y() = d[1];
//#if DIM == 3
//                            distance_vector.z()=d[2];
//#endif
//#ifdef DEEPTRACE
//                            DISTANCEVEC.push_back(distance_vector);
//#endif
//
//                            PosDist.push_back(std::sqrt(distMag));

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

void ForceCalc::WriteOptimalGPToFile() {
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
        std::string fPrefix = "optimalGP_ppo.vtk";

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

//
//class ForceCalc: public OptimalSurrogateGpLoop
//{
//    static constexpr  int numPoints = 1u << DIM;
//    static constexpr  int numFaces = 2*DIM;
//    ForceInfo localforce_info;
//    VecInfo ns_solution_;
//    std::vector<TALYFEMLIB::ZEROPTV> PosLocal;
//    std::vector<TALYFEMLIB::ZEROPTV> PosGlobal;
//    NSInputData *inputData_;
//
//public:
//        ForceCalc(DA * octDA, const std::vector<TREENODE> &treePart,  Vec nodalFalseElement, const DomainExtents &domain,
//                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, SubDomain subDomain,SubDomainBoundary *subDomainBoundary
//            , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees,const VecInfo & nsSolution);
//        void optimalSurrogatetraversaloperation(TALYFEMLIB::FEMElm & fe,const PetscScalar *values) override;
//        void calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe,const PetscScalar * values);
////        void getTotalForce(DENDRITE_REAL * globalForce, DENDRITE_REAL * projectedArea);
//
//        void getTotalForce(ForceInfo &globalForce);
//};
////ForceCalc::ForceCalc(DA *octDA, const std::vector<TREENODE> &treePart,  Vec nodalFalseElement, const DomainExtents &domain,
////                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> &eleMarkers, SubDomain subDomain,
////                     SubDomainBoundary *subDomainBoundary, NSInputData *inputData, const IMGA *imga,
////                     const std::vector<my_kd_tree_t *> kd_trees,const VecInfo & nsSolution): OptimalSurrogateGpLoop(octDA,treePart,{VecInfo(nodalFalseElement, 1, NSNodeData::NODE_ID)},domain,eleMarkers,&subDomain,subDomainBoundary,inputData,imga,kd_trees),ns_solution_(nsSolution),inputData_(inputData) {
////    this->set_postprocess();
////    this->traverse();
////}
//
//ForceCalc::ForceCalc(DA *octDA, const std::vector<TREENODE> &treePart,  Vec nodalFalseElement, const DomainExtents &domain,
//                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> &eleMarkers, SubDomain subDomain,
//                     SubDomainBoundary *subDomainBoundary, NSInputData *inputData, const IMGA *imga,
//                     const std::vector<my_kd_tree_t *> kd_trees,const VecInfo & nsSolution): OptimalSurrogateGpLoop(octDA,treePart,nsSolution,domain,eleMarkers,&subDomain,subDomainBoundary,inputData,imga,kd_trees),ns_solution_(nsSolution),inputData_(inputData) {
//    this->set_postprocess();
//    std::bitset<numFaces> mark;
//    mark.reset();
//
//    this->traverse();
//}
//
//void ForceCalc::optimalSurrogatetraversaloperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {
//    // Open a file for writing
//
////    I WANT TO CHECK WHAT VALUES ARE BEING PASSED HERE
//
//    while (fe.next_itg_pt()) {
//        bool OnWall = false;
//        static const double eps = 1e-10;
//        auto position = fe.position();
//        bool x_minus_wall = fabs(position.x() - inputData_->physDomain.min[0]) < eps;
//        bool y_minus_wall = fabs(position.y() - inputData_->physDomain.min[1]) < eps;
//        bool x_max_wall = fabs(position.x() - inputData_->physDomain.max[0]) < eps;
//        bool y_max_wall = fabs(position.y() - inputData_->physDomain.max[1]) < eps;
//        bool z_minus_wall = false;
//        bool z_max_wall = false;
//#if(DIM == 3)
//        z_minus_wall = fabs(position.z() - inputData_->physDomain.min[2]) < eps;
//        z_max_wall = fabs(position.z() - inputData_->physDomain.max[2]) < eps;
//#endif
//#if (DIM == 2)
//        if (x_minus_wall || y_minus_wall || x_max_wall || y_max_wall)
//#endif
//#if (DIM == 3)
//            if (x_minus_wall || y_minus_wall || x_max_wall || y_max_wall || z_minus_wall || z_max_wall)
//#endif
//        {
//            //TODO: fix hard code
//            OnWall = (inputData_->carved_out_geoms_def[0].outer_boundary == true) ? false : true;
//            //TODO: fix hard code
////            OnWall = false;
//        }
////                    std::cout << "OnWall = " << OnWall << "\n";
//
//        if (!OnWall) {
//            calcForceAtEachQuadraturePoint2D(fe, values);
//            PosLocal.push_back(fe.position());
//        }
//    }
//
//
//}
void ForceCalc::calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe, const PetscScalar *values,TALYFEMLIB::ZEROPTV distance_vector) {

//    let's loop over the petsc value to make sure that the values coming here as we expect
//    the sequence shoulbe u,p,v and mask
//    repeated over 4 nodes in 2D
//
//    for (int i = 0; i < 4; i++)
//    {
//        std::cout << "u = " << values[4*i] << "\n";
//        std::cout << "p = " << values[4*i+1] << "\n";
//        std::cout << "v = " << values[4*i+2] << "\n";
//        std::cout << "mask = " << values[4*i+3] << "\n";
//    }

//    using namespace TALYFEMLIB;
//    const int nsd = DIM;
//    double Coe_diff = 1.0 / inputData_->Re;
//    ZeroMatrix<double> du;
//    du.redim(nsd, nsd);
//    DENDRITE_REAL vald[(NSNodeData::NS_DOF+1)*DIM];
//    DENDRITE_REAL valc[NSNodeData::NS_DOF+1];
//    calcValueDerivativeFEM(fe,(NSNodeData::NS_DOF+1),values,vald);
//    calcValueFEM(fe,NSNodeData::NS_DOF+1,values,valc);
//
//    double  p = valc[NSNodeData::PRESSURE];
//
//    ZEROPTV Normal;
//
//    Normal.x()=distance_vector[0];
//    Normal.y()=distance_vector[1];
//    Normal.z()=0;
//#if (DIM==3)
//    Normal.z()=distance_vector[2];
//#endif
//    Normal.SafeNormalize();
////
//    for (int i = 0; i < nsd; i++) {
//        for (int j = 0; j < nsd; j++) {
//            du(i, j) = vald[i*DIM+j] + vald[j*DIM + i];
//        }
//    }
//    double GradPdotD;
//    for (int j = 0; j < DIM; j++){
//        GradPdotD += vald[DIM * DIM + j] * distance_vector[j];
//    }
//
//    DENDRITE_REAL valdd[NSNodeData::NS_DOF * DIM *DIM];
//    calcSecondDerivativeFEM(fe, NSNodeData::NS_DOF, values, valdd);
//    double GradTauDotD[DIM][DIM];
//
//    for (int k =0;k <DIM; k++)
//    {
//        for (int i = 0; i < DIM; i++)
//        {
//            for (int j = 0; j < DIM; j++)
//            {
//                GradTauDotD[i][j] +=  (valdd[i * DIM * DIM + k * DIM + j] + valdd[j * DIM * DIM + k * DIM + i])* distance_vector[k];
//            }
//        }
//    }
//    ZEROPTV diff;
//    ZEROPTV GradUdotD_NS;
//    ZEROPTV GradTauDotDDotN;
//    for (int i = 0; i < DIM; i++)
//    {
//        diff(i) = 0;
//        GradUdotD_NS(i) =0;
//        GradTauDotDDotN(i) = 0;
//        for (int j = 0; j < DIM; j++)
//        {
//            diff(i) += Coe_diff * du(i, j) * Normal(j);
//            GradTauDotDDotN(i) += Coe_diff * GradTauDotD[i][j] * Normal(j);
//            GradUdotD_NS(i) += vald[i * DIM + j] * distance_vector[j]; // from surrogate to true boundary
//        }
//    }
//
//    ZEROPTV Force_pressure;
//    ZEROPTV Force_viscous;
//    ZEROPTV Force_penalty;
//
//    for (int dim = 0; dim < DIM; dim++)
//    {
//        Force_pressure(dim) = (p + GradPdotD) * Normal(dim) * fe.surface()->normal().innerProduct(Normal) * fe.detJxW();
//        Force_viscous(dim) = -(diff(dim) + GradTauDotDDotN(dim)) *  fe.surface()->normal().innerProduct(Normal) * fe.detJxW();
////        Force_penalty(dim) = (Cb_f / hc) * (u[dim] + GradUdotD_NS(dim) - Ug(dim)) * SurrogateDotTrueNormal *
////                             fe.detJxW();
//    }
//    localforce_info.Force_pressure += Force_pressure;
//    localforce_info.Force_viscous += Force_viscous;
/*
 * This function calculates the force at each quadrature point below is implementation without distance vec
 */
    using namespace TALYFEMLIB;
    const int nsd = DIM;
    double Coe_diff = 1.0 /inputData_->Re ;
    ZeroMatrix<double> du;
    du.redim(nsd, nsd);
    DENDRITE_REAL vald[(NSNodeData::NS_DOF+1)* DIM];
    DENDRITE_REAL valc[NSNodeData::NS_DOF+1];
    calcValueDerivativeFEM(fe, NSNodeData::NS_DOF+1, values, vald);
    calcValueFEM(fe, NSNodeData::NS_DOF+1, values, valc);

    double p = valc[NSNodeData::PRESSURE];

//    std::cout<< "Pressure one node = " << p << "\n";
//    std::cout<< "Pressure second node = " << values[NSNodeData::PRESSURE+4] << "\n";

    for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
            du(i, j) = vald[i * DIM + j] + vald[j * DIM + i];
        }
    }
    ZEROPTV diff;
    ZEROPTV surface_normal = fe.surface()->normal();

    for (int i = 0; i < nsd; i++) {
        diff(i) = 0;
        for (int j = 0; j < nsd; j++) {
            diff(i) += Coe_diff * du(i, j) *surface_normal(j);
        }
    }
    ZEROPTV Force_pressure;
    ZEROPTV Force_viscous;

    for (int dim = 0; dim < DIM; dim++) {
        Force_pressure(dim) = (-p *surface_normal(dim)) * fe.detJxW();
        Force_viscous(dim) = diff(dim) * fe.detJxW();
    }
    localforce_info.Force_pressure += Force_pressure;
    localforce_info.Force_viscous += Force_viscous;
}

void ForceCalc::getTotalForce(ForceInfo &globalForce){
    MPI_Allreduce(&localforce_info.Force_pressure, &globalForce.Force_pressure, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
    MPI_Allreduce(&localforce_info.Force_viscous, &globalForce.Force_viscous, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());
    MPI_Allreduce(&localforce_info.Force_penalty, &globalForce.Force_penalty, 3, MPI_DOUBLE, MPI_SUM, m_octDA->getCommActive());


}

void ForceCalc::writeTotalForce()
{
    ForceInfo globalForce;
    getTotalForce(globalForce);

  PrintStatus(" DeepTrace force: ", globalForce.Force_all(true));
    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        std::string fname = "Force.dat";
        FILE *fp = fopen(fname.c_str(), "a");
        if (!fp)
        {
            throw TALYException() << "Cannot create file: " << fname;
        }
         PrintStatus("Deeptrace", " force: ", globalForce.Force_all(false));
#if (DIM == 3)
        fprintf(fp,
            "Timestep: %1d, Time: %.5f\n"
            "Force_pressure = % .10e, % .10e, % .10e\n"
            "Force_viscous = % .10e, % .10e, % .10e\n"
            "Force_penalty = % .10e, % .10e, % .10e\n",
                ti_.getTimeStepNumber(),
                ti_.getCurrentTime(),
            globalForce.Force_pressure.x(), globalForce.Force_pressure.y(), globalForce.Force_pressure.z(),
            globalForce.Force_viscous.x(), globalForce.Force_viscous.y(), globalForce.Force_viscous.z(),
            globalForce.Force_penalty.x(), globalForce.Force_penalty.y(), globalForce.Force_penalty.z());
#endif
#if (DIM == 2)
        fprintf(fp,
                "Force_pressure = % .10e, % .10e\n"
                "Force_viscous = % .10e, % .10e\n"
                "Force_penalty = % .10e, % .10e\n",
                globalForce.Force_pressure.x(), globalForce.Force_pressure.y(),
                globalForce.Force_viscous.x(), globalForce.Force_viscous.y(),
                globalForce.Force_penalty.x(), globalForce.Force_penalty.y());
#endif
        fclose(fp);
    }
}
//class ForceCalcOptimal: public OptimalSurrogateGpLoop{
//    ForceInfo globalForce;
//
//public:
//    ForceCalcOptimal(DA * octDA, const std::vector<TREENODE> &treePart, const VecInfo & v, const DomainExtents &domain,
//                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers, const SubDomain * subDomain,SubDomainBoundary *subDomainBoundary
//            , NSInputData *inputData,const IMGA *imga, const std::vector<my_kd_tree_t *> kd_trees);
//    void performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
//                                 const BoundarySurface &boundarySurface, const PetscScalar *values) override;
//
//};
//
//class ForceCalc: public SurfaceLoop{
//    DENDRITE_REAL Force_[2*DIM]{};
//    DENDRITE_REAL projectedArea_[DIM]{};
//    const DENDRITE_REAL Re_;
//    const OptimalSurrogateStruct optimalSurrogate_Data;
//
//public:
//    ForceCalc(DA * octDA,const std::vector<TREENODE> & treePart,const VecInfo & nsSolution, const DomainExtents & domain,
//              const DENDRITE_REAL Re, SubDomainBoundary *subDomainBoundary);
//
//    void performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
//                            const BoundarySurface &boundarySurface, const PetscScalar *values) override;
//
//    void calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe,const PetscScalar * values);
//
//    void getTotalForce(DENDRITE_REAL * globalForce, DENDRITE_REAL * projectedArea);
//
//
//};
//
//ForceCalc::ForceCalc(DA *octDA, const std::vector<TREENODE> &treePart, const VecInfo &nsSolution,
//                     const DomainExtents &domain, const DENDRITE_REAL Re, SubDomainBoundary *subDomainBoundary)
//        :SurfaceLoop(octDA,treePart,nsSolution,domain,subDomainBoundary),Re_(Re){
//    std::memset(Force_,0,sizeof(DENDRITE_REAL)*DIM*2);
//    std::memset(projectedArea_,0,sizeof(DENDRITE_REAL)*DIM);
//    this->traverse();
//}
//
//void ForceCalc::performSurfaceOperation(TALYFEMLIB::FEMElm & surfaceFe, const std::vector<TALYFEMLIB::ZEROPTV> &surfaceCoords,
//                                   const BoundarySurface &boundarySurface, const PetscScalar *values) {
////    TALYFEMLIB::PrintInfo("Inside the surface operation");
//    if((boundarySurface.boundaryType == BoundaryTypes::CIRCLE) or (boundarySurface.boundaryType == BoundaryTypes::SPHERE) or boundarySurface.boundaryType ==BoundaryTypes::FUNCTION){
//        while (surfaceFe.next_itg_pt()) {
//
//            TALYFEMLIB::PrintInfo("Going inside the surface loop");
//
//            calcForceAtEachQuadraturePoint2D(surfaceFe,values);
//        }
//    }
//
//}
//void ForceCalc::calcForceAtEachQuadraturePoint2D(const TALYFEMLIB::FEMElm &fe, const PetscScalar * values)
//{
//    using namespace TALYFEMLIB;
//    const int nsd = DIM;
//    double Coe_diff = 1.0 / Re_;
//    ZeroMatrix<double> du;
//    du.redim(nsd, nsd);
//    DENDRITE_REAL vald[NSNodeData::NS_DOF*DIM];
//    DENDRITE_REAL valc[NSNodeData::NS_DOF];
//    calcValueDerivativeFEM(fe,NSNodeData::NS_DOF,values,vald);
//    calcValueFEM(fe,NSNodeData::NS_DOF,values,valc);
//
//    for (int i = 0; i < nsd; i++) {
//        for (int j = 0; j < nsd; j++) {
//            du(i, j) = vald[i*DIM+j] + vald[j*DIM + i];
//        }
//    }
//    ZEROPTV diff;
//    ZEROPTV surface_normal = fe.surface()->normal();
//
//    for (int i = 0; i < nsd; i++) {
//        diff(i) = 0;
//        for (int j = 0; j < nsd; j++) {
//            //diff(i) += Coe_diff*du(i,
//            // j)*fe.normal(j);
//            diff(i) += Coe_diff * du(i, j) * surface_normal(j);
//        }
//    }
//    DENDRITE_REAL Force[DIM*2];
//    const DENDRITE_REAL & p = valc[NSNodeData::PRESSURE];
////    std::cout << fe.surface()->normal() << " " << p << "\n";
//    for (int j = 0; j < nsd; j++) {
//        //Force(j) = (-p*fe.normal(j) + diff(j))*detJxW;
//        Force[j*2 + 0] = (-p * surface_normal(j)) * fe.detJxW();
//        Force[j*2 + 1] = (diff(j)) * fe.detJxW();
//        projectedArea_[j]  += fabs(fe.detJxW()*surface_normal(j));
//
//    }
//    for (int d = 0; d < DIM; d++){
//        Force_[2*d + 0] += Force[2*d+0];
//        Force_[2*d + 1] += Force[2*d+1];
//    }
//}
//
//void ForceCalc::getTotalForce(DENDRITE_REAL * globalForce, DENDRITE_REAL * projectedArea){
//    MPI_Reduce(Force_,globalForce,DIM*2,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
//    MPI_Reduce(projectedArea_,projectedArea,DIM,MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
//
//}
//


//enum AnalyticType:short{
//    L2ERROR=0,
//    EVALFUNC=1
//};

//class Analytic: Traversal{
//    DENDRITE_REAL * L2error_ = nullptr;
//    DENDRITE_REAL *val_c_ = nullptr;
//    AnalyticFunction analFunction_;
//    AnalyticType analyticType;
//    DENDRITE_REAL time_;
//    DENDRITE_UINT sdof_;
//
//    void calculateL2error(DENDRITE_REAL * globalError);
//public:
//    /**
//     * @brief Constructor
//     * @param octDA
//     * @param v vector
//     * @param f analytic function
//     * @param time time
//     */
//    Analytic(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v,  const AnalyticFunction & f, const DomainExtents & domain, const DENDRITE_REAL time = 0);
//
//    /**
//     * @brief Overriding the travesal class operation
//     * @param fe
//     * @param values
//     */
//    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
//
//    /**
//     * @brief Prints the  L2 error
//     */
//    void getL2error();
//
//    /**
//     * @brief returns the L2 error
//     * @param error Retruns the L2 error. This must be allocated.
//     */
//    void getL2error(DENDRITE_REAL * error);
//
//    ~Analytic();
//};
//
//
//
//Analytic::Analytic(DA *octDA, const std::vector<TREENODE> & treePart,const VecInfo & v, const AnalyticFunction &f,  const DomainExtents & domain, const DENDRITE_REAL time)
//        :Traversal(octDA,treePart,v,domain){
//    analyticType = AnalyticType::L2ERROR;
//    analFunction_ = f;
//    L2error_ = new DENDRITE_REAL[v.ndof];
//    memset(L2error_,0, sizeof(DENDRITE_REAL)*v.ndof);
//    val_c_ = new DENDRITE_REAL[v.ndof];
//    time_ = time;
//    sdof_ = v.nodeDataIndex;
//    this->traverse();
//
//}
//
//void Analytic::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {
//
//    const DENDRITE_UINT ndof = this->getNdof();
//    while (fe.next_itg_pt()) {
//        calcValueFEM(fe, ndof, values, val_c_);
//        for (DENDRITE_UINT dof = 0; dof < ndof; dof++) {
//            DENDRITE_REAL val_a = analFunction_(fe.position(), dof + sdof_, time_);
//            L2error_[dof] += (val_c_[dof] - val_a) * (val_c_[dof] - val_a) * fe.detJxW();
//        }
//    }
//}
//
//void Analytic::getL2error() {
//    DENDRITE_REAL * globalL2Eror = new DENDRITE_REAL[this->getNdof()];
//    calculateL2error(globalL2Eror);
//    if (not(TALYFEMLIB::GetMPIRank())) {
//        std::cout << "L2Error : ";
//        for (int dof = 0; dof < this->getNdof(); dof++) {
//            std::cout << std::scientific << globalL2Eror[dof] << " ";
//        }
//        std::cout << std::defaultfloat << "\n";
//    }
//
//    delete [] globalL2Eror;
//
//}
//
//void Analytic::getL2error(DENDRITE_REAL *error) {
//    calculateL2error(error);
//
//}
//
//Analytic::~Analytic() {
//    if(analyticType == AnalyticType::L2ERROR){
//        delete [] L2error_;
//        delete [] val_c_;
//
//    }
//}
//
//void Analytic::calculateL2error(DENDRITE_REAL *globalL2Error) {
//    MPI_Reduce(L2error_,globalL2Error,this->getNdof(),MPI_DOUBLE,MPI_SUM,0,this->m_octDA->getCommActive());
//    if(TALYFEMLIB::GetMPIRank() == 0){
//        for(int dof = 0; dof < this->getNdof(); dof++) {
//            globalL2Error[dof] = sqrt(globalL2Error[dof]);
//        }
//    }
//}
//
//




#endif //LINNS_FORCECALC_H




//#endif //DENDRITEKT_NSPOST_H
