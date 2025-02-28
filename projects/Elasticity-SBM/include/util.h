#pragma once
#include "csv_reader.h"
#include <iostream>
#include <DendriteUtils.h>
#include <point.h>
#include <TalyEquation.h>
#include <PETSc/Solver/LinearSolver.h>
#include <PETSc/PetscUtils.h>
#include <IO/VTU.h>
#include <Traversal/Analytic.h>
#include <PETSc/IO/petscVTU.h>
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>
#include <Checkpoint/Checkpointer.h>
#include <talyfem/input_data/input_data.h>
#include <talyfem/talyfem.h>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
#include <map>
#include <vector>
#include <string>
#include <PETSc/VecBounds.h>
#include "LEInputData.h"
#include "LERefine.h"
#include <IMGA/IMGA.h>
#include <IMGA/Marker.h>
#include "nanoflann.hpp"
#include "BFS.h"
#include "ElementMarker.h"

/// SBM
#pragma mark OptSug
//#include "BFS.h"
//#include "ElementMarker.h"

// where QUAD is !
#define VTK_QUAD 9


struct TractionData {
    ZEROPTV pos;
    ZEROPTV traction;
};

// construct a kd-tree index:
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
PointCloud<double>, 3 /* dim */>;


namespace util_funcs
{
    void generateNeighborsOfFalseIntercepted(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain,
                                             std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,
                                             Vec & nodalFalseElement, LEInputData *idata,SubDomainBoundary *subDomainBoundary,const IMGA *imga
            ,  bool isAllocated = false){
        using CoordT = typename ot::DA<DIM>::C;
        using ot::RankI;
        OctToPhysical octToPhysical(domainExtents);
        const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
        const auto &treeNodes = distTree.getTreePartFiltered();
        std::vector<PetscInt> nodeIDs(treeNodes.size() * nPe, -1);

        // Get Global element ID
        {
            const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
            const size_t ghostedNodalSz = octDA->getTotalNodalSz();
            const TREENODE *odaCoords = octDA->getTNCoords();
            const bool visitEmpty = false;
            const unsigned int padLevel = 0;
            ot::MatvecBaseIn<DIM, RankI, false> treeLoopIn(ghostedNodalSz,
                                                           1,                // node id is scalar
                                                           octDA->getElementOrder(),
                                                           visitEmpty,
                                                           padLevel,
                                                           odaCoords,
                                                           &(*ghostedGlobalNodeId.cbegin()),
                                                           &(*treeNodes.cbegin()),
                                                           treeNodes.size(),
                                                           *octDA->getTreePartFront(),
                                                           *octDA->getTreePartBack());

            int eleCounter = 0;
            while (!treeLoopIn.isFinished()) {
                const ot::TreeNode<CoordT, DIM> subtree = treeLoopIn.getCurrentSubtree();
                const auto subtreeInfo = treeLoopIn.subtreeInfo();
                if (treeLoopIn.isPre() && subtreeInfo.isLeaf()) {
                    const RankI *nodeIdsFlat = subtreeInfo.readNodeValsIn();
                    const auto &octCoords = subtreeInfo.getNodeCoords();
                    const std::vector<bool> &nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
                    for (int i = 0; i < nPe; i++) {
                        //if (nodeNonhangingIn[i]) {
                        nodeIDs[eleCounter * nPe + i] = nodeIdsFlat[i];
                        //}
                    }
                    eleCounter++;
                }
                treeLoopIn.step();
            }
        }
        // Get Global node ID ends

        // Transfer cell information to nodes
        if(!isAllocated) {
            octDA->petscCreateVector(nodalFalseElement, false, false, 1);
        }
        VecSet(nodalFalseElement,0);
        for(int i = 0; i < treeNodes.size();i++){
            if(elementMarker[i].test(ElementMarker::SBM_FALSE_INTERCEPTED)){
                for(int n = 0; n < nPe; n++){
                    VecSetValue(nodalFalseElement,nodeIDs[nPe*i + n],NodeMarker::SBM_FALSE_INTERCEPTED_NODES,INSERT_VALUES);
                }
            }
        }
        VecAssemblyBegin(nodalFalseElement);
        VecAssemblyEnd(nodalFalseElement);
        // Transfer ends

        // Do one iteration of BFS to find the neighbors of false intercepted
        {
//            BFS(DA * octDA, const std::vector<TREENODE> &treePart, std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & eleMarkers,
//            const DomainExtents &domain, const VecInfo & v, const SubDomain * subDomain, const LEInputData *idata);
            BFS bfs(octDA, distTree.getTreePartFiltered(), elementMarker,domainExtents,VecInfo(nodalFalseElement, 1, 0), &subDomain,idata );

            std::vector<double> printMarker(elementMarker.size());
            for(int i = 0; i < elementMarker.size();i++){
                printMarker[i] = (double )(elementMarker[i].to_ulong());
            }
            const char *varname[] = {"marker"};
            IO::writeVecTopVtu(octDA,distTree.getTreePartFiltered(),printMarker.data(),"MarkerNeighbors","marker",varname,domainExtents,true);

        }
    }

    void generateNewMarkers(DA * octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain,
                            SubDomainBoundary *subDomainBoundary,IMGA *imga,
                            std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,
                            Vec &nodalFalseElement , LEInputData *idata){
        int RelativeOrder = 3; // TODO: hard code here
        SBMMarker marker(octDA, distTree.getTreePartFiltered(), domainExtents, elementMarker, &subDomain,
                         RelativeOrder, idata);
        octDA->petscCreateVector(nodalFalseElement, false, false, 1);

        generateNeighborsOfFalseIntercepted(octDA, distTree, domainExtents, subDomain, elementMarker, nodalFalseElement,
                                            idata, subDomainBoundary, imga,
                                            true);

        const char *varname[] = {"marker"};
        std::vector<double> printMarker(elementMarker.size());
        for (int i = 0; i < elementMarker.size(); i++) {
            printMarker[i] = static_cast<double>(elementMarker[i].to_ulong());
        }
        IO::writeVecTopVtu(octDA, distTree.getTreePartFiltered(), printMarker.data(), "SBMMarker", "marker", varname, domainExtents, true);

    }


    void GetTriangleCenter(const std::vector<GEOMETRY::STL *> &stls,PointCloud<double> &CenterPts)
    {
        for (int stlID = 0; stlID < stls.size(); stlID++)
        {
            const std::vector<GEOMETRY::Triangles> *m_triangles = &stls[stlID]->getTriangles();
            //std::cout<<"m_triangles.size() = " << m_triangles.size() << "\n";
            CenterPts.pts.resize(m_triangles->size());
            for (int i = 0;i<m_triangles->size();i++)
            {
                CenterPts.pts[i].x =(m_triangles->at(i).triangleCoord[0][0] + m_triangles->at(i).triangleCoord[1][0] + m_triangles->at(i).triangleCoord[2][0]) / 3;
                CenterPts.pts[i].y =(m_triangles->at(i).triangleCoord[0][1] + m_triangles->at(i).triangleCoord[1][1] + m_triangles->at(i).triangleCoord[2][1]) / 3;
                CenterPts.pts[i].z =(m_triangles->at(i).triangleCoord[0][2] + m_triangles->at(i).triangleCoord[1][2] + m_triangles->at(i).triangleCoord[2][2]) / 3;
            }
        }
    }

    double ReturnPenaltyParameters(LEInputData *idata_)
    {

#if (DIM == 2)
        if (idata_->caseType == CaseType::PLANESTRESS)
    {
      return idata_->planeStress.young/ (1- pow(idata_->planeStress.poisson,2));
    }
#endif

        if (idata_->caseType == CaseType::PLANESTRAIN)
        {
            return idata_->planeStrain.young/ (1- pow(idata_->planeStrain.poisson,2));
        }
        if (idata_->caseType == CaseType::LAME)
        {
            double lamda = idata_->lame.lamda;
            double mu = idata_->lame.mu;
            return mu * (3 * lamda + 2 * mu) / (mu + lamda) / (1- pow(lamda / 2 / (lamda + mu),2));
        }
    }


    void CalcCmatrix(LEInputData *idata_, std::vector<std::vector<double>> &Cmatrix)
    {
        for (auto &row : Cmatrix)
        {
            row.resize(3 * (DIM - 1));
        }

        /*
         * 3D do not have plane stress case
         */
#if (DIM == 2)
        if (idata_->caseType == CaseType::PLANESTRESS)
    {
      double young = idata_->planeStress.young;
      double poisson = idata_->planeStress.poisson;
      // C for plane stress
      Cmatrix[0][0] = young / (1 - pow(poisson, 2));
      Cmatrix[0][1] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 - pow(poisson, 2));
      Cmatrix[1][1] = young / (1 - pow(poisson, 2));
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (2 * (1 + poisson));
    }
#endif

        if (idata_->caseType == CaseType::PLANESTRAIN)
        {
            double young = idata_->planeStrain.young;
            double poisson = idata_->planeStrain.poisson;
            // C for plane strain
#if (DIM == 2)
            Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = young / (1 + poisson) / 2;
#endif
#if (DIM == 3)
            /*
             * this formulation is from FEM book page 241
             */
            Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][3] = 0;
            Cmatrix[0][4] = 0;
            Cmatrix[0][5] = 0;

            Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][3] = 0;
            Cmatrix[1][4] = 0;
            Cmatrix[1][5] = 0;

            Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][3] = 0;
            Cmatrix[2][4] = 0;
            Cmatrix[2][5] = 0;

            Cmatrix[3][0] = 0;
            Cmatrix[3][1] = 0;
            Cmatrix[3][2] = 0;
            Cmatrix[3][3] = young / (1 + poisson) / 2;
            Cmatrix[3][4] = 0;
            Cmatrix[3][5] = 0;

            // [previous bug here]
            Cmatrix[4][0] = 0;
            Cmatrix[4][1] = 0;
            Cmatrix[4][2] = 0;
            Cmatrix[4][3] = 0;
            Cmatrix[4][4] = young / (1 + poisson) / 2;
            Cmatrix[4][5] = 0;

            Cmatrix[5][0] = 0;
            Cmatrix[5][1] = 0;
            Cmatrix[5][2] = 0;
            Cmatrix[5][3] = 0;
            Cmatrix[5][4] = 0;
            Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
        }
        if (idata_->caseType == CaseType::LAME)
        {
            double lamda = idata_->lame.lamda;
            double mu = idata_->lame.mu;
            // C for lame parameters
#if (DIM == 2)
            Cmatrix[0][0] = lamda + 2 * mu;
      Cmatrix[0][1] = lamda;
      Cmatrix[0][2] = 0;
      Cmatrix[1][0] = lamda;
      Cmatrix[1][1] = lamda + 2 * mu;
      Cmatrix[1][2] = 0;
      Cmatrix[2][0] = 0;
      Cmatrix[2][1] = 0;
      Cmatrix[2][2] = mu;
#endif
#if (DIM == 3)
            double young = mu * (3 * lamda + 2 * mu) / (mu + lamda);
            double poisson = lamda / 2 / (lamda + mu);

            Cmatrix[0][0] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[0][3] = 0;
            Cmatrix[0][4] = 0;
            Cmatrix[0][5] = 0;

            Cmatrix[1][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][1] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][2] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[1][3] = 0;
            Cmatrix[1][4] = 0;
            Cmatrix[1][5] = 0;

            Cmatrix[2][0] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][1] = young * poisson / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][2] = young * (1 - poisson) / (1 + poisson) / (1 - 2 * poisson);
            Cmatrix[2][3] = 0;
            Cmatrix[2][4] = 0;
            Cmatrix[2][5] = 0;

            Cmatrix[3][0] = 0;
            Cmatrix[3][1] = 0;
            Cmatrix[3][2] = 0;
            Cmatrix[3][3] = young / (1 + poisson) / 2;
            Cmatrix[3][4] = 0;
            Cmatrix[3][5] = 0;

            // [previous bug here]
            Cmatrix[4][0] = 0;
            Cmatrix[4][1] = 0;
            Cmatrix[4][2] = 0;
            Cmatrix[4][3] = 0;
            Cmatrix[4][4] = young / (1 + poisson) / 2;
            Cmatrix[4][5] = 0;

            Cmatrix[5][0] = 0;
            Cmatrix[5][1] = 0;
            Cmatrix[5][2] = 0;
            Cmatrix[5][3] = 0;
            Cmatrix[5][4] = 0;
            Cmatrix[5][5] = young / (1 + poisson) / 2;
#endif
        }
    }

#pragma mark OptSug
//    void generateNeighborsOfFalseIntercepted(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain,
//                                             std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,
//                                             Vec & nodalFalseElement, const LEInputData *idata,  bool isAllocated = false){
//        using CoordT = typename ot::DA<DIM>::C;
//        using ot::RankI;
//        OctToPhysical octToPhysical(domainExtents);
//        const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
//        const auto &treeNodes = distTree.getTreePartFiltered();
//        std::vector<PetscInt> nodeIDs(treeNodes.size() * nPe, -1);
//
//        // Get Global element ID
//        {
//            const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
//            const size_t ghostedNodalSz = octDA->getTotalNodalSz();
//            const TREENODE *odaCoords = octDA->getTNCoords();
//            const bool visitEmpty = false;
//            const unsigned int padLevel = 0;
//            ot::MatvecBaseIn<DIM, RankI, false> treeLoopIn(ghostedNodalSz,
//                                                           1,                // node id is scalar
//                                                           octDA->getElementOrder(),
//                                                           visitEmpty,
//                                                           padLevel,
//                                                           odaCoords,
//                                                           &(*ghostedGlobalNodeId.cbegin()),
//                                                           &(*treeNodes.cbegin()),
//                                                           treeNodes.size(),
//                                                           *octDA->getTreePartFront(),
//                                                           *octDA->getTreePartBack());
//
//            int eleCounter = 0;
//            while (!treeLoopIn.isFinished()) {
//                const ot::TreeNode<CoordT, DIM> subtree = treeLoopIn.getCurrentSubtree();
//                const auto subtreeInfo = treeLoopIn.subtreeInfo();
//                if (treeLoopIn.isPre() && subtreeInfo.isLeaf()) {
//                    const RankI *nodeIdsFlat = subtreeInfo.readNodeValsIn();
//                    const auto &octCoords = subtreeInfo.getNodeCoords();
//                    const std::vector<bool> &nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
//                    for (int i = 0; i < nPe; i++) {
//                        //if (nodeNonhangingIn[i]) {
//                        nodeIDs[eleCounter * nPe + i] = nodeIdsFlat[i];
//                        //}
//                    }
//                    eleCounter++;
//                }
//                treeLoopIn.step();
//            }
//        }
//        // Get Global node ID ends
//
//        // Transfer cell information to nodes
//        if(!isAllocated) {
//            octDA->petscCreateVector(nodalFalseElement, false, false, 1);
//        }
//        VecSet(nodalFalseElement,0);
//        for(int i = 0; i < treeNodes.size();i++){
//            if(elementMarker[i].test(ElementMarker::SBM_FALSE_INTERCEPTED)){
//                for(int n = 0; n < nPe; n++){
//                    VecSetValue(nodalFalseElement,nodeIDs[nPe*i + n],NodeMarker::SBM_FALSE_INTERCEPTED_NODES,INSERT_VALUES);
//                }
//            }
//        }
//        VecAssemblyBegin(nodalFalseElement);
//        VecAssemblyEnd(nodalFalseElement);
//        // Transfer ends
//
//        // Do one iteration of BFS to find the neighbors of false intercepted
//        {
//            BFS bfs(octDA, distTree.getTreePartFiltered(), elementMarker,domainExtents,VecInfo(nodalFalseElement, 1, 0), &subDomain,idata);
//
//            std::vector<double> printMarker(elementMarker.size());
//            for(int i = 0; i < elementMarker.size();i++){
//                printMarker[i] = (double )(elementMarker[i].to_ulong());
//            }
//            const char *varname[] = {"marker"};
//            IO::writeVecTopVtu(octDA,distTree.getTreePartFiltered(),printMarker.data(),"MarkerNeighbors","marker",varname,domainExtents,true);
//
//        }
//    }
//
//    void generateNewMarkers(DA * octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain
//                            , std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker
//                            , std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarkerVisualize
//                            , Vec &nodalFalseElement
//                            , const LEInputData *idata){
//        int RelativeOrder = idata->RelOrderCheckActive;
//        SBMMarker marker(octDA,distTree.getTreePartFiltered(),domainExtents,elementMarker,&subDomain,RelativeOrder,idata);
//        marker.WriteDomainGPToFile();
//        const char *varname[] = {"marker"};
//        elementMarkerVisualize = elementMarker;
//
//        generateNeighborsOfFalseIntercepted(octDA, distTree, domainExtents,  subDomain,elementMarker,nodalFalseElement,idata, false);
//        PETSc::petscVectopvtu(octDA,distTree.getTreePartFiltered(),nodalFalseElement,"MarkerNeighbors","Nodes",varname,domainExtents,
//                              false);
//
//    }

#pragma mark OptSug end
  /**
 * Save the octree mesh to vtk binary file.
 * @param da
 * @param daScalingFactor
 * @param ti
 * @param ns
 * @param ht
 * @param prefix
 */
  PetscErrorCode save_timestep(DA *octDA,
                               const std::vector<TREENODE> &treePartition,
                               Vec vec,
                               unsigned int ndof,
                               const TimeInfo &ti,
                               const SubDomain &subDomain,
                               const std::string &prefix = "timestep",
                               const char **varName = nullptr /*,
                             bool saveGeometry = true,
                             bool writeData = true,
                             const std::vector<int> &dataIndex = {},
                             const std::function<double(double v_in, bool in_geom, unsigned int dof)> &f = nullptr*/
  )
  {
    // create directory for this timestep (if it doesnt already exist)
    char folder[PATH_MAX];
    snprintf(folder, sizeof(folder), "results_%05d", ti.getTimeStepNumber());
    int ierr = mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (ierr != 0 && errno != EEXIST)
    {
      PrintError("Could not create folder for storing results (", strerror(errno), ").");
      return 1;
    }

    char fname[PATH_MAX];
    snprintf(fname, sizeof(fname), "%s_%05d", prefix.c_str(), ti.getTimeStepNumber());
    PETSc::petscVectopvtu(octDA, treePartition, vec, folder, fname, varName,
                          subDomain.domainExtents(), false, false, ndof);
  }

    void writeVecTopVtu(DA *da,
                        const std::vector<TREENODE> & treePart,
                        Vec in,
                        const char *foldername,
                        const char *fprefix,
                        const char **varName,
                        const SubDomain &subDomain,
                        const bool isElemental,
                        const bool isGhosted ,
                        const unsigned int ndof){
        if(not(da->isActive())) {
            return;
        }
        if (not(TALYFEMLIB::GetMPIRank())) {
            int ierr = mkdir(foldername, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (ierr != 0 && errno != EEXIST) {
                TALYFEMLIB::PrintError("Could not create folder for storing results (", strerror(errno), ").");
                return;
            }
        }
        MPI_Barrier(da->getCommActive());
        char fname[PATH_MAX];
        snprintf(fname, sizeof(fname), "%s/%s", foldername, fprefix);
        PETSc::petscVectopvtu(da, treePart, in, foldername, fprefix, varName,
                              subDomain.domainExtents(), false, false, ndof);
    }


    int performRefinementSubDA(DA *&octDA, DomainExtents &domainExtents, DistTREE &dTree,
                                LEInputData &inputData, SubDomain *subdomain)
    {
        int no_refine = 0;
        while (true)
        {

            SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);

            LERefine refine(octDA, dTree.getTreePartFiltered(), domainExtents, &inputData, &subDomainBoundary);
            DA *newDA = refine.getRefineSubDA(dTree);
            if (newDA == NULL)
            {
                newDA = refine.getForceRefineSubDA(dTree);
                std::swap(newDA, octDA);
                break;
            }

            std::swap(newDA, octDA);

            delete newDA;

            subdomain->finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
            TALYFEMLIB::PrintStatus("Refinement ", no_refine, ", mesh count (node) = ", octDA->getGlobalNodeSz());
            no_refine++;
        }

        return no_refine;
    }


  /*
  void performRefinementSubDAIBM(DA *&octDA, const std::vector<TREENODE> &treeNode, DomainExtents &domainExtents, DistTREE &dTree,
                                 LEInputData &inputData, SubDomain *subdomain, std::vector<GEOMETRY::Geometry *> ibm_geoms)
  {
    int no_refine = 0;
    while (no_refine < 36)
    {
      SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);
      LERefine refine(octDA, treeNode, domainExtents, &inputData, &subDomainBoundary, ibm_geoms);
      DA *newDA = refine.getRefineSubDA(dTree);
      if (newDA == NULL)
      {
        newDA = refine.getForceRefineSubDA(dTree);
        std::swap(newDA, octDA);
        break;
      }

      std::swap(newDA, octDA);

      delete newDA;
      subdomain->finalize(octDA, dTree.getTreePartFiltered(), domainExtents);
      TALYFEMLIB::PrintStatus("Refinement ", no_refine, ", mesh count (node) = ", octDA->getGlobalNodeSz());
      no_refine++;
    }
  }
*/

  void read_msh(const std::string &fileName, std::vector<double> &node) // pass by reference
  {

    int nProc = TALYFEMLIB::GetMPISize();
    int rank = TALYFEMLIB::GetMPIRank();

    std::ifstream gmsh_line2d_file(fileName.c_str(), std::ios::in);
    if (!gmsh_line2d_file)
    {
      throw std::runtime_error("ERROR: COULD NOT READ GMSH 1D MESH FILE");
      assert(false);
    }

    node = {};
    LINE2D::Line2DData info("");
    std::string line;
    std::vector<std::vector<int>> connectivity;
    std::vector<Vector2d> line_nodes;

    while (std::getline(gmsh_line2d_file, line))
    {
      std::istringstream iss(line);
      std::vector<int> int_value_temp;
      int val_temp;
      if (iss.str() == "$Nodes")
      {
        std::getline(gmsh_line2d_file, line);
        iss.clear();
        iss.str(line);
        // read header
        while ((iss >> val_temp))
        {
          int_value_temp.push_back(val_temp);
        }
        int_value_temp.clear();
        while (std::getline(gmsh_line2d_file, line))
        {
          // read end of nodes
          iss.clear();
          iss.str(line);
          if (iss.str() == "$EndNodes")
          {
            break;
          }
          // read nodes in each section
          while ((iss >> val_temp))
          {
            int_value_temp.push_back(val_temp);
          }
          int no_nodes_per_section = int_value_temp.back();
          int_value_temp.clear();
          for (int i = 0; i < no_nodes_per_section; i++)
          {
            std::getline(gmsh_line2d_file, line);
          }
          for (int i = 0; i < no_nodes_per_section; i++)
          {
            std::getline(gmsh_line2d_file, line);
            iss.clear();
            iss.str(line);
            double val_coor;
            //std::vector<double> node;
            while ((iss >> val_coor))
            {
              node.push_back(val_coor);
              //std::cout<<"val_coor:"<<val_coor<<"\n";
            }
            //std::cout<<"\n";
            //assert(node.size() == 3 and fabs(node[2]) < 1e-15);
            line_nodes.emplace_back(node[0], node[1]);
          }
        }
        //for (int i=0;i<node.size();i++){
        //  std::cout<<"node:"<<node[i]<<"\n";
        //}
      }
    }
  }

#pragma mark Function for linear-elasticity
    Vec GetMag3D(Vec U_le, IS x_is, IS y_is, IS z_is)
    {

        Vec U_x, U_y, U_z;
        int n_node = 0;
        int first;
        VecGetLocalSize(U_le, &n_node);
        VecGetOwnershipRange(U_le, &first, PETSC_NULL);

        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first, 3, &x_is);
        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first + 1, 3, &y_is);
        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first + 2, 3, &z_is);

        // ISView(x_is,PETSC_VIEWER_STDOUT_WORLD);
        int nx = n_node / 3;
        VecGetSubVector(U_le, x_is, &U_x);
        VecGetSubVector(U_le, y_is, &U_y);
        VecGetSubVector(U_le, z_is, &U_z);

        PetscScalar min_val, max_val;
        PetscInt min_loc, max_loc;

//    PetscScalar scale_factor = 1e6;
//    PetscPrintf(PETSC_COMM_WORLD, "---------------scale_factor: %g ---------------\n", scale_factor);

        // Get min and max for U_x, then scale and print
        VecMin(U_x, &min_loc, &min_val);
        VecMax(U_x, &max_loc, &max_val);
//    PetscPrintf(PETSC_COMM_WORLD, "U_x: Min value: %f at location %D, Max value: %f at location %D\n", min_val * scale_factor, min_loc, max_val * scale_factor, max_loc);

        // Get min and max for U_y, then scale and print
        VecMin(U_y, &min_loc, &min_val);
        VecMax(U_y, &max_loc, &max_val);
//    PetscPrintf(PETSC_COMM_WORLD, "U_y: Min value: %f at location %D, Max value: %f at location %D\n", min_val * scale_factor, min_loc, max_val * scale_factor, max_loc);

        // Get min and max for U_z, then scale and print
        VecMin(U_z, &min_loc, &min_val);
        VecMax(U_z, &max_loc, &max_val);
//    PetscPrintf(PETSC_COMM_WORLD, "U_z: Min value: %f at location %D, Max value: %f at location %D\n", min_val * scale_factor, min_loc, max_val * scale_factor, max_loc);

        // Find the maximum absolute value among the min and max of U_x, U_y, and U_z
        PetscScalar max_abs_value = PetscMax(PetscAbsReal(max_val), PetscAbsReal(min_val));

        // Do the same for U_y and U_z
        VecMin(U_y, &min_loc, &min_val);
        VecMax(U_y, &max_loc, &max_val);
        max_abs_value = PetscMax(max_abs_value, PetscMax(PetscAbsReal(max_val), PetscAbsReal(min_val)));

        VecMin(U_z, &min_loc, &min_val);
        VecMax(U_z, &max_loc, &max_val);
        max_abs_value = PetscMax(max_abs_value, PetscMax(PetscAbsReal(max_val), PetscAbsReal(min_val)));

        // Calculate scale_factor so that the scaled values are in the range 0~10
        // if max_abs_value = 0, keep scale_factor as 1 to avoid division by zero
        PetscScalar scale_factor = (max_abs_value > 0) ? pow(10, floor(log10(max_abs_value) + 1) - 1) : 1;
        PetscPrintf(PETSC_COMM_WORLD, "---------------scale_factor: %g ---------------\n", scale_factor);

        // Now you can scale and print the vectors using the calculated scale_factor
        VecMin(U_x, &min_loc, &min_val);
        VecMax(U_x, &max_loc, &max_val);
        PetscPrintf(PETSC_COMM_WORLD, "U_x: Min value: %f at location %D, Max value: %f at location %D\n", min_val / scale_factor, min_loc, max_val / scale_factor, max_loc);

        VecMin(U_y, &min_loc, &min_val);
        VecMax(U_y, &max_loc, &max_val);
        PetscPrintf(PETSC_COMM_WORLD, "U_y: Min value: %f at location %D, Max value: %f at location %D\n", min_val / scale_factor, min_loc, max_val / scale_factor, max_loc);

        VecMin(U_z, &min_loc, &min_val);
        VecMax(U_z, &max_loc, &max_val);
        PetscPrintf(PETSC_COMM_WORLD, "U_z: Min value: %f at location %D, Max value: %f at location %D\n", min_val / scale_factor, min_loc, max_val / scale_factor, max_loc);


        VecPow(U_x, 2);
        VecPow(U_y, 2);
        VecPow(U_z, 2);
        VecAYPX(U_y, 1, U_x); // U_y=U_x+1*U_y
        VecAYPX(U_z, 1, U_y); // U_z=U_y+1*U_z
        VecSqrtAbs(U_z);

        ISDestroy(&x_is);
        ISDestroy(&y_is);
        ISDestroy(&z_is);
        VecDestroy(&U_x);
        VecDestroy(&U_y);

        return U_z;
    }

    void GetVec(const Vec &U_le, IS &x_is, IS &y_is, Vec &U_x, Vec &U_y)
    {
        int n_node = 0;
        int first;
        VecGetLocalSize(U_le, &n_node);
        VecGetOwnershipRange(U_le, &first, PETSC_NULL);

        ISCreateStride(PETSC_COMM_WORLD, n_node / 2, first, 2, &x_is);
        ISCreateStride(PETSC_COMM_WORLD, n_node / 2, first + 1, 2, &y_is); // PETSC_COMM_SELF=>XXX;PETSC_COMM_SELF

        // ISView(x_is,PETSC_VIEWER_STDOUT_WORLD);
        int nx = n_node / 2;
        VecGetSubVector(U_le, x_is, &U_x);
        VecGetSubVector(U_le, y_is, &U_y);
    }

    void GetVec(const Vec &U_le, IS &x_is, IS &y_is, IS &z_is, Vec &U_x, Vec &U_y, Vec &U_z)
    {
        int n_node = 0;
        int first;
        VecGetLocalSize(U_le, &n_node);
        VecGetOwnershipRange(U_le, &first, PETSC_NULL);

        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first, 3, &x_is);
        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first + 1, 3, &y_is); // PETSC_COMM_SELF=>XXX;PETSC_COMM_SELF
        ISCreateStride(PETSC_COMM_WORLD, n_node / 3, first + 2, 3, &z_is); // PETSC_COMM_SELF=>XXX;PETSC_COMM_SELF

        // ISView(x_is,PETSC_VIEWER_STDOUT_WORLD);
        int nx = n_node / 3;
        VecGetSubVector(U_le, x_is, &U_x);
        VecGetSubVector(U_le, y_is, &U_y);
        VecGetSubVector(U_le, z_is, &U_z);
    }

    Vec GetMag(Vec U_le, IS x_is, IS y_is)
    {

        Vec U_x, U_y;
        int n_node = 0;
        int first;
        VecGetLocalSize(U_le, &n_node);
        VecGetOwnershipRange(U_le, &first, PETSC_NULL);

        ISCreateStride(PETSC_COMM_WORLD, n_node / 2, first, 2, &x_is);
        ISCreateStride(PETSC_COMM_WORLD, n_node / 2, first + 1, 2, &y_is); // PETSC_COMM_SELF=>XXX;PETSC_COMM_SELF

        // ISView(x_is,PETSC_VIEWER_STDOUT_WORLD);
        int nx = n_node / 2;
        VecGetSubVector(U_le, x_is, &U_x);
        VecGetSubVector(U_le, y_is, &U_y);

        VecPow(U_x, 2);
        VecPow(U_y, 2);
        VecAYPX(U_y, 1, U_x); // U_y=U_x+1*U_y
        VecSqrtAbs(U_y);

        ISDestroy(&x_is);
        ISDestroy(&y_is);
        VecDestroy(&U_x);

        return U_y;
    }


    std::vector<TractionData> readCSV(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Could not open the file " << filename << "!" << std::endl;
            return {}; // return an empty vector
        }

        std::vector<TractionData> dataEntries;
        std::string line;

        // Read lines from the file
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            TractionData entry;

#if (DIM==2)
            ss >> entry.pos.x() >> entry.pos.y() >> entry.traction.x() >> entry.traction.y();
#endif

#if (DIM==3)
            ss >> entry.pos.x() >> entry.pos.y() >> entry.pos.z() >> entry.traction.x() >> entry.traction.y() >> entry.traction.z();
#endif

            dataEntries.push_back(entry);
        }

        return dataEntries;
    }

    std::vector<TractionData> readTXT(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Could not open the file " << filename << "!" << std::endl;
            return {}; // return an empty vector
        }

        std::vector<TractionData> dataEntries;
        std::string line;

        // Read lines from the file
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            TractionData entry;

#if (DIM==2)
            ss >> entry.pos.x() >> entry.pos.y() >> entry.traction.x() >> entry.traction.y();
#endif

#if (DIM==3)
            ss >> entry.pos.x() >> entry.pos.y() >> entry.pos.z() >> entry.traction.x() >> entry.traction.y() >> entry.traction.z();
#endif

            dataEntries.push_back(entry);
        }

        return dataEntries;
    }


    void GetLineEdgePt(const std::vector<GEOMETRY::MSH *> &mshs,PointCloud<double> &CenterPts)
    {

        for (int mshID = 0; mshID < mshs.size(); mshID++)
        {
            const std::vector<GEOMETRY::Lines> *m_lines = &mshs[mshID]->getLines();

            CenterPts.pts.resize(m_lines->size());

            for (int i = 0;i<m_lines->size();i++)
            {
                CenterPts.pts[i].x = m_lines->at(i).lineCoord[0][0];
                CenterPts.pts[i].y = m_lines->at(i).lineCoord[0][1];
                CenterPts.pts[i].z = 0;
            }
        }
    }

    void readTXT(const std::string& filename, PointCloud<double> & traction_pos , std::vector<ZEROPTV> &traction_vector) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Could not open the file " << filename << "!" << std::endl;
        }

        std::string line;

        // Read lines from the file
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            ZEROPTV pos;
            ZEROPTV vector;

#if (DIM==2)
            ss >> pos.x() >> pos.y() >> vector.x() >> vector.y();
            traction_pos.pts.push_back({pos[0],pos[1]});
#endif

#if (DIM==3)
            ss >> pos.x() >> pos.y() >> pos.z() >> vector.x() >> vector.y() >> vector.z();
            traction_pos.pts.push_back({pos[0],pos[1],pos[2]});
#endif

            traction_vector.push_back(vector);
        }
    }

    void removeDuplicates(PointCloud<double> &traction_pos, std::vector<ZEROPTV> &traction_vector) {
        for (size_t i = 0; i < traction_pos.pts.size(); ++i) {
            for (size_t j = i + 1; j < traction_pos.pts.size(); ) {
                ZEROPTV pt_i = {traction_pos.pts[i].x,traction_pos.pts[i].y,traction_pos.pts[i].z};
                ZEROPTV pt_j = {traction_pos.pts[j].x,traction_pos.pts[j].y,traction_pos.pts[j].z};
                if (pt_i.distanceTo(pt_j) == 0) {
#ifndef NDEBUG
//                    std::cout << "erase Duplicate Point\n";
#endif
                    traction_pos.pts.erase(traction_pos.pts.begin() + j);
                    traction_vector.erase(traction_vector.begin() + j);
                } else {
                    ++j;
                }
            }
        }
    }



    double linear_interpolation_extrapolation(ZEROPTV p1, ZEROPTV p2, double value1, double value2, ZEROPTV p3) {
        double dist_p1_p2 = p1.distanceTo(p2);
        double dist_p1_p3 = p1.distanceTo(p3);
        double dist_p2_p3 = p2.distanceTo(p3);

        // Check if p3 is between p1 and p2
        if (dist_p1_p3 + dist_p2_p3 == dist_p1_p2) {
            // Interpolate the value for p3
            double weight_p1 = dist_p2_p3 / dist_p1_p2;
            double weight_p2 = 1.0 - weight_p1;
            return weight_p1 * value1 + weight_p2 * value2;
        } else {
            // Otherwise, extrapolate as in the original function
            double slope = (value2 - value1) / ((dist_p1_p3 > dist_p2_p3) ? dist_p1_p2 : -dist_p1_p2);
            return slope * dist_p1_p3 + value1;
        }
    }

    ZEROPTV crossProduct(const ZEROPTV& u, const ZEROPTV& v)  {
        ZEROPTV result;

#ifdef ENABLE_4D
        if ((u(3) != 0.0) || (v(3) != 0.0)) {
        throw TALYException() << "Cross product does not work in 4D!";
    }
#endif

        result(0) = u(1) * v(2) - u(2) * v(1);
        result(1) = u(2) * v(0) - u(0) * v(2);
        result(2) = u(0) * v(1) - u(1) * v(0);

        return result;
    }

    ZEROPTV linear_interpolation_extrapolation(ZEROPTV p1, ZEROPTV p2, ZEROPTV value1, ZEROPTV value2, ZEROPTV query_pt) {
        double dist_p1_p2 = p1.distanceTo(p2);
        double dist_p1_p3 = p1.distanceTo(query_pt);
        double dist_p2_p3 = p2.distanceTo(query_pt);

        const double EPSILON = 1e-10;
        ZEROPTV result;

        if (std::fabs(dist_p1_p3 + dist_p2_p3 - dist_p1_p2) < EPSILON) {
            // Interpolation
            double weight_p1 = dist_p2_p3 / dist_p1_p2;
            double weight_p2 = 1.0 - weight_p1;
            result = value1 * weight_p1 + value2 * weight_p2;
        } else {
            ZEROPTV vec_p1_p2 = p2 - p1;
            ZEROPTV vec_p1_p3 = query_pt - p1;
            ZEROPTV cross_product = crossProduct(vec_p1_p2, vec_p1_p3);

            if (cross_product.distanceTo(ZEROPTV{0.0, 0.0, 0.0}) < EPSILON) {
                // Extrapolation
                ZEROPTV slope = (value2 - value1) * (1.0 / ((dist_p1_p3 > dist_p2_p3) ? dist_p1_p2 : -dist_p1_p2));
                result = slope * dist_p1_p3 + value1;
            } else {
                // Weighted Average
                double total_distance = dist_p1_p3 + dist_p2_p3;
                double weight_p1 = dist_p1_p3 / total_distance;
                double weight_p2 = dist_p2_p3 / total_distance;
                result = value1 * weight_p1 + value2 * weight_p2;
            }
        }

        // Check for inf or nan in the result
        if (std::isinf(result.x()) || std::isnan(result.x()) ||
            std::isinf(result.y()) || std::isnan(result.y()) ||
            std::isinf(result.z()) || std::isnan(result.z())) {
            throw std::runtime_error("Result contains inf or nan values.");
        }

        return result;
    }

}
