//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_NSUTILS_H
#define DENDRITEKT_NSUTILS_H

#include <DataTypes.h>
#include <NSRefine.h>
#include <SDARefine.h>
#include <PETSc/Solver/NonLinearSolver.h>
#include <NSInputData.h>
#include <NSNodeData.h>
#include <Traversal/Analytic.h>

#ifdef IBM
#include "SBMMarker.h"
#include "BFS.h"
#endif

#include "NSEquation.h"

using namespace PETSc;

//void performRefinement(DA *&octDA, const DomainExtents &domainInfo, std::vector<TREENODE> &treeNode, const NSInputData &inputData) {
//  while (true) {
//    NSRefine refine(octDA, treeNode,domainInfo, treeNode.size(), inputData.meshDef.refineLevelBoundary);
//    DA *newDA = refine.getRefineDA(treeNode);
//    if (newDA == NULL) {
//      break;
//    }
//    std::swap(newDA, octDA);
//    delete newDA;
//  }
//}


void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                       DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData){
    while(true) {
//        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel);
        SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, &inputData);

        DA *newDA = refine.getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            newDA = refine.getForceRefineSubDA(distTree);
            std::swap(newDA, octDA);
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }

}


void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
//    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel);
    SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, &inputData);
    DA * newDA = refine.getForceRefineSubDA(distTree,0.01);
    std::swap(octDA, newDA);
    delete newDA;
    subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);

}



static double AnalyticalSolution(const TALYFEMLIB::ZEROPTV &pt, int dof, const DENDRITE_REAL time) {
#if(DIM == 3)
    std::cout << "Analytical solution not supported \n";
    return 0;
#else

    if (dof == NSNodeData::VEL_X) {
        return (sin(M_PI * pt.x()) * cos(M_PI * pt.y()) * sin(2 * M_PI * (time)) + 2);
    }
    if (dof == NSNodeData::VEL_Y) {
        return (-cos(M_PI * pt.x()) * sin(M_PI * pt.y()) * sin(2 * M_PI * (time)) + 2);
    }
    if (dof == NSNodeData::PRESSURE) {
        return (sin(M_PI * pt.x()) * sin(M_PI * pt.y()) * cos(2 * M_PI * (time)) + 2);
    }
#endif
};
void setInitialCondition(DA *octDA, NSInputData &inputData, Vec &Solution) {
    std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
        if (not(inputData.ifMMS)) {
            var[NSNodeData::VEL_X] = 0;
            var[NSNodeData::VEL_Y] = 0;
#if(DIM == 3)
            var[NSNodeData::VEL_Z] = 0;
#endif
            var[NSNodeData::PRESSURE] = 0;
        } else {
            var[NSNodeData::VEL_X] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_X, 0.0);
            var[NSNodeData::VEL_Y] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::VEL_Y, 0.0);
            var[NSNodeData::PRESSURE] = AnalyticalSolution(TALYFEMLIB::ZEROPTV{x[0], x[1], 0.0}, NSNodeData::PRESSURE, 0.0);
        }
    };
    octDA->petscSetVectorByFunction(Solution, initial_condition, false, false, NSNodeData::NS_DOF);
}

#ifdef IBM

void generateNeighborsOfFalseIntercepted(DA *& octDA, DistTREE & distTree, DomainExtents & domainExtents, SubDomain & subDomain,
                                         std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>> & elementMarker,
                                         Vec & nodalFalseElement, NSInputData *idata,SubDomainBoundary *subDomainBoundary,const IMGA *imga
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
        BFS bfs(octDA, distTree.getTreePartFiltered(), elementMarker,domainExtents,VecInfo(nodalFalseElement, 1, 0), &subDomain,idata,
                subDomainBoundary, imga);

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
                        Vec &nodalFalseElement , NSInputData *idata){
    int RelativeOrder = 3; // TODO: hard code here
    SBMMarker marker(octDA, distTree.getTreePartFiltered(), domainExtents, elementMarker, &subDomain,
                     RelativeOrder, idata);
    octDA->petscCreateVector(nodalFalseElement, false, false, 1);

    if (imga->getIBMMethod() == SBM) {
        generateNeighborsOfFalseIntercepted(octDA, distTree, domainExtents, subDomain, elementMarker, nodalFalseElement,
                                            idata, subDomainBoundary, imga,
                                            true);
    }

}

#if (DIM ==2)
/**
 *
 * @param geom_def [IN]
 * @param mshs [OUT]
 * @param carved_geoms [OUT]
 * @param ibm_refinements [OUT]
 * @param inputData [IN]
 */
void addCarvedOutandIBMMSH2D(const CarvedOutGeom &geom_def, std::vector<GEOMETRY::MSH *> &mshs,
                             std::vector<GEOMETRY::Geometry *> &carved_geoms, std::vector<GeomRefinement> &ibm_refinements,
                             const NSInputData *inputData)
{
    if (geom_def.type == CarvedOutGeom::Type::MESHOBJECT_2D)
    {
        mshs.push_back(new GEOMETRY::MSH(geom_def.mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D));

        std::array<DENDRITE_REAL, DIM> shift{geom_def.InitialDisplacement[0],
                                             geom_def.InitialDisplacement[1]};
        auto geom_retain_side = RetainSide::OUT;
        if (geom_def.outer_boundary) {
            geom_retain_side = RetainSide::IN;
        }
        carved_geoms.push_back(new GEOMETRY::Geometry(mshs.back(), Point<DIM>(shift), geom_retain_side));
        ibm_refinements.emplace_back(geom_def.geomRefine);
    }
}
#endif

#if(DIM==3)
void addCarvedOutandIBMSTL3D(const CarvedOutGeom &geom_def, std::vector<GEOMETRY::STL *> &stls,
                             std::vector<GEOMETRY::Geometry *> &carved_geoms, std::vector<GeomRefinement> &ibm_refinements)
{
    if (geom_def.type == CarvedOutGeom::Type::MESHOBJECT)
    {
        stls.push_back(new GEOMETRY::STL(geom_def.mesh_path, GEOMETRY::InOutTest::RAY_TRACING));
        stls.back()->correctNormals();
        std::array<DENDRITE_REAL, DIM> shift{geom_def.InitialDisplacement[0],
                                             geom_def.InitialDisplacement[1],
                                             geom_def.InitialDisplacement[2]};

        // todo add multiple shift
        auto geom_retain_side = RetainSide::OUT;
        if (geom_def.outer_boundary) {
            geom_retain_side = RetainSide::IN;
        }
        carved_geoms.emplace_back(new GEOMETRY::Geometry(stls.back(), Point<DIM>(shift),geom_retain_side));
        ibm_refinements.emplace_back(geom_def.geomRefine);
    }
}
#endif

/**
 *
 * @param inputData [IN]
 * @param subDomain [OUT]
 * @param imga [OUT]
 * @param geomInfo [OUT]
 */
void InitializeGeometricObjects(NSInputData &inputData, SubDomain &subDomain, IMGA *imga, GeomInfo &geomInfo) {
#ifdef IBM
    std::vector<GeomRefinement> ibm_refinements;
#endif
    for (const auto &geom_def : inputData.carved_out_geoms_def) {
        switch (geom_def.type) {
#if (DIM == 2)
            case CarvedOutGeom::Type::MESHOBJECT_2D: {
#ifdef IBM
                addCarvedOutandIBMMSH2D(geom_def, geomInfo.mshs, geomInfo.carved_geoms, ibm_refinements, &inputData);
#else
                addCarvedOutMSH2D(geom_def, mshs, carved_geoms);
#endif
                break;
            }

            case CarvedOutGeom::Type::CIRCLE_2D_VOXEL: {
                const double coordsCenter[DIM]{geom_def.center_of_mass.x(), geom_def.center_of_mass.y()};
                subDomain.addObject(VOXEL::Circle(coordsCenter, geom_def.radius));
                TALYFEMLIB::PrintError("[NOTE] setting circle like this do not support IMGA staff!");
                break;
            }

            case CarvedOutGeom::Type::BOX_2D_VOXEL:{
                const double min[DIM]{geom_def.cube_dim[0].x(), geom_def.cube_dim[0].y()};
                const double max[DIM]{min[0] + geom_def.cube_dim[1].x(), min[1] + geom_def.cube_dim[1].y()};
                subDomain.addObject(VOXEL::Box(min, max));
                break;
            }
#endif

#if (DIM == 3)
                case CarvedOutGeom::Type::MESHOBJECT: {
#ifdef IBM
                    addCarvedOutandIBMSTL3D(geom_def, geomInfo.stls, geomInfo.carved_geoms, ibm_refinements);
#else
                    addCarvedOutSTL3D(geomInfo.geom_def, geomInfo.stls, carved_geoms);
#endif
                    break;
                }

                case CarvedOutGeom::Type::SPHERE_3D_VOXEL: {
                    const double coordsCenter[DIM]{geom_def.center_of_mass.x(), geom_def.center_of_mass.y(),
                                                   geom_def.center_of_mass.z()};
                    subDomain.addObject(VOXEL::Sphere(coordsCenter, geom_def.radius));
                    break;
                }

                case CarvedOutGeom::Type::CUBE_3D_VOXEL: {
                    const double min[DIM]{geom_def.cube_dim[0].x(), geom_def.cube_dim[0].y(), geom_def.cube_dim[0].z()};
                    const double max[DIM]{min[0] + geom_def.cube_dim[1].x(),
                                          min[1] + geom_def.cube_dim[1].y(),
                                          min[2] + geom_def.cube_dim[1].z()};
                    subDomain.addObject(VOXEL::Box(min, max));
                    break;
                }
#endif
        }
    }


    // Add carved geometries to the subdomain
#ifndef DEEPTRACE
    for (const auto &c : geomInfo.carved_geoms) {
        subDomain.addObject(c);
    }
#endif

#ifdef IBM
    // IBM specific geometry addition
    for (size_t i = 0; i < geomInfo.carved_geoms.size(); i++) {
        /// [todo] removing ibmGeom structure latter
        imga->addGeometry(geomInfo.carved_geoms.at(i), ibm_refinements.at(i));
    }
#endif

}

#endif

TalyEquation<NSEquation, NSNodeData>* InitializeEquations(
        NSInputData& inputData,
        TalyMesh<NSNodeData>& talyMesh,
        DA *octDA,
        SubDomainBoundary *subDomainBoundary,
        const DistTREE & dTree,
        const SubDomain& subDomain,
        const TimeInfo& ti,
        const IMGA* imga,
        const std::vector<my_kd_tree_ptr> &kd_trees,
        BoundaryDistanceData &boundaryDistanceData) {
    TalyEquation<NSEquation, NSNodeData> *nsEq;

    bool loopOverSurrogateBoundary = (imga->getIBMMethod() != NITSCHE);

    nsEq = new TalyEquation<NSEquation, NSNodeData>(&talyMesh, octDA, dTree.getTreePartFiltered(),
                                                    subDomain.domainExtents(), NSNodeData::NS_DOF,
                                                    &ti, loopOverSurrogateBoundary, subDomainBoundary, &inputData);

#ifdef IBM
    nsEq->equation()->setSbmCalc(imga, kd_trees);

    if (imga->getIBMMethod() == SBM) {
        nsEq->equation()->setDistanceArrays(boundaryDistanceData);
    }
#endif

    return nsEq;
}


void GetLineEdgePt(const std::vector<GEOMETRY::MSH *> &mshs,std::vector<PointCloud<double>> &CenterPts)
{

    CenterPts.resize(mshs.size());

    for (int mshID = 0; mshID < mshs.size(); mshID++)
    {
        const std::vector<GEOMETRY::Lines> *m_lines = &mshs[mshID]->getLines();

        CenterPts[mshID].pts.resize(m_lines->size());

        for (int i = 0;i<m_lines->size();i++)
        {
            CenterPts[mshID].pts[i].x = (m_lines->at(i).lineCoord[0][0]+m_lines->at(i).lineCoord[1][0])/2.0;
            CenterPts[mshID].pts[i].y = (m_lines->at(i).lineCoord[0][1]+m_lines->at(i).lineCoord[1][1])/2.0;
            CenterPts[mshID].pts[i].z = 0;
        }
    }
}

void GetTriangleCenter(const std::vector<GEOMETRY::STL *> &stls,std::vector<PointCloud<double>> &CenterPts)
{
    CenterPts.resize(stls.size());
    for (int stlID = 0; stlID < stls.size(); stlID++)
    {
        const std::vector<GEOMETRY::Triangles> *m_triangles = &stls[stlID]->getTriangles();
        //std::cout<<"m_triangles.size() = " << m_triangles.size() << "\n";
        CenterPts[stlID].pts.resize(m_triangles->size());
        for (int i = 0;i<m_triangles->size();i++)
        {
            CenterPts[stlID].pts[i].x =(m_triangles->at(i).triangleCoord[0][0] + m_triangles->at(i).triangleCoord[1][0] + m_triangles->at(i).triangleCoord[2][0]) / 3;
            CenterPts[stlID].pts[i].y =(m_triangles->at(i).triangleCoord[0][1] + m_triangles->at(i).triangleCoord[1][1] + m_triangles->at(i).triangleCoord[2][1]) / 3;
            CenterPts[stlID].pts[i].z =(m_triangles->at(i).triangleCoord[0][2] + m_triangles->at(i).triangleCoord[1][2] + m_triangles->at(i).triangleCoord[2][2]) / 3;
        }
    }
}

/**
 * generate points to contruct k-d tree
 * @param geomInfo
 * @return
 */
std::vector<PointCloud<double>>  SetupPointCloud( const GeomInfo &geomInfo) {

    std::vector<PointCloud<double>> CenterPts;

    // Setup for distance function calculation in SBM
#if (DIM == 3)
    GetTriangleCenter(geomInfo.stls, CenterPts);
#endif

#if (DIM == 2)
    GetLineEdgePt(geomInfo.mshs, CenterPts);
#endif
    return CenterPts;
}

#endif //DENDRITEKT_NSUTILS_H
