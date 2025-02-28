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
#include "CheckSurface.h"
#include <IMGA/Marker.h>
#include "IMGA/IMGASolverUtils.h"
#include "IMGALoop.h"
#include "OptimalSurrogateGpLoop.h"
#include "GetPoint2Interpolate.h"
#include "GetPoint2InterpolateIBM.h"
#include "BoundaryCalc_ExtraInter.h"
#include "Extrapolation.h"
#endif

#include "NSEquation.h"

#include "DARefine.h"

using namespace PETSc;

struct OptimalSurrogateStruct {
    std::vector<ZEROPTV> SurrogateGPpos, SurrogateGPnormal, TrueGPnormalGlobal,DISTANCEVEC;
    PointCloud<double> TrueGPpos_PointCloud, SurrogateGPpos_PointCloud;
    std::vector<double> TrueGPAreaGlobal;

    // Member function to clear all contents
    void clear() {
        SurrogateGPpos.clear();
        SurrogateGPnormal.clear();
        TrueGPnormalGlobal.clear();
        TrueGPAreaGlobal.clear();
        TrueGPpos_PointCloud.pts.clear();
        SurrogateGPpos_PointCloud.pts.clear();
    }
};


void performRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                       DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData,
                       const TimeInfo& ti){
    while(true) {
        Refinement *refine;
        SubDomainBoundary subDomainBoundary(&subDomain, octDA, domainInfo);
        if (inputData.InletBCType == NSInputData::LDC) {
            refine = new SDARefine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, &inputData, &subDomainBoundary, ti);
        } else {
            refine = new NSRefine(octDA, distTree.getTreePartFiltered(), domainInfo, &inputData, &subDomainBoundary,ti);
        }

        DA *newDA = refine->getRefineSubDA(distTree, 0.01,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateInByOut);
        if(newDA == nullptr){
            newDA = refine->getForceRefineSubDA(distTree);
            std::swap(newDA, octDA);
            break;
        }
        std::swap(octDA, newDA);
        delete newDA;
        subDomain.finalize(octDA, distTree.getTreePartFiltered(), domainInfo);
    }

}


void performNoChangeRefinement(DA *& octDA, DistTREE & distTree, DomainExtents &domainInfo,
                               DENDRITE_UINT maxLevel, SubDomain & subDomain, NSInputData &inputData,const TimeInfo& ti){
    // Use this function to perform no refinement but to remove the physical boundary octants in case of
    // retain inside case
//    SDARefine refine(octDA, distTree.getTreePartFiltered(),domainInfo,maxLevel);
    SubDomainBoundary subDomainBoundary(&subDomain, octDA, domainInfo);

    SDARefine refine(octDA, distTree.getTreePartFiltered(), domainInfo, maxLevel, &inputData, &subDomainBoundary,ti);
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

int checkElement(std::bitset<ElementMarker::MAX_ELMENT_TYPE>& elementMarker_alone) {
    if (elementMarker_alone.test(ElementMarker::SBM_FALSE_INTERCEPTED)) {
        return ElementMarker::SBM_FALSE_INTERCEPTED; // Returns 6
    } else if (elementMarker_alone.test(ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED)) {
        return ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED; // Returns 7
    } else if (elementMarker_alone.test(ElementMarker::IN_ELEMENT)) {
        return ElementMarker::IN_ELEMENT; // Returns 0
    } else if (elementMarker_alone.test(ElementMarker::OUT_ELEMENT)) {
        return ElementMarker::OUT_ELEMENT; // Returns 1
    } else if (elementMarker_alone.test(ElementMarker::INTERCEPTED_ELEMENT)) {
        return ElementMarker::INTERCEPTED_ELEMENT; // Returns 2
    } else if (elementMarker_alone.test(ElementMarker::IN_GP)) {
        return ElementMarker::IN_GP; // Returns 3
    } else if (elementMarker_alone.test(ElementMarker::OUT_GP)) {
        return ElementMarker::OUT_GP; // Returns 4
    } else if (elementMarker_alone.test(ElementMarker::INTERCEPTED_GP)) {
        return ElementMarker::INTERCEPTED_GP; // Returns 5
    }

    return -1; // No flag set or flag not checked
}

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
            printMarker[i] = static_cast<double>(/*elementMarker[i].to_ulong()*/checkElement(elementMarker[i]));
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

    const char *varname[] = {"marker"};
    std::vector<double> printMarker(elementMarker.size());
    for (int i = 0; i < elementMarker.size(); i++) {
//        std::cout << "checkElement(elementMarker[i]) = " << checkElement(elementMarker[i]) << "\n";
//        std::cout << "elementMarker[i].carved_geomsto_ulong() = " << static_cast<double>(elementMarker[i].to_ulong()) << "\n";
        printMarker[i] = static_cast<double>(/*elementMarker[i].to_ulong()*/checkElement(elementMarker[i]));
    }
    IO::writeVecTopVtu(octDA, distTree.getTreePartFiltered(), printMarker.data(), "Marker", "marker", varname, domainExtents, true);

    generateNeighborsOfFalseIntercepted(octDA, distTree, domainExtents, subDomain, elementMarker, nodalFalseElement,
                                            idata, subDomainBoundary, imga,
                                            true);

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
            case CarvedOutGeom::Type::DeepTrace:{
                TALYFEMLIB::PrintStatus("ALREADY ATTACHED");
            }
//
//                Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "test");
//
//                // Set up session options
//                Ort::SessionOptions session_options;
//                session_options.SetIntraOpNumThreads(1);
//
//
//
//                // Create a session with the loaded model
//                Ort::Session session(env, "model.onnx", session_options);
//
//                // Define names of the input and output tensors
//                const char *input_name = "input"; // Replace with the actual input name obtained from the model
//                const char *output_name = "output"; // Replace with the actual output name obtained from the model
//
//
//
//
//
//                const auto geomDef = [&](const double *coords) {
//                    float x = coords[0];
//                    float y = coords[1];
//                    float z = 0.0; // Initialize z with 0.0 by default
//#if (DIM == 3)
//                    z = coords[2];
//#endif
//
//                    // Prepare input tensor
//                    std::vector<float> input_tensor_values = {x, y, z};
//                    std::vector<int64_t> input_tensor_shape = {1, 3}; // Assuming the model expects a [1,3] shape tensor
//
//                    // Create memory info to describe the allocation info
//                    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
//
//                    // Create input tensor
//                    Ort::Value input_tensor = Ort::Value::CreateTensor<float>(
//                            memory_info, input_tensor_values.data(), input_tensor_values.size(), input_tensor_shape.data(),
//                            input_tensor_shape.size());
//
//                    // Perform the inference
//
//                    auto output_tensors = session.Run(Ort::RunOptions{nullptr}, &input_name, &input_tensor, 1, &output_name, 1);
//
//                    float *floatarr = output_tensors.front().GetTensorMutableData<float>();
//
//                    float dist = 0;
//
//                    if (inputData.stlRetainSide == RetainSide::IN) {
//                        dist = floatarr[0];
//                    } else {
//                        dist = -floatarr[0];
//                    }
//
//                    if (dist > 0) {
//                        return ibm::Partition::IN;
//                    }
//                    return ibm::Partition::OUT;
//                };
//
//
//
//                subDomain.addObject(geomDef);
//                break;
//            }
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
                case CarvedOutGeom::Type::DeepTrace:{
                        TALYFEMLIB::PrintStatus("ALREADY ATTACHED");
                        break;

                }
#endif
        }
    }

//#ifndef DEEPTRACE
    // Add carved geometries to the subdomain
    for (const auto &c : geomInfo.carved_geoms) {
        subDomain.addObject(c);
    }
//#endif
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
        BoundaryDistanceData &boundaryDistanceData,
        std::unordered_map<std::string, std::vector<double>> &queryToIndexMap) {
    TalyEquation<NSEquation, NSNodeData> *nsEq;

    bool loopOverSurrogateBoundary = (imga->getIBMMethod() == SBM);

    nsEq = new TalyEquation<NSEquation, NSNodeData>(&talyMesh, octDA, dTree.getTreePartFiltered(),
                                                    subDomain.domainExtents(), NSNodeData::NS_DOF,
                                                    &ti, loopOverSurrogateBoundary, subDomainBoundary, &inputData);

#ifdef IBM
#ifdef DEEPTRACE
    nsEq->equation()->setSbmCalc(imga);

    if (imga->getIBMMethod() == SBM) {
//        nsEq->equation()->setDistanceArrays(boundaryDistanceData);
        nsEq->equation()->setQueryToIndexMap(queryToIndexMap);
    }
#else
    nsEq->equation()->setSbmCalc(imga, kd_trees);

    if (imga->getIBMMethod() == SBM) {
//        nsEq->equation()->setDistanceArrays(boundaryDistanceData);
        nsEq->equation()->setQueryToIndexMap(queryToIndexMap);
    }
#endif
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

void SetElementMarkerForEqs(const NSInputData& inputData,
                            TalyEquation<NSEquation, NSNodeData> *nsEq,
                            const IMGA *imga,
                            DA *octDA,
                            const DistTREE & dTree,
                            std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>>& elementMarkers_withFalseIntercepted,
                            const Vec &nodalFalseElement,
                            const DomainExtents& domainExtents,
                            const SubDomain& subDomain,
                            SubDomainBoundary *subDomainBoundary) {
    if (inputData.ImmersedMethodType == NSInputData::ShiftedBoundary) {
        // Optimal SBM setup
        nsEq->assignIBMConstructs(imga, elementMarkers_withFalseIntercepted.data(), NSNodeData::NODE_ID);

        CheckSurface checkSurface(octDA, dTree.getTreePartFiltered(), {VecInfo(nodalFalseElement, 1, NSNodeData::NODE_ID)}, domainExtents, elementMarkers_withFalseIntercepted, &subDomain, subDomainBoundary);
        checkSurface.correctCycles();
    } else if (inputData.ImmersedMethodType == NSInputData::ImmersedBoundary) {
        // Immersed Boundary setup
        nsEq->assignIBMConstructs(imga, elementMarkers_withFalseIntercepted.data());
    }
}


/**
 * dirichlet nodes mean the nodal point that do not belong to any "active" elements (we do not infill volume GPs on that elements)
 * @param octDA
 * @param treePart
 * @param inputData
 * @param elementMarker
 * @return
 */
std::vector<PetscInt> GetIBMDirichletNodes(
        const DA* octDA,
        const DistTREE & treePart,
        const NSInputData& inputData,
        const Marker* elementMarker) {

    std::vector<PetscInt> dirichletNodes;

    if (inputData.ImmersedMethodType == NSInputData::ImmersedBoundary) {
        TALYFEMLIB::PrintStatus("We set Dirichlet nodes for IBM. Please note that!");
        getIBMDirichletNodes(octDA, treePart.getTreePartFiltered(), dirichletNodes, elementMarker);
    }

    return dirichletNodes;
}

void NSBoundaryErrorCalcIBM(const NSInputData *idata, const IMGA *imga, IMGALoop *imga_loop, const TimeInfo &ti)
{
    std::vector<double> NSBoundaryError; // resize in imga_loop->NSBoundaryError

    /// [todo] different geometries, different boundary errors
    imga_loop->NSBoundaryError(NSBoundaryError);
    // print to screen
    for (int dim = 0;dim<DIM;dim++) {
        PrintStatus("NS boundary error:", dim," direction = ", NSBoundaryError[dim]);
    }

    // write to file
    if (TALYFEMLIB::GetMPIRank() == 0)
    {
        std::string fname = "NSBoundaryError.dat";
        FILE *fp = fopen(fname.c_str(), "a");
        if (!fp)
        {
            throw TALYException() << "Cannot create file: " << fname;
        }
#if (DIM == 3)
        fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Ex = % .10e, Ey = % .10e, Ez = % .10e\n",
                ti.getTimeStepNumber(),
                ti.getCurrentTime(),
                NSBoundaryError[0], NSBoundaryError[1], NSBoundaryError[2]);
#endif
#if (DIM == 2)
        fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Ex = % .10e, Ey = % .10e\n",
                ti.getTimeStepNumber(),
                ti.getCurrentTime(),
                NSBoundaryError[0], NSBoundaryError[1]);
#endif
        fclose(fp);
    }

}

void ToEachProcessor(const std::vector<TALYFEMLIB::ZEROPTV> pts_position, std::vector<TALYFEMLIB::ZEROPTV> &pts_position_all)
{
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = pts_position.size();
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

    pts_position_all.resize(totalProcData);

    MPI_Datatype ZEROPTVtype;
    MPI_Type_contiguous(3, MPI_DOUBLE, &ZEROPTVtype);
    MPI_Type_commit(&ZEROPTVtype);
    MPI_Allgatherv(pts_position.data(), pts_position.size(), ZEROPTVtype, pts_position_all.data(), eachProcData.data(), disp.data(), ZEROPTVtype, MPI_COMM_WORLD);
}

void ToEachProcessor(const std::vector<double>& pts_position, std::vector<double>& pts_position_all)
{
    int nProc = TALYFEMLIB::GetMPISize(); // be careful
    int numNodes = pts_position.size();
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

    pts_position_all.resize(totalProcData);

    MPI_Allgatherv(pts_position.data(), pts_position.size(), MPI_DOUBLE, pts_position_all.data(), eachProcData.data(), disp.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
void DeepTrace2PointData(std::vector<ZEROPTV>&Surrogate,std::vector<ZEROPTV>&DISTVEC,std::vector<ZEROPTV> &location,std::vector<ZEROPTV> &normal, std::vector<double> &Area,DENDRITE_REAL area_)
{
//    location.resize(Surrogate.size());
    location.resize(Surrogate.size());
    normal.resize(Surrogate.size());
    Area.resize(Surrogate.size());

    for (int i=0;i<Surrogate.size();i++){
        location[i].x()  = Surrogate[i].x()-DISTVEC[i].x();
        location[i].y()  = Surrogate[i].y()-DISTVEC[i].y();
        normal[i].x()  = DISTVEC[i].x();
        normal[i].y()  = DISTVEC[i].y();
        Area[i] = area_;

#if(DIM==2)
        location[i].z()  = 0;
        normal[i].z()  = 0;
#endif
#if(DIM==2)
        location[i].z()  = Surrogate[i].z()-DISTVEC[i].z();
        normal[i].z()  = DISTVEC[i].z();
#endif
    }

}
void Imga2PointData(const IMGA *imga, std::vector<ZEROPTV> &location,std::vector<ZEROPTV> &normal, std::vector<double> &Area)
{
    location.resize(imga->getSurfaceGaussPoints().size());
    normal.resize(imga->getSurfaceGaussPoints().size());
    Area.resize(imga->getSurfaceGaussPoints().size());

    for (int i=0;i<imga->getSurfaceGaussPoints().size();i++){
        location[i].x()  = imga->getSurfaceGaussPoints()[i].location[0];
        location[i].y()  = imga->getSurfaceGaussPoints()[i].location[1];
        normal[i].x()  = imga->getSurfaceGaussPoints()[i].normal[0];
        normal[i].y()  = imga->getSurfaceGaussPoints()[i].normal[1];
        Area[i] = imga->getSurfaceGaussPoints()[i].elemArea;

#if(DIM==2)
        location[i].z()  = 0;
            normal[i].z()  = 0;
#endif
#if(DIM==2)
        location[i].z()  =  imga->getSurfaceGaussPoints()[i].location[2];
            normal[i].z()  = imga->getSurfaceGaussPoints()[i].normal[2];
#endif
    }
}

void ZEROPTV2PointCloud(const std::vector<ZEROPTV> *zeroptv, PointCloud<double> &pointcloud)
{
    pointcloud.pts.resize(zeroptv->size());
    for (int i = 0;i<zeroptv->size();i++)
    {
        pointcloud.pts[i].x = zeroptv->at(i).x();
        pointcloud.pts[i].y = zeroptv->at(i).y();
#if (DIM==2)
        pointcloud.pts[i].z = 0;
#endif
#if (DIM==3)
        pointcloud.pts[i].z = zeroptv->at(i).z();
#endif
    }
}


OptimalSurrogateStruct PreProcessSBM(NSInputData& inputData, DA* octDA, DistTREE& dTree,
                                     Vec nodalFalseElement, const DomainExtents& domainExtents,
                                     const SubDomain& subDomain, SubDomainBoundary* subDomainBoundary,
                                     IMGA* imga, IMGA *imga_TruePt_InterpPT, const std::vector<my_kd_tree_ptr> &kd_trees,
                                     std::vector<std::bitset<ElementMarker::MAX_ELMENT_TYPE>>  eleMarkers, const TimeInfo *ti){
    OptimalSurrogateStruct optimalSurrogateData;

    // Your logic here
    if (inputData.ImmersedMethodType == NSInputData::ShiftedBoundary) {


        std::vector<my_kd_tree_t*> kd_trees_;
        std::transform(kd_trees.begin(), kd_trees.end(), std::back_inserter(kd_trees_),
                       [](const auto& ptr) { return ptr.get(); });

        OptimalSurrogateGpLoop optimalSurrogateGpLoop(octDA, dTree.getTreePartFiltered(),
                                                      {VecInfo(nodalFalseElement, 1, NSNodeData::NODE_ID)},
                                                      domainExtents, eleMarkers, &subDomain,
                                                      subDomainBoundary, &inputData, imga, kd_trees_);

        /// TODO: DEBUG
//        optimalSurrogateGpLoop.WriteOptimalGPToFile();
#ifdef DEEPTRACE
        optimalSurrogateGpLoop.GetDistanceVector(optimalSurrogateData.DISTANCEVEC);
#endif

        optimalSurrogateGpLoop.GetPos(optimalSurrogateData.SurrogateGPpos);
        optimalSurrogateGpLoop.GetNormal(optimalSurrogateData.SurrogateGPnormal);
//        optimalSurrogateGpLoop.WriteOptimalGPDistToFile();
//        optimalSurrogateGpLoop.write_distance_vector();
            std::vector<ZEROPTV> SurrogateGPposGlobal;
            ToEachProcessor(optimalSurrogateData.SurrogateGPpos, SurrogateGPposGlobal);

            std::vector<ZEROPTV> TrueGPpos, TrueGPnormal, TrueGPposGlobal;
            std::vector<double> TrueGPArea;
            static DENDRITE_REAL area_ = computeOctantArea(domainExtents,inputData.meshDef.refineLevelBoundary);
//#ifdef DEEPTRACE
//            DeepTrace2PointData(optimalSurrogateData.SurrogateGPpos,optimalSurrogateData.DISTANCEVEC,TrueGPpos, TrueGPnormal, TrueGPArea,area_);
//#else
            Imga2PointData(imga, TrueGPpos, TrueGPnormal, TrueGPArea);
//#endif
            ToEachProcessor(TrueGPpos, TrueGPposGlobal);
            ToEachProcessor(TrueGPnormal, optimalSurrogateData.TrueGPnormalGlobal);
            ToEachProcessor(TrueGPArea, optimalSurrogateData.TrueGPAreaGlobal);

            ZEROPTV2PointCloud(&SurrogateGPposGlobal, optimalSurrogateData.SurrogateGPpos_PointCloud);

            my_kd_tree_t kd_tree_SurrogateGP_temp(3 /*dim*/, optimalSurrogateData.SurrogateGPpos_PointCloud, {10 /* max leaf */});

            /// TODO:DEBUG
//        std::ofstream fout("Point_step"+std::to_string(ti->getTimeStepNumber())+".txt", std::ios::app);
//        for (auto pt: optimalSurrogateData.SurrogateGPpos_PointCloud.pts) {
//#if (DIM == 2)
//            fout << pt.x << "," << pt.y
//                 << "\n";
//#endif
//        }
//
//        fout.close();


        // Using Surrogate GPs to move points a little bit inside the True Intercepted Elements
            std::vector<TALYFEMLIB::ZEROPTV> positionOfpointToInterpolateOn, positionOnSurface;
            GetPoint2Interpolate getPoint2Interpolate(octDA, imga, dTree.getTreePartFiltered(),
                                                      {VecInfo(nodalFalseElement, 1,
                                                               NSNodeData::NODE_ID)},
                                                      domainExtents,
                                                      &subDomain, &inputData, eleMarkers, subDomainBoundary,
                                                      SurrogateGPposGlobal, &kd_tree_SurrogateGP_temp);
            getPoint2Interpolate.WriteInterGPToFile();
            // NOTE: you need to call above before bottom
            getPoint2Interpolate.GetpositionOfpointToInterpolateOnAndpositionOnSurface(
                    positionOfpointToInterpolateOn,
                    positionOnSurface);

            /*
             * positionOfpointToInterpolateOn: the normal-based searched point inside the Active domain
             * positionOnSurface: the points on the true boundary
             */
            imga_TruePt_InterpPT->initIMGAComputation(octDA, dTree.getTreePartFiltered(), positionOfpointToInterpolateOn,
                                                      positionOnSurface);

            ZEROPTV2PointCloud(&TrueGPposGlobal, optimalSurrogateData.TrueGPpos_PointCloud);
    }
    else {
        std::vector<TALYFEMLIB::ZEROPTV> positionOfpointToInterpolateOn, positionOnSurface;

        GetPoint2InterpolateIBM getPoint2InterpolateIbm(octDA, imga, dTree.getTreePartFiltered(),
                                                  {VecInfo(nodalFalseElement, 1,
                                                           NSNodeData::NODE_ID)},
                                                  domainExtents,
                                                  &subDomain, &inputData, eleMarkers, subDomainBoundary);
        getPoint2InterpolateIbm.WriteInterGPToFile();
        // NOTE: you need to call above before bottom
        getPoint2InterpolateIbm.GetpositionOfpointToInterpolateOnAndpositionOnSurface(
                positionOfpointToInterpolateOn,
                positionOnSurface);

        for (int i= 0 ; i< positionOfpointToInterpolateOn.size();i++){
//            std::cout << "positionOfpointToInterpolateOn = " << positionOfpointToInterpolateOn[i].x() << " " << positionOfpointToInterpolateOn[i].y() << " " << positionOfpointToInterpolateOn[i].z() << "\n";
//            std::cout << "positionOnSurface = " << positionOnSurface[i].x() << " " << positionOnSurface[i].y() << " " << positionOnSurface[i].z() << "\n";
//
            // calculate the difference

//            if (fabs(positionOfpointToInterpolateOn[i].x() - positionOnSurface[i].x())>0) {
//                std::cout << "difference = " << positionOfpointToInterpolateOn[i].x() - positionOnSurface[i].x() << " "
//                          << positionOfpointToInterpolateOn[i].y() - positionOnSurface[i].y() << " " <<
//                          positionOfpointToInterpolateOn[i].z() - positionOnSurface[i].z() << "\n";
//            }
        }

        /*
         * positionOfpointToInterpolateOn: the normal-based searched point inside the Active domain
         * positionOnSurface: the points on the true boundary
         */
        imga_TruePt_InterpPT->initIMGAComputation(octDA, dTree.getTreePartFiltered(), positionOfpointToInterpolateOn,
                                                  positionOnSurface);

        std::vector<ZEROPTV> TrueGPpos, TrueGPnormal, TrueGPposGlobal;
        std::vector<double> TrueGPArea;
        Imga2PointData(imga, TrueGPpos, TrueGPnormal, TrueGPArea);
        ToEachProcessor(TrueGPpos, TrueGPposGlobal);
        ToEachProcessor(TrueGPnormal, optimalSurrogateData.TrueGPnormalGlobal);
        ToEachProcessor(TrueGPArea, optimalSurrogateData.TrueGPAreaGlobal);
        ZEROPTV2PointCloud(&TrueGPposGlobal, optimalSurrogateData.TrueGPpos_PointCloud);

    }


    return optimalSurrogateData;
}

//void ProcessBoundaryAndErrorCalculations(DA * octDA, IMGA *imga, DistTREE distTree,
//                                         PETSc::NonlinearSolver* nsSolver, DomainExtents domain, NSInputData &inputData,
//        TimeInfo ti, IMGA *imga_TruePt_InterpPT,
//        my_kd_tree_t *kd_tree_TrueGP, my_kd_tree_t *kd_tree_SurrogateGP,
//        OptimalSurrogateStruct optimalSurrogateData) {
//
//    if (inputData.ImmersedMethodType == NSInputData::ShiftedBoundary) {
//        if (fabs(inputData.RatioGPSBM - 1) < 1e-14) {
//            IMGALoop loop(octDA, imga, distTree.getTreePartFiltered(),
//                          {VecInfo(nsSolver->getCurrentSolution(), NSNodeData::NS_DOF, NSNodeData::VEL_X,
//                                   PLACEHOLDER_NONE)},
//                          domain, inputData, ti);
//            loop.WriteGPValue();
//            NSBoundaryErrorCalcIBM(&inputData, imga, &loop, ti);
//        } else {
//            BoundaryCalc_ExtraInter boundaryCalc(octDA, imga, imga_TruePt_InterpPT, distTree.getTreePartFiltered(),
//                                                 {VecInfo(nsSolver->getCurrentSolution(),
//                                                          NSNodeData::NS_DOF, NSNodeData::VEL_X, PLACEHOLDER_NONE)},
//                                                 domain, inputData, ti, /*point_map,*/
//                                                 kd_tree_TrueGP,
//                                                 kd_tree_SurrogateGP,
//                                                 optimalSurrogateData.TrueGPnormalGlobal,
//                                                 optimalSurrogateData.TrueGPAreaGlobal);
//            boundaryCalc.WriteGPBoundaryError();
//            boundaryCalc.WriteNSBoundaryError();
//        }
//    } else{
//
//        Extrapolation extrapolation(octDA, imga, imga_TruePt_InterpPT, distTree.getTreePartFiltered(),
//                                   {VecInfo(nsSolver->getCurrentSolution(), NSNodeData::NS_DOF, NSNodeData::VEL_X,
//                                            PLACEHOLDER_NONE)},
//                                   domain, ti,&inputData);
//
//        extrapolation.WriteGPBoundaryError();
//        extrapolation.WriteNSBoundaryError();
//        extrapolation.WriteTrueGPToFile();
//    }
//}

void forceCalcIBM(const NSInputData *idata, const IMGA *imga, IMGALoop *imga_loop, const TimeInfo &ti)
{
    std::vector<ForceInfo> globalForces;
    globalForces.resize(imga->getGeometries().size());
    imga_loop->computeForce(globalForces);
    // print to screen
    for (int i = 0; i < imga->getGeometries().size(); i++)
    {
        PrintStatus(idata->carved_out_geoms_def[i].name, " force: ", globalForces[i].Force_all(true));
        PrintStatus(idata->carved_out_geoms_def[i].name, " torque: ", globalForces[i].Torque);
    }
    // write to file
    for (int i = 0; i < imga->getGeometries().size(); i++)
    {
        if (TALYFEMLIB::GetMPIRank() == 0)
        {
            std::string fname = "Force_" + std::to_string(i) + "_" + idata->carved_out_geoms_def[i].name + ".dat";
            FILE *fp = fopen(fname.c_str(), "a");
            if (!fp)
            {
                throw TALYException() << "Cannot create file: " << fname;
            }
#if (DIM == 3)
            fprintf(fp,
                "Timestep: %1d, Time: %.5f\n"
                "Fx_pre = % .10e, Fy_pre = % .10e, Fz_pre = % .10e\n"
                "Fx_vis = % .10e, Fy_vis = % .10e, Fz_vis = % .10e\n"
                "Fx_pen = % .10e, Fy_pen = % .10e, Fz_pen = % .10e\n"
                "mass_conservation = % .10e\n"
                "Tx = % .10e, Ty = % .10e, Tz = % .10e\n",
                ti.getTimeStepNumber(),
                ti.getCurrentTime(),
                globalForces[i].Force_pressure.x(), globalForces[i].Force_pressure.y(), globalForces[i].Force_pressure.z(),
                globalForces[i].Force_viscous.x(), globalForces[i].Force_viscous.y(), globalForces[i].Force_viscous.z(),
                globalForces[i].Force_penalty.x(), globalForces[i].Force_penalty.y(), globalForces[i].Force_penalty.z(),
                globalForces[i].mass_conservation,
                globalForces[i].Torque.x(), globalForces[i].Torque.y(), globalForces[i].Torque.z());
#endif
#if (DIM == 2)
            fprintf(fp,
                    "Timestep: %1d, Time: %.5f\n"
                    "Fx_pre = % .10e, Fy_pre = % .10e\n"
                    "Fx_vis = % .10e, Fy_vis = % .10e\n"
                    "Fx_pen = % .10e, Fy_pen = % .10e\n"
                    "mass_conservation = % .10e\n"
                    "Tz = % .10e\n",
                    ti.getTimeStepNumber(),
                    ti.getCurrentTime(),
                    globalForces[i].Force_pressure.x(), globalForces[i].Force_pressure.y(),
                    globalForces[i].Force_viscous.x(), globalForces[i].Force_viscous.y(),
                    globalForces[i].Force_penalty.x(), globalForces[i].Force_penalty.y(),
                    globalForces[i].mass_conservation,
                    globalForces[i].Torque.z());
#endif
            fclose(fp);
        }
    }
}

void PostProcess_Cd_BoundaryError_OptSurrogate(NSInputData& inputData, DA* octDA, DistTREE& distTree,
                                                 PETSc::NonlinearSolver* nsSolver,
                                                 const DomainExtents& domain, const TimeInfo& ti,
                                                 IMGA* imga, IMGA *imga_TruePt_InterpPT,
                                                 my_kd_tree_t *kd_tree_TrueGP,
                                                 my_kd_tree_t *kd_tree_SurrogateGP,
                                                 const std::vector<my_kd_tree_ptr> &kd_trees,
                                                 OptimalSurrogateStruct &optimalSurrogateData,
                                                 std::unordered_map<std::string, uint32_t> & queryToIndexMap) {
#ifndef DEEPTRACE
    if (inputData.carved_out_geoms_def.empty()) {
        return; // Skip if no carved out geometries
    }
#endif

    auto vi = VecInfo(nsSolver->getCurrentSolution(),
                      NSNodeData::NS_DOF, NSNodeData::VEL_X, PLACEHOLDER_NONE);

    if (inputData.ImmersedMethodType == NSInputData::ShiftedBoundary) {
        if (fabs(inputData.RatioGPSBM - 1) < 1e-14) {
            IMGALoop loop(octDA, imga, distTree.getTreePartFiltered(),
                          vi,
                          domain, inputData, ti);
            loop.WriteGPValue();
            NSBoundaryErrorCalcIBM(&inputData, imga, &loop, ti);

            forceCalcIBM(&inputData, imga, &loop, ti);
        } else {

            BoundaryCalc_ExtraInter boundaryCalc(octDA, imga, imga_TruePt_InterpPT, distTree.getTreePartFiltered(), vi,
                                                 domain, inputData, ti, /*point_map,*/
                                                 kd_tree_TrueGP,
                                                 kd_tree_SurrogateGP,
                                                 kd_trees,
                                                 optimalSurrogateData.TrueGPnormalGlobal,
                                                 optimalSurrogateData.TrueGPAreaGlobal, queryToIndexMap);
//            boundaryCalc.WriteGPBoundaryError();
//            boundaryCalc.WriteNSBoundaryError();

            /// debug
//            boundaryCalc.WriteTrueGPToFile();

            if (queryToIndexMap.empty()) {
                boundaryCalc.gainQueryToIndexMap(queryToIndexMap);
            }
            boundaryCalc.forceCalcIBM(); // write Forces over time

        }
    } else { // IBM

        Extrapolation extrapolation(octDA, imga, imga_TruePt_InterpPT, distTree.getTreePartFiltered(), vi,
                                    domain, ti, &inputData,
                                    optimalSurrogateData.TrueGPnormalGlobal,
                                    kd_tree_TrueGP,
                                    queryToIndexMap);

        if (queryToIndexMap.empty()) {
            extrapolation.gainQueryToIndexMap(queryToIndexMap);
        }

//        extrapolation.WriteGPBoundaryError();
//        extrapolation.WriteNSBoundaryError();
//        extrapolation.WriteTrueGPToFile();
        extrapolation.forceCalc();
    }
}


void refineAndIntergrid(DA *&octDA,
                        DistTREE &dTree,
                        DomainExtents &domainExtents,
                        NSInputData &inputData,
                        Vec &vecNS1, Vec &vecNS2, Vec &vecNS3,
                        SubDomain &subDA,
                        bool doCoarsen,
                        DENDRITE_UINT maxLevel,
                        const TimeInfo &ti,
                        int &no_refine_return) {

    int no_refine = 0;
    while (true) {
        TALYFEMLIB::PrintStatus("Active comm = ",
                                octDA->getNpesActive(),
                                ", no_refine = ",
                                no_refine,
                                ", mesh count (node) = ",
                                octDA->getGlobalNodeSz());

        Refinement *refine;
        SubDomainBoundary subDomainBoundary(&subDA, octDA, domainExtents);

        if (inputData.InletBCType == NSInputData::LDC) {
            refine = new SDARefine(octDA, dTree.getTreePartFiltered(), domainExtents, maxLevel, &inputData, &subDomainBoundary, ti);
        } else {
            refine = new NSRefine(octDA, dTree.getTreePartFiltered(), domainExtents, &inputData, &subDomainBoundary, ti);
        }



        auto newDA = refine->getRefineSubDA(dTree, 0.01, RefinementStrategy::FULL_DOMAIN, ot::RemeshPartition::SurrogateInByOut);
//        DA *newDA = daRefine.getRefineSubDA(dTree);

//        std::cout << "here 1\n ";
        if (newDA == NULL) {
            break;
        }
        refine->initPetscIntergridTransfer();
//        std::cout << "here 2\n ";

        refine->petscIntergridTransfer(newDA, dTree, vecNS1, NSNodeData::NS_DOF, RefinementStage::REFINE);
        refine->petscIntergridTransfer(newDA, dTree, vecNS2, NSNodeData::NS_DOF, RefinementStage::REFINE);
        refine->petscIntergridTransfer(newDA, dTree, vecNS3, NSNodeData::NS_DOF, RefinementStage::REFINE);
//        std::cout << "here 3\n ";

        refine->finializeIntergridTransfer();

//        std::cout << "here 5\n ";

        std::swap(octDA, newDA);
        delete newDA;
        subDA.finalize (octDA, dTree.getTreePartFiltered (), domainExtents);
        no_refine++;
    }

    no_refine_return= no_refine;

#if (DIM == 3)
    static const char *ns_varname[]{"u", "v", "w", "p"};
#endif
#if (DIM == 2)
    static const char *ns_varname[]{"u", "v", "p"};
#endif

//    PETSc::petscVectopvtu(octDA,
//                          dTree.getTreePartFiltered(),
//                          vecNS1,
//                          "afterRefine",
//                          "afterRefineNS",
//                          ns_varname,
//                          domainExtents,
//                          false,
//                          false,
//                          NSNodeData::NS_DOF);

}

static void perfromNoChangeRefinement(DA *& octDA,DistTREE & dTree, DomainExtents & domainExtents,
                                      NSInputData * inputData, SubDomain *subdomain,  std::vector<VecInfo> & vecs,
                                      const TimeInfo &ti) {

    SubDomainBoundary subDomainBoundary(subdomain, octDA, domainExtents);

    SDARefine refine(octDA, dTree.getTreePartFiltered(), domainExtents, 8 /*we do not use this number inside refine*/, inputData,
                     &subDomainBoundary,ti);
    DA * newDA = refine.getForceRefineSubDA(dTree,0.1);

    refine.initPetscIntergridTransfer();
    for (int i =0;i<vecs.size();i++) {
        refine.petscIntergridTransfer(newDA, dTree, vecs[i].v, vecs[i].ndof);
    }
    refine.finializeIntergridTransfer();

    std::swap(newDA, octDA);
    delete newDA;
}


/**
 * From mehdish/samk just added it here
 * @brief Interleave two PETSc vectors based on a DOF (Degree of Freedom) pattern.
 *
 * This function interleaves two PETSc vectors, `vec1` and `vec2`, into a single interleaved vector.
 * The interleaving is based on the provided DOF pattern.
 *
 * @param vec1 The first PETSc vector to be interleaved.
 * @param vec2 The second PETSc vector to be interleaved.
 * @param interleavedVec Pointer to the resulting interleaved PETSc vector.
 * @param comm The MPI communicator.
 * @param dof1 The degree of freedom for the first vector.
 * @param dof2 The degree of freedom for the second vector.
 */
void interleaveVecs(Vec vec1, Vec vec2, Vec *interleavedVec, MPI_Comm comm, PetscInt dof1, PetscInt dof2) {
    PetscInt n1, n2, N1, N2;
    VecGetLocalSize(vec1, &n1);


    VecGetLocalSize(vec2, &n2);
    VecGetSize(vec1, &N1);
    VecGetSize(vec2, &N2);

    if (n1/dof1 != n2/dof2) {
        PetscPrintf(comm, "Error: Local sizes of vec1/dof1 and vec2/dof2 are not the same!\n");
        return;
    }

    PetscInt n_interleaved = n1 + n2;
    PetscInt N_interleaved = N1 + N2;

    // Create the interleaved vector
    VecCreateMPI(comm, n_interleaved, N_interleaved, interleavedVec);

    // Get pointers to the data arrays
    PetscScalar *array1, *array2, *array_interleaved;
    VecGetArray(vec1, &array1);
    VecGetArray(vec2, &array2);
    VecGetArray(*interleavedVec, &array_interleaved);

    // Interleave the data
    PetscInt i1 = 0, i2 = 0, i_interleaved = 0;
    while (i1 < n1 || i2 < n2) {
        for (PetscInt d = 0; d < dof1 && i1 < n1; ++d) {
            array_interleaved[i_interleaved++] = array1[i1++];
        }
        for (PetscInt d = 0; d < dof2 && i2 < n2; ++d) {
            array_interleaved[i_interleaved++] = array2[i2++];
        }
    }

    // Restore the arrays
    VecRestoreArray(vec1, &array1);
    VecRestoreArray(vec2, &array2);
    VecRestoreArray(*interleavedVec, &array_interleaved);
}
#endif //DENDRITEKT_NSUTILS_H
