//
// Created by maksbh on 5/20/20.
//

#ifndef DENDRITEKT_MATVECCOMMON_H
#define DENDRITEKT_MATVECCOMMON_H

#include <vector>
#include <sfcTreeLoop_matvec.h>
#include <PETSc/VecInfo.h>
#include "DataTypes.h"
#include "OctToPhysical.h"
#include <TalyDendroSync.h>
#include <SubDA/SubDomain.h>
#include <Boundary/SubDomainBoundary.h>
#include <DendriteUtils.h>
#ifdef IBM
#include <IMGA/IMGA.h>
#include <IMGA/IMGADataTypes.h>
#endif
/**
 * @brief Contains the functions that are needed for both matrix and vector assembly
 * @tparam NodeData
 */

enum SYNC_TYPE:bool{
  VECTOR_ONLY = false,
  ALL = true
};



template <typename NodeData>
class TalyMatVecCommon{
  const std::vector<TREENODE> & treePart_;
  ot::MatvecBase<DIM, PetscScalar> * treeloop = nullptr;

  TalyDendroSync sync_;
  DENDRITE_UINT eleOrder_;

  PetscScalar * ghostedVecs;
  bool finalized = false;
  const DomainInfo & fullDAdomain_;
  const DomainInfo & physicalDomain_;

  OctToPhysical octToPhysical;
  DENDRITE_UINT nPe_;

  SubDomainBoundary * domainBoundary_;

  DENDRITE_REAL * tempSyncVector_ = nullptr;
  const int * relativeOrder_ = nullptr;
  const int * InterceptedrelativeOrder_ = nullptr;
  DENDRITE_UINT lclElemID; /// This is not a data structure of Dendro. It is just for our convenience

#ifdef IBM
  const IMGA *imga_;
  const std::bitset<ElementMarker::MAX_ELMENT_TYPE> *elementMarker_;
  std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoIt_;
  std::vector<NodeAndValues<DENDRITE_REAL>>::const_iterator gpInfoEnd_;
  bool hasLocalIBMGaussPoints = true;

  int dofSBM_NodeID_ = -1;
  std::vector<PetscScalar> sbmElementNodes_;
#endif

  /**
   * @brief copies the value to ghosted array. Called only during vector assembly.
   * @param octDA
   * @param v
   */
  void copyToGhostedArray(const DA * octDA, std::vector<VecInfo> & v);


  DENDRITE_UINT getSurfaceFlags(const TALYFEMLIB::ZEROPTV & positon, const Point<DIM>& physMinimumDomain, const Point<DIM> & physMaximumDomain);

 protected:

  std::vector<DENDRITE_UINT > m_syncDof;
  DENDRITE_UINT m_totalDof = 0;
  SYNC_TYPE m_syncType = SYNC_TYPE::ALL;


 public:

  enum CALL:u_short {
    MAT_ASSEMBLY = 0,
    VEC_ASSEMBLY = 1,
    MATVEC_ASSEMBLY = 2
  };

  TalyMatVecCommon(const std::vector<TREENODE> & treePart,SubDomainBoundary * domainBoundary,const DomainInfo & fullDAdomain,const DomainInfo &physicalDomain);
  /**
   * @brief Initializes the class and allocates the vector. Must call finalize to complete the operation
   * @param octDA octDA
   * @param vec sync Vecs
   */
  void init(DA *octDA,const std::vector<VecInfo> & vec, const SYNC_TYPE & syncType);


  /**
   * @brief Initializes the matrix / vector /matVec operations. Called in preMat / preMatVec / preComputeVec
   * @param da octDA
   * @param vec vectors for synchronizing
   * @param domainMin  minimum of domain
   * @param domainMax maximum of domain
   */
  void initMatVecOperation(const DA * octDA, std::vector<VecInfo> & vec, const CALL & call);

  /**
   * Copies the values from array to gridfield.
   * @param gf gridfield
   */
  void mapToGridField(TALYFEMLIB::GridField<NodeData> * gf);

  const PetscScalar * mapToGridField();
  /**
   * @brief Destructor
   */
  ~TalyMatVecCommon();

  /**
   * @brief finalizes the matVecOperations. Called in postMat / postMatVec / postComputeVec
   */
  void finalizeMatVecoperation();

  /**
   * @brief Finalize the matVecCommon class
   */
  void finalize();

  /**
   * @brief Octants -> Physical coordinates
   * @param OctCoords Octant coordinates
   * @param PhysCoords Physical coordinates
   */
  void convertToPhys(const double * OctCoords, double * PhysCoords,std::array<DENDRITE_REAL,DIM> & gridDX);


  DENDRITE_UINT generateSurfaceFlags(const TALYFEMLIB::ZEROPTV & position);

  void getSurfaceID(const double *physCoords,std::array<BoundarySurface,2*DIM> & isBoundaryFace);

  inline void assignRelativeOrder(const int * relativeOrder) {
    relativeOrder_ = relativeOrder;
  }

    inline void assignInterceptedRelativeOrder(const int * InterceptedrelativeOrder) {
        InterceptedrelativeOrder_ = InterceptedrelativeOrder;
    }

  inline int getRelativeOrder() {
    int relativeOrder = 0;
    if(relativeOrder_) {
      relativeOrder = (*relativeOrder_);
    }
    return relativeOrder;
  }

  inline int getInterceptedRelativeOrder() {
      int relativeOrder = 5;
      if(InterceptedrelativeOrder_) {
          relativeOrder = (*InterceptedrelativeOrder_);
      }
      return relativeOrder;
  }

  inline void finishElementOperation() {
    lclElemID++;
  }

  inline int getLocalElemID() const {
    return lclElemID;
  }


#ifdef IBM
  void assignIMGAConstructs(const IMGA * imga,const std::bitset<ElementMarker::MAX_ELMENT_TYPE> * elementMarker);



  bool isInside(const TALYFEMLIB::ZEROPTV & position) const;

  bool getLocalSurfaceGaussPoint(TALYFEMLIB::FEMElm & fe, NodeAndValues<DENDRITE_REAL> & gaussPoint);

    std::bitset<ElementMarker::MAX_ELMENT_TYPE>  getElementMarkerType() const;

  inline IBM_METHOD getIBMMethod() const{
    return imga_->getIBMMethod();
  }

  inline void assignNodeID(const int SBM_NodeID) {
      dofSBM_NodeID_ = SBM_NodeID;

      sbmElementNodes_.resize(1u << DIM,-1);

  }
#endif

};

template <typename NodeData>
TalyMatVecCommon<NodeData>::TalyMatVecCommon(const std::vector<TREENODE> & treePart,SubDomainBoundary * domainBoundary,const DomainInfo & fullDAdomain,const DomainInfo &physicalDomain)
:treePart_(treePart),fullDAdomain_(fullDAdomain),physicalDomain_(physicalDomain){
  octToPhysical.init(fullDAdomain.min,fullDAdomain.max);
  domainBoundary_ = domainBoundary;
}

template <typename NodeData>
void TalyMatVecCommon<NodeData>::init(DA *octDA,const std::vector<VecInfo> & vec,const SYNC_TYPE & syncType) {
  m_totalDof = 0;
  m_syncType = syncType;
  m_syncDof.clear();
  for(const auto & v : vec){
    m_totalDof += v.ndof;
    for(DENDRITE_UINT dof = 0; dof < v.ndof; dof++) {
      m_syncDof.push_back(v.nodeDataIndex + dof);
    }
  }
  octDA->createVector(ghostedVecs, false,true,m_totalDof);
  octDA->createVector(tempSyncVector_, false,false,m_totalDof);
}



template <typename NodeData>
void TalyMatVecCommon<NodeData>::initMatVecOperation(const DA * octDA, std::vector<VecInfo> & vec, const CALL & call) {
  lclElemID = 0;
  nPe_ = octDA->getNumNodesPerElement();
  eleOrder_ = octDA->getElementOrder();
#ifdef IBM
  if(elementMarker_ == nullptr) {
    throw std::runtime_error("Need to add IBM construct for using IBM");
  }
  const auto & gpInfo = imga_->getSurfaceGaussPoints();
  if(gpInfo.empty()) {
    hasLocalIBMGaussPoints = false;
  }
  else {
    gpInfoIt_ = gpInfo.cbegin();
    gpInfoEnd_ = gpInfo.cend();
  }
#endif
  if(vec.size() == 0){
    finalized = true;
    return;
  }
  finalized = false;

  if(call == VEC_ASSEMBLY){
    copyToGhostedArray(octDA,vec);
  }

  const size_t totalNodalSz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  if(not(((call == MAT_ASSEMBLY) or (call == MATVEC_ASSEMBLY)) and (m_syncType == VECTOR_ONLY))) {
    treeloop = new ot::MatvecBase<DIM, PetscScalar>(totalNodalSz,
                                                    m_totalDof,
                                                    eleOrder_,
                                                    tnCoords,
                                                    ghostedVecs,
                                                    &(*(treePart_.cbegin())),
                                                    treePart_.size(),
                                                    *partFront,
                                                    *partBack);
  }
}



template < typename NodeData>
TalyMatVecCommon<NodeData>::~TalyMatVecCommon() {
  finalize();
}

template < typename NodeData>
void TalyMatVecCommon<NodeData>::finalize() {
  if(m_totalDof > 0) {
    delete[] ghostedVecs;
    delete [] tempSyncVector_;
  }
  m_totalDof = 0;
}

template < typename NodeData>
void TalyMatVecCommon<NodeData>::finalizeMatVecoperation() {

  if(not(finalized)) {
    delete treeloop;
    treeloop = nullptr;
    finalized = true;
  }

}

template < typename NodeData>
const PetscScalar * TalyMatVecCommon<NodeData>::mapToGridField() {
  if(treeloop == nullptr){
    return nullptr;
  }
  while (!(treeloop[0].isFinished())){
    if ((treeloop[0].isPre()) and (treeloop[0].subtreeInfo().isLeaf())) {
      const PetscScalar *nodeValsFlat = treeloop->subtreeInfo().readNodeValsIn();
//      for(int i = 0; i < m_totalDof*nPe_; i++) {
//        std::cout << nodeValsFlat[i] << " ";
//      }
//      std::cout << "\n";
      treeloop[0].next();
      return nodeValsFlat;
    }
    else
    {
      treeloop[0].step();
    }
  }
}

template < typename NodeData>
void TalyMatVecCommon<NodeData>::mapToGridField(TALYFEMLIB::GridField<NodeData> *gf) {
  if(treeloop == nullptr){
    return;
  }
  while (!(treeloop[0].isFinished())){
    if ((treeloop[0].isPre()) and (treeloop[0].subtreeInfo().isLeaf())) {
      const PetscScalar *nodeValsFlat = treeloop->subtreeInfo().readNodeValsIn();
      sync_.syncValues<NodeData>(nodeValsFlat, gf, m_totalDof, m_syncDof);
#ifdef IBM
      if(this->imga_->getIBMMethod() == SBM) {
          gf->NodeDataToArray(sbmElementNodes_.data(), dofSBM_NodeID_);
      }

#endif
      treeloop[0].next();

      return;
    }
    else
    {
      treeloop[0].step();
    }
  }
}
template <typename NodeData>
void TalyMatVecCommon<NodeData>::convertToPhys(const double *OctCoords, double *PhysCoords, std::array<DENDRITE_REAL,DIM> & gridDX) {
  std::memcpy(PhysCoords,OctCoords, sizeof(double)*DIM*nPe_);
  for(DENDRITE_UINT i = 0; i < nPe_; i++) {
    octToPhysical.convertCoordsToPhys(&PhysCoords[DIM*i]);
  }
  int oppositeElements[DIM];
  if(eleOrder_ == 1) {
    oppositeElements[0] = 1;
    oppositeElements[1] = 2;
#if (DIM >= 3)
    oppositeElements[2] = 4;
#endif
  }
  if(eleOrder_ == 2) {
    oppositeElements[0] = 2;
    oppositeElements[1] = 6;
#if (DIM >= 3)
    oppositeElements[2] = 9;
#endif
  }

 for(int dim = 0; dim < DIM; dim++) {
   gridDX[dim] = PhysCoords[oppositeElements[dim]*DIM + dim] - PhysCoords[dim];
 }
}
template <typename NodeData>
void TalyMatVecCommon<NodeData>::copyToGhostedArray(const DA * octDA,  std::vector<VecInfo> & vec)  {
  const size_t localNodalSz = octDA->getLocalNodalSz();

  for(DENDRITE_UINT i = 0; i < vec.size(); i++){
    VecGetArrayRead(vec[i].v,&vec[i].val);
  }

  DENDRITE_UINT  currentDof(0) ;
  for(size_t node = 0; node < localNodalSz; node++) {
    currentDof = 0;
    for (DENDRITE_UINT i = 0; i < vec.size(); i++) {
      std::memcpy(&tempSyncVector_[node*m_totalDof + currentDof],&vec[i].val[node*vec[i].ndof], sizeof(PetscScalar)*vec[i].ndof);
      currentDof += vec[i].ndof;
    }
  }
  for(DENDRITE_UINT i = 0; i < vec.size(); i++){
    VecRestoreArrayRead(vec[i].v,&vec[i].val);
  }

  octDA->nodalVecToGhostedNodal(tempSyncVector_,ghostedVecs, true,m_totalDof);
  octDA->readFromGhostBegin(ghostedVecs,m_totalDof);
  octDA->readFromGhostEnd(ghostedVecs,m_totalDof);
}

template<typename NodeData>
DENDRITE_UINT TalyMatVecCommon<NodeData>::generateSurfaceFlags(const TALYFEMLIB::ZEROPTV &position) {
  return(getSurfaceFlags(position,physicalDomain_.min,physicalDomain_.max));
}

template<typename NodeData>
void TalyMatVecCommon<NodeData>::getSurfaceID(const double *physCoords, std::array<BoundarySurface,2*DIM> & isBoundaryFace){
  static constexpr unsigned numPoints = 1u << DIM;
  if(eleOrder_ > 2){
    throw std::runtime_error("Not implmeneted");
  }
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


  std::array<Point<DIM>,numPoints>arrayOfPoints;
  for(int i = 0; i < numPoints; i++){
    arrayOfPoints[i] = Point<DIM>(&physCoords[cornerMap[eleOrder_ - 1][i]*DIM]);
  }

  std::array<DENDRITE_UINT ,numPoints> ids{};
  std::array<std::bitset<MAX_BOUNDARY_TYPES> ,numPoints> boundaryBits;
  for(int i = 0; i < numPoints; i++){
    domainBoundary_->generateBoundaryFlags(arrayOfPoints[i],ids[i]);
    boundaryBits[i] = domainBoundary_->getBoundary();
  }
  BoundarySurface surface{static_cast<DENDRITE_UINT>(-1),static_cast<DENDRITE_UINT>(-1)};
  isBoundaryFace.fill(surface);



  for(DENDRITE_UINT id = 0; id < 2*DIM; id++) {
      bool didSurfaceIntegration = false;
    for (DENDRITE_UINT i = 0; i < MAX_BOUNDARY_TYPES; i++) {

      bool check = boundaryBits[faceID[id][0]].test(i);
      DENDRITE_UINT  objectID  = ids[faceID[id][0]];

      for (DENDRITE_UINT face = 1; (face < numNodesPerFace) and (check == true); face++) {
        check = boundaryBits[faceID[id][face]].test(i) and (ids[faceID[id][face]] == objectID);
      }
      if (check == true) {
        isBoundaryFace[id].boundaryType = i;
        isBoundaryFace[id].id = objectID;
        didSurfaceIntegration = true;
        break;
      }
    }

#ifdef IBM
    // In case of SBM this will only be executed for Intercepted or neighbors of False Intercepted
    // If it is intercepted we can proceed as usual
    // If it is neighbors of False intercepted, we need additional check
    // Note for neighbors of freshly intercepted which are not intercepted the above check will be false
    if(this->imga_->getIBMMethod() == IBM_METHOD::SBM) {
        // Check:
        auto marker = this->getElementMarkerType();
        if(!didSurfaceIntegration and (marker.test(ElementMarker::SBM_NEIGHBORS_FALSE_INTERCEPTED))) {

              // This is not a solid logic but a one night discussion with Dr. Baskar and Dr. Adarsh.
              // If the boundary nodes belong to following combination, then only integrate over the surface:
              // Combination 1: Any one of the nodes is false intercepted and other nodes are marked In (belong to inactive region)
              // Combination 2: All nodes are false intercepted
              // I (Kumar) am not sure who is going to test it thoroughly. Mostly not me.

              // Testing Combination 1

              bool SBMCheck = true;
              const auto subDomain = domainBoundary_->getSubDomain();
              for (DENDRITE_UINT face = 0; (face < numNodesPerFace); face++) {
                  SBMCheck = SBMCheck and ( (subDomain->functionToRetainPhysNodes(&physCoords[cornerMap[eleOrder_ - 1][faceID[id][face]]*DIM]) == ibm::IN)// This node belong to inactive region
                                             or ((sbmElementNodes_[cornerMap[eleOrder_ - 1][faceID[id][face]]] == NodeMarker::SBM_FALSE_INTERCEPTED_NODES)));

              }
            // Testing Combination 2 : I do not think it is required. It falls within the condition of Combination 1
            if (SBMCheck == true) {

                isBoundaryFace[id].boundaryType = 0;
                isBoundaryFace[id].id = -1;
            }
          }


      }
#endif
  }






}

template<typename NodeData>
DENDRITE_UINT TalyMatVecCommon<NodeData>::getSurfaceFlags(const TALYFEMLIB::ZEROPTV & position, const Point<DIM> & minPhysDomain, const Point<DIM> & maxPhysDomain){

  DENDRITE_UINT flags = 0;

  if (FEQUALS(position.x(),minPhysDomain.x(0))) {

    flags |= (1u << 1u);
  }
  if (FEQUALS(position.x(),maxPhysDomain.x(0))) {
    flags |= (1u << 2u);
  }
  if (FEQUALS(position.y(),minPhysDomain.x(1))) {
    flags |= (1u << 3u);
  }
  if (FEQUALS(position.y(),maxPhysDomain.x(1))) {
    flags |= (1u << 4u);
  }
#if(DIM == 3) /// Not sure if holds for 4D too
  if (FEQUALS(position.z(),minPhysDomain.x(2))) {
    flags |= (1u << 5u);
  }
  if (FEQUALS(position.z(),maxPhysDomain.x(2))) {
    flags |= (1u << 6u);
  }
#endif
  return flags;
}

#ifdef IBM
template<typename NodeData>
void TalyMatVecCommon<NodeData>::assignIMGAConstructs(const IMGA * imga,const std::bitset<ElementMarker::MAX_ELMENT_TYPE> * marker){
  imga_ = imga;
  elementMarker_ = marker;
}
template<typename NodeData>
bool TalyMatVecCommon<NodeData>::isInside(const TALYFEMLIB::ZEROPTV & position) const{
#pragma message "Add normal based in - out test"
  return (imga_->ifInside(position.data()));
}

template<typename NodeData>
bool TalyMatVecCommon<NodeData>::getLocalSurfaceGaussPoint(TALYFEMLIB::FEMElm & fe, NodeAndValues<DENDRITE_REAL> & gaussPoint){
  if(hasLocalIBMGaussPoints) {
    if ((gpInfoIt_->localElemID==lclElemID) and (gpInfoIt_ < gpInfoEnd_)) {
      gaussPoint = *gpInfoIt_;
      /// TODO Fix it. Seems stupid;
      TALYFEMLIB::ZEROPTV ptvg, ptvl;
      std::memcpy(ptvg.data(), gaussPoint.location, sizeof(double)*DIM);
      GetLocalPtv(fe, ptvg, ptvl);
      fe.calc_at(ptvl);

      static constexpr DENDRITE_REAL weight = (DIM==2) ? 2.0 : 3.0;
      fe.set_jacc_x_w(gaussPoint.elemArea/weight);
      gpInfoIt_ = std::next(gpInfoIt_);
      return true;
    }
  }
  return false;
}

template<typename NodeData>
std::bitset<ElementMarker::MAX_ELMENT_TYPE> TalyMatVecCommon<NodeData>::getElementMarkerType() const {
  const auto & marker = elementMarker_[lclElemID];

  return (marker);
}


#endif


#endif //DENDRITEKT_MATVECCOMMON_H
