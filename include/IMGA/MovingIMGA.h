//
// Created by maksbh on 7/31/21.
//

#ifndef DENDRITEKT_MOVINGIMGA_H
#define DENDRITEKT_MOVINGIMGA_H

#include <DataTypes.h>
#include <IMGA/Marker.h>


static void getInNodes(DA *octDA, const std::vector<TREENODE> &treeNodes, Vec &inNodes, const IMGA *imga,
                       const Marker *marker, const DomainExtents &domain, bool isAllocated = false) {
  if (octDA->isActive()) {
    const auto &eleMarker = marker->getMarkers();
    if (not isAllocated) {
      octDA->petscCreateVector(inNodes, false, false, 1);
    }
    VecSet(inNodes, 0.0);

    std::vector<PetscInt> inNodeList;
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;
    OctToPhysical octToPhysical(domain);
    const size_t ghostedNodalSz = octDA->getTotalNodalSz();
    const TREENODE *odaCoords = octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
    const bool visitEmpty = false;
    const unsigned int padLevel = 0;
    const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
    DENDRITE_REAL *physCoords = new DENDRITE_REAL[DIM * nPe];
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
        if ((eleMarker[eleCounter] == OUT_GP) or eleMarker[eleCounter] == OUT_ELEMENT) { // Out of the geometry
//          continue;
        } else if ((eleMarker[eleCounter] == IN_GP) or eleMarker[eleCounter] == IN_ELEMENT) { // In of the geometry
          for (int i = 0; i < nPe; i++) {
            if (nodeNonhangingIn[i]) {
              inNodeList.push_back(nodeIdsFlat[i]);
            }
          }
        } else if ((eleMarker[eleCounter] == INTERCEPTED_GP) or
                   eleMarker[eleCounter] == INTERCEPTED_ELEMENT) { // Intercepted
          std::memcpy(physCoords, octCoords, sizeof(DENDRITE_REAL) * nPe * DIM);
          octToPhysical.convertCoordsToPhys(physCoords, nPe);
          for (int i = 0; i < nPe; i++) {
            if (nodeNonhangingIn[i]) {
              if (imga->ifInside(&physCoords[i * DIM])) {
                inNodeList.push_back(nodeIdsFlat[i]);
              }
            }
          }
        }
        eleCounter++;
      }
      treeLoopIn.step();
    }
    delete[] physCoords;
    std::sort(inNodeList.begin(), inNodeList.end());
    auto last = std::unique(inNodeList.begin(), inNodeList.end());
    inNodeList.erase(last, inNodeList.end());
    std::vector<PetscScalar> values(inNodeList.size(), 1.0);
    VecSetValues(inNodes, inNodeList.size(), inNodeList.data(), values.data(), ADD_VALUES);
    VecAssemblyBegin(inNodes);
    VecAssemblyEnd(inNodes);
    PetscScalar *array;
    VecGetArray(inNodes, &array);
    for (int i = 0; i < octDA->getLocalNodalSz(); i++) {
      if (array[i] > 0) {
        array[i] = 1.0;
      }
    }
    VecRestoreArray(inNodes, &array);
    // All IN Nodes must get the value of 1 and out nodes the value of 0 after this point .
  }
}

// Same as previous function except that it does not require marker.
static void getInNodes(DA *octDA, const std::vector<TREENODE> &treeNodes, Vec &inNodes, const IMGA *imga,
                       const DomainExtents &domain, bool isAllocated = false) {
  if (octDA->isActive()) {

    if (not isAllocated) {
      octDA->petscCreateVector(inNodes, false, false, 1);
    }
    VecSet(inNodes, 0.0);

    std::vector<PetscInt> inNodeList;
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;
    OctToPhysical octToPhysical(domain);
    const size_t ghostedNodalSz = octDA->getTotalNodalSz();
    const TREENODE *odaCoords = octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
    const bool visitEmpty = false;
    const unsigned int padLevel = 0;
    const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();
    DENDRITE_REAL *physCoords = new DENDRITE_REAL[DIM * nPe];
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


    while (!treeLoopIn.isFinished()) {
      const ot::TreeNode<CoordT, DIM> subtree = treeLoopIn.getCurrentSubtree();
      const auto subtreeInfo = treeLoopIn.subtreeInfo();
      if (treeLoopIn.isPre() && subtreeInfo.isLeaf()) {
        const RankI *nodeIdsFlat = subtreeInfo.readNodeValsIn();
        const auto &octCoords = subtreeInfo.getNodeCoords();
        const std::vector<bool> &nodeNonhangingIn = subtreeInfo.readNodeNonhangingIn();
        std::memcpy(physCoords, octCoords, sizeof(DENDRITE_REAL) * nPe * DIM);
        octToPhysical.convertCoordsToPhys(physCoords, nPe);
        for (int i = 0; i < nPe; i++) {
          if (nodeNonhangingIn[i]) {
            if (imga->ifInside(&physCoords[i * DIM])) {
              inNodeList.push_back(nodeIdsFlat[i]);
            }
          }
        }
      }
      treeLoopIn.step();
    }
    delete[] physCoords;
    std::sort(inNodeList.begin(), inNodeList.end());
    auto last = std::unique(inNodeList.begin(), inNodeList.end());
    inNodeList.erase(last, inNodeList.end());
    std::vector<PetscScalar> values(inNodeList.size(), 1.0);
    VecSetValues(inNodes, inNodeList.size(), inNodeList.data(), values.data(), ADD_VALUES);
    VecAssemblyBegin(inNodes);
    VecAssemblyEnd(inNodes);
    PetscScalar *array;
    VecGetArray(inNodes, &array);
    for (int i = 0; i < octDA->getLocalNodalSz(); i++) {
      if (array[i] > 0) {
        array[i] = 1.0;
      }
    }
    VecRestoreArray(inNodes, &array);
    // All IN Nodes must get the value of 1 and out nodes the value of 0 after this point .
  }
}

struct FreshlyClearedNodesPosition{
  PetscInt nodeID;
  DENDRITE_REAL coords[DIM];

  bool operator==(FreshlyClearedNodesPosition const &other) const {
    return ((this->nodeID) == other.nodeID);
  }//end

  bool operator<(FreshlyClearedNodesPosition const &other) const {
    return ((this->nodeID) < other.nodeID);
  }//end
  bool operator>(FreshlyClearedNodesPosition const &other) const {
    return not(*this < other);
  }//end
  static MPI_Datatype dataType() {
    static bool first = true;
    static MPI_Datatype _datatype;
    if (first) {
      first = false;
      MPI_Type_contiguous(sizeof(FreshlyClearedNodesPosition), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }

};
struct FreshlyClearedNodes {
  PetscInt nodeID;
  DENDRITE_UINT rank;

  bool operator==(FreshlyClearedNodes const &other) const {
    return ((this->nodeID) == other.nodeID);
  }//end

  bool operator<(FreshlyClearedNodes const &other) const {
    return ((this->nodeID) < other.nodeID);
  }//end
  bool operator>(FreshlyClearedNodes const &other) const {
    return not(*this < other);
  }//end
  static MPI_Datatype dataType() {
    static bool first = true;
    static MPI_Datatype _datatype;
    if (first) {
      first = false;
      MPI_Type_contiguous(sizeof(FreshlyClearedNodes), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }
};

template<int dof>
struct NeighborList {
  FreshlyClearedNodes clearedNodes;
  DENDRITE_REAL coords[DIM];
  PetscInt neighborID;

  DENDRITE_REAL values[dof];

  bool operator==(NeighborList const &other) const {
    return ((this->clearedNodes) == (other.clearedNodes) and (this->neighborID == other.neighborID));
  }//end

  bool operator<(NeighborList const &other) const {
//    return ((this->clearedNodes) < other.clearedNodes);
     if (this->clearedNodes.nodeID == other.clearedNodes.nodeID) {
       return ((this->neighborID) < other.neighborID);
     }
     return ((this->clearedNodes) < other.clearedNodes);
  }//end

  bool operator>(NeighborList const &other) const {
    return not(*this < other);
  }//end
  static bool compareWithRank(const NeighborList<dof> &first,const NeighborList<dof> & other){
    if(first.clearedNodes.rank < other.clearedNodes.rank){
      return true;
    }
    else if(first.clearedNodes.rank > other.clearedNodes.rank){
      return false;
    }
    else{
      return first < other;
    }
  }
  static MPI_Datatype dataType() {
    static bool first = true;
    static MPI_Datatype _datatype;
    if (first) {
      first = false;
//      MPI_Type_contiguous(sizeof(NeighborList), MPI_BYTE, &_datatype);
      MPI_Type_contiguous(sizeof(NeighborList<dof>), MPI_BYTE, &_datatype);
      MPI_Type_commit(&_datatype);
    }
    return _datatype;
  }
};

template<int dof>
static void
getFreshlyClearedNodes(DA *octDA, const std::vector<TREENODE> &treeNodes, const Vec &oldInNodes, const Vec &newInNodes,
                       VecInfo &valuesVec,
                       const DomainExtents &domain, bool printFreshlyCleared = false) {
  assert(valuesVec.ndof == dof);
  // TODO: Fix in case of inactive comms
  if (octDA->isActive()) {
    using CoordT = typename ot::DA<DIM>::C;
    using ot::RankI;
    OctToPhysical octToPhysical(domain);
    const size_t ghostedNodalSz = octDA->getTotalNodalSz();
    const TREENODE *odaCoords = octDA->getTNCoords();
    const std::vector<RankI> &ghostedGlobalNodeId = octDA->getNodeLocalToGlobalMap();
    const bool visitEmpty = false;
    const unsigned int padLevel = 0;
    const DENDRITE_UINT nPe = octDA->getNumNodesPerElement();

    Vec freshlyClearedVec;
    octDA->petscCreateVector(freshlyClearedVec, false, false, 1);
    /// fresh = old - new
    VecCopy(oldInNodes, freshlyClearedVec);
    VecAXPY(freshlyClearedVec, -1.0, newInNodes);
    // Freshly cleared nodes are those which have old = 1 and new = 0. So, old - new = 1;
#ifndef NDEBUG
    if (printFreshlyCleared) {
      static const char *varname[]{"fresh"};
      PETSc::petscVectopvtu(octDA, treeNodes, freshlyClearedVec, "freshlyCleared", "cleared", varname, domain, false,
                            false, 1);
    }
#endif
    const double *array;
    VecGetArrayRead(freshlyClearedVec, &array);
    PetscInt startVecOwnershipRange, endVecOwnershipRange;
    VecGetOwnershipRange(freshlyClearedVec, &startVecOwnershipRange, &endVecOwnershipRange);
    std::vector<PetscInt> freshlyClearedListnodeIds;
    for (int i = startVecOwnershipRange; i < endVecOwnershipRange; i++) {
      if (array[i - startVecOwnershipRange] == 1) {
        freshlyClearedListnodeIds.push_back(i);
      }
    }
    VecRestoreArrayRead(freshlyClearedVec, &array);

    DENDRITE_REAL *physCoords = new DENDRITE_REAL[DIM * nPe];
    std::vector<FreshlyClearedNodes> freshlyClearedList;
    std::vector<FreshlyClearedNodesPosition> freshlyClearedCoordsPosition;


    FreshlyClearedNodes temp;
    FreshlyClearedNodesPosition tempPosition;
    // nodeIDs store the ids of each element. -1 means hanging nodes
    std::vector<PetscInt> nodeIDs(treeNodes.size() * nPe, -1);

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
        std::memcpy(physCoords, octCoords, sizeof(DENDRITE_REAL) * nPe * DIM);
        octToPhysical.convertCoordsToPhys(physCoords, nPe);
        for (int i = 0; i < nPe; i++) {
          if (nodeNonhangingIn[i]) {
            nodeIDs[eleCounter * nPe + i] = nodeIdsFlat[i];
            if ((nodeIdsFlat[i] >= startVecOwnershipRange) and (nodeIdsFlat[i] < endVecOwnershipRange)) {
              if (std::binary_search(freshlyClearedListnodeIds.begin(), freshlyClearedListnodeIds.end(),
                                     nodeIdsFlat[i])) {
                tempPosition.nodeID = nodeIdsFlat[i];
                std::memcpy(tempPosition.coords,&physCoords[i*DIM], sizeof(DENDRITE_REAL)*DIM);
                temp.nodeID = nodeIdsFlat[i];
                temp.rank = octDA->getRankActive();
                freshlyClearedList.push_back(temp);
                freshlyClearedCoordsPosition.push_back(tempPosition);
              }
            }
          }
        }
        eleCounter++;
      }
      treeLoopIn.step();
    }

    VecDestroy(&freshlyClearedVec);

    // This sort does not need rank information.
    {
      std::sort(freshlyClearedList.begin(), freshlyClearedList.end());
      auto last = std::unique(freshlyClearedList.begin(), freshlyClearedList.end());
      freshlyClearedList.erase(last, freshlyClearedList.end());
    }
    {
      std::sort(freshlyClearedCoordsPosition.begin(), freshlyClearedCoordsPosition.end());
      auto last = std::unique(freshlyClearedCoordsPosition.begin(), freshlyClearedCoordsPosition.end());
      freshlyClearedCoordsPosition.erase(last, freshlyClearedCoordsPosition.end());
    }
    assert(freshlyClearedCoordsPosition.size()==freshlyClearedList.size());
    // Communicating the freshly node ids to all the processors
    std::vector<int> numberOfLocalFreshlyCleared(octDA->getNpesAll());
    DENDRITE_UINT numLocalCleared = freshlyClearedList.size();
    MPI_Allgather(&numLocalCleared, 1, MPI_INT, numberOfLocalFreshlyCleared.data(), 1, MPI_INT, octDA->getCommActive());
    DENDRITE_UINT totalFreshlyCleared = std::accumulate(numberOfLocalFreshlyCleared.begin(),
                                                        numberOfLocalFreshlyCleared.end(), 0);
    std::vector<FreshlyClearedNodes> globalFreshlyCleared(totalFreshlyCleared);
    std::vector<int> displacement(numberOfLocalFreshlyCleared.size(), 0);
    for (int i = 1; i < displacement.size(); i++) {
      displacement[i] = displacement[i - 1] + numberOfLocalFreshlyCleared[i - 1];
    }
    MPI_Allgatherv(freshlyClearedList.data(), freshlyClearedList.size(),
                   FreshlyClearedNodes::dataType(), globalFreshlyCleared.data(),
                   numberOfLocalFreshlyCleared.data(), displacement.data(),
                   FreshlyClearedNodes::dataType(), octDA->getCommActive());

    // For this sort every elements should be unique. We cannot take rank into account
    std::sort(globalFreshlyCleared.begin(), globalFreshlyCleared.end());
    // globalFreshlyCleared contains the data for all the freshly cleared nodes.

    // Now we loop over all the elements to find the nodes that is shared with the freshly cleared nodes.

    std::vector<NeighborList<dof>> neighborList;
    NeighborList<dof> tempNeighbors;

    // We need to interpolate the values only from those nodes that were not inside the geometry
    Vec valuesVecOldDAInNodes;
    VecInfo combinedVec(valuesVecOldDAInNodes, dof + 1, 0);
    {
      std::vector<VecInfo> inVecs(2);
      inVecs[0] = VecInfo(valuesVec.v, valuesVec.ndof, 0);
      inVecs[1] = VecInfo(oldInNodes, 1, 0);
      PETSc::recombineVec(octDA, inVecs, combinedVec, false);
    }

    VecGetArrayRead(combinedVec.v, &combinedVec.val);
    PetscScalar *ghostedArray;
    octDA->nodalVecToGhostedNodal(combinedVec.val, ghostedArray, false, combinedVec.ndof);
    octDA->readFromGhostBegin(ghostedArray, combinedVec.ndof);
    octDA->readFromGhostEnd(ghostedArray, combinedVec.ndof);
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    const int &ndof = combinedVec.ndof;
    ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, octDA->getElementOrder(), tnCoords, ghostedArray,
                                              &(*treeNodes.cbegin()),
                                              treeNodes.size(), *partFront, *partBack);
    eleCounter = 0;
    while (!treeloop.isFinished()) {
      if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf()) {
        const double *octCoords = treeloop.subtreeInfo().getNodeCoords();
        const PetscScalar *nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
        const unsigned int level = treeloop.getCurrentSubtree().getLevel();
        for (int k = 0; k < nPe; k++) {
          temp.nodeID = nodeIDs[eleCounter * nPe + k];
          temp.rank = octDA->getRankAll();
          if (std::binary_search(globalFreshlyCleared.begin(), globalFreshlyCleared.end(), temp)) {
            // If this node is in the list of freshly cleared nodes, then add the element nodes that were out and non-hanging into the neighbor list
            auto lowerBound = std::lower_bound(globalFreshlyCleared.begin(), globalFreshlyCleared.end(), temp);
            int dist = std::distance(globalFreshlyCleared.begin(), lowerBound);
            temp.rank = globalFreshlyCleared[dist].rank;
            // Node k exists. Now we need to push all the neighbors that lie within this element
            std::memcpy(physCoords, octCoords, sizeof(DENDRITE_REAL) * DIM * nPe);
            octToPhysical.convertCoordsToPhys(physCoords, nPe);
            for (int n = 0; n < nPe; n++) {
              if ((nodeIDs[eleCounter * nPe + k] != -1) /* Ignore hanging nodes */ and
              (nodeValsFlat[n * (dof + 1) + dof] == 0) /*Only OUT nodes*/  and (n != k)) {

                tempNeighbors.clearedNodes = temp;
                tempNeighbors.neighborID = nodeIDs[eleCounter * nPe + n];
                std::memcpy(tempNeighbors.coords, &physCoords[n * DIM], sizeof(DENDRITE_REAL) * DIM);
                std::memcpy(tempNeighbors.values, &nodeValsFlat[n * (dof + 1)],
                            sizeof(DENDRITE_REAL) * dof); /* We have dof + 1 entries */
                neighborList.push_back(tempNeighbors);
              }
            }

          }
        }
        eleCounter++;
        treeloop.next();
      } else
        treeloop.step();
    }

    size_t writtenSz = treeloop.finalize(ghostedArray);

    if (sz > 0 && writtenSz == 0)
      std::cerr << "Warning: matvec() did not write any data! Loop misconfigured?\n";

    VecRestoreArrayRead(combinedVec.v, &combinedVec.val);
    VecDestroy(&combinedVec.v);
    delete[] ghostedArray;
    delete[] physCoords;


    // First remove duplicates
    std::sort(neighborList.begin(), neighborList.end());
    auto lastNeighborList = std::unique(neighborList.begin(), neighborList.end());
    neighborList.erase(lastNeighborList, neighborList.end());

    // Now sort it such that all rank are grouped together. Can we combine both of them?
    std::sort(neighborList.begin(), neighborList.end(),NeighborList<dof>::compareWithRank);

    // Now lets start communicating to each processor via P2P communication

    // Computing send counts : How much I am going to send to other processors
    std::vector<int> sendCounts(octDA->getNpesAll(), 0);
    for (const auto &neighbor: neighborList) {
      sendCounts[neighbor.clearedNodes.rank]++;
    }


    std::vector<int> receiveCount(octDA->getNpesAll());
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, receiveCount.data(), 1, MPI_INT, octDA->getCommActive());

    std::vector<int> recvDisplacements(receiveCount.size(), 0);
    std::vector<int> sendDisplacements(sendCounts.size(), 0);

    for (int i = 1; i < sendDisplacements.size(); i++) {
      sendDisplacements[i] = sendDisplacements[i - 1] + sendCounts[i - 1];
      recvDisplacements[i] = recvDisplacements[i - 1] + receiveCount[i - 1];
    }

    std::vector<MPI_Request *> recieveRequests;
    int totalNeighborListSize = std::accumulate(receiveCount.begin(), receiveCount.end(), 0);
    std::vector<NeighborList<dof>> totalNeighborList(totalNeighborListSize);

    // Start receiving communication
    for (int i = 0; i < receiveCount.size(); i++) {
      if ((receiveCount[i] > 0) and (i != octDA->getRankAll())) {
        MPI_Request *req = new MPI_Request();
        recieveRequests.push_back(req);
        MPI_Irecv(&totalNeighborList[recvDisplacements[i]], receiveCount[i], NeighborList<dof>::dataType(), i, MPI_ANY_TAG,
                  octDA->getCommActive(), req);

      }
    }
    // Start sending communication
    for (int i = 0; i < receiveCount.size(); i++) {
      if ((sendCounts[i] > 0) and (i != octDA->getRankAll())) {
        MPI_Send(&neighborList[sendDisplacements[i]], sendCounts[i], NeighborList<dof>::dataType(), i,
                  octDA->getRankAll(), octDA->getCommActive());
      }
    }


    // Copy own data that
    int rank = octDA->getRankAll();
    assert(sendCounts[rank] == receiveCount[rank]);

    std::memcpy(&totalNeighborList[recvDisplacements[rank]], &neighborList[sendDisplacements[rank]],
                sizeof(NeighborList<dof>) * sendCounts[rank]);

    for(int i = 0; i < recieveRequests.size(); i++){
      MPI_Wait(recieveRequests[i],MPI_STATUS_IGNORE);
    }


    for(int i = 0; i < recieveRequests.size(); i++){
      delete recieveRequests[i];
    }

    // Do we need to free individual requests also?

    // Here it is important not to consider the rank and sort without the rank. As then all equal nodeIDs will be grouped together.
    // However, the rank of all the processor must be same here.
    std::sort(totalNeighborList.begin(), totalNeighborList.end());
    auto lasttotalNeighborList = std::unique(totalNeighborList.begin(), totalNeighborList.end());
    totalNeighborList.erase(lasttotalNeighborList, totalNeighborList.end());


    // TODO: Interpolate from neighbours to freshly cleared.
  }

}

#endif //DENDRITEKT_MOVINGIMGA_H
