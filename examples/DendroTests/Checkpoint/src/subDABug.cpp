//
// Created by maksbh on 3/3/22.
//

#include <oda.h>
#include <checkPoint.h>
#include <sfcTreeLoop_matvec_io.h>

// Is giving different answer on 142 vs 143 processor"
static constexpr int DIM = 2;

typedef ot::TreeNode<unsigned int, DIM> TREENODE;
using DistTree = ot::DistTree<uint, DIM>;
using DA = ot::DA<DIM>;

void readDACheckpoint(const std::string &foldername, DistTree &dtree, const int rank) {
  std::vector<TREENODE> reloadTreeNode;
  std::string fname = foldername + "/oct_rank_" + std::to_string(rank);
  io::checkpoint::readOctFromFile(fname.c_str(), reloadTreeNode);

  DistTree dtree_(reloadTreeNode, MPI_COMM_WORLD);
  std::swap(dtree, dtree_);
}

void getfineNumLevel(DA * octDA, DistTree & dtree, std::vector<int> & numRefineLevel){
  const size_t sz = octDA->getTotalNodalSz();
  auto partFront = octDA->getTreePartFront();
  auto partBack = octDA->getTreePartBack();
  const auto tnCoords = octDA->getTNCoords();
  int elemID=0;

  const int nPe = octDA->getNumNodesPerElement();
  double coords[4*DIM]; // only for linear in 2D
  int eleOrder = 1;
  double scalingFactor = 4;
  ot::MatvecBaseCoords <DIM> loop(sz,eleOrder, false,0,tnCoords,
                                  &(*dtree.getTreePartFiltered().cbegin()),
                                  dtree.getTreePartFiltered().size(),*partFront,*partBack);
  while(!loop.isFinished()){
    if (loop.isPre() && loop.subtreeInfo().isLeaf()) {
      const double *nodeCoordsFlat = loop.subtreeInfo().getNodeCoords();

      int level = loop.getCurrentSubtree().getLevel();
      std::memcpy(coords,nodeCoordsFlat, sizeof(double)*nPe*DIM);
      for(int i = 0; i < nPe*DIM; i++){
        coords[i] *= scalingFactor; // Converting domain from [0.25,1] -> [1,4]
      }

      { // Some condition to refine.
        double amplitude = 0.1;

        /// This is interface
        double delta = 0.00125;

        double locationInt = 2.0;

        double phiCapillaryWave[4];
        /// Gaussian
        for (int i = 0; i < nPe; i++) {
          double x = coords[i * DIM + 0];
          double y = coords[i * DIM + 1];
          phiCapillaryWave[i] = tanh((locationInt - y + (amplitude * cos(2 * M_PI * x))) / (sqrt(2) * delta));
        }

        int numRefineLevelElement = 0;
        for (int i = 0; i < nPe; i++) {
          if (fabs(phiCapillaryWave[i]) < 0.99999) {
            numRefineLevelElement = 1;
            break;
          }
        }
        numRefineLevel[elemID] = numRefineLevelElement;
      }
      elemID++;
      loop.next();

    }
    else{
      loop.step();

    }
  }
}


void DACreate(bool enableMultiLevelRefine) {


  _InitializeHcurve(DIM);
  m_uiMaxDepth = 30;
  std::array<double,DIM> physMin, physMax;
  physMin.fill(0.0);
  physMax.fill(1.0);
  physMax[0] = 0.25;
  DistTree::BoxDecider boxDecider(physMin,physMax);
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  double sfcTol = 0.1;
  int elemOrder = 1;
  DistTree regDistTree= DistTree::constructSubdomainDistTree(7,boxDecider,MPI_COMM_WORLD,sfcTol);
  const std::vector<TREENODE>&treePart = regDistTree.getTreePartFiltered();
  DA * reg_subDA = new DA(regDistTree, MPI_COMM_WORLD, elemOrder,100,sfcTol);
  DA *newDA = nullptr;
  std::vector<int> refineLevel(regDistTree.getTreePartFiltered().size(),0);
  getfineNumLevel(reg_subDA,regDistTree,refineLevel);
  DistTree fineDistTree;
  if(enableMultiLevelRefine) {
    if(!comm_rank){
      std::cout << "Using Multi-Level Refine\n";
    }
    DistTree::distRefine(regDistTree, std::move(refineLevel), fineDistTree, sfcTol);
  }
  else{
    // Applicable only when numLevel is 1/0
    if(!comm_rank){
      std::cout << "Using Single Level Refine\n";
    }
    std::vector<ot::OCT_FLAGS::Refine> refineFlags(refineLevel.size());
    for(int i = 0; i < refineFlags.size();i++){
      refineFlags[i] = (ot::OCT_FLAGS::Refine) refineLevel[i];
    }
    DistTree surrDistTree;
    DistTree::distRemeshSubdomain(regDistTree, refineFlags, fineDistTree, surrDistTree,
                                  ot::RemeshPartition::SurrogateInByOut, sfcTol);

  }

  newDA = new DA(fineDistTree, MPI_COMM_WORLD, elemOrder, 100, sfcTol);
//  DomainInfo physDomain, cubeDomain;
//  cubeDomain.min.fill(0.0);
//  cubeDomain.max.fill(4.0);
//  physDomain.min.fill(0.0);
//  physDomain.max.fill(1.0);
//  physDomain.max[1] = 4.0;
//  DomainExtents domainExtents(cubeDomain, physDomain);
//  IO::writeBoundaryElements(newDA, fineDistTree.getTreePartFiltered(), "Boundary", "bnd", domainExtents);
  if(!comm_rank){
    std::cout << "Number of global nodes = " << newDA->getGlobalNodeSz() << "\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);

}

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin() ;

  if(argc < 2){
    std::cout << "Usage [enable_multiLevel]\n" ;

    MPI_Abort(MPI_COMM_WORLD,0);
  }
  bool enableMultiLevel = static_cast<bool>(std::atoi(argv[1]));
  DACreate(enableMultiLevel);

  DendroScopeEnd();
  PetscFinalize();
}