//
// Created by maksbh on 3/21/22.
//

#include <oda.h>
#include <checkPoint.h>
#include <sfcTreeLoop_matvec_io.h>
#include <coarseToFine.hpp>

static constexpr int DIM = 2;
typedef ot::TreeNode<unsigned int, DIM> TREENODE;
using DistTree = ot::DistTree<uint, DIM>;
using DA = ot::DA<DIM>;

int main(int argc, char *argv[]){
  PetscInitialize(&argc, &argv, NULL, NULL);
  _InitializeHcurve(DIM);
  m_uiMaxDepth = 30;
  std::array<double,DIM> physMin, physMax;
  physMin.fill(0.0);
  physMax.fill(1.0);
  physMax[0] = 0.5;
  DistTree::BoxDecider boxDecider(physMin,physMax);
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  double sfcTol = 0.1;
  int elemOrder = 1;
  DistTree regDistTree= DistTree::constructSubdomainDistTree(4,boxDecider,MPI_COMM_WORLD,sfcTol);
  const std::vector<TREENODE>&treePart = regDistTree.getTreePartFiltered();
  DA * da = new DA(regDistTree, MPI_COMM_WORLD, elemOrder,100,sfcTol);
  std::vector<std::size_t > bdy_index;
  da->getBoundaryNodeIndices(bdy_index);
  double coords[DIM];
  for (int i = 0; i < bdy_index.size(); i++) {
    ot::treeNode2Physical(da->getTNCoords()[bdy_index[i] + da->getLocalNodeBegin()], da->getElementOrder(), coords);
    std::cout << coords[0]*4 << " " << coords[1]*4 << "\n";
  }
  PetscFinalize();

}