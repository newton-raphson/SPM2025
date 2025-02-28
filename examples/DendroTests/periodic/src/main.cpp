//
// Created by maksbh on 8/12/21.
//
#include <DendriteUtils.h>
int main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);
  _InitializeHcurve(DIM);
  MPI_Comm comm = MPI_COMM_WORLD;

  int comm_size, comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  const int fineLevel = 3;
  periodic::PCoord<uint, DIM>::periods({(1u<<m_uiMaxDepth), (1u<<m_uiMaxDepth)/2});
  std::array<double, DIM> physDims = {1.0, 0.5};

  ot::DistTree<uint, DIM> distTree =
    ot::DistTree<uint, DIM>::constructSubdomainDistTree(fineLevel, comm);
  distTree.filterTree((typename ot::DistTree<unsigned int, DIM>::BoxDecider)(physDims));
  ot::DA<DIM> * octDA = new ot::DA<DIM>(distTree, comm, 1);


  std::cout << octDA->getLocalElementSz() << " " << octDA->getGlobalNodeSz() << "\n";
  std::cout << octDA->getNpesActive() << "\n";
  MPI_Finalize();
}