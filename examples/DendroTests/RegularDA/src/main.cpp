//
// Created by maksbh on 9/16/19.
//

#include <iostream>
#include <oda.h>
#include <point.h>

#define DIM 2

inline ibm::Partition domainDecider(const double *elemPhysCoords, double elemPhysSize){
  return ibm::Partition::OUT;
}
int main(int argc, char * argv[]){

  PetscInitialize(&argc, &argv, NULL, NULL);
  _InitializeHcurve(DIM);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int eleOrder = 1;
  m_uiMaxDepth = 30;
  int level = 7;
  auto distTree=ot::DistTree<uint,DIM>::constructSubdomainDistTree(level,domainDecider,MPI_COMM_WORLD,0.1);
  const auto & treePart = distTree.getTreePartFiltered();
  auto subDA = new ot::DA<DIM>(distTree, MPI_COMM_WORLD, eleOrder,100,0.1);
  if(!rank) {
   std::cout <<  "Num Active processor = " <<  subDA->getNpesActive() <<"\n";
  }

  std::cout << subDA->getLocalNodalSz() << "\n";

  PetscFinalize();
}
