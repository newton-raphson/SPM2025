//
// Created by maksbh on 2/28/22.
//

#include <oda.h>
#include <checkPoint.h>
#include <sfcTreeLoop_matvec_io.h>
#include <coarseToFine.hpp>
#include "PETSc/IO/petscVTU.h"
static const char *ch_varname[]{"phi", "mu"};

//static constexpr int DIM = 2;

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

void readPetscVec(const std::string & foldername,DA *octDA, Vec & vec, int rank){
  octDA->petscCreateVector(vec,false,false,2);
  double *array;
  VecGetArray(vec,&array);
  PetscInt lsz(0),lstart(0),lend(0);
  PetscInt flsz(0),flstart(0),flend(0);
  if(octDA->isActive()) {
    VecGetLocalSize(vec, &lsz);
    VecGetOwnershipRange(vec, &lstart, &lend);
    std::string fname = foldername+ "/val_rank_" + std::to_string(rank) + "_6";
    FILE *infile = fopen(fname.c_str(), "r");

    fread(&flsz, sizeof(PetscInt), 1, infile);
    fread(&flstart, sizeof(PetscInt), 1, infile);
    fread(&flend, sizeof(PetscInt), 1, infile);
    if (infile == NULL)
    {
      std::cout << fname << " file open failed " << std::endl;
      return  ;
    }
    if(lsz != flsz) {
      std::cout << "Local node Mismatch " << fname << "\n";
      return;
    }
    if((lstart != flstart) or (lend != flend)) {
      std::cout << " ownership range Mismatch " << fname << "\n";
      return;
    }
    fread(array, sizeof(PetscScalar),lsz, infile);
    fclose(infile);
  }
  VecRestoreArray(vec,&array);
}

void getVectorfromPetscVec(const DA * octDA, const Vec & vec,std::vector<PetscScalar> & stdVec){
  if(octDA->isActive()){
    const PetscScalar *array = nullptr;
    PetscInt lsz ;
    VecGetArrayRead(vec,&array);
    VecGetLocalSize(vec,&lsz);
    stdVec.resize(lsz);
    std::memcpy(stdVec.data(),array, sizeof(PetscScalar)*lsz);
    VecRestoreArrayRead(vec,&array);
  }
}

void loadFromCheckpoint() {
  std::string folderName = "CheckPoint";
  _InitializeHcurve(DIM);
  m_uiMaxDepth = 30;
  std::array<double,DIM> physMin, physMax;
  physMin.fill(0.0);
  physMax.fill(1.0);
  DistTree::BoxDecider boxDecider(physMin,physMax);
  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank, comm_size;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);
  std::array<uint,DIM> per;
  for (int dim = 0; dim < DIM; dim++) {
    per[dim] = ((1u << m_uiMaxDepth));
  }
  periodic::PCoord<uint , DIM>::periods(per);
  DistTree coarse_dtree, fine_dtree;
  readDACheckpoint(folderName, coarse_dtree, comm_rank);
  coarse_dtree.filterTree(boxDecider);

  DA *coarseDA = new DA(coarse_dtree, MPI_COMM_WORLD, 1, 100, 0.1);

  if (!comm_rank) {
    std::cout << " Global nodal size = " << coarseDA->getGlobalNodeSz() << "\n";

  }
  Vec coarseDAVec;
  readPetscVec(folderName,coarseDA,coarseDAVec,comm_rank);
  int ndof = 2;
  bool isBalanced = ot::is2to1Balanced(coarse_dtree.getTreePartFiltered(), MPI_COMM_WORLD);

  if(not isBalanced) {
    std::cout << isBalanced;
  }

  DomainInfo physDomain;
  physDomain.min.fill(0.0);
  physDomain.max.fill(4.0);
  PETSc::petscVectopvtu(coarseDA,coarse_dtree.getTreePartFiltered(),coarseDAVec,
                        "BeforeIntergrid","ch",ch_varname,physDomain,false,false,ndof);

  if (!comm_rank) {
    std::cout << "Checkpoint Loaded Successfully" << "\n";

  }
  {
    std::vector<int> numLevel(coarse_dtree.getTreePartFiltered().size(), 0);
    std::ifstream fin("meshAdaption_0_" + std::to_string(comm_rank));
    int numEntry;
    fin >> numEntry;
    if(numEntry != coarseDA->getLocalElementSz()){
      std::cout << comm_rank << " " << numEntry << " " << coarseDA->getLocalElementSz() << "\n";
    }
    assert(numEntry == coarseDA->getLocalElementSz());
    for (int i = 0; i < numEntry; i++) {
      fin >> numLevel[i];
    }
    fin.close();
    DistTree::distRefine(coarse_dtree, std::move(numLevel), fine_dtree, 0.1);
  }


  if (!comm_rank) {
    std::cout << "Refine Flags Loaded Successfully" << "\n";

  }

  DA *DAAfterRefine = new DA(fine_dtree, MPI_COMM_WORLD, 1, 100, 0.1);
  isBalanced = ot::is2to1Balanced(fine_dtree.getTreePartFiltered(), MPI_COMM_WORLD);

  if(not comm_rank) {
    std::cout << "Fine dtree is Balanced = " << isBalanced << "\n";
  }

  std::vector<double> stdInVec, stdOutVec;
  std::vector<double> coarseCellInDof,coarseCellOutDof;
  DAAfterRefine->createVector(stdOutVec,false,false,ndof);

  getVectorfromPetscVec(DAAfterRefine,coarseDAVec,stdInVec);

  fem::coarse_to_fine(coarse_dtree, coarseDA, ndof, 0,stdInVec , coarseCellInDof,fine_dtree, DAAfterRefine, stdOutVec, coarseCellOutDof);



  IO::writeVecTopVtu(DAAfterRefine,fine_dtree.getTreePartFiltered(),stdOutVec.data(),"AfterIntergrid","ch",ch_varname,physDomain,false,
                     false,ndof);





  MPI_Barrier(MPI_COMM_WORLD);

}

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, NULL, NULL);
  DendroScopeBegin() ;
  loadFromCheckpoint();

  DendroScopeEnd();
  PetscFinalize();
}