//
// Created by maksbh on 6/13/20.
//

#include <Traversal/Refinement.h>
#include <DendriteUtils.h>
#include <intergridTransfer.h>
#include <PETSc/PetscUtils.h>
#include <coarseToFine.hpp>
#include <surrogate_cell_transfer.hpp>

#ifdef PROFILING
#include "Profiling/Timer.h"
#endif
Refinement::Refinement(DA *octDA, const std::vector<TREENODE> &treePart, const DomainExtents &domainInfo,const bool enableMultiRefine)
  : Traversal(octDA, treePart, domainInfo),enableMultiRefine_(enableMultiRefine) {
  assert(octDA->getLocalElementSz() == treePart.size());
  if(enableMultiRefine_){
    numRefineLevel_.resize(treePart.size(), 0);
  }
  else{
    refineFlags_.resize(treePart.size(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  }
  m_refinementType = RefinementType::POSITION_BASED;
  coords_.resize(octDA->getNumNodesPerElement());

}
Refinement::Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v, const DomainExtents & domainInfo, const bool enableMultiRefine)
  :Traversal(octDA,treePart,v,domainInfo),enableMultiRefine_(enableMultiRefine){
  assert(octDA->getLocalElementSz() == treePart.size());
  if(enableMultiRefine_){
    numRefineLevel_.resize(treePart.size(), 0);
  }
  else{
    refineFlags_.resize(treePart.size(), ot::OCT_FLAGS::Refine::OCT_NO_CHANGE);
  }

  m_refinementType = RefinementType::VALUE_BASED;
  coords_.resize(octDA->getNumNodesPerElement());
}

void Refinement::traverseOperation(TALYFEMLIB::FEMElm &fe) {

  coordsToZeroptv(m_coords, coords_, this->m_octDA->getElementOrder(), true);
  if(enableMultiRefine_){
    assert(counter_ < numRefineLevel_.size());
    numRefineLevel_[counter_] = getNumLevelRefine(fe,coords_);
  }
  else {
    assert(counter_ < refineFlags_.size());
    refineFlags_[counter_] = getRefineFlags(fe, coords_);
  }
  counter_++;
}

void Refinement::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar * values) {

  coordsToZeroptv(m_coords, coords_, this->m_octDA->getElementOrder(), true);
  if(enableMultiRefine_){
    assert(counter_ < numRefineLevel_.size());
    numRefineLevel_[counter_] = getNumLevelRefine(fe,coords_,values);
  }
  else{
    assert(counter_ < refineFlags_.size());
    refineFlags_[counter_] = getRefineFlags(fe, coords_,values);
  }

  counter_++;
}

DA *Refinement::getRefineDA(std::vector<ot::TreeNode<DENDRITE_UINT, DIM>> &oldDAtreeNode,const double sfcTol) {
  throw  std::logic_error("Use getRefineSubDA");
}

DA *Refinement::getRefineSubDASingleLevel(DistTREE &oldDistTree,const double sfcTol,const RefinementStrategy strategy, ot::RemeshPartition partition) {
  bool doRefine = true;
  if (std::all_of(refineFlags_.begin(), refineFlags_.end(),
                  [](ot::OCT_FLAGS::Refine v) { return v == ot::OCT_FLAGS::Refine::OCT_NO_CHANGE; })) {
    doRefine = false;
  }

  bool gRefine = false;
  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

  if (not(gRefine)) {
    TALYFEMLIB::PrintStatus("Active Processors [No Remesh] = ", m_octDA->getNpesActive());
    return NULL;
  }
  {

    if (strategy == RefinementStrategy::SUBDOMAIN) {
      DOLLAR("distRemesh-subdomain")
      DistTREE::distRemeshSubdomain(oldDistTree, refineFlags_, newDistTree_, surrDistTree_,
                                    partition, sfcTol);
    } else {
      DOLLAR("distRemesh-fullDomain")
      TALYFEMLIB::PrintWarning("Using deprecated routine");
      DistTREE::distRemeshSubdomainViaWhole(oldDistTree, refineFlags_, newDistTree_, surrDistTree_,
                                            partition, sfcTol);
    }
  }
  DENDRITE_UINT oldDASz = m_octDA->getLocalElementSz();
  DENDRITE_UINT newDASz = newDistTree_.getFilteredTreePartSz();
  doRefine = (oldDASz!=newDASz);
  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
  if (not(gRefine)) {
    TALYFEMLIB::PrintStatus("Active Processors [No Remesh] = ", m_octDA->getNpesActive());
    return NULL;
  }
  DA *newDA = nullptr;
  {
    DOLLAR("DACreation-singleLevel")
    newDA = new DA(newDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
  }

  std::swap(newDistTree_, oldDistTree);
  TALYFEMLIB::PrintStatus("Active Processors [Remesh] = ", newDA->getNpesActive());
  return newDA;


}

DA *Refinement::getRefineSubDAMultiLevel(DistTREE &oldDistTree,const double sfcTol,const RefinementStrategy strategy) {
#ifdef PROFILING
  Timer timer("getRefineSubDAMultiLevel");
#endif
  bool doRefine = true;
  if (std::all_of(numRefineLevel_.begin(), numRefineLevel_.end(),
                  [](int v) { return v == 0; })) {
    doRefine = false;
  }
  bool gRefine = false;
  MPI_Allreduce(&doRefine, &gRefine, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

  if (not(gRefine)) {
    TALYFEMLIB::PrintStatus("Active Processors [No Remesh] = ", m_octDA->getNpesActive());
    return NULL;
  }
  {
    DOLLAR("distRefine-multilevel")
    DistTREE::distRefine(oldDistTree, std::move(numRefineLevel_), newDistTree_, sfcTol,true);
  }
  DA *newDA = nullptr;
  {
    DOLLAR("NewDA-multilevel")
    newDA = new DA(newDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
  }
  std::swap(newDistTree_,oldDistTree);
  TALYFEMLIB::PrintStatus("Active Processors [Remesh] = ", newDA->getNpesActive());
  return newDA;

}

DA *Refinement::getRefineSubDA(DistTREE &oldDistTree,const double sfcTol,const RefinementStrategy strategy, ot::RemeshPartition partition) {
  if(this->enableMultiRefine_){
#ifdef PROFILING
    DOLLAR("Multi-level")
    double start = MPI_Wtime();
#endif
    return getRefineSubDAMultiLevel(oldDistTree,sfcTol,strategy);
#ifdef PROFILING
    double end = MPI_Wtime();
    TALYFEMLIB::PrintInfo("[Time][Refine] = ",end - start);
#endif
  }
  else{
    DOLLAR("Single-level")
    return getRefineSubDASingleLevel(oldDistTree,sfcTol,strategy,partition);
  }
}

DA *Refinement::getForceRefineSubDA(DistTREE &oldDistTree, const double sfcTol) {
  if(TALYFEMLIB::GetMPISize() != m_octDA->getNpesActive()) {
    TALYFEMLIB::PrintStatus("Using refine logic for re-partitioning");
    enableMultiRefine_ = true;
    numRefineLevel_.resize(oldDistTree.getFilteredTreePartSz(), 0);
    DistTREE::distRefine(oldDistTree, std::move(numRefineLevel_), newDistTree_, sfcTol, true);
  }
  else{
    TALYFEMLIB::PrintStatus("Using repartition logic");
    newDistTree_ = oldDistTree;
    newDistTree_ = std::move(newDistTree_).repartitioned(sfcTol);

  }

  DA *newDA = new DA(newDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder(), 100, sfcTol);
  TALYFEMLIB::PrintStatus("Active Processors [Forced] = ", newDA->getNpesActive());
  std::swap(newDistTree_,oldDistTree);

  return newDA;
}

void Refinement::initRefinement() {
  DOLLAR("RefinementFlags")
#ifdef PROFILING
  double start = MPI_Wtime();
#endif
  this->traverse();
#ifdef PROFILING
  double end = MPI_Wtime();
  TALYFEMLIB::PrintInfo("[Time][flags] = ",end - start);
#endif


}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & vec){
  if(vec.size() == 1){
    petscIntergridTransfer(newDA,newDistTree,vec[0].v,vec[0].ndof);
    return;
  }
  int ndof = 0;
  for(const auto & v: vec){
    ndof+= v.ndof;
  }
  Vec outVec;
  VecInfo vecInfoOutVec(outVec,ndof,0);

  if(m_octDA->isActive()) {
    PETSc::recombineVec(m_octDA, vec, vecInfoOutVec, false);
  }
  petscIntergridTransfer(newDA,newDistTree,vecInfoOutVec.v,vecInfoOutVec.ndof);

  if(m_octDA->isActive()) {
    for(auto & v: vec){
      VecDestroy(&v.v);
    }
  }
  for(auto & v: vec){
    newDA->petscCreateVector(v.v, false, false,v.ndof);
  }
  if(newDA->isActive()){
    PETSc::scatterVec(newDA,vec,vecInfoOutVec,true);
    VecDestroy(&vecInfoOutVec.v);
  }
}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & nodalVec,
                                        std::vector<PetscScalar> & cellVec, const int cellDof, bool isCoarsenStage,
                                        const fem::CellCoarsen coarsen){
#ifdef PROFILING
  double start = MPI_Wtime();
#endif
  int ndof = 0;
  for(const auto & v: nodalVec){
    ndof+= v.ndof;
  }

  std::vector<PetscScalar> combinedNodalVec;
  {
    DOLLAR("RecombineVec")
    if (m_octDA->isActive()) {
      PETSc::recombineVec(m_octDA, nodalVec, combinedNodalVec);
    }
  }
  {
    DOLLAR("intergrid-transfer")
    nodalAndCellIntergridTransfer(newDA, newDistTree, combinedNodalVec, cellVec, ndof, cellDof, isCoarsenStage,coarsen);
  }

  {
    DOLLAR("DestroyVec")
    if (m_octDA->isActive()) {
      for (auto &v: nodalVec) {
        VecDestroy(&v.v);
      }
    }
  }
  {
    DOLLAR("allocateVec")
    for (auto &v: nodalVec) {
      newDA->petscCreateVector(v.v, false, false, v.ndof);
    }
  }
  {
    DOLLAR("ScatterVec")
    if (newDA->isActive()) {
      PETSc::scatterVec(newDA, nodalVec, combinedNodalVec, true);
    }
  }
#ifdef PROFILING
  double end = MPI_Wtime();
  TALYFEMLIB::PrintInfo("[Time][IntergridTransfer] = ", end - start);
#endif
}
void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecUtils> & vec){
  if(vec.size() == 1){
    petscIntergridTransfer(newDA,newDistTree,vec[0].petscVec,vec[0].ndof);
    return;
  }
  int ndof = 0;
  for(const auto & v: vec){
    ndof+= v.ndof;
  }
  Vec outVec;
  VecUtils vecInfoOutVec(outVec,ndof);

  if(m_octDA->isActive()) {
    PETSc::recombineVec(m_octDA, vec, vecInfoOutVec, false);
  }

  petscIntergridTransfer(newDA,newDistTree,vecInfoOutVec.petscVec,vecInfoOutVec.ndof);

  if(m_octDA->isActive()) {
    for(auto & v: vec){
      VecDestroy(&v.petscVec);
    }
  }
  for(auto & v: vec){
    newDA->petscCreateVector(v.petscVec, false, false,v.ndof);
  }
  if(newDA->isActive()){
    PETSc::scatterVec(newDA,vec,vecInfoOutVec, true);
    VecDestroy(&vecInfoOutVec.petscVec);
  }


}
void Refinement::getVectorfromPetscVec(const Vec & vec,std::vector<PetscScalar> & stdVec){
  if(m_octDA->isActive()){
    const PetscScalar *array = nullptr;
    PetscInt lsz ;
    VecGetArrayRead(vec,&array);
    VecGetLocalSize(vec,&lsz);
    stdVec.resize(lsz);
    std::memcpy(stdVec.data(),array, sizeof(PetscScalar)*lsz);
    VecRestoreArrayRead(vec,&array);
  }
}
void Refinement::setPetscVecFromVector(DA* newDA, Vec & vec, const std::vector<PetscScalar> & stdVec){
  double *newDAVecArray;
  if(newDA->isActive()) {
    VecGetArray(vec, &newDAVecArray);
    std::memcpy(newDAVecArray,stdVec.data(), sizeof(PetscScalar)*stdVec.size());
    VecRestoreArray(vec,&newDAVecArray);
  }

}

void Refinement::petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, Vec &inVec, const DENDRITE_UINT ndof, const RefinementStage refinementStage) {
  if(enableMultiRefine_){
      if(refinementStage == RefinementStage::REFINE) {
          this->petscIntergridTransferMultiLevel(newDA, newDistTree, inVec, ndof);
      }
      else{
          throw std::logic_error("Function not yet implemented") ;
      }
  }
  else{
    if(refinementStage == RefinementStage::REFINE) {
        this->petscIntergridTransferSingleLevelRefine(newDA, newDistTree, inVec, ndof);
    }
    else{
        this->petscIntergridTransferSingleLevelCoarsen(newDA, newDistTree, inVec, ndof);
    }


  }
}

void Refinement::nodalAndCellIntergridTransfer(DA *newDA, const DistTREE &newDistTree,std::vector<PetscScalar> & nodalVec, std::vector<PetscScalar> & cellVec, const int nodalDof, const int cellDof, bool isCoarsenStage, const fem::CellCoarsen coarsen) {
  if(enableMultiRefine_ and not isCoarsenStage){
    this->nodalAndCellIntergridTransferMultiLevel(newDA,newDistTree,nodalVec,cellVec,nodalDof,cellDof);
  }
  else{
      if(isCoarsenStage){
          this->nodalIntergridTransferSingleLevelCoarsen(newDA,newDistTree,nodalVec,nodalDof);
      }
      else{
          this->nodalIntergridTransferSingleLevelRefine(newDA,newDistTree,nodalVec,nodalDof);
      }

    this->cellIntergridTransfer(newDA,newDistTree,cellVec,cellDof,isCoarsenStage,coarsen);
  }
}



void Refinement::cellIntergridTransfer(DA * newDA, const  DistTREE & newDistTree, std::vector<double> & cellValues,
                                       const  DENDRITE_UINT ndof, bool isCoarsenStage,const fem::CellCoarsen coarsen) {
  std::vector<double> newDACellValues(newDA->getLocalElementSz()*ndof);
  if(enableMultiRefine_){
    std::vector<double> temp;
    fem::coarse_to_fine(newDistTree_, this->m_octDA,0,ndof,std::vector<double>{},cellValues,newDistTree,newDA,temp,newDACellValues);

  }
  else{
    if(isCoarsenStage){
      fem::cell_transfer_coarsen(newDistTree_,this->m_octDA,ndof,cellValues.data(),surrDistTree_,surrDA_,newDA,newDACellValues.data(),coarsen);
    }
    else{
      fem::cell_transfer_refine(this->m_octDA,ndof,cellValues.data(),surrDistTree_,surrDA_,newDistTree,newDA,newDACellValues.data());
    }

  }
  std::swap(cellValues,newDACellValues);
}

void Refinement::petscIntergridTransferMultiLevel(DA *newDA, const DistTREE &newDistTree, Vec &inVec, const DENDRITE_UINT ndof) {
  std::vector<PetscScalar> stdInVec, stdOutVec;
  std::vector<PetscScalar> coarseCellInDof,coarseCellOutDof;
  newDA->createVector(stdOutVec,false,false,ndof);
  this->getVectorfromPetscVec(inVec,stdInVec);

  fem::coarse_to_fine(newDistTree_, this->m_octDA, ndof, 0,stdInVec , coarseCellInDof,newDistTree, newDA, stdOutVec, coarseCellOutDof);
  if(this->m_octDA->isActive()){
    VecDestroy(&inVec);
  }
  newDA->petscCreateVector(inVec, false, false, ndof);
  this->setPetscVecFromVector(newDA,inVec,stdOutVec);
}

void Refinement::nodalAndCellIntergridTransferMultiLevel(DA *newDA, const DistTREE &newDistTree, std::vector<PetscScalar>& nodalVec, std::vector<PetscScalar> & cellVec, const int nodaldof, const int cellDof) {
  std::vector<PetscScalar>  newDAnodalVec;
  std::vector<PetscScalar>  newDACellVec;
  newDA->createVector(newDAnodalVec,false,false,nodaldof);
  newDA->createVector(newDACellVec,true,false,cellDof);

  fem::coarse_to_fine(newDistTree_, this->m_octDA, nodaldof, cellDof,nodalVec , cellVec,newDistTree, newDA, newDAnodalVec, newDACellVec);

  std::swap(cellVec,newDACellVec);
  std::swap(nodalVec,newDAnodalVec);
}
void Refinement::petscIntergridTransferSingleLevelCoarsen(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof){
    // fine -> intergrid surrogate -> shift coarse

    // m_octDA = fineDA newDA = coarseDA
    static std::vector<VECType> fineGhosted;
    static std::vector<VECType> surrGhosted;
    m_octDA->template createVector(fineGhosted, false, true, ndof);
    surrDA_->template createVector<VECType>(surrGhosted, false, true, ndof);
    std::fill(surrGhosted.begin(), surrGhosted.end(), 0);
    PetscScalar *array = nullptr;
    if(m_octDA->isActive()) {
        VecGetArray(inVec, &array);
    }
    m_octDA->nodalVecToGhostedNodal(array, fineGhosted.data(), ndof);
    m_octDA->readFromGhostBegin(fineGhosted.data(), ndof);
    m_octDA->readFromGhostEnd(fineGhosted.data(), ndof);

    fem::MeshFreeInputContext<VECType, TREENODE>
            inctx = fem::mesh_free_input_context(fineGhosted.data(), m_octDA, newDistTree_);

    fem::MeshFreeOutputContext<VECType, TREENODE>
            outctx = fem::mesh_free_output_context(surrGhosted.data(), surrDA_, surrDistTree_);

    // Hack: Only to writeToGhosts, when not all owned nodes touch owned cells.
    // (When fixed, i.e. owned nodes touch owned cells, can use readFromGhosts).
    static std::vector<char> outDirty;
    const char zero = 0;
    surrDA_->template createVector<char>(outDirty, false, true, 1);
    surrDA_->setVectorByScalar(outDirty.data(), &zero, false, true, 1);

    const RefElement *refel = newDA->getReferenceElement();
    fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

    surrDA_->template writeToGhostsBegin<VECType>(surrGhosted.data(), ndof, &(*outDirty.cbegin()));
    surrDA_->template writeToGhostsEnd<VECType>(surrGhosted.data(), ndof, false, &(*outDirty.cbegin()));

    if(m_octDA->isActive()) {
        VecRestoreArray(inVec, &array);
        VecDestroy(&inVec);
    }

    newDA->petscCreateVector(inVec, false, false, ndof);
    double *newDAVecArray;
    if(newDA->isActive()) {
        VecGetArray(inVec, &newDAVecArray);
    }
    ot::distShiftNodes(*surrDA_, surrGhosted.data() + ndof * surrDA_->getLocalNodeBegin(),
                       *newDA, newDAVecArray,
                       ndof);
    if(newDA->isActive()) {
        VecRestoreArray(inVec, &newDAVecArray);
    }
}

void Refinement::petscIntergridTransferSingleLevelRefine(DA *newDA, const DistTREE &newDistTree, Vec &inVec, const DENDRITE_UINT ndof) {
  static std::vector<VECType> fineGhosted, surrGhosted;
  newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
  surrDA_->template createVector<VECType>(surrGhosted, false, true, ndof);
  std::fill(fineGhosted.begin(), fineGhosted.end(), 0);

  VECType *fineGhostedPtr = fineGhosted.data();
  VECType *surrGhostedPtr = surrGhosted.data();
  PetscScalar *array = nullptr;
  if(m_octDA->isActive()) {
    VecGetArray(inVec, &array);
  }

  ot::distShiftNodes(*this->m_octDA, array,
                     *surrDA_, surrGhostedPtr + ndof * surrDA_->getLocalNodeBegin(),
                     ndof);

  surrDA_->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
  surrDA_->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);

  fem::MeshFreeInputContext<VECType, TREENODE>
    inctx{surrGhostedPtr,
          surrDA_->getTNCoords(),
          (unsigned) surrDA_->getTotalNodalSz(),
          &(*surrDistTree_.getTreePartFiltered().cbegin()),
          surrDistTree_.getTreePartFiltered().size(),
          *surrDA_->getTreePartFront(),
          *surrDA_->getTreePartBack()};

  fem::MeshFreeOutputContext<VECType, TREENODE>
    outctx{fineGhostedPtr,
           newDA->getTNCoords(),
           (unsigned) newDA->getTotalNodalSz(),
           &(*newDistTree.getTreePartFiltered().cbegin()),
           newDistTree.getTreePartFiltered().size(),
           *newDA->getTreePartFront(),
           *newDA->getTreePartBack()};

  // Hack for the current version
  std::vector<char> outDirty;
  newDA->template createVector<char>(outDirty, false, true, 1);
  newDA->setVectorByScalar(&(*outDirty.begin()), (std::array<char, 1>{0}).data(), false, true, 1);

  const RefElement *refel = newDA->getReferenceElement();
  fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

  newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof, &(*outDirty.cbegin()));
  newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof, false, &(*outDirty.cbegin()));
  if(m_octDA->isActive()) {
    VecRestoreArray(inVec, &array);
    VecDestroy(&inVec);
  }

  newDA->petscCreateVector(inVec, false, false, ndof);
  double *newDAVecArray;
  if(newDA->isActive()) {
    VecGetArray(inVec, &newDAVecArray);
  }
  newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, newDAVecArray, true, ndof);
  if(newDA->isActive()) {
    VecRestoreArray(inVec, &newDAVecArray);
  }
}

void Refinement::nodalIntergridTransferSingleLevelRefine(DA *newDA, const DistTREE &newDistTree, std::vector<double> & nodalVec, const DENDRITE_UINT ndof) {
  static std::vector<VECType> fineGhosted, surrGhosted;
  newDA->template createVector<VECType>(fineGhosted, false, true, ndof);
  surrDA_->template createVector<VECType>(surrGhosted, false, true, ndof);
  std::fill(fineGhosted.begin(), fineGhosted.end(), 0);

  VECType *fineGhostedPtr = fineGhosted.data();
  VECType *surrGhostedPtr = surrGhosted.data();


  ot::distShiftNodes(*this->m_octDA, nodalVec.data(),
                     *surrDA_, surrGhostedPtr + ndof * surrDA_->getLocalNodeBegin(),
                     ndof);

  surrDA_->template readFromGhostBegin<VECType>(surrGhostedPtr, ndof);
  surrDA_->template readFromGhostEnd<VECType>(surrGhostedPtr, ndof);

  fem::MeshFreeInputContext<VECType, TREENODE>
      inctx{surrGhostedPtr,
            surrDA_->getTNCoords(),
            (unsigned) surrDA_->getTotalNodalSz(),
            &(*surrDistTree_.getTreePartFiltered().cbegin()),
            surrDistTree_.getTreePartFiltered().size(),
            *surrDA_->getTreePartFront(),
            *surrDA_->getTreePartBack()};

  fem::MeshFreeOutputContext<VECType, TREENODE>
      outctx{fineGhostedPtr,
             newDA->getTNCoords(),
             (unsigned) newDA->getTotalNodalSz(),
             &(*newDistTree.getTreePartFiltered().cbegin()),
             newDistTree.getTreePartFiltered().size(),
             *newDA->getTreePartFront(),
             *newDA->getTreePartBack()};

  // Hack for the current version
  std::vector<char> outDirty;
  newDA->template createVector<char>(outDirty, false, true, 1);
  newDA->setVectorByScalar(&(*outDirty.begin()), (std::array<char, 1>{0}).data(), false, true, 1);

  const RefElement *refel = newDA->getReferenceElement();
  fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

  newDA->template writeToGhostsBegin<VECType>(fineGhostedPtr, ndof, &(*outDirty.cbegin()));
  newDA->template writeToGhostsEnd<VECType>(fineGhostedPtr, ndof, false, &(*outDirty.cbegin()));


  newDA->createVector(nodalVec, false, false, ndof);
  double *newDAVecArray = nodalVec.data();

  newDA->template ghostedNodalToNodalVec<VECType>(fineGhostedPtr, newDAVecArray, true, ndof);
}


void Refinement::nodalIntergridTransferSingleLevelCoarsen(DA *newDA, const DistTREE &newDistTree, std::vector<double> & nodalVec, const DENDRITE_UINT ndof) {
    static std::vector<VECType> fineGhosted, surrGhosted;

    m_octDA->template createVector(fineGhosted, false, true, ndof);
    surrDA_->template createVector<VECType>(surrGhosted, false, true, ndof);
    std::fill(surrGhosted.begin(), surrGhosted.end(), 0);

    m_octDA->nodalVecToGhostedNodal(nodalVec.data(), fineGhosted.data(), ndof);
    m_octDA->readFromGhostBegin(fineGhosted.data(), ndof);
    m_octDA->readFromGhostEnd(fineGhosted.data(), ndof);

    fem::MeshFreeInputContext<VECType, TREENODE>
            inctx = fem::mesh_free_input_context(fineGhosted.data(), m_octDA, newDistTree_);

    fem::MeshFreeOutputContext<VECType, TREENODE>
            outctx = fem::mesh_free_output_context(surrGhosted.data(), surrDA_, surrDistTree_);

    // Hack: Only to writeToGhosts, when not all owned nodes touch owned cells.
    // (When fixed, i.e. owned nodes touch owned cells, can use readFromGhosts).
    static std::vector<char> outDirty;
    const char zero = 0;
    surrDA_->template createVector<char>(outDirty, false, true, 1);
    surrDA_->setVectorByScalar(outDirty.data(), &zero, false, true, 1);

    const RefElement *refel = newDA->getReferenceElement();
    fem::locIntergridTransfer(inctx, outctx, ndof, refel, &(*outDirty.begin()));

    surrDA_->template writeToGhostsBegin<VECType>(surrGhosted.data(), ndof, &(*outDirty.cbegin()));
    surrDA_->template writeToGhostsEnd<VECType>(surrGhosted.data(), ndof, false, &(*outDirty.cbegin()));


    newDA->createVector(nodalVec, false, false, ndof);

    ot::distShiftNodes(*surrDA_, surrGhosted.data() + ndof * surrDA_->getLocalNodeBegin(),
                       *newDA,  nodalVec.data(),
                       ndof);

}

void Refinement::initPetscIntergridTransfer() {
  if(not(enableMultiRefine_)) {
    DOLLAR("SurrogateDA")
    surrDA_ = new DA(surrDistTree_, MPI_COMM_WORLD, this->m_octDA->getElementOrder());
  }
}

void Refinement::finializeIntergridTransfer(){
  if(not(enableMultiRefine_)) {
    delete surrDA_;
  }
}