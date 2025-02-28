//
// Created by maksbh on 6/13/20.
//

#ifndef DENDRITEKT_REFINEMENT_H
#define DENDRITEKT_REFINEMENT_H

#include "Traversal.h"
#include <IO/VTU.h>
#include <PETSc/VecUtils.h>
#include "surrogate_cell_transfer.hpp"
enum RefinementType : short{
  POSITION_BASED = 0,
  VALUE_BASED = 1
};



enum RefinementStrategy : bool{
  SUBDOMAIN = true,
  FULL_DOMAIN = false
};

enum RefinementStage : int{
    COARSEN = 0,
    REFINE = 1,
    NO_CHANGE = 2
};

class Refinement : public Traversal{
  /// The refinement data that you want to store: REFINE, NO_CHANGE, COARSEN
  std::vector<ot::OCT_FLAGS::Refine> refineFlags_;

  std::vector<int> numRefineLevel_;
  /// counter
  size_t counter_ = 0;

  bool enableMultiRefine_;
  /// surrogate Tree Node for DA
  std::vector<TREENODE> surrTreeNode_;
  std::vector<TREENODE> newDATreeNode_;

  /// surrogate Tree Node for subDA
  DistTREE newDistTree_;
  DistTREE surrDistTree_;

  std::vector<TALYFEMLIB::ZEROPTV> coords_;
  DA * surrDA_ = nullptr;


  void getVectorfromPetscVec(const Vec & vec, std::vector<PetscScalar> & stdVec);
  void setPetscVecFromVector(DA * newDA, Vec & vec, const std::vector<PetscScalar> & stdVec);
protected:
  RefinementType m_refinementType;
  using Traversal::m_coords;

public:
  /***
   * @brief Constructor to do Position Based Refinement
   * @param octDA
   * @param treePart treePartition
   */
  Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const DomainExtents & domainInfo, const bool enableMultiRefine = false);

  /**
   * @brief Constructor to do value Based Refinement
   * @param octDA
   * @param treePart treePartition
   * @param v vecInfo for syncing vector
   * @param domainInfo domainInfo
   */
  Refinement(DA * octDA, const std::vector<TREENODE> & treePart, const VecInfo & v, const DomainExtents & domainInfo,const bool enableMultiRefine = false);



  /**
   * @ brief This is override of the traverse operation
   * @param fe
   */
  virtual void traverseOperation(TALYFEMLIB::FEMElm & fe) override;
  virtual void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;
  /**
   * @brief The user needs to override this function to provide the refinement strategy
   * depending on the position of the elements.
   * @param coords vector of coordinates
   * @param fe
   * @return Refinement flags
   */
  virtual ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide refinement flags");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }

  virtual int getNumLevelRefine(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide number of Refinement level");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }
/**
   * @brief The user needs to override this function to provide the refinement strategy
   * depending on the position of the elements.
   * @param coords vector of coordinates
   * @param fe
   * @param values the values of vector that is used to sync
   * @return Refinement flags
   */
  virtual ot::OCT_FLAGS::Refine getRefineFlags(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords, const PetscScalar * values){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide refinement flags");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }

  virtual int getNumLevelRefine(TALYFEMLIB::FEMElm & fe, const std::vector<TALYFEMLIB::ZEROPTV>&coords, const PetscScalar * values){
    TALYFEMLIB::PrintInfo("This function needs to be overridden to provide number of Refinement level");
    throw  TALYFEMLIB::TALYException() << " You need to override " << __func__ << "\n";
  }
  /**
   * @brief Returns the newDA
   * @param oldDAtreeNode The treeNode corresponding to oldDA node.
   * @return New DA
   */
  DA * getRefineDA(std::vector<TREENODE>& oldDAtreeNode,const double sfcTol = 0.001);



  /**
   *
   * @param oldDAtreeNode
   * @return
   */
  DA * getRefineSubDASingleLevel(DistTREE & olddistTree,const double sfcTol = 0.001, const RefinementStrategy strategy = RefinementStrategy::SUBDOMAIN,ot::RemeshPartition partition = ot::RemeshPartition::SurrogateOutByIn);
  DA * getRefineSubDAMultiLevel(DistTREE & olddistTree,const double sfcTol = 0.001, const RefinementStrategy strategy = RefinementStrategy::SUBDOMAIN);
  DA * getRefineSubDA(DistTREE & olddistTree,const double sfcTol = 0.001, const RefinementStrategy strategy = RefinementStrategy::SUBDOMAIN,ot::RemeshPartition partition = ot::RemeshPartition::SurrogateOutByIn);

  /**
   * @brief init the traversal operation
   */
  void initRefinement();


  /**
   * @brief This is a temporary fix in order to initiate the refinement.
   * This solves the serial Vs Parallel issue for now.
   * @param olddistTree
   * @return
   */
  DA * getForceRefineSubDA(DistTREE & olddistTree,const double sfcTol = 0.001);

  /**
   *
   * @param [in] newDA new DA
   * @param [in] newDistTree new DistTree
   * @param [in,out] inVec petsc Vec
   * @param [in] ndof  ndof
   */
  void petscIntergridTransferSingleLevelRefine(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof) ;
  void petscIntergridTransferSingleLevelCoarsen(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof) ;
  void nodalIntergridTransfer(DA *newDA, const DistTREE &newDistTree,  std::vector<PetscScalar> & nodalVec, std::vector<PetscScalar> & cellVec, const int nodalDof, const int cellDof, bool isCoarsenStage = false);
  void petscIntergridTransferMultiLevel(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof) ;
  void nodalAndCellIntergridTransferMultiLevel(DA * newDA, const  DistTREE & newDistTree, std::vector<PetscScalar> & nodalVec, std::vector<PetscScalar> & cellVec,const  int nodaldof, const int cellDof) ;

  void petscIntergridTransfer(DA * newDA, const  DistTREE & newDistTree, Vec &inVec,const  DENDRITE_UINT ndof, const RefinementStage refinementStage = RefinementStage::REFINE);
  void petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & nodalVec, std::vector<PetscScalar> & cellVec, const int cellDof,bool isCoarsenStage,const fem::CellCoarsen coarsen);
  void cellIntergridTransfer(DA * newDA, const  DistTREE & newDistTree, std::vector<double> & cellValues,const  DENDRITE_UINT ndof,
                             bool isCoarsenStage = true, fem::CellCoarsen coarsen = fem::CellCoarsen::Sum ) ;
  void petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecInfo> & vec);
  void petscIntergridTransfer(DA *newDA, const DistTREE &newDistTree, std::vector<VecUtils> & vec);
  void nodalAndCellIntergridTransfer(DA *newDA, const DistTREE &newDistTree,  std::vector<PetscScalar> & nodalVec, std::vector<PetscScalar> & cellVec, const int nodalDof, const int cellDof,
                                     bool isCoarsenStage = false, const fem::CellCoarsen coarsen = fem::CellCoarsen::Sum);

  void nodalIntergridTransferSingleLevelRefine(DA *newDA, const DistTREE &newDistTree, std::vector<double> & nodalVec, const DENDRITE_UINT ndof) ;
  void nodalIntergridTransferSingleLevelCoarsen(DA *newDA, const DistTREE &newDistTree, std::vector<double> & nodalVec, const DENDRITE_UINT ndof) ;



  void initPetscIntergridTransfer();
  void finializeIntergridTransfer();
  void printRefineFlags(const std::string & folderName, const std::string & filename, const DomainExtents & domainExtents) const{
    if(enableMultiRefine_){
      std::ofstream fout(filename + std::to_string(TALYFEMLIB::GetMPIRank()));

      std::vector<double> refineFlags(numRefineLevel_.size());
      for (int i = 0; i < refineFlags.size(); i++) {
        refineFlags[i] = numRefineLevel_[i];
      }
      static const char *varname[]{"flags"};
      IO::writeVecTopVtu(this->m_octDA, this->m_treePart, refineFlags.data(), folderName.c_str(), filename.c_str(),
                         varname, domainExtents,
                         true, false, 1);
      fout << refineFlags.size() << "\n";
      for (int i = 0; i < refineFlags.size(); i++) {
        fout << refineFlags[i] << "\n";
      }
      fout.close();
      MPI_Barrier(MPI_COMM_WORLD);
    }
    else {
      std::ofstream fout(filename + std::to_string(TALYFEMLIB::GetMPIRank()));
      std::vector<double> refineFlags(refineFlags_.size());

      for (int i = 0; i < refineFlags.size(); i++) {
        refineFlags[i] = refineFlags_[i];
      }
      fout << refineFlags.size() << "\n";
      for (int i = 0; i < refineFlags.size(); i++) {
        fout << refineFlags[i] << "\n";
      }
      fout.close();
      static const char *varname[]{"flags"};
      IO::writeVecTopVtu(this->m_octDA, this->m_treePart, refineFlags.data(), folderName.c_str(), filename.c_str(),
                         varname, domainExtents,
                         true, false, 1);
    }
  }
};
#endif //DENDRITEKT_REFINEMENT_H
