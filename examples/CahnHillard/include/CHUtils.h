//
// Created by maksbh on 2/11/23.
//

#ifndef DENDRITEKT_CHUTILS_H
#define DENDRITEKT_CHUTILS_H
#include <CHRefine.h>
#include "Checkpoint/Checkpointer.h"
#include "MergeChildrenMass.h"

void checkIntergridTransfer(const Vec & vec, ot::DA<DIM> * octDA,const ot::DistTree<unsigned int, DIM> &distTree, const unsigned int ndof){
    const double * array;
    VecGetArrayRead(vec,&array);
    double *ghostedArray;
    octDA->nodalVecToGhostedNodal(array,ghostedArray, false,ndof);
    octDA->readFromGhostBegin(ghostedArray,ndof);
    octDA->readFromGhostEnd(ghostedArray,ndof);
    const size_t sz = octDA->getTotalNodalSz();
    auto partFront = octDA->getTreePartFront();
    auto partBack = octDA->getTreePartBack();
    const auto tnCoords = octDA->getTNCoords();
    const unsigned int nPe = octDA->getNumNodesPerElement();
    ot::MatvecBase<DIM, PetscScalar> treeloop(sz, ndof, octDA->getElementOrder(), tnCoords, ghostedArray, &(*distTree.getTreePartFiltered().cbegin()), distTree.getTreePartFiltered().size(), *partFront, *partBack);
    bool testPassed = true;
    while (!treeloop.isFinished())
    {
        if (treeloop.isPre() && treeloop.subtreeInfo().isLeaf())
        {
            const double * nodeCoordsFlat = treeloop.subtreeInfo().getNodeCoords();
            const PetscScalar * nodeValsFlat = treeloop.subtreeInfo().readNodeValsIn();
            for(int i = 0; i < nPe; i++){
                double correctValue = 0;
                for(int dim = 0; dim < DIM; dim++){
                    correctValue += nodeCoordsFlat[DIM*i+dim];
                }
                for(int dof = 0; dof < ndof; dof++) {
                    double interpolatedValue = nodeValsFlat[i * ndof + dof];
                    if (fabs(interpolatedValue - correctValue) > 1E-6) {
                        std::cout << "Value at (" << nodeCoordsFlat[DIM * i + 0] << " ," << nodeCoordsFlat[DIM * i + 1]
                                  << ") = " << interpolatedValue << "\n";
                        testPassed = false;
                    }
                }

            }
            treeloop.next();
        }
        else
            treeloop.step();
    }
    bool gtestPassed;
    MPI_Allreduce(&testPassed,&gtestPassed,1,MPI_CXX_BOOL,MPI_LAND,MPI_COMM_WORLD);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(!rank) {
        if (!gtestPassed) {
            std::cout << RED << "TEST failed" << NRM << "\n";
        }
    }
    delete [] ghostedArray;
    VecRestoreArrayRead(vec,&array);

    MPI_Barrier(MPI_COMM_WORLD);
    if(!gtestPassed){
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
}

template<RefinementStage stage>
bool performRefinement(DA *& octDA, DistTREE & dtree, VecInfo & vec, DomainExtents domainExtents, const CHInputData * chInputData,
                       std::vector<double> & cellmass,std::vector<double> & cellmassPerChild,std::vector<double> & getCellCoarsened,
                       bool doCellIntergrid = false){

    bool DAChanged = false;
    int counter = 0;
    getCellCoarsened.resize(octDA->getLocalElementSz(),1);
    while (true) {
        if(stage == RefinementStage::REFINE){
            TALYFEMLIB::PrintStatus("Refining");
        }
        else{
            TALYFEMLIB::PrintStatus("Coarsening");
        }
//        std::string foldername;
//        if(stage == RefinementStage::REFINE){
//            foldername = "Refine";
//        }
//        else{
//            foldername = "Coarse";
//        }
        std::function<void(const double *, double *)> linearFunction = [&](const double *x, double *var) {
          var[0] = 0.0;
          for (int d = 0; d < DIM; d++){
            var[0] += x[d];
          }


        };

//        Checkpointer checkpointBefore(1,foldername+"_Before");
//        Checkpointer checkpointAfter(1,foldername+"_After");


        Vec linearVec;
        octDA->petscCreateVector(linearVec, false, false,1);
        octDA->petscSetVectorByFunction(linearVec,linearFunction,false,false,1);
        std::vector<VecInfo> vecsBefore = {VecInfo(linearVec,1,0)};
//        checkpointBefore.storeCheckpoint(octDA,&dtree,vecsBefore,domainExtents, nullptr);
        CHRefine<stage> chRefine(octDA, dtree, vec, domainExtents, chInputData,false);
//        chRefine.printRefineFlags("RefineFlags",foldername,domainExtents);
        ot::RemeshPartition partition = stage == RefinementStage::REFINE ? ot::SurrogateInByOut : ot::SurrogateOutByIn;
        const bool isCoarsenStage = (stage == RefinementStage::COARSEN);


        auto newDA = chRefine.getRefineSubDA(dtree, chInputData->mesh_def.sfc_tol, RefinementStrategy::FULL_DOMAIN, partition);
        if (newDA == nullptr) {
            return DAChanged;
        }
        DAChanged = true;
        chRefine.initPetscIntergridTransfer();
        chRefine.petscIntergridTransfer(newDA, dtree, vec.v, vec.ndof, stage);
        chRefine.petscIntergridTransfer(newDA, dtree, linearVec, 1, stage);


        if(doCellIntergrid) {


            static constexpr int numItg = 1u << DIM;
            static constexpr int numChild = 1u << DIM;

            std::vector<double> tempCell(octDA->getLocalElementSz(),1);



            chRefine.cellIntergridTransfer(newDA,dtree,cellmass,CHNodeData::CH_DOF*numItg,isCoarsenStage,fem::CellCoarsen::Sum);
            chRefine.cellIntergridTransfer(newDA,dtree,cellmassPerChild,CHNodeData::CH_DOF*numChild,isCoarsenStage,fem::CellCoarsen::Sum);
            chRefine.cellIntergridTransfer(newDA,dtree,tempCell,1,isCoarsenStage,fem::CellCoarsen::Sum);
            chRefine.cellIntergridTransfer(newDA,dtree,getCellCoarsened,1,isCoarsenStage,fem::CellCoarsen::Sum);

            for(int i = 0; i < tempCell.size(); i++){
                if(FEQUALS(tempCell[i],numChild)){
                    for(int child = 0; child < numChild; child++){
                        // This is with assumption numChild = numDof
                        cellmass[i*numChild*CHNodeData::CH_DOF + child*CHNodeData::CH_DOF + CHNodeData::PHI_DOF]=
                                cellmassPerChild[i*numChild*CHNodeData::CH_DOF + child*CHNodeData::CH_DOF + CHNodeData::PHI_DOF];
                        cellmass[i*numChild*CHNodeData::CH_DOF + child*CHNodeData::CH_DOF + CHNodeData::MU_DOF]=
                                cellmassPerChild[i*numChild*CHNodeData::CH_DOF + child*CHNodeData::CH_DOF + CHNodeData::MU_DOF];

                    }
                }
            }

            std::vector<double> vecOut;
            MergeChildrenMass mergeChildrenMass(newDA,dtree.getTreePartFiltered(),domainExtents,vecOut,cellmass,tempCell,CHNodeData::CH_DOF);

            std::swap(vecOut,cellmassPerChild);
        }

        chRefine.finializeIntergridTransfer();
        std::swap(octDA, newDA);
        counter++;
        delete newDA;
        TALYFEMLIB::PrintStatus("Total number of nodes = ", octDA->getGlobalNodeSz());
        std::vector<VecInfo> vecsAfter = {VecInfo(linearVec,1,0)};
//        checkpointAfter.storeCheckpoint(octDA,&dtree,vecsAfter,domainExtents, nullptr);
        checkIntergridTransfer(linearVec, octDA,dtree, 1);
        VecDestroy(&linearVec);

    }

}

template<RefinementStage stage>
bool performCoarsening(DA *& octDA, DistTREE & dtree, VecInfo & vec, DomainExtents domainExtents, const CHInputData * chInputData,
                       std::vector<double> & cellmass,bool doCellIntergrid = false){
    bool DAChanged = false;
    int counter = 0;

    std::function<void(const double *, double *)> linearFunction = [&](const double *x, double *var) {
        var[0] = x[0] + x[1];

    };


    Vec linearVec;
    octDA->petscCreateVector(linearVec, false, false, 1);
    octDA->petscSetVectorByFunction(linearVec, linearFunction, false, false, 1);
    std::vector<VecInfo> vecsBefore = {VecInfo(linearVec, 1, 0)};
    CHRefine<stage> chRefine(octDA, dtree, vec, domainExtents, chInputData, true);
    ot::RemeshPartition partition = stage == RefinementStage::REFINE ? ot::SurrogateInByOut : ot::SurrogateOutByIn;
    const bool isCoarsenStage = (stage == RefinementStage::COARSEN);
    auto newDA = chRefine.getRefineSubDA(dtree, chInputData->mesh_def.sfc_tol, RefinementStrategy::FULL_DOMAIN,
                                         partition);
    DAChanged = true;
    chRefine.initPetscIntergridTransfer();
    chRefine.petscIntergridTransfer(newDA, dtree, vec.v, vec.ndof, stage);
    chRefine.petscIntergridTransfer(newDA, dtree, linearVec, 1, stage);

    if (doCellIntergrid) {
        std::vector<double> tempCell(octDA->getLocalElementSz(), 1);
        double beforeMass = 0;
        for(int i =0 ; i < cellmass.size(); i+=2){
            beforeMass += cellmass[i];
        }
        std::cout << "Before local mass = " << beforeMass << "\n";
        chRefine.cellIntergridTransfer(newDA, dtree, cellmass, CHNodeData::CH_DOF, isCoarsenStage, fem::CellCoarsen::Sum);
        chRefine.cellIntergridTransfer(newDA, dtree, tempCell, 1, isCoarsenStage, fem::CellCoarsen::Sum);
        double afterMass = 0;
        for(int i =0 ; i < cellmass.size(); i+=2){
            afterMass += cellmass[i];
        }
        std::cout << "Before local mass = " << afterMass << "\n";
        for (int i = 0; i < tempCell.size(); i++) {
            cellmass[i * 2 + 0] = cellmass[i * 2 + 0] / tempCell[i];
            cellmass[i * 2 + 1] = cellmass[i * 2 + 1] / tempCell[i];

        }
    }
    chRefine.finializeIntergridTransfer();
    std::swap(octDA, newDA);
    counter++;
    delete newDA;
    TALYFEMLIB::PrintStatus("Total number of nodes = ", octDA->getGlobalNodeSz());
    std::vector<VecInfo> vecsAfter = {VecInfo(linearVec, 1, 0)};
    checkIntergridTransfer(linearVec, octDA, dtree, 1);
    VecDestroy(&linearVec);
    return DAChanged;
}
#endif //DENDRITEKT_CHUTILS_H
