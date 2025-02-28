//
// Created by maksbh on 1/27/23.
//

#ifndef DENDRITEKT_UTILS_H
#define DENDRITEKT_UTILS_H

#include "Traversal/Refinement.h"
#include "SubDA/SubDomain.h"
#include "Refinement.h"
static void performRefinement(DA *& octDA,DistTREE & distTree, VecInfo & ch, DomainExtents & domain, SubDomain & subDA, std::vector<double> & localMass,
                              std::vector<double> & localMass_perChild){
    double massSumBefore = std::accumulate(localMass.begin(),localMass.end(),0.0);
    std::vector<double> coarseFlag;
    MRefine refine(octDA,distTree.getTreePartFiltered(),domain,coarseFlag);
    DA *newDA = refine.getRefineSubDA (distTree,0.001,RefinementStrategy::FULL_DOMAIN,ot::RemeshPartition::SurrogateOutByIn);
    if(newDA == nullptr){
        return;
    }
    static constexpr int numItg = 1u << DIM;
    static constexpr int numChild = 1u << DIM;
    refine.initPetscIntergridTransfer();
    refine.petscIntergridTransfer(newDA,distTree,ch.v,1,RefinementStage::COARSEN);
    refine.cellIntergridTransfer(newDA,distTree,localMass,numItg,true,fem::CellCoarsen::Sum);
//    std::cout << "Before -----------\n";
//    for(int i = 0; i < octDA->getLocalElementSz(); i++){
//        std::cout << "( ";
//        for(int j = 0; j < numChild; j++){
//            std::cout << localMass_perChild[i*numChild + j] << " ";
//        }
//        std::cout << ")\n";
//    }
    refine.cellIntergridTransfer(newDA,distTree,localMass_perChild,numChild,true,fem::CellCoarsen::Sum);
//    std::cout << "After -----------\n";
//    for(int i = 0; i < newDA->getLocalElementSz(); i++){
//        std::cout << "( ";
//        for(int j = 0; j < numChild; j++){
//            std::cout << localMass_perChild[i*numChild + j] << " ";
//        }
//        std::cout << ")\n";
//    }
    refine.cellIntergridTransfer(newDA,distTree,coarseFlag,1,true,fem::CellCoarsen::Sum);
    refine.finializeIntergridTransfer();
    std::swap(octDA,newDA);
    delete newDA;
    subDA.finalize(octDA,distTree.getTreePartFiltered(),domain);
    double massSumAfter = std::accumulate(localMass.begin(),localMass.end(),0.0);
    TALYFEMLIB::PrintStatus("Local mass sum = ", massSumBefore, " ", massSumAfter);
    assert(localMass.size() == coarseFlag.size()*numItg);
    for(int i = 0; i < coarseFlag.size(); i++){
        if(coarseFlag[i] == numChild){

            double sum = 0;
            for(int j = 0; j < numItg; j++) {
                sum += localMass[i*numItg + j];
            }
            for (int j = 0; j < numItg; j++) {
                localMass[i * numItg + j] = sum/(numItg*numChild);
            }
        }
        else {
            if(coarseFlag[i] != 1){
                throw std::runtime_error("Not possible");
            }
            for (int j = 0; j < numItg; j++) {
                localMass[i * numItg + j] /= coarseFlag[i];
            }
        }


    }
}

std::tuple<double,double> findDifference(std::vector<double> & vec1,std::vector<double> & vec2, MPI_Comm comm){
    assert(vec1.size() == vec2.size());
    double maxDifference = -std::numeric_limits<double>::infinity();
    double totalDifference = 0;
    for(int i =0; i < vec1.size();i++){
        double diff = std::abs(vec1[i] - vec2[i]);
        maxDifference = std::max(diff,maxDifference);
        totalDifference += diff;
    }
    double globalMaxDifference,globalTotalDifference;
    MPI_Reduce(&maxDifference,&globalMaxDifference,1,MPI_DOUBLE,MPI_MAX,0,comm);
    MPI_Reduce(&totalDifference,&globalTotalDifference,1,MPI_DOUBLE,MPI_SUM,0,comm);
    return std::make_tuple(globalMaxDifference,globalTotalDifference);
}
#endif //DENDRITEKT_UTILS_H
