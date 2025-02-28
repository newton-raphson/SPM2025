//
// Created by maksbh on 2/11/23.
//

#ifndef DENDRITEKT_CHPOSTPROCESSING_H
#define DENDRITEKT_CHPOSTPROCESSING_H
#include "CHNodeData.h"
#include "Traversal/Traversal.h"
class CHPostProcessing:Traversal{
    enum EnergyAndMass:int{
        MASS = 0,
        FREE_ENERGY = 1,
        INTERFACIAL_ENRGY = 2,
        MAX = 3
    };
    double localQuantity_[EnergyAndMass::MAX];
public:
    CHPostProcessing(DA * octDA, const DistTREE & distTree, const DomainExtents & domainExtents, const VecInfo & chVec);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    void getGlobalQuantity(double * globalQuantity, MPI_Comm comm);

};

CHPostProcessing::CHPostProcessing(DA *octDA, const DistTREE & distTree, const DomainExtents & domainExtents,
                                   const VecInfo & chVec)
    : Traversal(octDA,distTree.getTreePartFiltered(),chVec,domainExtents){
    std::memset(localQuantity_,0, sizeof(double )*EnergyAndMass::MAX);
    this->traverse();

}

void CHPostProcessing::traverseOperation(TALYFEMLIB::FEMElm &fe, const PetscScalar *values) {

    double valc[CHNodeData::CH_DOF];
    double vald[CHNodeData::CH_DOF*DIM];
    while(fe.next_itg_pt()){
        calcValueFEM(fe,CHNodeData::CH_DOF,values,valc);
        calcValueDerivativeFEM(fe,CHNodeData::CH_DOF,values,vald);
        const double phi = valc[CHNodeData::PHI_DOF];
        localQuantity_[EnergyAndMass::MASS] += phi*fe.detJxW();
        localQuantity_[EnergyAndMass::FREE_ENERGY] += 0.25 * (phi *phi - 1) * (phi*phi - 1)*fe.detJxW();
        for(int d = 0; d < DIM; d++){
            localQuantity_[EnergyAndMass::INTERFACIAL_ENRGY] += 0.5 * (vald[DIM*CHNodeData::PHI_DOF + d] * vald[DIM*CHNodeData::PHI_DOF + d])*fe.detJxW();
        }

    }
}

void CHPostProcessing::getGlobalQuantity(double * globalQuantity,MPI_Comm comm){
    MPI_Reduce(localQuantity_,globalQuantity,EnergyAndMass::MAX,MPI_DOUBLE,MPI_SUM,0,comm);

}



#endif //DENDRITEKT_CHPOSTPROCESSING_H
