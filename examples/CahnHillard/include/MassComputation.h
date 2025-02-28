//
// Created by maksbh on 1/27/23.
//

#ifndef DENDRITEKT_MASSCOMPUTATION_H
#define DENDRITEKT_MASSCOMPUTATION_H

#include "Traversal/Traversal.h"
#include "DendriteUtils.h"
#include <DataTypes.h>

class MassCalculator: Traversal {

    std::vector<double> & elemental_mass_;
    std::vector<double> & elemental_massLocal_;
    int elemID = 0;
    double mass_local_ = 0;
    static constexpr int numItg = 1u << DIM;
    static constexpr int numChild = 1u << DIM;
public:
    MassCalculator(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch, const DomainExtents &domain,
                   std::vector<double> & elemental_mass,std::vector<double> & elemental_massLocal);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    double getTotalMass()const;

};

MassCalculator::MassCalculator(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch,
                               const DomainExtents &domain,std::vector<double> & elemental_mass,std::vector<double> & elementalMassLocal)
        :Traversal(octDA,treePart,ch,domain),elemental_mass_(elemental_mass),elemental_massLocal_(elementalMassLocal){

    elemental_mass_.resize(octDA->getLocalElementSz()*CHNodeData::CH_DOF*numItg);
    elemental_massLocal_.resize(octDA->getLocalElementSz()*CHNodeData::CH_DOF*numChild);
    this->traverse();
}

void MassCalculator::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values){
    double val_c[CHNodeData::CH_DOF];
    double lphi = 0;
    double lmu = 0;

    while (fe.next_itg_pt()) {
        calcValueFEM(fe, CHNodeData::CH_DOF, values, val_c);
        elemental_mass_[elemID*numItg*CHNodeData::CH_DOF + fe.cur_itg_pt_num()*CHNodeData::CH_DOF + (CHNodeData::PHI_DOF)] = val_c[0];
        elemental_mass_[elemID*numItg*CHNodeData::CH_DOF + fe.cur_itg_pt_num()*CHNodeData::CH_DOF + (CHNodeData::MU_DOF)] = val_c[1];
        lphi += val_c[CHNodeData::PHI_DOF];
        lmu += val_c[CHNodeData::MU_DOF];
        mass_local_ += val_c[CHNodeData::PHI_DOF]*fe.detJxW();
    }
    elemental_massLocal_[elemID*numChild*CHNodeData::CH_DOF + m_childNum*CHNodeData::CH_DOF + (CHNodeData::PHI_DOF)] = lphi/fe.n_itg_pts();
    elemental_massLocal_[elemID*numChild*CHNodeData::CH_DOF + m_childNum*CHNodeData::CH_DOF + (CHNodeData::MU_DOF)] = lmu/fe.n_itg_pts();


    elemID++;

}
double MassCalculator::getTotalMass() const{
    double mass_total;
    MPI_Reduce(&mass_local_, &mass_total, 1, MPI_DOUBLE, MPI_SUM, 0, this->m_octDA->getCommActive());

    return mass_total;
}

#endif //DENDRITEKT_MASSCOMPUTATION_H
