//
// Created by maksbh on 2/28/23.
//

#ifndef DENDRITEKT_MASSCALCULATORTOTAL_H
#define DENDRITEKT_MASSCALCULATORTOTAL_H
#include "Traversal/Traversal.h"
#include <DataTypes.h>

class MassCalculatorTotal: Traversal {

    double mass_local_ = 0;
    int elemID = 0;
    static constexpr int numItg = 1u << DIM;
    static constexpr int numChild = 1u << DIM;

public:
    MassCalculatorTotal(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch, const DomainExtents &domain);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

    double getTotalMass();

};

MassCalculatorTotal::MassCalculatorTotal(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch,
                                         const DomainExtents &domain)
        :Traversal(octDA,treePart,ch,domain){

    this->traverse();
}

void MassCalculatorTotal::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values){
    double val_c;
    double lmass = 0;


    while (fe.next_itg_pt()) {
        calcValueFEM(fe, 1, values, &val_c);
        lmass += val_c;
        mass_local_ += val_c*fe.detJxW();
    }

    elemID++;

}

double MassCalculatorTotal::getTotalMass(){
    double mass_total;
    MPI_Reduce(&mass_local_, &mass_total, 1, MPI_DOUBLE, MPI_SUM, 0, this->m_octDA->getCommActive());

    return mass_total;
}
#endif //DENDRITEKT_MASSCALCULATORTOTAL_H
