//
// Created by maksbh on 1/27/23.
//

#ifndef DENDRITEKT_MASSCOMPUTATION_H
#define DENDRITEKT_MASSCOMPUTATION_H

#include "Traversal/Traversal.h"
#include <DataTypes.h>

class MassCalculator: Traversal {

    std::vector<double> & elemental_mass_;
    std::vector<double> & elemental_massLocal_;
    int elemID = 0;
    static constexpr int numItg = 1u << DIM;
    static constexpr int numChild = 1u << DIM;

public:
    MassCalculator(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch, const DomainExtents &domain,
                   std::vector<double> & elemental_mass,std::vector<double> & elemental_massLocal);

    void traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values) override;

};

MassCalculator::MassCalculator(DA *octDA, const std::vector<TREENODE> & treePart, const VecInfo & ch,
                               const DomainExtents &domain,std::vector<double> & elemental_mass,std::vector<double> & elemental_massLocal)
        :Traversal(octDA,treePart,ch,domain),elemental_mass_(elemental_mass),elemental_massLocal_(elemental_massLocal){

    elemental_mass_.resize(octDA->getLocalElementSz()*numItg,0);
    elemental_massLocal_.resize(octDA->getLocalElementSz()*numChild,0);

    this->traverse();
}

void MassCalculator::traverseOperation(TALYFEMLIB::FEMElm & fe, const PetscScalar * values){
    double val_c;
    double lmass = 0;
//    fe.refill(0,-1);
//    fe.refill(0,-1);
//    while(fe.next_itg_pt()){
//
//        std::cout <<fe.nbf() << "\n";
//        for(int a = 0; a < fe.nbf(); a++){
//            std::cout << fe.N(a) << " " << fe.detJxW() << " " << fe.detJxW()/fe.jacc() << " ";
//        }
//        std::cout << "\n";
//    }

    while (fe.next_itg_pt()) {
        calcValueFEM(fe, 1, values, &val_c);
        elemental_mass_[elemID*numItg + fe.cur_itg_pt_num()] = val_c;
        lmass += val_c;
    }
    elemental_massLocal_[elemID*numChild + m_childNum] = lmass/fe.n_itg_pts();
    elemID++;

}

#endif //DENDRITEKT_MASSCOMPUTATION_H
