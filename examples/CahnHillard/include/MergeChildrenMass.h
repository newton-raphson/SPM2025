//
// Created by maksbh on 2/28/23.
//

#ifndef DENDRITEKT_MERGECHILDRENMASS_H
#define DENDRITEKT_MERGECHILDRENMASS_H
#include "Traversal/Traversal.h"
#include "DendriteUtils.h"
#include <DataTypes.h>
#include "CHNodeData.h"

class MergeChildrenMass: Traversal {

    const std::vector<double> & vecIn_;
    std::vector<double> & vecOut_;
    const std::vector<double> & coarsenInformation_;
    int elemID = 0;

    static constexpr int numItg = 1u << DIM;
    static constexpr int numChild = 1u << DIM;
    const int vecdof_;
public:
    MergeChildrenMass(DA *octDA, const std::vector<TREENODE> & treePart, const DomainExtents &domain,
                      std::vector<double> & vecOut,const std::vector<double> & vecIn,
                      const std::vector<double> & coarsenInformation,const int vecDof);

    void traverseOperation(TALYFEMLIB::FEMElm & fe) override;



};

MergeChildrenMass::MergeChildrenMass(DA *octDA, const std::vector<TREENODE> & treePart,
                               const DomainExtents &domain,std::vector<double> & vecOut,const std::vector<double> & vecIn,
                               const std::vector<double> & coarsenInformation,const int vecDof)
        :Traversal(octDA,treePart,domain),vecIn_(vecIn),vecOut_(vecOut),coarsenInformation_(coarsenInformation),vecdof_(vecDof){

    vecOut.resize(vecIn.size(),0.0);
    assert(coarsenInformation_.size()*numChild*vecDof == vecOut.size());
    this->traverse();
}

void MergeChildrenMass::traverseOperation(TALYFEMLIB::FEMElm & fe){

    double sum[CHNodeData::CH_DOF] = {0.,0.};
    for(int i = 0; i < numChild; i++){
        sum[CHNodeData::PHI_DOF] += vecIn_[elemID*numChild*CHNodeData::CH_DOF + i*CHNodeData::CH_DOF + CHNodeData::PHI_DOF];
        sum[CHNodeData::MU_DOF] += vecIn_[elemID*numChild*CHNodeData::CH_DOF + i*CHNodeData::CH_DOF + CHNodeData::MU_DOF];
    }
    vecOut_[elemID*numChild*CHNodeData::CH_DOF + m_childNum*CHNodeData::CH_DOF + CHNodeData::PHI_DOF] = sum[CHNodeData::PHI_DOF]/numChild;
    vecOut_[elemID*numChild*CHNodeData::CH_DOF + m_childNum*CHNodeData::CH_DOF + CHNodeData::MU_DOF] = sum[CHNodeData::MU_DOF]/numChild;


    elemID++;

}

#endif //DENDRITEKT_MERGECHILDRENMASS_H
