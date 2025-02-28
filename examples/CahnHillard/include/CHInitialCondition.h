//
// Created by maksbh on 2/11/23.
//

#ifndef DENDRITEKT_CHINITIALCONDITION_H
#define DENDRITEKT_CHINITIALCONDITION_H

#include "CHInputData.h"
#include "CHNodeData.h"
enum CHInitType:int{
    RANDOM = 0,
    SINOSOIDAL = 1
};
class CHInitialCondition{
    const CHInitType chInitType_;

    void randomInitialCondition(const double * x,  double * var);

public:
    CHInitialCondition(const CHInitType chInitType);

    void setInitialCondition(DA * octDA, Vec initialConditionVec);

};

CHInitialCondition::CHInitialCondition(const CHInitType chInitType)
:chInitType_(chInitType){

}
void CHInitialCondition::randomInitialCondition(const double * x, double * var) {
    double r = rand()/(RAND_MAX*1.0);
    var[CHNodeData::PHI_DOF] = (2*r - 1)*0.1;
//    var[CHNodeData::PHI_DOF] = abs(cos(M_PI*x[0]));
    var[CHNodeData::MU_DOF] = 0.0;
}

void CHInitialCondition::setInitialCondition(DA * octDA, Vec initialConditionVec) {
    std::function<void(const double *, double *)> initial_condition = [&](const double *x, double *var) {
        if(chInitType_ == CHInitType::RANDOM){
            this->randomInitialCondition(x,var);
        }
        else{
            throw std::runtime_error("Inital Condition Not supported");
        }

    };
    octDA->petscSetVectorByFunction(initialConditionVec, initial_condition, false, false, CHNodeData::CH_DOF);


}

#endif //DENDRITEKT_CHINITIALCONDITION_H
