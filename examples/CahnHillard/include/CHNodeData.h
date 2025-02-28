//
// Created by maksbh on 1/21/19.
//
#pragma once
#include <stdexcept>
class CHNodeData {
public:
    double phi;
    double phi_prev;
    double mu;
    double mu_prev;
    static constexpr int PHI_DOF = 0;
    static constexpr int MU_DOF = 1;
    static constexpr int CH_DOF = 2;


    double &value(int index) {
        switch (index) {
            case 0: return phi;
            case 1: return mu;
            case 2: return phi_prev;
            case 3: return mu_prev;
            default: throw std::runtime_error("Invalid CHNodeData index");
        }
    }

    inline double value(int index) const {
        return const_cast<CHNodeData *>(this)->value(index);
    }

    static const char *name(int index) {
        switch (index) {
            case 0: return "phi";
            case 1: return "mu";
            case 2: return "phi_prev";
            case 3: return "mu_prev";
            default: throw std::runtime_error("Invalid CHNodeData index");
        }
    }

    static int valueno() {
        return 4;
    }
};
