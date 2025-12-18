#pragma once


// copy from Thermoelasticity branch
#include <exception>
#include <assert.h>
#include <DataTypes.h>

class LENodeData {
public:
    /// LE degrees of freedom
    static constexpr unsigned int LE_DOF = DIM;
    /// Store the values of the degrees of freedom at the current timestep
    DENDRITE_REAL a[LE_DOF];
    DENDRITE_REAL v[LE_DOF];
    DENDRITE_REAL u[LE_DOF];
    /// Store the values at the degrees of freedom at the previous timestep (n-1)
    DENDRITE_REAL a_pre1[LE_DOF];
    DENDRITE_REAL v_pre1[LE_DOF];
    DENDRITE_REAL u_pre1[LE_DOF];
    /// Store the values at the degrees of freedom at the second previous timestep (n-2)
    DENDRITE_REAL a_pre2[LE_DOF];
    DENDRITE_REAL v_pre2[LE_DOF];
    DENDRITE_REAL u_pre2[LE_DOF];

    static constexpr unsigned int HT_DOF = 1;
    /// number of  variables in the LENodeData
    static constexpr unsigned int NUM_VARS = LE_DOF * 3 + HT_DOF; // we count for a, v, u, T

    /// SBM
    DENDRITE_REAL nodeid;

    enum Vars : DENDRITE_UINT {
        /// Current time step
        AX = 0,
        AY = 1,
#if (DIM ==3)
        AZ = 2,
#endif
        VX = DIM,
        VY = DIM+1,
#if (DIM ==3)
        VZ = DIM+2,
#endif
        UX = 2*DIM,
        UY = 2*DIM+1,
#if (DIM ==3)
        UZ = 2*DIM+2,
#endif
        TEMPERATURE = 3*DIM,

        /// n-1 timestep
        AX_PRE1 = NUM_VARS,
        AY_PRE1 = NUM_VARS+1,
#if (DIM ==3)
        AZ_PRE1 = NUM_VARS+2,
#endif
        VX_PRE1 = NUM_VARS+DIM,
        VY_PRE1 = NUM_VARS+DIM+1,
#if (DIM ==3)
        VZ_PRE1 = NUM_VARS+DIM+2,
#endif
        UX_PRE1 = NUM_VARS+2*DIM,
        UY_PRE1 = NUM_VARS+2*DIM+1,
#if (DIM ==3)
        UZ_PRE1 = NUM_VARS+2*DIM+2,
#endif
        TEMPERATURE_PRE1 = NUM_VARS+3*DIM,

        /// n-2 timestep
        AX_PRE2 = 2*NUM_VARS,
        AY_PRE2 = 2*NUM_VARS+1,
#if (DIM ==3)
        AZ_PRE2 = 2*NUM_VARS+2,
#endif
        VX_PRE2 = 2*NUM_VARS+DIM,
        VY_PRE2 = 2*NUM_VARS+DIM+1,
#if (DIM ==3)
        VZ_PRE2 = 2*NUM_VARS+DIM+2,
#endif
        UX_PRE2 = 2*NUM_VARS+2*DIM,
        UY_PRE2 = 2*NUM_VARS+2*DIM+1,
#if (DIM ==3)
        UZ_PRE2 = 2*NUM_VARS+2*DIM+2,
#endif
        TEMPERATURE_PRE2 = 2*NUM_VARS+3*DIM,

        TEMPERATURE_SUM =  2*NUM_VARS+3*DIM+1,
        /// SBM
        NODE_ID = TEMPERATURE_PRE2 +1

    };

    LENodeData() {
        std::memset(a, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(v, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(u, 0, sizeof(DENDRITE_REAL) * (LE_DOF+HT_DOF));
        std::memset(a_pre1, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(v_pre1, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(u_pre1, 0, sizeof(DENDRITE_REAL) * (LE_DOF+HT_DOF));
        std::memset(a_pre2, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(v_pre2, 0, sizeof(DENDRITE_REAL) * LE_DOF);
        std::memset(u_pre2, 0, sizeof(DENDRITE_REAL) * (LE_DOF+HT_DOF));
    }



    /**

     * Returns reference to the given value in the object
     *
     * @param index the index of the desired item
     * @return reference to the desired data item
     */
    inline double &value(int index) {

        if (index >= 0 && index < LE_DOF) {
            return a[index];
        }
        if (index >= LE_DOF && index < LE_DOF * 2) {
            return v[index - LE_DOF];
        }
        if (index >= LE_DOF * 2 && index < NUM_VARS) {
//        std::cout<< "index - 2 * LE_DOF = " << index - 2 * LE_DOF << "\n";
            return u[index - 2 * LE_DOF];
        }

        if (index >= NUM_VARS  && index < NUM_VARS +LE_DOF) {
            return a_pre1[index - NUM_VARS];
        }
        if (index >= NUM_VARS + LE_DOF && index < NUM_VARS + LE_DOF*2) {
            return v_pre1[index - NUM_VARS-LE_DOF];
        }
        if (index >= NUM_VARS + LE_DOF * 2 && index < 2*NUM_VARS) {
            return u_pre1[index - NUM_VARS- 2*LE_DOF];
        }

        if (index >= 2*NUM_VARS && index < 2*NUM_VARS +LE_DOF) {
            return a_pre2[index - 2 * NUM_VARS];
        }
        if (index >= 2*NUM_VARS + LE_DOF && index < 2*NUM_VARS + LE_DOF*2) {
            return v_pre2[index - 2 * NUM_VARS - LE_DOF];
        }
        if (index >= 2*NUM_VARS + LE_DOF *2 && index < 3*NUM_VARS) {
            return u_pre2[index - 2 * NUM_VARS - 2*LE_DOF];
        }

        if (index ==NODE_ID){
            //std::cout<< "nodeid = " << nodeid << "\n";
            return nodeid;
        }

        TALYFEMLIB::TALYException() << "Invalid variable index!";
    }

    inline double value(int index) const {
        return const_cast<LENodeData *>(this)->value(index);
    }

    /**
     * Returns the name of the given data value in the object
     * @param index the index of the desired item
     * @return name of the specified data item
     */
    static const char *name(int index) {
        switch (index) {
            case AX: return "AX";
            case AY: return "AY";
            case VX: return "VX";
            case VY: return "VY";
            case UX: return "UX";
            case UY: return "UY";
#if(DIM == 3)
            case AZ: return "AZ";
            case VZ: return "VZ";
            case UZ: return "UZ";
#endif

            case NODE_ID: return "nodeID";
            default: throw std::runtime_error("Invalid LENodeData index");
        }
        return nullptr;
    }

    /**
     * Returns the number of the data items in the object
     * @return number of the data items in the object
     */
    static int valueno() {
        return NUM_VARS;
    }
};
