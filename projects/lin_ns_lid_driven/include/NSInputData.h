//
// Created by maksbh on 6/14/20.
//

#ifndef DENDRITEKT_NSINPUTDATA_H
#define DENDRITEKT_NSINPUTDATA_H

#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
#include <point.h>

/**
 * Parameters for background mesh
 * You can either pass a scaling Factor directly
 * or pass a max and get the scalingFactor.
 */
struct MeshDef : public DomainInfo {

    /// refinement level
    DENDRITE_UINT baseLevel = 0;
    DENDRITE_UINT refineLevelBoundary = 0;


    /// Usage: dimensions which are non-periodic are set to zero.  Periodic dimensions are set to ratio of their
    // length to largest length (in some direction).  Periodic dimension is always set to a power of 2.
    /// e.g for a 2d mesh with dimensions [4,2],
    /// 1) setting periodicScale to [0.0,0.5] would set Y periodicity
    /// 2) setting periodicScale to [1.0, 0.0] would set X periodicity
    /// 3) setting periodicScale to [1.0, 0.5] would set both X and Y periodicity
    /// 4) setting periodicScale to [0.0, 0.0] would no periodicity
    /// Note that the values passed are the ratios to the largest dimension, and powers of 2.
    TALYFEMLIB::ZEROPTV periodicScale; ///< Variable for specifying periodicity, iniliazed to zero for non-periodicity


    /**
     * read channel mesh from config
     * @param root
     */
    void read_from_config(const libconfig::Setting &root) {
        baseLevel = (DENDRITE_UINT) root["baseLevel"];
        if (root.exists("refineLevelBoundary")) {
            refineLevelBoundary = (DENDRITE_UINT) root["refineLevelBoundary"];
        }
        for (DENDRITE_UINT dim = 0; dim < DIM; dim++) {
            min[dim] = (DENDRITE_REAL) root["min"][dim];
            max[dim] = (DENDRITE_REAL) root["max"][dim];
        }


        if (root.exists("periodicBoundariesAndScale")) {
            periodicScale = TALYFEMLIB::ZEROPTV ((double) root["periodicBoundariesAndScale"][0], (double) root["periodicBoundariesAndScale"][1],
                                                 (double) root["periodicBoundariesAndScale"][2]);
        } else {
            periodicScale = {0.0, 0.0, 0.0};
        }

    }
};

class NSInputData : public TALYFEMLIB::InputData {
public:

    enum TIME_STEPPING : u_short {
        BACKWARD_EULER = 0,
        CRANK_NICHOLSON = 1,
        BDF2 = 2,
    };



    /// Mesh definition
    MeshDef meshDef;
    /// Basis function
    std::string bfStr;
    /// Matrix  free
    bool mfree = false;
    /// Time stepper
    TIME_STEPPING timeStepping = TIME_STEPPING::CRANK_NICHOLSON;
    std::vector<double> dt;
    std::vector<double> totalT;

    /// ns parameters

    DENDRITE_REAL Re;
    DENDRITE_REAL Re_tmp;
    bool DiffFineTerms = false;
    bool SecondViscosity = false;
    bool DoReRamping = false;
    DENDRITE_REAL RampingTime;
    DENDRITE_REAL RampingRe;

    SolverOptions solverOptionsNS;
    bool ifMMS = false;
    bool dump_vec = false;

    DomainInfo physDomain;

    std::string stlFileName;
    RetainSide stlRetainSide;
    bool tempRetainSide = false;


    ///< Checkpointing options
    DENDRITE_UINT checkpointFrequency;
    DENDRITE_UINT numberOfBackups;

    DENDRITE_UINT OutputSpan; // number of timesteps between each data write (NOT time)


    /// region refinement
    DENDRITE_REAL rectLowerLeft[DIM];
    DENDRITE_REAL rectUpperRight[DIM];
    DENDRITE_UINT rect_lvl;


    DENDRITE_UINT VelocityExtrapolationOrder = 1;


    static TIME_STEPPING readTimeStepper(libconfig::Setting &root,
                                         const char *name) {
        std::string str;
        /// If nothing specified stays stabilizedNS
        if (root.lookupValue(name, str)) {
            if (str == "CN") {
                return TIME_STEPPING::CRANK_NICHOLSON;
            } else if (str == "BE") {
                return TIME_STEPPING::BACKWARD_EULER;
            } else if (str == "BDF2") {
                return TIME_STEPPING::BDF2;
            } else {
                throw TALYFEMLIB::TALYException() << "Unknown time stepping for NS: " << name << str;
            }
        } else {
            throw TALYFEMLIB::TALYException() << "Must specify timeStepping: CN, BE or BDF2";
        }

    }

    bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
        ReadConfigFile(filename);
        /// mesh size and level
        meshDef.read_from_config(cfg.getRoot()["background_mesh"]);

        if (!ReadValue("stlFileName", stlFileName)) {return false;}
        if (!ReadValue("stlRetainInside", tempRetainSide)) {return false;}
        stlRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;

        // Checkpointing options
        checkpointFrequency = 1.0;
        numberOfBackups = 2.0;
        ReadValue("checkpointFrequency", checkpointFrequency);
        ReadValue("checkpointNumberOfBackups", numberOfBackups);

        /// File writing
        if (ReadValue("OutputSpan", OutputSpan)) {}


        for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
            physDomain.min[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMin"][dim];
            physDomain.max[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMax"][dim];
        }


        for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
            rectLowerLeft[dim] = (DENDRITE_REAL) cfg.getRoot()["rectLowerLeft"][dim];
            rectUpperRight[dim] = (DENDRITE_REAL) cfg.getRoot()["rectUpperRight"][dim];
        }
        if (ReadValue("rect_lvl", rect_lvl)) {}


        /// basis function order
        if (ReadValue("basisFunction", bfStr)) {}
        basisFunction = bfStr.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bfStr);
        /// timestep control
        timeStepping = readTimeStepper(cfg.getRoot(), "TimeStepper");
        {
            if (cfg.exists("dt_V")) {
                const libconfig::Setting &dt_ = cfg.getRoot()["dt_V"];
                for (int i = 0; i < dt_.getLength(); ++i) {
                    dt.push_back(dt_[i]);
                }
            } else {
                double dt_const;
                ReadValueRequired("dt", dt_const);
                dt.push_back(dt_const);
            }

            if (cfg.exists("totalT_V")) {
                const libconfig::Setting &totalT_ = cfg.getRoot()["totalT_V"];
                for (int i = 0; i < totalT_.getLength(); ++i) {
                    totalT.push_back(totalT_[i]);
                }
            } else {
                double totalT_const;
                ReadValueRequired("totalT", totalT_const);
                totalT.push_back(totalT_const);
            }
        }

        /// NS parameter
        ReadValueRequired("Re", Re);
        Re_tmp = Re;
        if (ReadValue("DoReRamping", DoReRamping)) {}
        if (ReadValue("RampingRe", RampingRe)) {}
        if (ReadValue("RampingTime", RampingTime)) {}
        if (ReadValue("mfree", mfree)) {}
        if (ReadValue("DiffFineTerms", DiffFineTerms)) {}
        if (ReadValue("SecondViscosity", SecondViscosity)) {}
        if (ReadValue("MMS", ifMMS)) {}
        if (ReadValue("dump_vec", dump_vec)) {}
        if (ReadValue("VelocityExtrapolationOrder", VelocityExtrapolationOrder)) {}
        solverOptionsNS = read_solver_options(cfg, "solver_options");
        return true;
    }

    /// check if the input are valid
    bool CheckInputData() {
        /// Matrix free version of the code cannot have pre-conditioner
        if (mfree) {
            if (solverOptionsNS.vals.count("pc_type") == 1) {
                solverOptionsNS.vals.at("pc_type") = "none";
                TALYFEMLIB::PrintWarning("mfree = True, changing pc_type to 'none' automatically!");
            }
        }
        return true;
    }
};

#endif //DENDRITEKT_NSINPUTDATA_H
