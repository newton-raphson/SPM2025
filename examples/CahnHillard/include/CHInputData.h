#pragma once

#include <talyfem/input_data/input_data.h>

using TALYFEMLIB::ZEROPTV;
using TALYFEMLIB::PrintStatusStream;
using TALYFEMLIB::PrintWarning;
using TALYFEMLIB::PrintInfo;
using TALYFEMLIB::PrintStatus;
using TALYFEMLIB::PrintError;
static const char *ch_varname[]{"phi", "mu"};

/**
 * Parameters for background mesh
 */
struct MeshDef {
    /// two corners of domain
    ZEROPTV channel_min;
    ZEROPTV channel_max;
    /// refinement level
    int refine_lvl_base;
    int refine_lvl_interface;
    double refine_tol;
    DomainInfo fullDADomain; /// The domain from which its carved out.
    DomainInfo physDomain;/// The actual SubDAdomain
    double sfc_tol;
    /**
     * read channel mesh from config
     * @param root
     */
    void read_from_config(const libconfig::Setting &root) {
        refine_lvl_base = (int) root["refine_lvl_base"];
        refine_lvl_interface = (int) root["refine_lvl_interface"];
        refine_tol = (double) root["refine_tol"];
        sfc_tol = (double) root["sfc_tol"];
        channel_min = ZEROPTV((double) root["min"][0], (double) root["min"][1], (double) root["min"][2]);
        channel_max = ZEROPTV((double) root["max"][0], (double) root["max"][1], (double) root["max"][2]);
        fullDADomain.min.fill(0.0);
        fullDADomain.max.fill(*std::max_element(channel_max.data(),channel_max.data()+DIM));

        physDomain.min.fill(0.0);
        std::memcpy(physDomain.max.data(),channel_max.data(),sizeof(double)*DIM);
    }

};

struct CHInputData : public TALYFEMLIB::InputData {
    /// Mesh free
    bool mfree = false;

    /// Time stepper
    double dt;
    double totalT;
    double OutputStartTime = 0.0;
    int OutputInterval = 1;

    int elementOrder = 1;
    double epsilon;
    /// Mesh definition
    MeshDef mesh_def;
    bool useMassConservingInterpolation = false;

    /// PETSC options (in config.txt)
    SolverOptions solverOptionsCH;
    SolverOptions solverOptionsMass;

    CHInputData() : InputData() {}

    ~CHInputData() = default;

    /// read configure from file
    bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
        ReadConfigFile(filename);
        mesh_def.read_from_config(cfg.getRoot()["background_mesh"]);

        solverOptionsCH = read_solver_options(cfg, "solver_options_ch");
        /// Output control
        if (ReadValue("OutputStartTime", OutputStartTime)) {}
        if (ReadValue("OutputInterval", OutputInterval)) {}
        if (ReadValue("useMassConservingInterpolation", useMassConservingInterpolation)) {
            if(useMassConservingInterpolation){
                TALYFEMLIB::PrintInfo("Using mass conserving interpolation");
                solverOptionsMass = read_solver_options(cfg, "solver_options_mass");
                solverOptionsMass.apply_to_petsc_options ("-mass_");
            }
        }
        ReadValueRequired("dt",dt);
        ReadValueRequired("totalT",totalT);

        ReadValueRequired("epsilon",epsilon);
        return true;
    }

    /// check if the input are valid
    bool CheckInputData() {
        /// Matrix free version of the code cannot have pre-conditioner
        if (mfree) {
            if (solverOptionsCH.vals.count("pc_type") == 1) {
                solverOptionsCH.vals.at("pc_type") = "none";
                PrintWarning("mfree = True, changing pc_type to 'none' automatically!");
            }
        }
        return true;
    }



};
