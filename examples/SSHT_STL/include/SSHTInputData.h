#pragma once

#include <iomanip>
#include <talyfem/input_data/input_data.h>
#include <DataTypes.h>
using TALYFEMLIB::ZEROPTV;
using TALYFEMLIB::PrintStatusStream;
using TALYFEMLIB::PrintWarning;
using TALYFEMLIB::PrintInfo;
using TALYFEMLIB::PrintStatus;
using TALYFEMLIB::PrintError;

/**
 * Parameters for background mesh
 */
struct MeshDef: public DomainInfo {
  DENDRITE_UINT refine_lvl;
  void read_from_config(const libconfig::Setting &root){
    refine_lvl = (int) root["refine_lvl"];
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
      min[dim] = (DENDRITE_REAL) root["min"][dim];
      max[dim] = (DENDRITE_REAL) root["max"][dim];
    }
  }
};


struct SSHTInputData : public TALYFEMLIB::InputData {

  std::string bf_str;

  /// Mesh definition
  MeshDef mesh_def;

  /// PETSC options (in config.txt)
  SolverOptions solverOptionsSSHT;

  /// dump_vec (for regression test)
  bool dump_vec = false;

  /// mfree
  bool mfree = false;


  DENDRITE_UINT refine_lvl;


  SSHTInputData() : InputData() {}

  ~SSHTInputData() = default;

    // Cube Domain min and max coordinates
    DENDRITE_REAL cubeDomainMin[DIM];
    DENDRITE_REAL cubeDomainMax[DIM];

    // Physical Domain min and max coordinates
    DENDRITE_REAL physDomainMin[DIM];
    DENDRITE_REAL physDomainMax[DIM];


#if(DIM == 2)
    // Ex-3 .msh file name and retain side (outside: True, inside: False)
    std::string mshFileName;
    RetainSide mshRetainSide;
#endif

#if(DIM == 3)
    // Ex-3 .stl file name and retain side (outside: True, inside: False)
  std::string stlFileName;
  RetainSide stlRetainSide;
#endif



    /// read configure from file
  bool ReadFromFile(const std::string &filename = std::string("config.txt")) {

      bool tempRetainSide = false;
    /// fill cfg, but don't initialize default fields
    ReadConfigFile(filename);

//    mesh_def.read_from_config(cfg.getRoot()["background_mesh"]);




    solverOptionsSSHT = read_solver_options(cfg, "solver_options_ssht");

    if (ReadValue("dump_vec", dump_vec)) {}

    if (ReadValue("mfree", mfree)) {}

    if (ReadValue("basisFunction", bf_str)) {}
    basisFunction = bf_str.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bf_str);

    if (ReadValue("refine_lvl", refine_lvl)) {}

      for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
          cubeDomainMin[dim] = (DENDRITE_REAL) cfg.getRoot()["cubeDomainMin"][dim];
          cubeDomainMax[dim] = (DENDRITE_REAL) cfg.getRoot()["cubeDomainMax"][dim];
          physDomainMin[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMin"][dim];
          physDomainMax[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMax"][dim];
      }


#if(DIM == 2)
        // Read Ex-3 .msh file name and retain side (inside: True, outside: False)
        if (!ReadValue("mshFileName", mshFileName)) {return false;}
        if (!ReadValue("mshRetainInside", tempRetainSide)) {return false;}
        mshRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;
#endif


#if(DIM == 3)
        // Read Ex-3 .stl file name and retain side (inside: True, outside: False)
    if (!ReadValue("stlFileName", stlFileName)) {return false;}
    if (!ReadValue("stlRetainInside", tempRetainSide)) {return false;}
    stlRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;
#endif


    return true;
  }

  /// check if the input are valid
  bool CheckInputData() {
    return true;
  }

};
