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


struct UniMeshInputData : public TALYFEMLIB::InputData {

  std::string bf_str;

  DENDRITE_UINT refine_lvl;

  /// Mesh definition
//  MeshDef mesh_def;

//  /// PETSC options (in config.txt)
//  SolverOptions solverOptionsSSHT;

  /// dump_vec (for regression test)
  bool dump_vec = false;

  /// mfree
//  bool mfree = false;

  /// Geometry parameters
  // Cube Domain min and max coordinates
  DENDRITE_REAL cubeDomainMin[DIM];
  DENDRITE_REAL cubeDomainMax[DIM];

  // Physical Domain min and max coordinates
  DENDRITE_REAL physDomainMin[DIM];
  DENDRITE_REAL physDomainMax[DIM];

  // 2 Dim example parameters
#if(DIM == 2)
  // Ex-1 Circle center, radius and retain side (outside: True, inside: False)
  DENDRITE_REAL circleCenter[DIM];
  DENDRITE_REAL circleRadius;
  RetainSide circleRetainSide;


  // Ex-2 Rectangle lower left, upper right corners and retain side (outside: True, inside: False)
  DENDRITE_REAL rectLowerLeft[DIM];
  DENDRITE_REAL rectUpperRight[DIM];
  RetainSide rectRetainSide;

  // Ex-3 .msh file name and retain side (outside: True, inside: False)
  std::string mshFileName;
  RetainSide mshRetainSide;
#endif
#if(DIM == 3)
  // 3 Dim example parameters
  // Ex-1 Sphere center, radius and retain side (outside: True, inside: False)
  DENDRITE_REAL sphereCenter[DIM];
  DENDRITE_REAL sphereRadius;
  RetainSide sphereRetainSide;

  // Ex-2 Box lower left, upper right corners and retain side (outside: True, inside: False)
  DENDRITE_REAL boxLowerLeft[DIM];
  DENDRITE_REAL boxUpperRight[DIM];
  RetainSide boxRetainSide;

  // Ex-3 .stl file name and retain side (outside: True, inside: False)
  std::string stlFileName;
  RetainSide stlRetainSide;
#endif

  UniMeshInputData() : InputData() {}

  ~UniMeshInputData() = default;

  /// read configure from file
  bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
    /// fill cfg, but don't initialize default fields
    ReadConfigFile(filename);

//    mesh_def.read_from_config(cfg.getRoot()["background_mesh"]);

//    solverOptionsSSHT = read_solver_options(cfg, "solver_options_ssht");

    bool tempRetainSide = false;
    if (ReadValue("dump_vec", dump_vec)) {}

//    if (ReadValue("mfree", mfree)) {}

    if (ReadValue("refine_lvl", refine_lvl)) {}

    if (ReadValue("basisFunction", bf_str)) {}
    basisFunction = bf_str.empty() ? basisFunction : TALYFEMLIB::basis_string_to_enum(bf_str);

    // Read Geometry parameters
    // Read Cube Domain min and max coordinates using ReadArray template function

    for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
        cubeDomainMin[dim] = (DENDRITE_REAL) cfg.getRoot()["cubeDomainMin"][dim];
        cubeDomainMax[dim] = (DENDRITE_REAL) cfg.getRoot()["cubeDomainMax"][dim];
        physDomainMin[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMin"][dim];
        physDomainMax[dim] = (DENDRITE_REAL) cfg.getRoot()["physDomainMax"][dim];
    }
//    if (!ReadArray("cubeDomainMin", cubeDomainMin, DIM)) {return false;}
//    if (!ReadArray("cubeDomainMax", cubeDomainMax, DIM)) {return false;}

    // Read Physical Domain min and max coordinates using ReadArray template function
//    if (!ReadArray("physDomainMin", physDomainMin, DIM)) {return false;}
//    if (!ReadArray("physDomainMax", physDomainMax, DIM)) {return false;}

#if(DIM == 2)
    // Read 2 Dim example parameters
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
        circleCenter[dim] = (DENDRITE_REAL) cfg.getRoot()["circleCenter"][dim];
        rectLowerLeft[dim] = (DENDRITE_REAL) cfg.getRoot()["rectLowerLeft"][dim];
        rectUpperRight[dim] = (DENDRITE_REAL) cfg.getRoot()["rectUpperRight"][dim];
    }
    // Read Ex-1 Circle center, radius and retain side (inside: True, outside: False)
    if (!ReadValue("circleRadius", circleRadius)) {return false;}
    if (!ReadValue("circleRetainInside", tempRetainSide)) { return false;}
      circleRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;

    // Read Ex-2 Rectangle lower left, upper right corners and retain side (inside: True, outside: False)
    if (!ReadValue("rectRetainInside", tempRetainSide)) {return false;}
      rectRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;

    // Read Ex-3 .msh file name and retain side (inside: True, outside: False)
    if (!ReadValue("mshFileName", mshFileName)) {return false;}
    if (!ReadValue("mshRetainInside", tempRetainSide)) {return false;}
      mshRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;
#endif
#if(DIM == 3)
    // Read 3 Dim example parameters
    for (DENDRITE_UINT dim = 0; dim < DIM; dim++){
        sphereCenter[dim] = (DENDRITE_REAL) cfg.getRoot()["sphereCenter"][dim];
        boxLowerLeft[dim] = (DENDRITE_REAL) cfg.getRoot()["boxLowerLeft"][dim];
        boxUpperRight[dim] = (DENDRITE_REAL) cfg.getRoot()["boxUpperRight"][dim];
    }
    // Read Ex-1 Sphere center, radius and retain side (inside: True, outside: False)
    if (!ReadValue("sphereRadius", sphereRadius)) {return false;}
    if (!ReadValue("sphereRetainInside", tempRetainSide)) {return false;}
      sphereRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;

    // Read Ex-2 Box lower left, upper right corners and retain side (inside: True, outside: False)
    if (!ReadValue("boxRetainInside", tempRetainSide)) {return false;}
      boxRetainSide = tempRetainSide ? RetainSide::IN : RetainSide::OUT;

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
