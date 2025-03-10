#pragma once
#include "LEInputDataStructs.h"
#include <time.h>
#include <talyfem/input_data/input_data.h>

// construct a kd-tree index:
using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
PointCloud<double>, 3 /* dim */>;

struct Planest
{
  double young;
  double poisson;

  void read_from_config(const libconfig::Setting &root)
  {
    young = root["young"];
    poisson = root["poisson"];
  }
};

struct Lame
{
  double lamda;
  double mu;

  void read_from_config(const libconfig::Setting &root)
  {
    lamda = root["lamda"];
    mu = root["mu"];
  }
};

enum CaseDir : DENDRITE_UINT
{
  RIGHT = 0,
  TOP = 1,

  MAX_DIR_CASE_TYPE = 2
};

struct TractionBC
{
  double traction;
  CaseDir direction;

  void read_from_config(const libconfig::Setting &root)
  {
    traction = root["traction"];
    direction = read_dir_type(root["direction"]);
  }

private:
  static CaseDir read_dir_type(const std::string &str)
  {
    if (str == "right")
    {
      return CaseDir::RIGHT;
    }
    if (str == "top")
    {
      return CaseDir::TOP;
    }
  }
};

struct BottomTractionBC
{
    double traction;

    bool NeumannFromSBM = false;

    void read_from_config(const libconfig::Setting &root)
    {
        traction = root["traction"];
        if (root.exists("NeumannFromSBM")) {
            NeumannFromSBM = root["NeumannFromSBM"];
            if (NeumannFromSBM){
                PrintStatus("Setting NeumannFromSBM = true");
            }
        } else {
            PrintStatus("Default setting NeumannFromSBM = false");
        }
    }

};

struct TractionTopBC
{
  double traction;

  void read_from_config(const libconfig::Setting &root)
  {
    traction = root["traction"];
  }
};

struct Traction4Side
{
  std::vector<double> traction;

  void read_from_config(const libconfig::Setting &root)
  {
    ReadVectorRoot(root, "traction", traction);
  }
};

struct DisplacmentBC
{
  double displacement;

  void read_from_config(const libconfig::Setting &root)
  {
    displacement = root["displacement"];
  }
};

struct CSV_TractionBC
{
    enum Side {
        INVALID = -1,
        X_MINUS = 0,
        X_PLUS = 1,
        Y_MINUS = 2,
        Y_PLUS = 3,
        Z_MINUS = 4,
        Z_PLUS = 5,
    };

    std::string filename;
    Side fixside;
    std::vector<double> shift_pos;

    void read_from_config(const libconfig::Setting &root)
    {
        // Check if the setting exists and is a string before assigning
        if(root.exists("filename") && root["filename"].getType() == libconfig::Setting::TypeString)
        {
            filename = root["filename"].c_str();
        }

        // Assuming FixSide is always a string
        fixside = read_side_type(root["FixSide"].c_str());


        ReadVectorRoot(root, "shift_pos", shift_pos);

    }

private:
    static Side read_side_type(const std::string &str)
    {
        if (str == "x+")
        {
            return X_PLUS;
        }
        else if (str == "x-")
        {
            return X_MINUS;
        }
        else if (str == "y+")
        {
            return Y_PLUS;
        }
        else if (str == "y-")
        {
            return Y_MINUS;
        }
        else if (str == "z+")
        {
            return Z_PLUS;
        }
        else if (str == "z-")
        {
            return Z_MINUS;
        }
        return INVALID; // Return INVALID if no matches
    }
};


/// Declare enum to store the type of LE cases
enum CaseType : DENDRITE_UINT
{
  PLANESTRESS = 0,
  PLANESTRAIN = 1,
  LAME = 2,

  MAX_CASE_TYPE = 3
};

/// Declare enum to store the type of LE BC cases
enum BCCaseType : DENDRITE_UINT
{
  NORMAL_TRACTION = 0,
  DISPLACEMENT_BOTH_SIDE = 1,
  FIXED_AT_WALL = 2,
  ZERO_TRACTION = 3,
  HALF_BEAM = 4,
  TRACT4SIDE = 5,
  BOTTOM_FORCE = 6,
  CSV_FORCE = 7,
  POSITION_DISPLACEMENT = 8,

  MAX_BCCASE_TYPE = 9
};

struct RadialBodyForce
{
  int br_pow;
  double br_v;

  void read_from_config(const libconfig::Setting &root)
  {
    br_pow = root["BR_pow"];
    br_v = root["BR_v"];
  }
};

static const char *caseTypeName[]{"PLANESTRESS", "PLANESTRAIN", "LAME"};

class LEInputData : public TALYFEMLIB::InputData
{
public: // need to put the variable need to use in the other subroutine here!
  static constexpr int nsd = DIM;

  DENDRITE_UINT elemOrder = 1;
  bool ifMatrixFree = false;
  bool ifHessian = false;
  bool BaselvlFromArgument = false;

  CaseType caseType;
  BCCaseType bccaseType;
  Planest planeStress;
  Planest planeStrain;
  TractionBC NormalTraction;
  TractionTopBC HalfBeam;
  BottomTractionBC BottomTract;
  DisplacmentBC DisplacementBothSide;
  Lame lame;
  RadialBodyForce radialbodyforce;
  Traction4Side traction4side;
  CSV_TractionBC CsvForce;

  std::vector<std::pair<double, double>> minmax_cantilever;

#ifdef DEEPTRACE
    const char *model_path;
    double scale=1.0;
#endif
    /// this is just a special note for the onnx based mode for deeptrace
    /// this is not a general purpose variable
    /// region to perform inference
#ifdef DEEPTRACE
    DENDRITE_REAL inference_rectLowerLeft[DIM];
    DENDRITE_REAL inference_rectUpperRight[DIM];
#endif

  /// For passing data to other class
  my_kd_tree_t* kdTree_ = nullptr; // Pointer to the KD-Tree
    PointCloud<double> traction_position_;
    std::vector<ZEROPTV> traction_vector_;
    ZEROPTV shift_;


    /// Linear Elasticity
  TALYFEMLIB::ZEROPTV DomainMax;
  TALYFEMLIB::ZEROPTV DomainMin;

  TALYFEMLIB::ZEROPTV BodyForce;
  double scaleFactor = 1;
  double rho = 1;
  std::vector<std::vector<double>> Cmatrix;


  /// Time stepper
  std::vector<double> dt;
  std::vector<double> totalT;
  double OutputStartTime;
  int OutputInterval;
  int CheckpointInterval = 1;
  int CheckpointNumbackup = 5;

  bool CalcUxError = true;

  /// weakBC parameters
  ParameterDef Cb_e;
  ParameterDef TauN;
  bool weakBCglobal = false;
  double orderOfhb = 1.0;

  std::vector<BoundaryDef> boundary_def;

  /// IBM geometries
  std::vector<IBMGeomDef> ibm_geom_def;

  /// SBM geomtry type
  enum SBMGeo
  {
    ////// all this geometries are for SPM2025 Paper
    BUNNY = 1,
    GYROID = 2,
    PLANE = 3,
    PLANT = 4,
    EIFFEL = 5,
    NONE = 5 // fix bug without the geometries
  };
  /// Declare the geo type of SBM
  SBMGeo SbmGeo = NONE;

  /// Dist Calc in 3D
  enum typeDistCalc
  {
    NORMAL_BASED = 0,
    GP_BASED = 1,
    KD_TREE = 2
  };

  /// Declare the dist calc type
  typeDistCalc DistCalcType;

  /// SBM
  double RatioGPSBM = 1;
  int RelOrderCheckActive = 3;
  int RelOrderInterceptedElement = 3;



    /// [SBM]
  bool InsideSBM = false;
  bool IfAdjointConsistency = false;

    /// Solver options for PETSc are handled in these structures
  SolverOptions solverOptionsLE;

  /// Setup the meshDef object for subDA parameters
  MeshDef mesh_def;
  std::vector<RegionalRefine> region_refine;

  ~LEInputData() = default;

  bool ReadFromFile(const std::string &filename = std::string("config.txt")) // call this in main
  {
    ReadConfigFile(filename);
    ReadValue("elemOrder", elemOrder);
    ReadValue("ifMatrixFree", ifMatrixFree);
    ReadValue("ifHessian", ifHessian);


#ifdef DEEPTRACE
      if (!ReadValue("model_path", model_path)) {return false;}
      if (!ReadValue("scale", scale)) {return false;}
#endif



      if (ReadValue("BaselvlFromArgument", BaselvlFromArgument))
    {
    }

    if (ReadValue("InsideSBM", InsideSBM))
    {
    }

    if (ReadValue("IfAdjointConsistency", IfAdjointConsistency))
    {
    }

      /// always have dim*2 boundary_def in the order of x-, x+, y-, y+, z-, z+
    boundary_def.resize(DIM * 2);
    boundary_def[0].side = BoundaryDef::Side::X_MINUS;
    boundary_def[1].side = BoundaryDef::Side::X_PLUS;
    boundary_def[2].side = BoundaryDef::Side::Y_MINUS;
    boundary_def[3].side = BoundaryDef::Side::Y_PLUS;
#if (DIM == 3)
    boundary_def[4].side = BoundaryDef::Side::Z_MINUS;
    boundary_def[5].side = BoundaryDef::Side::Z_PLUS;
#endif

    caseType = read_LEcase(cfg.getRoot(), "LEcaseType");
    //ReadValueRequired("caseType",str);
    //caseType = static_cast<CaseType>(convertEnumToStrings<caseTypeName,CaseType::MAX_CASE_TYPE>(str.c_str()));

    if (caseType == CaseType::PLANESTRESS)
    {
      planeStress.read_from_config(cfg.getRoot()["planestress"]);
    }
    if (caseType == CaseType::PLANESTRAIN)
    {
      planeStrain.read_from_config(cfg.getRoot()["planestrain"]); // [fix bug]
    }
    if (caseType == CaseType::LAME)
    {
      lame.read_from_config(cfg.getRoot()["lame"]);
    }

    bccaseType = read_LEBCcase(cfg.getRoot(), "LEBCcaseType");
    if (bccaseType == BCCaseType::NORMAL_TRACTION)
    {
      NormalTraction.read_from_config(cfg.getRoot()["NormalTraction"]);

      bool x_minus_wall = false;
      bool y_minus_wall = false;
      bool x_max_wall = false;
      bool y_max_wall = false;
      bool z_minus_wall = false;
      bool z_max_wall = false;

      if (NormalTraction.direction == CaseDir::RIGHT)
      {
        x_max_wall = true;
      }
      if (NormalTraction.direction == CaseDir::TOP)
      {
        y_max_wall = true;
      }

#if (DIM == 2)
      std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall};
#endif
#if (DIM == 3)
      std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall, z_minus_wall, z_max_wall};
#endif
      std::vector<int> traction_dir = {0, 0, 1, 1, 2, 2};

      for (int wallID = 0; wallID < DIM * 2; wallID++)
      {
        if (walls[wallID])
        {
          boundary_def[wallID].Traction(traction_dir[wallID]) = NormalTraction.traction;
          boundary_def[wallID].disp_type = BoundaryDef::Disp_Type::NEUMANN;
        }
      }
    }

      if (bccaseType == BCCaseType::CSV_FORCE){
          CsvForce.read_from_config(cfg.getRoot()["CSV_FORCE"]);

          bool x_minus_wall = true;
          bool y_minus_wall = true;
          bool x_max_wall = true;
          bool y_max_wall = true;
          bool z_minus_wall = true;
          bool z_max_wall = true;


#if (DIM == 2)
          std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall};
#endif
#if (DIM == 3)
          std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall, z_minus_wall, z_max_wall};
#endif

          walls[CsvForce.fixside] = false;

          std::vector<int> traction_dir = {0, 0, 1, 1, 2, 2};

          for (int wallID = 0; wallID < DIM * 2; wallID++)
          {
              if (walls[wallID])
              {
                  boundary_def[wallID].disp_type = BoundaryDef::Disp_Type::NEUMANN;
              }
          }
      }

    if (bccaseType == BCCaseType::HALF_BEAM)
    {
      HalfBeam.read_from_config(cfg.getRoot()["HalfBeam"]);

      bool x_minus_wall = false;
      bool y_minus_wall = false;
      bool x_max_wall = false;
      bool y_max_wall = false;
      bool z_minus_wall = false;
      bool z_max_wall = false;

      y_max_wall = true;

#if (DIM == 2)
      std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall};
#endif
#if (DIM == 3)
      std::vector<bool> walls = {x_minus_wall, x_max_wall, y_minus_wall, y_max_wall, z_minus_wall, z_max_wall};
#endif
      std::vector<int> traction_dir = {0, 0, 1, 1, 2, 2};

      for (int wallID = 0; wallID < DIM * 2; wallID++)
      {
        if (walls[wallID])
        {
          boundary_def[wallID].Traction(traction_dir[wallID]) = HalfBeam.traction;
          boundary_def[wallID].disp_type = BoundaryDef::Disp_Type::NEUMANN;
        }
      }
    }

    if (bccaseType == BCCaseType::TRACT4SIDE)
    {
      traction4side.read_from_config(cfg.getRoot()["Tract4Side"]);

      std::vector<int> traction_dir = {0, 0, 1, 1, 2, 2};

      for (int wallID = 0; wallID < DIM * 2; wallID++)
      {
        boundary_def[wallID].Traction(traction_dir[wallID]) = traction4side.traction[wallID];
        boundary_def[wallID].disp_type = BoundaryDef::Disp_Type::NEUMANN;
      }
    }

    if (bccaseType == BCCaseType::DISPLACEMENT_BOTH_SIDE)
    {
      DisplacementBothSide.read_from_config(cfg.getRoot()["DisplacementBothSide"]);
    }

    /// SubDA (channel parameters)
    mesh_def.read_from_config(cfg.getRoot()["channel_mesh"]); //  some config file in KT is "background_mesh"
    if (cfg.exists("region_refine"))
    {
      const auto &cfg_refine = cfg.getRoot()["region_refine"];
      region_refine.resize(cfg_refine.getLength());
      for (unsigned int i = 0; i < region_refine.size(); i++)
      {
        region_refine[i].read_from_config(cfg_refine[i]);
      }
    }

    /// Linear Elasticity
    ReadVectorOrValue("bodyforce", bodyforce);
    for (int i = 0; i < DIM; i++)
    {
      BodyForce(i) = bodyforce[i];
    }

    radialbodyforce.read_from_config(cfg.getRoot()["radialbodyforce"]);
    ReadValue("rho", rho);
    ReadValue("scaleFactor", scaleFactor);

      if (cfg.exists("geometries_ibm"))
      {
          const auto &geometries = cfg.getRoot()["geometries_ibm"];
          ibm_geom_def.resize(geometries.getLength());
          for (unsigned int i = 0; i < ibm_geom_def.size(); i++)
          {
              ibm_geom_def[i].read_from_config(geometries[i]);
          }
      }

    /// Read Geo for SBM
    if (ibm_geom_def.size() != 0) {
        SbmGeo = read_SbmGeo(cfg.getRoot(), "SBMGeo");

    }


    /// timestep control
    ReadVectorOrValue("dt", dt);
    ReadVectorOrValue("totalT", totalT);

    /// Output control
    if (ReadValue("OutputStartTime", OutputStartTime))
    {
    }
    if (ReadValue("OutputInterval", OutputInterval))
    {
    }
    CheckpointInterval = OutputInterval;
    if (ReadValue("CheckpointInterval", CheckpointInterval))
    {
    }
    if (ReadValue("CheckpointNumbackup", CheckpointNumbackup))
    {
    }

    /// WeakBC parameters
    Cb_e.read_from_config(cfg, "Cb_e");
    TauN.read_from_config(cfg, "TauN");
    PrintStatus("Cb_e = ",Cb_e.value_at(100000));
    PrintStatus("TauN = ",TauN.value_at(100000));

    /// SBM
    if (ReadValue("RatioGPSBM",RatioGPSBM) || ReadValue("lambda",RatioGPSBM))
    {
    }

    if (ReadValue("RelOrderCheckActive",RelOrderCheckActive))
    {
    }
    if (ReadValue("RelOrderInterceptedElement",RelOrderInterceptedElement))
    {
    }

    /// Solver Options
    solverOptionsLE = read_solver_options(cfg, "solver_options_le");

    return true;
  }

  /// Function for reading a vector or a single value (stored in vector)
  template <typename T>
  void ReadVectorOrValue(const std::string &key_name, std::vector<T> &value)
  {
    if (cfg.exists(key_name + "_V"))
    {
      InputData::ReadVector(cfg, key_name + "_V", value);
    }
    else
    {
      double value_const;
      ReadValueRequired(key_name, value_const);
      value.push_back(value_const);
    }
  }

  /**
   * Printout every item of inputdata for debug purpose.
   */
  void PrintInputData()
  {
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank)
    {
      std::ofstream fout("InputDataOutput.txt", std::ios::app);
      time_t my_time = time(NULL);
      fout << "##############################"
           << "\n";
      fout << ctime(&my_time);
      fout << "Total number of processor = " << TALYFEMLIB::GetMPISize() << "\n";
      fout << "size of DendroInt " << sizeof(DendroIntL) << "\n";
      fout << "size of PetscInt " << sizeof(PetscInt) << "\n";

      fout << "Dimension: " << nsd << "\n";
      fout << "basisFunctionOrder: " << elemOrder << "\n";
      fout << "mfree: " << ifMatrixFree << "\n";
      fout << "hessian: " << ifHessian << "\n";

      fout << "====== mesh_def ======"
           << "\n";
      mesh_def.PrintMeshDef(fout);
      fout << "====================="
           << "\n\n";

      fout << "====== timestepper ======"
           << "\n";
      PrintVector(fout, "dt", dt);
      PrintVector(fout, "totalT", totalT);
      fout << "====================="
           << "\n\n";

      fout << "===== parameters ======="
           << "\n";
      Cb_e.PrintParameterDef(fout, "Cb_e");
      fout << "\n";
      TauN.PrintParameterDef(fout, "TauN");

      fout << "====================="
           << "\n\n";

      fout << "\nregion_refine: {\n";
      for (const auto &r : region_refine)
      {
        r.PrintRegionRefineDef(fout);
        fout << "}\n{\n";
      }
      fout << "========== solver settings ==========="
           << "\n";
      fout << "]"
           << "\n";
      fout << "solverOptionsLE: ["
           << "\n";
      for (auto &val : solverOptionsLE.vals)
      {
        fout << val.first << " --> " << val.second << "\n";
      }
      fout << "]"
           << "\n";

      fout.close();
    }
  }

private:
  std::string str;
  std::vector<double> bodyforce;

  /// function reading geo of SBM
  static SBMGeo read_SbmGeo(libconfig::Setting &root, const char *name)
  {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str))
    {
        if(str == "PLANE")
        {
            return PLANE;
        }
        else if (str == "PLANT")
        {
            return PLANT;
        }
        else if (str == "EIFFEL")
        {
            return EIFFEL;
        }
        else if (str == "GYROID")
        {
            return GYROID;
        }
      else if (str == "BUNNY")
      {
        return BUNNY;
      }
      else if (str == "NONE")
      {
          return NONE;
      }
      else
      {
        throw TALYFEMLIB::TALYException() << "Unknown SBM geo-- " << name << str;
      }
    }
    else
    {
//      std::cout
//          << "Must specify SBMGeo \n";
        PrintStatus("-------------------------------------");
        PrintStatus("User do not set any SBMGeo, But this code is designed to run with SBMGeo for SPM25");
        PrintStatus("We set it as NONE");
        PrintStatus("-------------------------------------");
        return NONE;

    }
  }

  /// Function for reading type of distance function calculation in 3D
  static typeDistCalc read_DistCalc(libconfig::Setting &root, const char *name)
  {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str))
    {
      if (str == "NORMAL_BASED")
      {
        return NORMAL_BASED;
      }
      else if (str == "GP_BASED")
      {
        return GP_BASED;
      }
      else
      {
        throw TALYFEMLIB::TALYException() << "Unknown solver name for DistCalc: " << name << str;
      }
    }
    else
    {
        PrintStatus("User do not set any DistCalcType, we set it as KD_TREE");
        return KD_TREE;
    }
  }

  static CaseType read_LEcase(libconfig::Setting &root, const char *name)
  {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str))
    {
      if (str == "planestress")
      {
#if (DIM == 3)
        PrintError("3D do not have this type, please use planestrain");
        exit(EXIT_FAILURE);
#endif
#if (DIM == 2)
        PrintStatus("[LE case] PLANESTRESS");
#endif
        return PLANESTRESS;
      }
      else if (str == "planestrain")
      {
        PrintStatus("[LE case] PLANESTRAIN");
        return PLANESTRAIN;
      }
      else if (str == "lame")
      {
        PrintStatus("[LE case] LAME");
        return LAME;
      }
      else
      {
        throw TALYFEMLIB::TALYException() << "Unknown case name for LE: " << name << str;
      }
    }
    else
    {
      throw TALYFEMLIB::TALYException() << "Must specify case: planestress, planestrain, or lame";
    }
  }

  static BCCaseType read_LEBCcase(libconfig::Setting &root, const char *name)
  {
    std::string str;
    /// If nothing specified stays stabilizedNS
    if (root.lookupValue(name, str))
    {
      if (str == "NormalTraction")
      {
        PrintStatus("[LE BC case] NORMAL_TRACTION");
        return NORMAL_TRACTION;
      }
      else if (str == "DisplacementBothSide")
      {
        PrintStatus("[LE BC case] DISPLACEMENT_BOTH_SIDE");
        return DISPLACEMENT_BOTH_SIDE;
      }
      else if (str == "FixedAtWall")
      {
        PrintStatus("[LE BC case] FIXED_AT_WALL");
        return FIXED_AT_WALL;
      }
      else if (str == "ZeroTraction")
      {
        PrintStatus("[LE BC case] ZERO_TRACTION");
        return ZERO_TRACTION;
      }
      else if (str == "HalfBeam")
      {
        PrintStatus("[LE BC case] HALF_BEAM");
        return HALF_BEAM;
      }
      else if (str == "Tract4Side")
      {
        return TRACT4SIDE;
      }
      else if (str == "BOTTOM_FORCE")
      {
          return  BOTTOM_FORCE;
      }
      else if (str == "CSV_FORCE"){
          PrintStatus("[LE BC case] CSV_FORCE");
          return CSV_FORCE;
      }
      else if (str == "POSITION_DISPLACEMENT")
      {
          PrintStatus("[LE BC case] POSITION_DISPLACEMENT; Everything is Setup at Code");
          return POSITION_DISPLACEMENT;
      }
      else
      {
        throw TALYFEMLIB::TALYException() << "Unknown BC case name for LE: " << name << str;
      }
    }
    else
    {
      throw TALYFEMLIB::TALYException() << "Must specify BC case: NormalTraction, DisplacementBothSide, FixedAtWall, HALF_BEAM, or ZERO_TRACTION";
    }
  }
};
