#pragma once
#include "LEInputDataStructs.h"
#include <time.h>
#include <talyfem/input_data/input_data.h>

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
struct PlantProperty
{
  double young_vascular,young_pith,young_rind;
  double poisson_vascular,poisson_pith,poisson_rind;
  void read_from_config(const libconfig::Setting &root)
  {
    young_vascular = root["young_vascular"];
    young_pith = root["young_pith"];
    young_rind = root["young_rind"];
    poisson_vascular = root["poisson_vascular"];
    poisson_pith = root["poisson_pith"];
    poisson_rind = root["poisson_rind"];
  }
};
struct PlantGeometry
{
  ///// Geometry parameters (not initialized)
  double ellipse_cx;
  double ellipse_cy;
  double ellipse_rx;
  double ellipse_ry;
  double major_circle_radius;
  double height;

  //// vascular bundle pattern

  double circle_radius;
  double offset_threshold;
  double radius_max;
  int num_rings;
  int num_sectors;
  double radius_tolerance;
  double angle_jitter_min;
  double angle_jitter_max;
  double radius_jitter_min;
  double radius_jitter_max;
  unsigned rng_seed;
  /// vascular bundle refinement
  unsigned bundle_refinement;
  //// this is the center for circles
  using Point = std::pair<double, double>;
  std::vector<Point> circle_centers;


  // Load from config file
  void read_from_config(const libconfig::Setting& root) {
    ellipse_cx = root["ellipse_cx"];
    ellipse_cy = root["ellipse_cy"];
    ellipse_rx = root["ellipse_rx"];
    ellipse_ry = root["ellipse_ry"];
    major_circle_radius = root["major_circle_radius"];
    height = root["height"];

    circle_radius = root["circle_radius"];
    offset_threshold = root["offset_threshold"];
    radius_max = root["radius_max"];
    num_rings = root["num_rings"];
    num_sectors = root["num_sectors"];
    radius_tolerance = root["radius_tolerance"];
    angle_jitter_min = root["angle_jitter_min"];
    angle_jitter_max = root["angle_jitter_max"];
    radius_jitter_min = root["radius_jitter_min"];
    radius_jitter_max = root["radius_jitter_max"];
    rng_seed = root["rng_seed"];

    generateRadialPattern();  // Automatically compute and cache the circle centers
  }

private:
  void generateRadialPattern() {
    std::default_random_engine rng(rng_seed);
    std::uniform_real_distribution<double> angle_jitter(angle_jitter_min, angle_jitter_max);
    std::uniform_real_distribution<double> radius_jitter(radius_jitter_min, radius_jitter_max);

    double angular_step = 2.0 * M_PI / num_sectors;
    double base_radius = radius_max / num_rings;

    circle_centers.clear();
    for (int r = 1; r <= num_rings; ++r) {
      double radial_distance = r * base_radius;
      for (int a = 0; a < num_sectors; ++a) {
        double angle = a * angular_step + angle_jitter(rng);
        double radius = radial_distance + radius_jitter(rng);
        double x = radius * std::cos(angle);
        double y = radius * std::sin(angle);

        bool overlaps = false;
        for (const auto& c : circle_centers) {
          double dx = x - c.first;
          double dy = y - c.second;
          if (std::sqrt(dx * dx + dy * dy) < 2.1 * radius_tolerance) {
            overlaps = true;
            break;
          }
        }

        if (!overlaps) {
          circle_centers.emplace_back(x, y);
        }
      }
    }
  }
  void read_RadialPattern()
  {
    /// read a csv file names bundle_center.csv and
    throw std::runtime_error("Not implemented");

  }

public:
  // Signed distance to ellipse ⊖ circle
  double signedDistance(double x, double y) const {
    double dx = (x - ellipse_cx) / ellipse_rx;
    double dy = (y - ellipse_cy) / ellipse_ry;
    double ellipse_phi = (std::sqrt(dx * dx + dy * dy) - 1.0) * std::min(ellipse_rx, ellipse_ry);
    double circle_phi = std::sqrt(x * x + y * y) - circle_radius;
    return std::max(circle_phi, -ellipse_phi);
  }

  bool isInOffsetRegion(const std::vector<ZEROPTV>& coords) const {
    for (const auto& p : coords) {
      if (signedDistance(p.x(), p.y()) > offset_threshold) {
        return true;
      }
    }
    return false;
  }

  bool isInsideRadialCircle(const std::vector<ZEROPTV>& coords) const {
    for (const auto& pt : coords) {
      for (const auto& center : circle_centers) {
        double dx = pt.x() - center.first;
        double dy = pt.x() - center.second;
        if (std::sqrt(dx * dx + dy * dy) <= radius_tolerance) {
          return true;
        }
      }
    }
    return false;
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

  MAX_CASE_TYPE = 3,
  PLANTPROPERTY = 4,
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


struct PlantFiberProp
{
  double MatrixE;
  double Matrixmu;
  double FiberE;
  double Fibermu;
  bool hardThreshold;
  double hardThresholdValue;
  /// read string for the image path
  std::string Image_Path;
  void read_from_config(const libconfig::Setting &root)
  {
    MatrixE = root["MatrixE"];
    Matrixmu = root["Matrixmu"];
    FiberE = root["FiberE"];
    Fibermu = root["Fibermu"];
    hardThreshold = root["hardormixture"];
    Image_Path = static_cast<const char *>(root["Image_Path"]);

    if (hardThreshold)
    {
      hardThresholdValue = root["hardThresholdValue"];
    }

  }

};

enum BCTYPE {e100,e010,e001};



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

  PlantFiberProp planeFiberProp;



  std::vector<std::pair<double, double>> minmax_cantilever;

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

    }

    planeFiberProp.read_from_config(cfg.getRoot()["planeFiberProp"]);


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
      else if (str == "PLANTPROPERTY")
      {
        PrintStatus("[LE case] PLANTPROPERTY");
        return PLANTPROPERTY;
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
