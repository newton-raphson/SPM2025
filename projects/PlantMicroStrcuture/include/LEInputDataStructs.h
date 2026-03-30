#pragma once
#include "util.h"

/// read size of 3 vector as ZEROPTV
void ReadZEROPTV(const libconfig::Setting &root, const std::string &key_name, ZEROPTV &value, bool required = true) {
  if (root.exists(key_name.c_str())) {
    if (root[key_name.c_str()].getLength() >= DIM) {
      value(0) = (double) root[key_name.c_str()][0];
      value(1) = (double) root[key_name.c_str()][1];
#if (DIM == 2)
      value(2) = 0.0;
#endif
#if (DIM == 3)
      value(2) = (double) root[key_name.c_str()][2];
#endif
    } else {
      throw TALYException() << key_name + " have size of " + std::to_string(root[key_name.c_str()].getLength());
    }
  } else if (required) {
    throw TALYException() << key_name + " doesn't exist!";
  }
}

/// Function for reading a vector from root
template<typename T>
void ReadVectorRoot(const libconfig::Setting &root, const std::string &key_name, std::vector<T> &value) {
  const libconfig::Setting &config_v = root[key_name.c_str()];
  if (config_v.isNumber()) {
    value.push_back(config_v);
  } else {
    value.reserve(config_v.getLength());
    for (int i = 0; i < config_v.getLength(); ++i) {
      value.emplace_back(config_v[i]);
    }
  }
}

/*
 * maximum value of an array
 */
double maximum(std::array<DENDRITE_REAL, DIM> &array) {
  DENDRITE_REAL max = array[0];
  for (int i = 1; i < DIM; i++) {
    if (array[i] > max) {
      max = array[i];
    }
  }
  return max;
}
double minimum(std::array<DENDRITE_REAL, DIM> &array)
{
  DENDRITE_REAL min = array[0];
  for (int i = 1; i < DIM; i++)
  {
    if (array[i] < min)
        min = array[i];
  }
  return min;
}

template<typename T>
void PrintVector(std::ofstream &fstream, const std::string &name, const std::vector<T> &vec) {
  int rank = TALYFEMLIB::GetMPIRank();
  if (!rank and fstream.is_open()) {
    fstream << name << ": [";
    for (const auto &v : vec) {
      fstream << v << ", ";
    }
    fstream << "]\n";
  }
}

/**
 * Time-dependent parameter
 */
struct ParameterDef {
  enum RampingMethod {
    CONSTANT = 0,
    LINEAR = 1,
    SPLINE = 2,
  };
  std::vector<double> value;
  std::vector<double> ramping_vec;
  RampingMethod ramping_method = RampingMethod::CONSTANT;

  void read_from_config(libconfig::Config &cfg, const std::string &key_name, double default_value = 0.0) {
    if (cfg.exists(key_name + "_V")) {
      InputData::ReadVector(cfg, key_name + "_V", value);
      InputData::ReadVector(cfg, key_name + "_ramping", ramping_vec);
      if (cfg.exists(key_name + "_method")) {
        ramping_method = read_ramping_method_from_config(cfg.getRoot(), key_name + "_method");
      }
    } else if (cfg.exists(key_name)) {
      double value_const;
      bool result = cfg.lookupValue(key_name, value_const);
      if (!result) {
        throw TALYException() << "Can not find" + key_name;
      }
      value.push_back(value_const);
      ramping_vec.push_back(0.0);
    } else {
      value.push_back(default_value);
      ramping_vec.push_back(0.0);
    }
  }

  static RampingMethod read_ramping_method_from_config(const libconfig::Setting &root, const std::string &key_name) {
    return str_to_ramping_method(root[key_name.c_str()]);
  }

  double value_at(double time) {
    if (ramping_method == RampingMethod::CONSTANT) {
      unsigned int interval = ramping_vec.size() - 1;
      if (fabs(ramping_vec[0]) > 1e-15) {
        throw TALYException() << "First element of ramping have to be 0";
      }
      if (value.size() != ramping_vec.size()) {
        throw TALYException() << "Ramping vector size must be same with value";
      }
      for (unsigned int i = 0; i < ramping_vec.size(); i++) {
        if (time < ramping_vec[i]) {
          interval = i - 1;
          break;
        }
      }
      double current_variable = value[interval];
      return current_variable;
    } else if (ramping_method == RampingMethod::LINEAR) {
      unsigned int interval = ramping_vec.size();
      for (unsigned int i = 0; i < ramping_vec.size(); i++) {
        if (time < ramping_vec[i]) {
          interval = i;
          break;
        }
      }
      double current_variable;
      if (interval == 0) {
        current_variable = value[0] * (time / ramping_vec[0]);
      } else if (interval == ramping_vec.size()) {
        current_variable = value[interval - 1];
      } else {
        double per1 = (time - ramping_vec[interval - 1]) / (ramping_vec[interval] - ramping_vec[interval - 1]);
        double per2 = (ramping_vec[interval] - time) / (ramping_vec[interval] - ramping_vec[interval - 1]);
        current_variable = value[interval] * per1 + value[interval - 1] * per2;
      }
      return current_variable;
    } else if (ramping_method == RampingMethod::SPLINE) {
      throw TALYException() << "not implentmentd yet";
    }
  }

  void PrintParameterDef(std::ofstream &fstream, const std::string &name) {
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()) {
      fstream << name << ": {\n";
      fstream << "value: [";
      for (const auto &v : value) {
        fstream << v << ", ";
      }
      fstream << "]\n";
      fstream << "ramping_vec: [";
      for (const auto &v : ramping_vec) {
        fstream << v << ", ";
      }
      fstream << "]\n";
      if (ramping_method == RampingMethod::CONSTANT) {
        fstream << "ramping_method: CONSTANT\n";
      } else if (ramping_method == RampingMethod::LINEAR) {
        fstream << "ramping_method: LINEAR\n";
      } else if (ramping_method == RampingMethod::SPLINE) {
        fstream << "ramping_method: SPLINE\n";
      }
      fstream << "}\n";
    }
  }

 private:
  static RampingMethod str_to_ramping_method(const std::string &str) {
    if (str == "constant") {
      return RampingMethod::CONSTANT;
    } else if (str == "linear") {
      return RampingMethod::LINEAR;
    } else if (str == "spline") {
      return RampingMethod::SPLINE;
    } else {
      throw TALYFEMLIB::TALYException() << "Invalid Ramping method (constant, linear, spline)";
    }
  }

};


/**
 * Domain definition
 */
struct MeshDef {
  DomainInfo fullDADomain; /// The domain from which its carved out.
  DomainInfo physDomain;   /// The actual SubDAdomain

  std::array<DENDRITE_REAL, DIM> channel_max; ///< end of domain
  std::array<DENDRITE_REAL, DIM> channel_min;

  ///< Dendro options
  DENDRITE_UINT refine_lvl_base;
  DENDRITE_UINT refine_lvl_channel_wall;
  /// refine_lvl_interface

  bool refine_walls = false; ///< whether or not fx_refine refines at channel_min/channel_max

  void read_from_config(const libconfig::Setting &root) {

    refine_lvl_base = (unsigned int) root["refine_lvl_base"];

    channel_max[0] = (DENDRITE_REAL) root["max"][0];
    channel_max[1] = (DENDRITE_REAL) root["max"][1];
#if (DIM == 3)
    channel_max[2] = (DENDRITE_REAL) root["max"][2];
#endif

    /// channel_min
    channel_min[0] = (DENDRITE_REAL) root["min"][0];
    channel_min[1] = (DENDRITE_REAL) root["min"][1];
#if (DIM == 3)
    channel_min[2] = (DENDRITE_REAL) root["min"][2];
#endif


    if (!root.lookupValue("refine_walls", refine_walls)) {
      refine_walls = false;
    } else {
      refine_walls = root["refine_walls"];
    }

    if (refine_walls) {
      refine_lvl_channel_wall = (unsigned int) root["refine_lvl_channel_wall"];
      // refine_lvl_pillar_wall = (unsigned int)root["refine_lvl_pillar_wall"];
      if (refine_lvl_base < 1 || refine_lvl_base >= 31 ||
          refine_lvl_channel_wall < 1 || refine_lvl_channel_wall >= 31) {
        PrintWarning("Invalid refine_ level - should be 1 < refine_lvl <= 31.");
      }
    }
    if (!refine_walls) {
      refine_lvl_channel_wall = refine_lvl_base;
    }

    if (refine_lvl_base > refine_lvl_channel_wall) {
      PrintWarning("refine_lvl_base > refine_lvl_channel_wall");
    }

    PrintStatus("refine_lvl_base: ", refine_lvl_base);
    PrintStatus("refine_lvl_channel_wall: ", refine_lvl_channel_wall);

    double sizeOfBox = maximum(channel_max);
    double minimum_cube = minimum(channel_min);

    fullDADomain.min.fill(minimum_cube);
    fullDADomain.max.fill(sizeOfBox);


    physDomain.min[0] = channel_min[0];
    physDomain.min[1] = channel_min[1];

    physDomain.max[0] = channel_max[0];
    physDomain.max[1] = channel_max[1];
#if (DIM == 3)
    physDomain.min[2] = channel_min[2];
    physDomain.max[2] = channel_max[2];
#endif

    PrintStatus("fullDADomain: ", fullDADomain.max[0], ", ", fullDADomain.max[1], ", ",
                fullDADomain.max[2]);
    PrintStatus("physicalDomain: ", physDomain.max[0], ", ", physDomain.max[1], ", ",
                physDomain.max[2]);
  }

  void PrintMeshDef(std::ofstream &fstream) {
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()) {
      fstream << "channel_max: [";
      for (int i = 0; i < DIM; i++) {
        fstream << channel_max[i] << ", ";
      }
      fstream << "]\n";

      fstream << "fullDADomain-min: [";
      for (int i = 0; i < DIM; i++) {
        fstream << fullDADomain.min[i] << ", ";
      }
      fstream << "]\n";
      fstream << "fullDADomain-max: [";
      for (int i = 0; i < DIM; i++) {
        fstream << fullDADomain.max[i] << ", ";
      }
      fstream << "]\n";

      fstream << "physDomain-min: [";
      for (int i = 0; i < DIM; i++) {
        fstream << physDomain.min[i] << ", ";
      }
      fstream << "]\n";
      fstream << "physDomain-max: [";
      for (int i = 0; i < DIM; i++) {
        fstream << physDomain.max[i] << ", ";
      }
      fstream << "]\n";

      fstream << "refine_lvl_base: " << refine_lvl_base << "\n";
      fstream << "refine_lvl_channel_wall: " << refine_lvl_channel_wall << "\n";
      fstream << "refine_walls: " << refine_walls << "\n";
    }
  }
};

/**
 * Parameters for boundary condition
 */
struct BoundaryDef {
  enum Side {
    INVALID = -1,
    X_MINUS = 0,
    X_PLUS = 1,
    Y_MINUS = 2,
    Y_PLUS = 3,
    Z_MINUS = 4,
    Z_PLUS = 5,
  };

  enum Disp_Type {
    DIRICHLET_DX = 0,
    NEUMANN_DX = 1,
    WEAK_DX = 2,
    DIRICHLET_DY = 3,
    NEUMANN_DY = 4,
    WEAK_DY = 5,
    DIRICHLET = 6,
    NEUMANN = 7,
    WEAK = 8,
    DIRICHLET_DZ = 9,
    NEUMANN_DZ = 10,
    WEAK_DZ = 10,
    ROTATION_D=11,
    SBM_NEUMANN_WALL = 12,
    SBM_ZERO_DIRICHLET = 13
  };
  
  Side side = INVALID;
  //Linear
  Disp_Type disp_type;
  
  /// for rotation BC (Linear Elasticity)
  DENDRITE_REAL theta;

  /// for Dirichlet BC (Linear Elasticity)
  DENDRITE_REAL UX;
  DENDRITE_REAL UY;
  DENDRITE_REAL UZ;
  
  /// for Neumann BC (Linear Elasticity)
  ZEROPTV Traction;
  
  bool ifRefine = true;
  
 private:


};

/**
 * Parameters for Regional refine
 */
struct RegionalRefine {
  enum RetainType {
    OUTSIDE = 0,
    INSIDE = 1,
    ON_BOUNDARY = 2,
  };
  enum Type {
    INVALID = -1,
    CUBE = 0,
    SPHERE = 1,
    CYLINDER = 2,
    MESHOBJECT = 10,
    MESHOBJECT_2D = 11,
    DeepTrace = 12,
  };
  bool forRetain = false;
  // max refinement level
  int refine_region_lvl = 0;
  int refine_region_lvl_boundary = 0;
  std::string mesh_path;
  ZEROPTV shift;
  GEOMETRY::Geometry *geometry;
  GEOMETRY::STL *stl;
  GEOMETRY::MSH *msh;

  void read_from_config(const libconfig::Setting &root) {
    if (root.exists("forRetain")) {
      forRetain = (bool) root["forRetain"];
    }
    refine_type = read_type(root["type"]);
    if (forRetain) {
      if (refine_type != CUBE and refine_type != MESHOBJECT) {
        throw TALYException() << "Not supported retain type!";
      }
    } else {
      refine_region_lvl = (int) root["refine_region_lvl"];
      refine_region_lvl_boundary = refine_region_lvl;
      if (root.exists("refine_region_lvl_boundary")) {
        refine_region_lvl_boundary = (int) root["refine_region_lvl_boundary"];
      }
    }
    if (refine_type == CUBE) {
      ReadZEROPTV(root, "min_c", min_corner);
      ReadZEROPTV(root, "max_c", max_corner);
    }
#if (DIM == 3)
    if (refine_type == MESHOBJECT) {
      mesh_path = (const char *) root["mesh_path"];
      ReadZEROPTV(root, "shift", shift);

      stl = new  GEOMETRY::STL(mesh_path, GEOMETRY::InOutTest::RAY_TRACING);
      std::array<DENDRITE_REAL, DIM> point;
      point[0] = shift[0];
      point[1] = shift[1];
      point[2] = shift[2];
      geometry = new GEOMETRY::Geometry(stl, Point<DIM>(point));
    }
#endif

#if (DIM == 2)
    if (refine_type == MESHOBJECT_2D) {
      mesh_path = (const char *) root["mesh_path"];
      ReadZEROPTV(root, "shift", shift);

      msh = new GEOMETRY::MSH(mesh_path, GEOMETRY::InOutTest2D::RAY_TRACING_2D);
      std::array<DENDRITE_REAL, DIM> point;
      point[0] = shift[0];
      point[1] = shift[1];
      geometry = new GEOMETRY::Geometry(msh, Point<DIM>(point));
    }
#endif
    if (refine_type == SPHERE) {
      ReadZEROPTV(root, "center", center_sphere);
      radius_sphere = (double) root["radius"];
      if (root.exists("radius_in")) {
        radius_sphere_in = (double) root["radius_in"];
      }
    }
    if (refine_type == CYLINDER) {
      ReadZEROPTV(root, "c1", min_corner);
      ReadZEROPTV(root, "c2", max_corner);
      radius_cylinder = (double) root["radius"];
      if (root.exists("radius_in")) {
        radius_cylinder_in = (double) root["radius_in"];
      }
    }

  }

 public:
  /**
   * @return true for outside retain cube (boundary nodes included (for apply BC)), false for inside retain cube
   */
  RetainType out_retain(ZEROPTV p) {
    if (not(forRetain)) {
      throw TALYException() << "Calling function with forRetain = false";
    }
    switch (refine_type) {
      case CUBE: {
        const double eps = 1e-8;
        if (p.x() < min_corner.x() - eps ||
            p.y() < min_corner.y() - eps ||
            p.z() < min_corner.z() - eps ||
            p.x() > max_corner.x() + eps ||
            p.y() > max_corner.y() + eps ||
            p.z() > max_corner.z() + eps) {
          return OUTSIDE;
        } else if (p.x() > min_corner.x() + eps &&
            p.y() > min_corner.y() + eps &&
            p.z() > min_corner.z() + eps &&
            p.x() < max_corner.x() - eps &&
            p.y() < max_corner.y() - eps &&
            p.z() < max_corner.z() - eps) {
          return INSIDE;
        } else {
          return ON_BOUNDARY;
        }
      }
        break;
      case MESHOBJECT: {
        bool ifInside = geometry->ifInside(p.data());
        if (ifInside) {
          return INSIDE;
        } else {
          return OUTSIDE;
        }
      }
        break;
      default:throw TALYFEMLIB::TALYException() << "Wrong type with in/out test for retain!";
    }

  }

  /**
   * @return false for outside region, true for inside region
   */
  bool in_region(ZEROPTV p) {
    if (forRetain) {
      throw TALYException() << "Calling function with forRetain = true";
    }
    switch (refine_type) {
      case CUBE: {
        return !(p.x() < min_corner.x() ||
            p.y() < min_corner.y() ||
            p.z() < min_corner.z() ||
            p.x() > max_corner.x() ||
            p.y() > max_corner.y() ||
            p.z() > max_corner.z());
      }
        break;
      case SPHERE: {
        double distance = (p - center_sphere).norm();
        return (distance < radius_sphere and distance > radius_sphere_in);
      }
        break;
      case CYLINDER: {
        ZEROPTV AB = min_corner - max_corner;
        ZEROPTV temp = AB * (1.0 / AB.innerProduct(AB));
        ZEROPTV AP = p - max_corner;
        ZEROPTV proj_point = max_corner + temp * AP.innerProduct(AB);
        double distance_s = (proj_point - p).innerProduct(proj_point - p);
        if (distance_s < radius_cylinder * radius_cylinder
            and (distance_s - radius_cylinder_in * radius_cylinder_in) > -1e-6) {
          double max_x = std::max(max_corner.x(), min_corner.x());
          double max_y = std::max(max_corner.y(), min_corner.y());
          double max_z = std::max(max_corner.z(), min_corner.z());
          double min_x = std::min(max_corner.x(), min_corner.x());
          double min_y = std::min(max_corner.y(), min_corner.y());
          double min_z = std::min(max_corner.z(), min_corner.z());
          return !(proj_point.x() > max_x ||
              proj_point.y() > max_y ||
              proj_point.z() > max_z ||
              proj_point.x() < min_x ||
              proj_point.y() < min_y ||
              proj_point.z() < min_z);
        } else {
          return false;
        }
      }
        break;
      case MESHOBJECT: {
        bool ifInside = geometry->ifInside(p.data());
        return ifInside;
      }
      case MESHOBJECT_2D: {
        bool ifInside = geometry->ifInside(p.data());
        return ifInside;
      }
        break;
      default:throw TALYFEMLIB::TALYException() << "Wrong type!";
    }
  }

  Type GetRefineType() const {
    return refine_type;
  }

  void PrintRegionRefineDef(std::ofstream &fstream) const {
    int rank = TALYFEMLIB::GetMPIRank();
    if (!rank and fstream.is_open()) {

      if (refine_type == Type::INVALID) {
        fstream << "refine_type: INVALID\n";
      } else if (refine_type == Type::CUBE) {
        fstream << "refine_type: CUBE\n";
        fstream << "min_corner: [";
        for (int i = 0; i < 3; i++) {
          fstream << min_corner[i] << ", ";
        }
        fstream << "]\n";
        fstream << "max_corner: [";
        for (int i = 0; i < 3; i++) {
          fstream << max_corner[i] << ", ";
        }
        fstream << "]\n";
      } else if (refine_type == Type::SPHERE) {
        fstream << "refine_type: SPHERE\n";
        fstream << "center_sphere: [";
        for (int i = 0; i < 3; i++) {
          fstream << center_sphere[i] << ", ";
        }
        fstream << "]\n";
        fstream << "radius_sphere: " << radius_sphere << "\n";
        fstream << "radius_sphere_in: " << radius_sphere_in << "\n";
      } else if (refine_type == Type::CYLINDER) {
        fstream << "refine_type: CYLINDER\n";
        fstream << "min_corner: [";
        for (int i = 0; i < 3; i++) {
          fstream << min_corner[i] << ", ";
        }
        fstream << "]\n";
        fstream << "max_corner: [";
        for (int i = 0; i < 3; i++) {
          fstream << max_corner[i] << ", ";
        }
        fstream << "]\n";

        fstream << "radius_cylinder: " << radius_cylinder << "\n";
        fstream << "radius_cylinder_in: " << radius_cylinder_in << "\n";
      } else if (refine_type == Type::MESHOBJECT) {
        fstream << "refine_type: MESHOBJECT\n";
      } else if (refine_type == Type::MESHOBJECT_2D) {
        fstream << "refine_type: MESHOBJECT_2D\n";
      }

      fstream << "shift: [";
      for (int i = 0; i < 3; i++) {
        fstream << shift[i] << ", ";
      }
      fstream << "]\n";

      fstream << "forRetain: " << forRetain << "\n";
      fstream << "refine_region_lvl: " << refine_region_lvl << "\n";
      fstream << "refine_region_lvl_boundary: " << refine_region_lvl_boundary << "\n";
      fstream << "mesh_path: " << mesh_path << "\n";
    }
  }

 protected:
  Type refine_type = INVALID;
  // for sphere
  ZEROPTV center_sphere;
  double radius_sphere = 0.0;
  double radius_sphere_in = 0.0;
  // for cube
  ZEROPTV min_corner;
  ZEROPTV max_corner;
  // for cylinder
//  ZEROPTV min_corner;
//  ZEROPTV max_corner;
  double radius_cylinder = 0.0;
  double radius_cylinder_in = 0.0;
  static Type read_type(const std::string &str) {
    if (str == "sphere") {
      return SPHERE;
    } else if (str == "cube") {
      return CUBE;
    } else if (str == "cylinder") {
      return CYLINDER;
    } else if (str == "meshobject") {
      return MESHOBJECT;
    } else if (str == "meshobject_2d") {
      return MESHOBJECT_2D;
    }
    else if( str == "deeptrace")
    {
        return DeepTrace;
    }
    else
      throw TALYFEMLIB::TALYException() << "Invalid Regional refinement type";
  }
};


