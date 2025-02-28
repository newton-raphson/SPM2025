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
    enum Vel_Type {
        NO_SLIP = 0,
        NO_GRADIENT_VEL = 1,
        DIRICHLET_VEL = 2,
        WEAK_VEL = 3,
        PARABOLIC = 4,
        TOP_HAT = 5,
        DISTURBANCE = 6
    };
    enum Pressure_Type {
        NO_GRADIENT_P = 0,
        DIRICHLET_P = 1,
        BACKFLOW_STAB = 2
    };

    Side side = INVALID;
    Vel_Type vel_type = NO_GRADIENT_VEL;
    Pressure_Type pressure_type = NO_GRADIENT_P;
    std::vector<TALYFEMLIB::ZEROPTV> wall_velocity;
    DENDRITE_REAL pressure = 0.0;
    DENDRITE_REAL BetaNS = 0.5;
    std::vector<DENDRITE_REAL> ramping;
    bool ifRefine = false;


    static Side read_side_from_config(const libconfig::Setting &root) {
        return str_to_side(root["side"]);
    }

    void read_from_config(const libconfig::Setting &root) {
        if (root.exists("ifRefine")) {
            ifRefine = (bool) root["ifRefine"];
        }
        if (root.exists("vel_type")) {
            vel_type = read_vel_type(root["vel_type"]);
        }
        if (root.exists("pressure_type")) {
            pressure_type = read_pressure_type(root["pressure_type"]);
        }
        if (root.exists("ramping_vec")) {
            const libconfig::Setting &ramping_time_vec = root["ramping_vec"];
            for (int i = 0; i < ramping_time_vec.getLength(); ++i) {
                ramping.push_back(ramping_time_vec[i]);
            }
        } else {
            ramping.push_back(0.0);
        }

        if (vel_type == Vel_Type::DIRICHLET_VEL || vel_type == Vel_Type::DISTURBANCE || vel_type == Vel_Type::WEAK_VEL || vel_type == Vel_Type::PARABOLIC) {
            const libconfig::Setting &wall_vel_vec = root["wall_vel"];
            for (int i = 0; i < wall_vel_vec.getLength() / 3; ++i) {
                TALYFEMLIB::ZEROPTV temp_vel = TALYFEMLIB::ZEROPTV((DENDRITE_REAL) root["wall_vel"][3 * i],
                                           (DENDRITE_REAL) root["wall_vel"][3 * i + 1],
                                           (DENDRITE_REAL) root["wall_vel"][3 * i + 2]);
                wall_velocity.push_back(temp_vel);
            }
            if (wall_velocity.size() != ramping.size()) {
                throw TALYFEMLIB::TALYException() << "Wall_vel vector not matching with ramping";
            }
        }

        if (pressure_type == Pressure_Type::DIRICHLET_P) {
            pressure = (DENDRITE_REAL) root["pressure"];
        }

        if (pressure_type == Pressure_Type::BACKFLOW_STAB) {
            if (root.exists("BetaNS")) {
                BetaNS = (DENDRITE_REAL) root["BetaNS"];
                TALYFEMLIB::PrintStatus("[Set] BetaNS =", BetaNS);
            } else {
                TALYFEMLIB::PrintStatus("[Default] BetaNS =", BetaNS);
            }
        }
    }

    TALYFEMLIB::ZEROPTV interpolate_csv_vel(const TALYFEMLIB::ZEROPTV &p) const {
        DENDRITE_REAL coord_x = p(csv_vel_coord[0]);
        DENDRITE_REAL interpolated_value = 0.0;
        /*
         * vel_csv_data[0] = coordinate
         * vel_csv_data[1] = velocity
         */
        DENDRITE_REAL min_x = *(std::min_element(vel_csv_data[0].begin(), vel_csv_data[0].end()));
        DENDRITE_REAL max_x = *(std::max_element(vel_csv_data[0].begin(), vel_csv_data[0].end()));
        if (coord_x < min_x or coord_x > max_x) { // fall out off the range
            interpolated_value = 0.0;
        } else {
            for (int i = 0; i < vel_csv_data[0].size(); i++) {
                if (vel_csv_data[0][i] > coord_x) {
                    interpolated_value =
                            (coord_x - vel_csv_data[0][i - 1]) / (vel_csv_data[0][i] - vel_csv_data[0][i - 1]) * (vel_csv_data[1][i] - vel_csv_data[1][i - 1]) + (vel_csv_data[1][i - 1]);
                    break;
                }
            }
        }
        TALYFEMLIB::ZEROPTV return_ptv;
        return_ptv(csv_vel_direction[0]) = interpolated_value;
        return return_ptv;
    }

    DENDRITE_REAL interpolate_csv_t(const TALYFEMLIB::ZEROPTV &p) const {
        DENDRITE_REAL coord_x = p(csv_t_coord[0]);
        DENDRITE_REAL interpolated_value = 0.0;
        DENDRITE_REAL min_x = *(std::min_element(t_csv_data[0].begin(), t_csv_data[0].end()));
        DENDRITE_REAL max_x = *(std::max_element(t_csv_data[0].begin(), t_csv_data[0].end()));
        if (coord_x < min_x or coord_x > max_x) {
            interpolated_value = 0.0;
        } else {
            for (int i = 0; i < t_csv_data[0].size(); i++) {
                if (t_csv_data[0][i] > coord_x) {
                    interpolated_value =
                            (coord_x - t_csv_data[0][i - 1]) / (t_csv_data[0][i] - t_csv_data[0][i - 1]) * (t_csv_data[1][i] - t_csv_data[1][i - 1]) + (t_csv_data[1][i - 1]);
                    break;
                }
            }
        }
        return interpolated_value;
    }
    template<class T>
    T getRamping(const DENDRITE_REAL time, const std::vector<T> &v) const {
        DENDRITE_UINT interval = ramping.size();
        for (DENDRITE_UINT i = 0; i < ramping.size(); i++) {
            if (time < ramping[i]) {
                interval = i;
                break;
            }
        }
        T current_variable;
        if (interval == 0) {
            current_variable = v[0] * (time / ramping[0]);
        } else if (interval == ramping.size()) {
            current_variable = v[interval - 1];
        } else {
            DENDRITE_REAL per1 = (time - ramping[interval - 1]) / (ramping[interval] - ramping[interval - 1]);
            DENDRITE_REAL per2 = (ramping[interval] - time) / (ramping[interval] - ramping[interval - 1]);
            current_variable = v[interval] * per1 + v[interval - 1] * per2;
        }
        return current_variable;
    }

    /*
     * returns the parabolic profile for x- locaiton (assert position[0] == 0)
     * */
    TALYFEMLIB::ZEROPTV parabolic_vel(const std::array<DENDRITE_REAL,DIM> &channel_min,
                          const std::array<DENDRITE_REAL,DIM> &channel_max,
                          const TALYFEMLIB::ZEROPTV &position,
                          const TALYFEMLIB::ZEROPTV &wall_vel,
                          TALYFEMLIB::ZEROPTV &inlet_vel) {
        if (fabs(position.x()) > 1e-6) {
            throw TALYFEMLIB::TALYException() << "parabolic profile should be x- position only";
        }
#if(DIM==3)
        /// this ensures that the total mass flow is the same with uniform velocity inlet
    double scaling_coeff = 36.0/pow(channel_max[1] - channel_min[1], 2)
        /pow(channel_max[2] - channel_min[2], 2);
    /// parabolic profile
    double parabolic = (position.y() - channel_min[1])*(channel_max[1] - position.y())
        *(position.z() - channel_min[2])*(channel_max[2] - position.z());
    inlet_vel = {scaling_coeff*wall_vel.norm()*parabolic, 0, 0};
#endif
#if(DIM==2)
        /// this ensures that the total mass flow is the same with uniform velocity inlet
        double scaling_coeff = 6.0/pow(channel_max[1] - channel_min[1], 2);
        /// parabolic profile
        double parabolic = (position.y() - channel_min[1])*(channel_max[1] - position.y());
        inlet_vel = {scaling_coeff*wall_vel.norm()*parabolic, 0, 0};
#endif
    }

    void PrintBoundaryDef(std::ofstream &fstream) const {
        int rank = TALYFEMLIB::GetMPIRank();
        if (!rank and fstream.is_open()) {

            if (side == Side::INVALID) {
                fstream << "side: INVALID\n";
            } else if (side == Side::X_MINUS) {
                fstream << "side: X_MINUS\n";
            } else if (side == Side::X_PLUS) {
                fstream << "side: X_PLUS\n";
            } else if (side == Side::Y_MINUS) {
                fstream << "side: Y_MINUS\n";
            } else if (side == Side::Y_PLUS) {
                fstream << "side: Y_PLUS\n";
            } else if (side == Side::Z_MINUS) {
                fstream << "side: Z_MINUS\n";
            } else if (side == Side::Z_PLUS) {
                fstream << "side: Z_PLUS\n";
            }

            if (vel_type == Vel_Type::NO_SLIP) {
                fstream << "vel_type: NO_SLIP\n";
            } else if (vel_type == Vel_Type::NO_GRADIENT_VEL) {
                fstream << "vel_type: NO_GRADIENT_VEL\n";
            } else if (vel_type == Vel_Type::DIRICHLET_VEL) {
                fstream << "vel_type: DIRICHLET_VEL\n";
            } else if (vel_type == Vel_Type::WEAK_VEL) {
                fstream << "vel_type: WEAK_VEL\n";
            }

            if (pressure_type == Pressure_Type::NO_GRADIENT_P) {
                fstream << "pressure_type: NO_GRADIENT_P\n";
            } else if (pressure_type == Pressure_Type::DIRICHLET_P) {
                fstream << "pressure_type: DIRICHLET_P\n";
            } else if (pressure_type == Pressure_Type::BACKFLOW_STAB) {
                fstream << "pressure_type: BACKFLOW_STAB\n";
            }

            fstream << "wall_velocity: [";
            for (const auto &v : wall_velocity) {
                for (int i = 0; i < 3; i++) {
                    fstream << v[i] << ", ";
                }
            }
            fstream << "]\n";
            fstream << "pressure: " << pressure << "\n";
            fstream << "ifRefine: " << ifRefine << "\n";
        }
    }

private:
    std::string vel_csv_path;
    std::vector<int> csv_vel_coord;
    std::vector<int> csv_vel_direction;
    std::vector<std::vector<DENDRITE_REAL>> vel_csv_data;

    std::string t_csv_path;
    std::vector<int> csv_t_coord;
    std::vector<int> csv_t_direction;
    std::vector<std::vector<DENDRITE_REAL>> t_csv_data;

    static Side str_to_side(const std::string &str) {
        if (str == "x-") {
            return Side::X_MINUS;
        } else if (str == "x+") {
            return Side::X_PLUS;
        } else if (str == "y-") {
            return Side::Y_MINUS;
        } else if (str == "y+") {
            return Side::Y_PLUS;
        } else if (str == "z-") {
            return Side::Z_MINUS;
        } else if (str == "z+") {
            return Side::Z_PLUS;
        } else {
            throw TALYFEMLIB::TALYException() << "Invalid BC side";
        }
    }

    static Vel_Type read_vel_type(const std::string &str) {
        if (str == "no_slip") {
            return Vel_Type::NO_SLIP;
        } else if (str == "no_gradient") {
            return Vel_Type::NO_GRADIENT_VEL;
        } else if (str == "dirichlet") {
            return Vel_Type::DIRICHLET_VEL;
        } else if (str == "weak") {
            return Vel_Type::WEAK_VEL;
        } else if (str == "parabolic") {
            return Vel_Type::PARABOLIC;
        } else if (str == "TOP_HAT") {
            return Vel_Type::TOP_HAT;
        } else if (str == "disturbance") {
            return Vel_Type::DISTURBANCE;
        } else {
            throw TALYFEMLIB::TALYException() << "Invalid BC type for vel";
        }
    }

    static Pressure_Type read_pressure_type(const std::string &str) {
        if (str == "no_gradient") {
            return Pressure_Type::NO_GRADIENT_P;
        } else if (str == "dirichlet") {
            return Pressure_Type::DIRICHLET_P;
        } else if (str == "BACKFLOW_STAB"){
            return Pressure_Type::BACKFLOW_STAB;
        } else {
            throw TALYFEMLIB::TALYException() << "Invalid BC type for pressure";
        }
    }

};


/**
 * Parameters for Carved out geometry
 */
struct CarvedOutGeom
{
    enum Type
    {
        INVALID = 0,
        SPHERE = 1,
        PILLAR = 2,
        MESHOBJECT = 3,
        MESHOBJECT_2D = 4,
        CYLINDER = 5,
        CUBE = 6,
        CIRCLE_2D_VOXEL = 7,
        SPHERE_3D_VOXEL = 8,
        BOX_2D_VOXEL = 9,
        CUBE_3D_VOXEL = 10,
        D
    };
    enum FileFormat
    {
        NOTKNOWN = 0,
        GMSH = 1,
        STL = 2,
        ML = 3
    };
    enum BCType
    {
        INVALID_BC = 0,
        WEAK = 1,
        SBM = 2,
        PENALTY_SBM = 3
    };
    enum InoutTYPE
    {
        UNKNOWN = -1,
        PRIMITIVE = 0,
        RT = 1,
    };

    Type type = INVALID;
    ///< path to surface mesh (currently must be triangular)
    std::string mesh_path;
    ///< name of the geometry (used for printing file and etc..)
    std::string name;
    FileFormat fileformat = NOTKNOWN;
    ///< refinement level for future use
    unsigned int refine_lvl = 0;
    ///< initial displacement for the geometry (used only once to move geometry to desired location)
    TALYFEMLIB::ZEROPTV InitialDisplacement;
    ///< retain inside or outside (default is false == retain outside)
    bool outer_boundary = false;

    /// variables depending on type...
    TALYFEMLIB::ZEROPTV center_of_mass;
    ///< radius for sphere or cylinder types
    DENDRITE_REAL radius = 0.0;
    ///< cylinder orientation, used for in_out test
    int cylinderDir = 2;
    double height = -1.0;
    TALYFEMLIB::ZEROPTV bottom_center;
    ///< first one is bottom left back corner, the second one is dimension in x-y-z
    std::vector<TALYFEMLIB::ZEROPTV> cube_dim{2};

    /// auto refine for IMGA loop
#ifdef IBM
    GeomRefinement geomRefine;
#endif



    /// boundary conditions for geometry
    std::vector<CarvedOutGeom::BCType> bc_type_V = {};
    std::vector<DENDRITE_REAL> dirichlet_V = {};

    std::vector<DENDRITE_REAL> getBC(int dof) const
    {
        if (bc_type_V[dof] == CarvedOutGeom::BCType::WEAK
            or bc_type_V[dof] == CarvedOutGeom::BCType::SBM
            or bc_type_V[dof] == CarvedOutGeom::BCType::PENALTY_SBM)
        {
            return std::vector<DENDRITE_REAL>{dirichlet_V.at(dof)};
        }
    }

    void read_from_config(const libconfig::Setting &root)
    {
        if (root.exists("refine_lvl"))
        {
            refine_lvl = (unsigned int)root["refine_lvl"];
        }

        if (root.exists("outer_boundary"))
        {
            outer_boundary = (bool)root["outer_boundary"];
        }

        type = str_to_type(root["type"]);
        if (type != Type::CIRCLE_2D_VOXEL and type != Type::SPHERE_3D_VOXEL and
            type != Type::BOX_2D_VOXEL and type != Type::CUBE_3D_VOXEL)
        {
            mesh_path = (const char *)root["mesh_path"];
            if (root.exists("name"))
            {
                name = (const char *)root["name"];
            }
            else
            {
                if (mesh_path.find_last_of('/'))
                {
                    name = mesh_path.substr(mesh_path.find_last_of('/') + 1,
                                            mesh_path.find_last_of('.') - mesh_path.find_last_of('/') - 1);
                }
                else
                {
                    name = mesh_path.substr(0, mesh_path.find_last_of('.') - 1);
                }
            }
            fileformat = mesh_path_to_file_format(mesh_path);
        }
        else
        {
            name = (const char *)root["name"];
        }
        /// extra parameter for each primitive type
        switch (type)
        {
            case Type::SPHERE:
                radius = (DENDRITE_REAL)root["radius"];
                ReadZEROPTV(root, "com", center_of_mass);
                break;
            case Type::PILLAR:
                radius = (DENDRITE_REAL)root["radius"];
                ReadZEROPTV(root, "com", center_of_mass);
                if (root.exists("cylinderDir"))
                {
                    cylinderDir = (int)root["cylinderDir"];
                }
                break;
            case Type::CYLINDER:
                radius = (DENDRITE_REAL)root["radius"];
                ReadZEROPTV(root, "com", center_of_mass);
                height = (DENDRITE_REAL)root["height"];
                if (root.exists("cylinderDir"))
                {
                    cylinderDir = (int)root["cylinderDir"];
                }
                ReadZEROPTV(root, "bottom_center", bottom_center);
                bottom_center += InitialDisplacement;
                break;
            case Type::CUBE:
                ReadZEROPTV(root, "min_corner", cube_dim[0]);
                ReadZEROPTV(root, "dim", cube_dim[1]);
                break;
            case Type::CIRCLE_2D_VOXEL:
                radius = (DENDRITE_REAL)root["radius"];
                ReadZEROPTV(root, "com", center_of_mass);
                break;
            case Type::SPHERE_3D_VOXEL:
                radius = (DENDRITE_REAL)root["radius"];
                ReadZEROPTV(root, "com", center_of_mass);
                break;
            case Type::BOX_2D_VOXEL:
                ReadZEROPTV(root, "min_corner", cube_dim[0]);
                ReadZEROPTV(root, "dim", cube_dim[1]);
                break;
            case Type::CUBE_3D_VOXEL:
                ReadZEROPTV(root, "min_corner", cube_dim[0]);
                ReadZEROPTV(root, "dim", cube_dim[1]);
                break;
            default:
                assert(type == CarvedOutGeom::MESHOBJECT || type == CarvedOutGeom::MESHOBJECT_2D);
        }

        geomRefine.maxSplitIteration = 100;
        geomRefine.octantLevel = refine_lvl;
        geomRefine.ratioArea = 0.25;

        const libconfig::Setting &bc_dof = root["bc_type_V"];
        for (int i = 0; i < bc_dof.getLength(); ++i)
        {
            std::string a = bc_dof[i];
            CarvedOutGeom::BCType temp = str_to_bctype(bc_dof[i]);
            bc_type_V.push_back(temp);
        }
        for (int i = 1; i < 3; i++)
        {
            if (bc_type_V[i] != bc_type_V[0])
            {
                throw TALYFEMLIB::TALYException() << "weakBC for velocity have to be consistent!";
            }
        }
        int d_iter = 0, n_iter = 0, r_iter = 0;
        for (int i = 0; i < bc_type_V.size(); i++)
        {
            if ( bc_type_V.at(i) == CarvedOutGeom::BCType::WEAK
                or bc_type_V.at(i) == CarvedOutGeom::BCType::SBM
                or bc_type_V.at(i) == CarvedOutGeom::BCType::PENALTY_SBM)
            {
                const libconfig::Setting &temp1 = root["dirichlet_V"];
                dirichlet_V.push_back(DENDRITE_REAL(temp1[d_iter++]));
            }
        }

        if (type != Type::CIRCLE_2D_VOXEL and type != Type::SPHERE_3D_VOXEL and
            type != Type::BOX_2D_VOXEL and type != Type::CUBE_3D_VOXEL)
        {
            ReadZEROPTV(root, "position", InitialDisplacement);
        }

    }

protected:

    /// read size of 3 vector as ZEROPTV
    void ReadZEROPTV(const libconfig::Setting &root, const std::string &key_name, TALYFEMLIB::ZEROPTV &value, bool required = true) {
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
                throw TALYFEMLIB::TALYException() << key_name + " have size of " + std::to_string(root[key_name.c_str()].getLength());
            }
        } else if (required) {
            throw TALYFEMLIB::TALYException() << key_name + " doesn't exist!";
        }
    }

    /// read mesh type
    Type str_to_type(const std::string &str) const
    {
        if (str == "sphere")
        {
            return Type::SPHERE;
        }
        else if (str == "pillar")
        {
            return Type::PILLAR;
        }
        else if (str == "cylinder")
        {
            return Type::CYLINDER;
        }
        else if (str == "cube")
        {
            return Type::CUBE;
        }
        else if (str == "meshobject")
        {
            return Type::MESHOBJECT;
        }
        else if (str == "meshobject_2d")
        {
            return Type::MESHOBJECT_2D;
        }
        else if (str == "circle_2d")
        {
            return Type::CIRCLE_2D_VOXEL;
        }
        else if (str == "sphere_3d")
        {
            return Type::SPHERE_3D_VOXEL;
        }
        else if (str == "box_2d")
        {
            return Type::BOX_2D_VOXEL;
        }
        else if (str == "cube_3d")
        {
            return Type::CUBE_3D_VOXEL;
        }
        else if (str == "deeptrace")
        {
            return Type::
        }
        else
        {
            throw TALYFEMLIB::TALYException() << "Invalid geometry type '" << str
                                              << "' (expected sphere, pillar, cylinder cube or meshobject(2D) )"
                                              << "\n for voxel, circle_2d or sphere_3d";
        }
    }

    FileFormat mesh_path_to_file_format(const std::string &str) const
    {
        if (str.substr(str.find_last_of('.') + 1) == "msh")
        {
            return FileFormat::GMSH;
        }
        else if (str.substr(str.find_last_of('.') + 1) == "stl")
        {
            return FileFormat::ML;
        }
        else
        {
            throw TALYFEMLIB::TALYException() << "Invalid file extension '" << str.substr(str.find_last_of('.') + 1)
                                              << "' (support .msh, .stl)";
        }
    }

    /// read bc type
    BCType str_to_bctype(const std::string &str) const
    {
# ifdef IBM
        if (str == "sbm")
        {
            return BCType::SBM;
        }
        else if (str == "penalty_sbm")
        {
            return BCType::PENALTY_SBM;
        }
        else if (str == "weak")
        {
            return BCType::WEAK;
        }
        else
        {
            throw TALYFEMLIB::TALYException() << "Invalid bounday condition type '" << str
                                              << "' (expected weak, sbm or penalty_sbm)";
        }
#endif
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

    /// Boundary conditions
    std::vector<BoundaryDef> boundary_def;

    /// Curved-out geometries
    std::vector<CarvedOutGeom> carved_out_geoms_def;



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
        if (ReadValue("MMS", ifMMS)) {}
        if (ReadValue("dump_vec", dump_vec)) {}
        if (ReadValue("VelocityExtrapolationOrder", VelocityExtrapolationOrder)) {}
        solverOptionsNS = read_solver_options(cfg, "solver_options");


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
        if (cfg.exists("boundary")) {
            const auto &bc = cfg.getRoot()["boundary"];
            for (auto &bc_def : boundary_def) {
                for (int j = 0; j < bc.getLength(); j++) {
                    /// If the side of bc in config matches with the preset bc
                    if (bc_def.side == BoundaryDef::read_side_from_config(bc[j])) {
                        bc_def.read_from_config(bc[j]);
                    }
                }
            }
        }

        if (cfg.exists("geometries")) {
            const auto &geometries = cfg.getRoot()["geometries"];
            carved_out_geoms_def.resize(geometries.getLength());
            for (unsigned int i = 0; i < carved_out_geoms_def.size(); i++) {
                carved_out_geoms_def[i].read_from_config(geometries[i]);
            }
        }

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

    int elemOrder = 1; // Hard coded for now
};

#endif //DENDRITEKT_NSINPUTDATA_H
