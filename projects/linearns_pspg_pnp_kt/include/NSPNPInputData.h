//
// Created by kumar on 7/16/20.
//

#ifndef DENDRITEKT_NSINPUTDATA_H
#define DENDRITEKT_NSINPUTDATA_H
#include <talyfem/input_data/input_data.h>
#include <talyfem/talyfem.h>
#include <DataTypes.h>
#include <point.h>
#include <DendriteUtils.h>
#include <map>
#include <vector>
#include <string>



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


class NSPNPInputData: public TALYFEMLIB::InputData {
 public:
    static  constexpr int nsd = DIM;
	enum NS_FORMULATION: u_short {
			PSPG = 0,
			PROJECTION = 1
	};

	enum TIME_STEPPING: u_short {
			BDF = 0,
			THETA = 1
	};

	/// Integrands type
	enum IntegrandsType: unsigned int {
			NS = 0,
			Stokes = 1,
			Heat = 2,
	};

	DENDRITE_UINT nsOrder;
	bool ifUseRotationalForm = false;
	DomainInfo physDomain;

	DENDRITE_REAL Re;
	DENDRITE_REAL Fr;

	/// VMS model parameters
	DENDRITE_REAL Ci_f;
	DENDRITE_REAL scale;

	DENDRITE_REAL timeStep;
	DENDRITE_REAL TotalTime;


	SolverOptions solverOptionsMomentum, solverOptionsPNP;
	bool ifMMS = false;
	DENDRITE_UINT pExtrapOrder = 2;

	/// Order of extrapolation
	int advectionExtrapolationOrder;

    /// Basic parameters of the case
    DENDRITE_UINT elemOrder;
    bool ifMatrixFree = false;

	/// Integrands type
	unsigned int integrandType;
	/// Time scheme
	TIME_STEPPING timescheme;

	/// Whether to use linearNS
	bool ifLinearNS = false;

	/// PNP Stuff
	std::vector<DENDRITE_REAL> A_,B_,C_,D_; ///<coefficients
	DENDRITE_REAL E_; /// Coefficient for Poisson
	DENDRITE_REAL Ci_f_pnp = 4.0; ///< for vms

	/// PNP BC
	/// Initial and boundary conditions
	std::vector<DENDRITE_REAL> initU;     /// initial velocity [m/s]
	DENDRITE_REAL initP;             /// initial Pressure [Pa]
	std::vector<DENDRITE_REAL> initC;     /// initial concentration of each ion [mol/m^3]
	DENDRITE_REAL initPhi;           /// initial potential [V]
	int numDBCUx;             /// number of Ux Dirichlet boundaries
	std::vector<DENDRITE_UINT> DBCIdxUx;     /// index of Ux Dirichlet BC
	std::vector<Expression> DBCUx;     /// value of Ux Dirichlet BC
	int numDBCUy;             /// number of Uy Dirichlet boundaries
	std::vector<DENDRITE_UINT> DBCIdxUy;     /// index of Uy Dirichlet BC
	std::vector<Expression> DBCUy;     /// value of Uy Dirichlet BC
	std::vector<DENDRITE_REAL> pointPreBC; /// position for point pressure BC [x, y, z (if exist), value]
	int numDBCP;              /// number of pressure Dirichlet boundaries
	std::vector<DENDRITE_UINT> DBCIdxP;      /// index of pressure Dirichlet BC
	std::vector<Expression> DBCP;      /// value of pressure Dirichlet BC
	int numDBCPhi;             /// number of potential Dirichlet boundaries
	std::vector<DENDRITE_UINT> DBCIdxPhi;    /// index of potential Dirichlet BC
	std::vector<Expression> DBCPhi;    /// value of potential Dirichlet BC
	int numDBCC1;              /// number of species 1 Dirichlet BC
	std::vector<DENDRITE_UINT> DBCIdxC1;     /// index of species 1 Dirichlet BC
	std::vector<Expression> DBCC1;     /// value of species 1 Dirichlet BC
	int numDBCC2;              /// number of species 2 Dirichlet BC
	std::vector<DENDRITE_UINT> DBCIdxC2;     /// index of species 2 Dirichlet BC
	std::vector<Expression> DBCC2;     /// value of species 2 Dirichlet BC
	int numWCBC;              /// number of wall charge BC
	std::vector<DENDRITE_UINT> WCBCIdx;     /// index of wall charge BC
	std::vector<Expression> WCBC;     /// value of wall charge BC  [C/m^2]
	int numWeakBCPhi;              /// number of weak BC
	std::vector<DENDRITE_UINT> weakBCPhiIdx;     /// index of weak BC
	std::vector<Expression> weakBCPhi;     /// value of weak BC  [C/m^2]
	int numWeakBCC1;              /// number of weak BC
	std::vector<DENDRITE_UINT> weakBCC1Idx;     /// index of weak BC
	std::vector<Expression> weakBCC1;     /// value of weak BC  [C/m^2]
	int numWeakBCNS;  /// Number of NS weak bc boundaries
	std::vector<int> weakBCUIdx;
	std::vector<Expression> weakBCU; /// weakBCU = [u_ij] where i is idx for coordinate, j is weak boundary idx

	DENDRITE_REAL gamma_NS;
	DENDRITE_REAL Cb_NS;
	DENDRITE_REAL Cb_Poisson;  /// weak BC penalty term coeff Cb for Poisson equation
	DENDRITE_REAL Cb_NP;       /// weak BC penalty term coeff Cb for Nernst-Planck equation
	DENDRITE_REAL gamma_Poisson;  /// weak BC adjoint term coeff gamma for Poisson equation (either 1 or -1)
	DENDRITE_REAL gamma_NP;       /// weak BC adjoint term coeff gamma for Nernst-Planck equation (either 1 or -1)

	DENDRITE_UINT numSp4FLX; /// number of species for boundary flx calculation
	int numB4FLX = false; /// number of boundaries to be calculated
	std::vector<DENDRITE_UINT> BFLXIdx;  /// physical boundary index for boundary flux calculation
	std::vector<DENDRITE_UINT> BFLXSpecies;  /// species to calculate boundary flux

    /// Usage: dimensions which are non-periodic are set to zero.  Periodic dimensions are set to ratio of their
    // length to largest length (in some direction).  Periodic dimension is always set to a power of 2.
    /// e.g for a 2d mesh with dimensions [4,2],
    /// 1) setting periodicScale to [0.0,0.5] would set Y periodicity
    /// 2) setting periodicScale to [1.0, 0.0] would set X periodicity
    /// 3) setting periodicScale to [1.0, 0.5] would set both X and Y periodicity
    /// 4) setting periodicScale to [0.0, 0.0] would no periodicity
    /// Note that the values passed are the ratios to the largest dimension, and powers of 2.
    std::vector<DENDRITE_REAL> periodicScale; ///< Variable for specifying periodicity, iniliazed to zero for non-periodicity


    ///< Checkpointing options
    DENDRITE_UINT checkpointFrequency;
    DENDRITE_UINT numberOfBackups;

    /// File writing options
    DENDRITE_UINT OutputSpan; // number of timesteps between each data write (NOT time)

    /// Mesh definition
    MeshDef meshDef;

    // Physical Domain min and max coordinates
    std::array<DENDRITE_REAL,DIM> physDomainMin;
    std::array<DENDRITE_REAL,DIM> physDomainMax;


    /// NS Non-dimensionalized parameters
    DENDRITE_REAL ndcf_time_;
    DENDRITE_REAL ndcf_conv_;
    DENDRITE_REAL ndcf_diff_;
    DENDRITE_REAL ndcf_pres_;
    DENDRITE_REAL ndcf_pnp_coupling_;

    DENDRITE_REAL rectLowerLeft[DIM];
    DENDRITE_REAL rectUpperRight[DIM];
    DENDRITE_UINT rect_lvl;
    DENDRITE_REAL restartTimeEof = 0.0;


    /// Parameters for EOF Test case /// initialized to numbers for L = 8 case
    DENDRITE_REAL eofPNPonlybottomPotential = -2.707;
    DENDRITE_REAL eofInletPotential = 0.30944;
    DENDRITE_REAL eofbottomPotetialSlope = -0.03868;
    DENDRITE_REAL eofbottomPotetialIntercept = -2.39756;


	bool ifPurePNP;


	static NS_FORMULATION readFormulationType(libconfig::Setting &root,
	                                          const char *name) {
		std::string str;
		/// If nothing specified stays stabilizedNS
		if (root.lookupValue (name, str)) {
			if (str == "PSPG") {
				return PSPG;
			} else if (str == "PROJECTION") {
				return PROJECTION;
			} else {
				throw TALYFEMLIB::TALYException () << "Unknown solver name for NS: " << name << str;
			}
		} else {
			throw TALYFEMLIB::TALYException () << "Must specify Formulation: PSPG or PROJECTION";
		}

	}

	static TIME_STEPPING readTimeStepper(libconfig::Setting &root,
	                                     const char *name) {
		std::string str;
		/// If nothing specified
		if (root.lookupValue (name, str)) {
			if (str == "CN") {
				return TIME_STEPPING::THETA;
			} else if (str == "BE") {
				return TIME_STEPPING::BDF;
			} else {
				throw TALYFEMLIB::TALYException () << "Unknown time stepping for NS: " << name << str;
			}
		} else {
			throw TALYFEMLIB::TALYException () << "Must specify timeStepping: CN or BE";
		}

	}
	bool ReadFromFile(const std::string &filename = std::string ("config.txt")) {
		ReadConfigFile (filename);

        /// mesh size and level
        meshDef.read_from_config(cfg.getRoot()["background_mesh"]);


        /// Basic parameters of the case
        elemOrder = 1;
        ReadValueRequired ("elemOrder" , elemOrder);
        ifMatrixFree = false;
        ReadValueRequired ("ifMatrixFree", ifMatrixFree);


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

		ReadValueRequired ("Re", Re);
		ReadValueRequired ("Fr", Fr);


        /// NS Non-dimensionalized parameters
        ReadValueRequired ("ndcf_time", ndcf_time_);
        ReadValueRequired ("ndcf_conv", ndcf_conv_);
        ReadValueRequired ("ndcf_diff", ndcf_diff_);
        ReadValueRequired ("ndcf_pres", ndcf_pres_);
        ReadValueRequired ("ndcf_pnp_coupling", ndcf_pnp_coupling_);


        if (ReadValue("restartTimeEof", restartTimeEof)) {
        }


        ReadValueRequired("eofPNPonlybottomPotential",eofPNPonlybottomPotential);
        ReadValueRequired("eofInletPotential",eofInletPotential);
        ReadValueRequired("eofbottomPotetialSlope",eofbottomPotetialSlope);
        ReadValueRequired("eofbottomPotetialIntercept",eofbottomPotetialIntercept);



		/// VMS model parameters
		ReadValueRequired ("Ci_f", Ci_f);
		ReadValueRequired ("scale", scale);



		ReadValueRequired ("nsOrder", nsOrder);
		solverOptionsMomentum = read_solver_options (cfg, "solver_options_momentum");
		solverOptionsPNP = read_solver_options(cfg, "solver_options_pnp");

		if (ReadValue ("MMS", ifMMS)) {}

		if (ReadValue ("ifLinearNS", ifLinearNS)) {}

		ifUseRotationalForm = false;
		if (ReadValue ("ifUseRotationalForm", ifUseRotationalForm)) {}
		pExtrapOrder = 1;
		if (ReadValue ("PressureExtrapolationOrder", pExtrapOrder)) {}
		advectionExtrapolationOrder = 2;
		ReadValue("advectionExtrapolationOrder", advectionExtrapolationOrder);

		timescheme = BDF; // default method
		std::string time_scheme_name ("bdf");
		if (ReadValue ("timeScheme", time_scheme_name)) {
			if (time_scheme_name == "theta") {
				timescheme = THETA;
			}
		}
		ReadValue ("dt",timeStep);
		ReadValue ("totalTime",TotalTime);
		integrandType = IntegrandsType::NS; // by default NS
		if (ReadValue ("integrandType", integrandType)) {}




		//// PNP stuff
		if (ReadArray("A_", A_, 2)) {
		}
		if (ReadArray("B_", B_, 2)) {
		}
		if (ReadArray("C_", C_, 2)) {
		}
		if (ReadArray("D_", D_, 2)) {
		}
		if (ReadValue("E_", E_)) {
		}
		if (ReadValue("Ci_f_pnp", Ci_f_pnp)) {
		}
		if (ReadValue("ifPurePNP", ifPurePNP)) {
		}

		/// Initial and boundary conditions
		if (ReadArray("initU", initU, nsd)) {
		}
		if (ReadValue("initP", initP)) {
		}
		if (ReadValue("initPhi", initPhi)) {
		}
		if (ReadValue("numDBCUx", numDBCUx)) {
			/// receive Dirichlet BC if exist
			if (numDBCUx != 0)  {
				if (ReadArray("DBCIdxUx", DBCIdxUx, numDBCUx)) {
				}
				if (ReadArray("DBCUx", DBCUx, numDBCUx)) {
				}
			}
		}
		if (ReadValue("numDBCUy", numDBCUy)) {
			/// receive Dirichlet BC if exist
			if (numDBCUy != 0)  {
				if (ReadArray("DBCIdxUy", DBCIdxUy, numDBCUy)) {
				}
				if (ReadArray("DBCUy", DBCUy, numDBCUy)) {
				}
			}
		}
		if (ReadArray("pointPreBC", pointPreBC, nsd+1))  { /// if pressure is defined at a single point
			numDBCP = 0;
		} else {
			if (ReadValue("numDBCP", numDBCP)) {
				/// receive Dirichlet BC if exist
				if (numDBCP != 0) {
					if (ReadArray("DBCIdxP", DBCIdxP, numDBCP)) {
					}
					if (ReadArray("DBCP", DBCP, numDBCP)) {
					}
				}
			}
		}
		if (ReadValue("numDBCPhi", numDBCPhi)) {
			/// receive Dirichlet BC if exist
			if (numDBCPhi != 0) {
				if (ReadArray("DBCIdxPhi", DBCIdxPhi, numDBCPhi)) {
				}
				if (ReadArray("DBCPhi", DBCPhi, numDBCPhi)) {
				}
			}
		}
		if (ReadValue("numDBCC1", numDBCC1)) {
			/// receive Dirichlet BC if exist
			if (numDBCC1 != 0)  {
				if (ReadArray("DBCIdxC1", DBCIdxC1, numDBCC1)) {
				}
				if (ReadArray("DBCC1", DBCC1, numDBCC1)) {
				}
			}
		}
		if (ReadValue("numDBCC2", numDBCC2)) {
			/// receive Dirichlet BC if exist
			if (numDBCC2 != 0)  {
				if (ReadArray("DBCIdxC2", DBCIdxC2, numDBCC2)) {
				}
				if (ReadArray("DBCC2", DBCC2, numDBCC2)) {
				}
			}
		}
		if (ReadValue("numWeakBCNS",numWeakBCNS)) {
			if (ReadArray("weakBCUIdx",weakBCUIdx,numWeakBCNS)) {
			}
			if (ReadArray("weakBCU",weakBCU,nsd*numWeakBCNS)) {
			}
		}
		if (ReadValue("numWCBC", numWCBC)) {
			/// receive Dirichlet BC if exist
			if (numWCBC != 0)  {
				if (ReadArray("WCBCIdx", WCBCIdx, numWCBC)) {
				}
				if (ReadArray("WCBC", WCBC, numWCBC)) {
				}
			}
		}
		if (ReadValue("numWeakBCPhi", numWeakBCPhi)) {
			/// receive Dirichlet BC if exist
			if (numWeakBCPhi != 0)  {
				if (ReadArray("weakBCPhiIdx", weakBCPhiIdx, numWeakBCPhi)) {
				}
				if (ReadArray("weakBCPhi", weakBCPhi, numWeakBCPhi)) {
				}
			}
		}
		if (ReadValue("numWeakBCC1", numWeakBCC1)) {
			/// receive Dirichlet BC if exist
			if (numWeakBCC1 != 0)  {
				if (ReadArray("weakBCC1Idx", weakBCC1Idx, numWeakBCC1)) {
				}
				if (ReadArray("weakBCC1", weakBCC1, numWeakBCC1)) {
				}
			}
		}

		if (ReadValue("Cb_Poisson", Cb_Poisson)) {}
		if (ReadValue("Cb_NP", Cb_NP)) {}
		if (ReadValue("gamma_Poisson", gamma_Poisson)) {}
		if (ReadValue("gammam_NP", gamma_NP)) {}

		if (ReadValue("numSp4FLX", numSp4FLX)) {
			/// receive Dirichlet BC if exist
			if (numSp4FLX != 0) {
				if (ReadValue("numB4FLX", numB4FLX)) {}
				if (ReadArray("BFLXIdx", BFLXIdx, numB4FLX)) {
				}
				if (ReadArray("BFLXSpecies", BFLXSpecies, numSp4FLX)) {
				}
			}
		}

		return true;
	}
};
#endif //DENDRITEKT_NSINPUTDATA_H
