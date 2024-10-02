//
// Created by gkluhana on 25/03/24.
//
#ifndef IBIMPLICIT_UTILITIES_H
#define IBIMPLICIT_UTILITIES_H
#include <petscksp.h>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <petscdm.h>
#include <petsclog.h>
#include <debug.h>
#include <fstream>
namespace ibImplicit {
    /// Type of Boundary Condition
    PetscErrorCode getTruncatedMGLevels(PetscInt powerOfTwo, PetscInt *result);
    PetscInt RoundUpToNextPowerOf2(PetscInt n);
    PetscInt RoundDownToPowerOf2(PetscInt n);
    inline PetscBool stringToPetscBool(const std::string& str) {
        // Check if the string is "1"
        if (str == "1") {
            return PETSC_TRUE;
        }
            // Check if the string is "0"
        else if (str == "0") {
            return PETSC_FALSE;
        }
            // If the string is neither "1" nor "0", handle the error case or return a default value
        else {
            // Handle error or return a default value
            return PETSC_FALSE; // Default value
        }
    }

    typedef enum {
        DIRICHLET,
        NEUMANN,
        TRUNCATED_BC
    } BCType;


    typedef enum {
        EXACT,
        LSC,
        TRUNCATED_SC
    } SCType;

    typedef struct BCIndices_ {
        //Boundary Indices
        std::vector<int> u_e;
        std::vector<int> u_w;
        std::vector<int> v_s;
        std::vector<int> v_n;
        //Near Boundary Indices
        std::vector<int> u_s_near;
        std::vector<int> u_n_near;
        std::vector<int> v_e_near;
        std::vector<int> v_w_near;
        // Pressure Boundary Indices
        std::vector<int> p_s_near;
        std::vector<int> p_e_near;
        std::vector<int> p_n_near;
        std::vector<int> p_w_near;
    }BCIndicesData;
    typedef struct BCStaggeredType_ {
        BCType u_n;
        BCType u_s;
        BCType u_e;
        BCType u_w;
        BCType p_n;
        BCType p_s;
        BCType p_e;
        BCType p_w;
    }BCStaggeredTypeData;
    typedef BCStaggeredTypeData *BCStaggeredType;
    typedef struct BCFunctions2D_{
        void (*fe)(double x, double y,double mu, double hx, double hy, double *f);
        void (*fw)(double x, double y,double mu, double hx, double hy, double *f);
        void (*fs)(double x, double y,double mu, double hx, double hy, double *f);
        void (*fn)(double x, double y,double mu, double hx, double hy, double *f);
        void (*ge)(double x, double y,double mu, double hx, double hy, double *g);
        void (*gw)(double x, double y,double mu, double hx, double hy, double *g);
        void (*gs)(double x, double y,double mu, double hx, double hy, double *g);
        void (*gn)(double x, double y,double mu, double hx, double hy, double *g);
    }BCFunctions2DData;
    typedef BCFunctions2DData *BCFunctions2D;
    typedef BCIndicesData *BCIndices;
    typedef struct BCSpec_ {
        BCIndices bc_idx;
        BCStaggeredType bc_stag_type;
        BCFunctions2D bc_funcs;
    }BCSpecData;
    typedef BCSpecData *BCSpec;
    struct SpringParametersData{
        // TODO: Convert these to Petsc Vecs
        std::vector<int> master_idx;
        std::vector<int> slave_idx;
        std::vector<float> kappa;
        std::vector<float> rest_length;
        PetscScalar uniform_spring_stiffness;
        PetscInt num_edges;
    };
    typedef struct SpringParametersData SpringParametersData;
    typedef SpringParametersData *SpringParameters;

    struct LevelCtxData{
        DM          dm_level, dm_faces, dm_center;
        PetscReal   hx, hy, hz;
    };
    typedef struct LevelCtxData LevelCtxData;

    typedef LevelCtxData *LevelCtx;
    typedef struct MGParameters_{
        PetscInt n_levels;
        PetscInt level_curr;
        PetscBool faces_only;
        PetscInt RefRatio;
        LevelCtx *levels;
    }MGParametersData;
    typedef MGParametersData *MGParameters;

    struct SystemParametersData{
        PetscReal lx,Lx; // x dim is [lx, Lx]
        PetscReal ly,Ly; // y dim is [ly, Ly]
        PetscInt Nx, Ny;
        PetscReal hx, hy;
        PetscInt nfluid, np, nu, nv;
        PetscReal dt;
        PetscReal mu; // viscosity
        PetscReal rho;
        PetscBool matrix_free;
        BCSpec bcSpec; //boundary specification
        PetscInt ns, *nS;

        PetscInt n_structs;
        std::string *structName; // structure File Name
        std::string *structNames;
        Vec* X0; //initial configuration of structure
        Vec* Xs; //inital configuration of all structures
        SpringParameters springData;
        SpringParameters *springDatas;
        MGParameters mg_params;
        PetscInt pre_block_type;
        PetscInt truncated_pre;
        PetscInt truncated_sys;
        BCType truncated_BCtype;
        SCType s1_type;
        SCType s2_type;

        void recomputeGridParameters();
        void recomputeGridParametersTruncated(PetscReal lx_f, PetscReal ly_f);
        PetscErrorCode clean();
        ~SystemParametersData(){
            if (X0)VecDestroy(X0);
            if (Xs)VecDestroy(Xs);
        }
    };
    typedef struct SystemParametersData SystemParametersData;
    typedef SystemParametersData *SystemParameters;
    PetscErrorCode createSystemParameters(SystemParameters *sys, std::unordered_map<std::string,std::string>* data);
    PetscErrorCode dumpISToFile(ISLocalToGlobalMapping is, std::string str);
    PetscErrorCode dumpVecToFile(Vec v, std::string str);
    PetscErrorCode dumpMatToFile(Mat A, std::string str);
    PetscErrorCode dumpMatNestToFile(Mat M, std::string str);
    /// Error Function
    class ibError : public std::runtime_error {
    public:
        ibError() : std::runtime_error("No error") {} // Default constructor
        ibError(const std::string &message) : std::runtime_error(message) {}
    };

}
#endif
