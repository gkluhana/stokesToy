//
// Created by gkluhana on 25/03/24.
//
#ifndef IBIMPLICIT_STAGGERED_H
#define IBIMPLICIT_STAGGERED_H
#include <iostream>
#include <boundary.h>
#include <petscvec.h>
#include <utilities.h>
#include <lagrangian.h>
#include <petscdm.h>
#include <petscdmstag.h>
#include <petscdmda.h>
/* Shorter, more convenient names for DMStagStencilLocation entries */
#define DOWN_LEFT  DMSTAG_DOWN_LEFT
#define DOWN       DMSTAG_DOWN
#define DOWN_RIGHT DMSTAG_DOWN_RIGHT
#define LEFT       DMSTAG_LEFT
#define ELEMENT    DMSTAG_ELEMENT
#define RIGHT      DMSTAG_RIGHT
#define UP_LEFT    DMSTAG_UP_LEFT
#define UP         DMSTAG_UP
#define UP_RIGHT   DMSTAG_UP_RIGHT

namespace ibImplicit {

    PetscErrorCode DumpStokesData(DM dm, Vec x, std::string var_name);
    PetscErrorCode DumpPressureData(DM dm, Vec x, std::string var_name);
    PetscErrorCode DumpVelocityData(DM dm, Vec x, std::string var_name);
    PetscErrorCode GenerateFluidGridMatrices(SystemParameters sys);

    /// Staggered Grid Class
    class StaggeredGrid {
    public:
        StaggeredGrid(){
            // Default Constructor
        }
        StaggeredGrid(std::unordered_map<std::string,std::string>* data){
            PetscErrorCode ierr= createSystemParameters(&d_sys,data);

            if (ierr!=0){
             throw ibError("Unable To Create System Parameters for Staggered Grid");
            }
        }
        StaggeredGrid(StaggeredGrid *sgrid){
            StaggeredGrid(sgrid->d_sys);
        }
        StaggeredGrid(SystemParameters sys): d_sys(sys){
            MatCreate(PETSC_COMM_WORLD,&d_A);
            MatCreate(PETSC_COMM_WORLD,&d_B);
            MatCreate(PETSC_COMM_WORLD,&d_Bt);
            //TODO: Avoid Making These Matrices if Matrix-Free Implementation
            AssembleStokes2D();
            AssembleCellCenteredGrad2D();
            AssembleStaggeredDiv2D();
        }
        ~StaggeredGrid(){
            if(d_p) (VecDestroy(&d_p));
            if(d_rhsp) VecDestroy(&d_rhsp);
            if(d_rhsu) VecDestroy(&d_rhsu);
            if(d_A) MatDestroy(&d_A);
            if(d_B) MatDestroy(&d_B);
            if(d_Bt) MatDestroy(&d_Bt);
        }
        //Matrices
        PetscErrorCode AssembleStokes2D();
        PetscErrorCode AssembleCellCenteredGrad2D();
        PetscErrorCode AssembleStaggeredDiv2D();
        //MatrixFree
        PetscErrorCode ApplyStokes2D(Vec* x, Vec* y);
        PetscErrorCode ApplyStaggeredDiv2D(Vec*x , Vec* y);
        PetscErrorCode ApplyCellCenteredGrad2D(Vec*x , Vec* y);
        //RHS
        PetscErrorCode AssembleVelocityRHS(Vec*x , Vec* y);
        PetscErrorCode AssemblePressureRHS( Vec* y);
        PetscErrorCode AssembleFluidRHS(Vec *prhsf, Vec* prhsp);
        //getters
        Mat& getA(){
            return d_A;
        }
        Mat& getB(){
            return d_B;
        }
        Mat& getBt(){
            return d_Bt;
        }
        SystemParameters getSys(){
            return d_sys;
        }
        void resetData(SystemParameters sys){
            d_sys = sys;
            MatDestroy(&d_A);
            MatDestroy(&d_B);
            MatDestroy(&d_Bt);
            MatCreate(PETSC_COMM_WORLD,&d_A);
            MatCreate(PETSC_COMM_WORLD,&d_B);
            MatCreate(PETSC_COMM_WORLD,&d_Bt);
            AssembleStokes2D();
            AssembleStaggeredDiv2D();
            AssembleCellCenteredGrad2D();
        }
    private:
        PetscErrorCode CheckFluidSystem(SystemParameters sys);
        SystemParameters d_sys;
        Mat d_A, d_B, d_Bt;
        Vec d_u, d_p, d_rhsu, d_rhsp;
    };

    PetscErrorCode CheckFluidSystem(StaggeredGrid *grid);
    PetscErrorCode CheckRhs(Vec *prhsf, Vec* prhsp,SystemParameters sys);
    PetscErrorCode CheckMatrixFreeOps(StaggeredGrid *grid);
    PetscErrorCode AttachStokesNullspace(DM dmStokes, Mat *A_stokes);
    PetscErrorCode AttachPressureNullspace(DM dmPressure, Mat *A_stokes);
    PetscErrorCode CheckDMStag(SystemParameters sys);
    PetscErrorCode CreateSystem(DM dmSol, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys);
    PetscErrorCode CreateGradientCenterOperator(DM dm_stokes, Mat *pB, SystemParameters sys);
    PetscErrorCode CreateUnSteadyStokesSystem(DM dm_stokes, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys);
    PetscErrorCode CreateSteadyStokesSystem(DM dm_stokes, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys);
    PetscErrorCode CreateVelocityFacesPoissonSystem(DM dmSol, Mat *pA, Vec *pRhs, SystemParameters sys);
    PetscErrorCode CreateVelocityFacesUnsteadySystem(DM dmSol, Mat *pA, Vec *pRhs, SystemParameters sys);
    PetscErrorCode CreatePressureCenterPoissonSystem(DM dm_p, Mat *pA, Vec *pRhs, SystemParameters sys);
    PetscErrorCode CreatePressureCenterExactSolution(DM dm_p, Vec *pRhs, SystemParameters sys);

    PetscScalar uxRef(PetscScalar x, PetscScalar y,PetscBool is_top, PetscInt truncated_sys);
    PetscScalar uxRef(PetscScalar x, PetscScalar y,PetscInt truncated_sys);
    PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscInt truncated_sys);
    PetscScalar pRef(PetscScalar x, PetscScalar y, PetscInt truncated_sys);

    PetscScalar fx(PetscScalar x, PetscScalar y, PetscInt truncated_sys);
    PetscScalar fy(PetscScalar x, PetscScalar y, PetscInt truncated_sys);

    PetscScalar g(PetscScalar x, PetscScalar y, PetscInt truncated_sys);
    PetscScalar g(PetscScalar x, PetscScalar y, PetscInt ex, PetscInt ey, PetscInt truncated_sys);
}
#endif