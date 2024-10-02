//
// Created by gkluhana on 15/07/24.
//

#ifndef IBIMPLICIT_STOKES_H
#define IBIMPLICIT_STOKES_H
#include <staggered.h>
namespace ibImplicit{

    /**
     * @brief Stokes Solver Class
     */
    class StokesSolver{
    public:
        StokesSolver(){
            //Default Constructor
        }
        StokesSolver(SystemParameters sys) {
            d_sys = sys;
        }

        ~StokesSolver(){
        }
        /**
         * @brief Update System Parameters
         * @param sys
         * @return
         */
        PetscErrorCode updateSys(SystemParameters sys);
        /**
         * @brief Initialize the Stokes Solver
         * @return
         */
        PetscErrorCode initialize();
        /**
         * @brief Solve the Stokes System
         * @return
         */
        PetscErrorCode solve();
        /**
         * @brief Solve the Stokes System with given vectors
         * @param x : Right Hand Side
         * @param y : Solution Vector
         * @return
         */
        PetscErrorCode solve(Vec x, Vec y);
        /**
         * @brief Clean the Stokes Solver
         * @return
         */
        PetscErrorCode clean();
        PetscErrorCode linearInterpGhostCells(Vec x, Vec y);
        /**
         * @brief Apply the Upper Triangular Preconditioner
         * @param x_eul : Input Vector
         * @param f_s : External force from immersed boundary
         * @param y : Output Vector
         * @return
         */
        PetscErrorCode UpperPreconditionerApply(Vec x_eul, Vec f_s, Vec y_eul);

        PetscErrorCode solveVelocity(Vec x,Vec y);

        DM getDMStokes(){
            return dm_stokes;
        }
        SystemParameters getSys(){
            return d_sys;
        }
        Mat getA_stokes(){
            return A_s;
        }
        Mat getB_stokes(){
            return B;
        }
        Mat getBt_stokes(){
            return Bt;
        }
        Vec getRhs_stokes(){
            return rhs_stokes;
        }
        PetscInt get_its(){
            return its;
        }
        void reset_its(){
            its = 0;
        }
        void avg_its(PetscReal total, PetscReal *avg_its){
            if(total==0){
               *avg_its = 0;
            }
            else{
                *avg_its = its/total;
            }
        }

        void set_inner(){
            inner_solve = PETSC_TRUE;
        }
        void reset_inner(){
            inner_solve = PETSC_FALSE;
        }
        PetscInt get_total_its(){
            return total_its;
        }
        PetscInt get_total_inner_its(){
            return total_inner_its;
        }
        DM getDMVelocity(){
            return ctx_dm_uv;
        }
        DM getDMPressure(){
            return ctx_dm_p;
        }
    private:
        SystemParameters d_sys;
        DM dm_stokes;
        Mat A_s;
        Vec rhs_stokes;
        PetscBool pinPressure, rhs_assembled=PETSC_FALSE;
        PetscBool inner_solve=PETSC_FALSE;
        PetscReal dt,rho,mu;
        DM ctx_dm_p,ctx_dm_uv;
        Vec x_p, x_uv, y_uv, u_tilde, f_uv, f_p, Bu_tilde, y_p, y_S1_p, Bty_p, y_S1_uv, y_eul_p;
        Mat A_uv, A_p, A_c;
        Mat B, Bt;
        KSP ksp_uv,ksp_p, ksp_dummy, ksp_S1;
        PC S1_pc;
        Vec u_star,Bu_star,Btphi,phi;
        MatNullSpace matNullSpacePressure, matNullSpaceStokes;
        PetscInt its=0, inner_its = 0, total_its=0, total_inner_its = 0;
        PetscErrorCode fillBoundary(Mat *pA, Vec *pRhs);
        PetscErrorCode CreateStokesSystem(Mat *pA, Vec *pRhs);
        PetscErrorCode CreateVelocitySystem(Mat *pA, Vec *pRhs);
        PetscErrorCode CreateTruncatedUnSteadyStokesSystem_NEUMANN(Mat *pA);
        PetscErrorCode CreateTruncatedUnSteadyStokesSystem_TRUNCATED(Mat *pA);
        PetscErrorCode CreateTruncatedVelocityFacesUnsteadySystem_NEUMANN(DM dm_uv,
                                                                          Mat *pA);
        PetscErrorCode CreateTruncatedVelocityFacesUnsteadySystem_TRUNCATED(DM dm_uv,
                                                                            Mat *pA);
        PetscErrorCode CreateUnSteadyStokesSystem(Mat *pA,
                                                  Vec *pRhs);
        PetscErrorCode CreateVelocityFacesUnsteadySystem(DM dm_uv,
                                                         Mat *pA,
                                                         Vec *pRhs);
        PetscErrorCode CreatePressureCenterPoissonSystem(DM dm_p, Mat *pA, Vec *pRhs);
        PetscErrorCode AttachPressureNullspace(DM dmPressure, Mat *A);
        PetscErrorCode CreatePressureCenterUnsteadySystem(DM dm_p, Mat *pA, Vec *pRhs);
        PetscErrorCode PCApply_Projection(Vec x, Vec y);
        PetscErrorCode S1_Apply(Vec x, Vec y);
        PetscErrorCode S1_Solve_LSC(Vec x, Vec y);
        static PetscErrorCode S1_Apply_Wrapper(Mat A, Vec x, Vec y){
            PetscFunctionBeginUser;
            void *ctx;
            PetscCall(MatShellGetContext(A,&ctx));
            StokesSolver *self = static_cast<StokesSolver *>(ctx);
            self->S1_Apply(x,y);
            PetscFunctionReturn(PETSC_SUCCESS);
        }
        static PetscErrorCode PCApply_Projection_Wrapper(PC pc,Vec x, Vec y){
            PetscFunctionBeginUser;
            void* ctx;
            PetscCall(PCShellGetContext(pc,&ctx));
            StokesSolver *self = static_cast<StokesSolver*>(ctx);
            PetscCall(self->PCApply_Projection(x,y));
            PetscFunctionReturn(PETSC_SUCCESS);
        }
    };


}

#endif //IBIMPLICIT_STOKES_H
