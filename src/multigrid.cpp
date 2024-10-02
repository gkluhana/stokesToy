//
// Created by gkluhana on 29/04/24.
//
#include "../include/multigrid.h"
namespace ibImplicit{
    PetscErrorCode IsMatrixSPD(Mat A, PetscBool *isSPD) {
        PetscErrorCode ierr;
        PetscBool isSymmetric = PETSC_FALSE;
        MatFactorInfo info;

        // Initialize SPD flag to false
        *isSPD = PETSC_FALSE;

        // Check if the matrix is symmetric
        ierr = MatIsSymmetric(A, 1e-10, &isSymmetric); CHKERRQ(ierr);
        if (!isSymmetric) {
            return 0; // Matrix is not symmetric, hence not SPD
        }

        // Try Cholesky factorization
        ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);
        Mat F;
        ierr = MatGetFactor(A, MATSOLVERPETSC, MAT_FACTOR_CHOLESKY, &F); CHKERRQ(ierr);
        ierr = MatCholeskyFactorSymbolic(F, A, NULL, &info); CHKERRQ(ierr);
        ierr = MatCholeskyFactorNumeric(F, A, &info);
        if (ierr == 0) {
            *isSPD = PETSC_TRUE; // Cholesky factorization succeeded, matrix is SPD
        }
        ierr = MatDestroy(&F); CHKERRQ(ierr);

        return 0;
    }
    PetscErrorCode MultigridOnVelPoissonDM(SystemParameters sys){
        PetscFunctionBeginUser;
        /*
        This is a testing/experiments function that applies MG-CG on the velocity sub-block
         */
        DM dm_main = sys->mg_params->levels[sys->mg_params->n_levels - 1]->dm_level;
        Mat A, A_level, *A_faces=NULL;
        Vec rhs, *rhs_level=NULL;
        Vec sol_faces;
        PetscBool pinPressure;
        KSP ksp;
        Vec sol;
        MGParameters mgCtx = sys->mg_params;
        PetscInt finest_level = sys->mg_params->n_levels-1;
        const PetscInt dof0 = 0, dof1 = 1, dof2 = 0; /* 1 dof on each edge and element center */
        const PetscInt stencilWidth = 1;

        // Create DMStag
        PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, sys->Nx, sys->Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm_main));
        PetscCall(DMSetFromOptions(dm_main));
        PetscCall(DMSetUp(dm_main));
        PetscCall(DMStagSetUniformCoordinatesProduct(dm_main, sys->lx, sys->Lx, sys->ly , sys->Ly, 0.0, 0.0));
        PetscCall(DMCreateGlobalVector(dm_main,&rhs));
        pinPressure = PETSC_FALSE; // If this is set to False we need to attach a constant nullspace with the DM
        PetscMalloc1(mgCtx->n_levels,&A_faces);
        PetscMalloc1(mgCtx->n_levels,&rhs_level);
        PetscCall(CreateVelocityFacesPoissonSystem(dm_main, &A_faces[finest_level], &rhs,  sys));
        //dumpMatToFile(&A_faces[finest_level],"A_Vel_Poisson");
        //dumpVecToFile(&rhs,"rhs_Vel_Poisson");
        PetscInt startx, starty, startz, hx, hy, hz, nExtrax, nExtray, nExtraz;
        PetscInt Nx, Ny;
        DMStagGetCorners(dm_main, &startx, &starty, &startz, &hx, &hy, &hz, &nExtrax, &nExtray, &nExtraz);
        DMStagGetGlobalSizes(dm_main, &Nx, &Ny, PETSC_NULLPTR);

        // Create DMs for all levels
        sys->mg_params->levels[finest_level]->dm_faces= dm_main;
        DM *dm_coarse, *dm_fine;
        PetscInt NxC, NyC;
        for (PetscInt level=finest_level; level > 0; --level){
            dm_fine  = &sys->mg_params->levels[level]->dm_faces;
            dm_coarse= &sys->mg_params->levels[level-1]->dm_faces;
            PetscCall(DMCoarsen(*dm_fine, PETSC_COMM_WORLD, dm_coarse));
            PetscCall(DMStagSetUniformCoordinatesProduct(*dm_coarse, sys->lx, sys->Lx, sys->ly, sys->Ly, 0.0, 0.0));
            DMStagGetGlobalSizes(*dm_coarse, &NxC,&NyC, PETSC_NULLPTR);
            sys->mg_params->levels[level-1]->hx = (sys->Lx - sys->lx)/(NxC);
            sys->mg_params->levels[level-1]->hy = (sys->Ly - sys->ly)/(NyC);
        }


        //KSP Solve
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetType(ksp,KSPCG);
        KSPSetOperators(ksp,A_faces[finest_level],A_faces[finest_level]);
        KSPSetDM(ksp, dm_main);
        KSPSetDMActive(ksp,PETSC_FALSE);
        KSPSetFromOptions(ksp);
        PetscCall(KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED));
        PC pc;
        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCMG);
        PetscCall(PCMGSetLevels(pc, sys->mg_params->n_levels, NULL));
        PetscCall(PCMGSetCycleType(pc, PC_MG_CYCLE_V));
        PetscCall(PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH));

       //mg settings
       for (PetscInt level=0; level< sys->mg_params->n_levels; ++level){
           KSP ksp_level;
           PC pc_level;

            // Smoothers
            PetscCall(PCMGGetSmoother(pc,level,&ksp_level));
            PetscCall(KSPGetPC(ksp_level,&pc_level));
            //PetscCall(KSPSetOperators(ksp_level,A_faces[level],A_faces[level]));
            //dumpMatToFile(&A_faces[level],"A_faces_"+ std::to_string(level));

            if (level>0) {
                PetscCall(KSPSetType(ksp_level, KSPRICHARDSON));
                PetscCall(PCSetType(pc_level, PCSOR));
                PetscCall(PCSORSetOmega(pc_level,2.0/3.0));
            }

            // Transfer Operators
            if(level>0){
                Mat restriction, interpolation;
                DM dm_level = sys->mg_params->levels[level]->dm_faces;
                DM dm_coarser = sys->mg_params->levels[level-1] ->dm_faces;
                PetscCall(DMCreateInterpolation(dm_coarser, dm_level, &interpolation, NULL));
                PetscCall(PCMGSetInterpolation(pc,level,interpolation));
                if(debug) dumpMatToFile(interpolation,"interp_"+ std::to_string(level));
                PetscCall(MatDestroy(&interpolation));
                PetscCall(DMCreateRestriction(dm_coarser, dm_level, &restriction));
                PetscCall(PCMGSetRestriction(pc, level, restriction));
                if(debug) dumpMatToFile(restriction,"restrict_"+ std::to_string(level));
                PetscCall(MatDestroy(&restriction));
            }
       }
        // Coarse Grid Solver
        KSP ksp_coarse;
        PC pc_coarse;
        //PCMGGetCoarseSolve(pc,&ksp_coarse);
        //KSPSetType(ksp_coarse,KSPCG);
        //KSPGetPC(ksp_coarse,&pc_coarse);
        //PCSetType(pc_coarse,PCILU);

        //Solve
        KSPSetUp(ksp);
        PetscCall(DMCreateGlobalVector(dm_main,&sol_faces));
        PetscCall(KSPSolve(ksp,rhs,sol_faces));
        DumpVelocityData(dm_main, sol_faces, "u");
        if(debug){
            KSPConvergedReason reason;
            PetscCall(KSPGetConvergedReason(ksp, &reason));
            PetscCall(PetscInfo(reason,"KSP Converged Reason: %d",reason));
            PetscCheck(reason >= 0, PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED, "Linear solve failed on Velocity");
        }
        KSPDestroy(&ksp);
        VecDestroy(&sol_faces);
        DMDestroy(&dm_main);
        PetscFunctionReturn(0);

    }
    PetscErrorCode MultigridOnPressurePoissonDM(SystemParameters sys){
        PetscFunctionBeginUser;
        /*
        This is a testing/experiments function that applies MG-CG on the pressure sub-block
         */
        DM dm_p = sys->mg_params->levels[sys->mg_params->n_levels - 1]->dm_level;
        Mat A_p;
        Vec rhs, *rhs_level=NULL;
        KSP ksp;
        Vec sol_center;
        MGParameters mgCtx = sys->mg_params;
        PetscInt finest_level = sys->mg_params->n_levels-1;
        const PetscInt dof0 = 0, dof1 = 0, dof2 =1; /* 1 dof on each edge and element center */
        const PetscInt stencilWidth = 1;

        // Create DMStag
        PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, sys->Nx, sys->Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm_p));
        PetscCall(DMSetFromOptions(dm_p));
        PetscCall(DMSetUp(dm_p));
        PetscCall(DMStagSetUniformCoordinatesProduct(dm_p, sys->lx, sys->Lx, sys->ly , sys->Ly, 0.0, 0.0));
        PetscCall(DMCreateGlobalVector(dm_p, &rhs));
        sys->mg_params->levels[finest_level]->hx = sys->hx;
        sys->mg_params->levels[finest_level]->hy = sys->hy;
        PetscCall(CreatePressureCenterPoissonSystem(dm_p, &A_p,&rhs, sys));
        MatNullSpace nullsp;
        Vec constantPressure;
        PetscReal nrm;
        PetscCall(DMGetGlobalVector(dm_p, &constantPressure));
        PetscCall(VecSet(constantPressure, 1.0));
        PetscCall(VecNorm(constantPressure, NORM_2, &nrm));
        PetscCall(VecScale(constantPressure, 1.0 / nrm));
        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,1,&constantPressure,&nullsp);
        MatSetNullSpace(A_p,nullsp);
        MatNullSpaceDestroy(&nullsp);
        dumpMatToFile(A_p, "A_P_Poisson");
        dumpVecToFile(rhs,"rhs_P_Poisson");
        PetscInt startx, starty, startz, hx, hy, hz, nExtrax, nExtray, nExtraz;
        PetscInt Nx, Ny;
        DMStagGetCorners(dm_p, &startx, &starty, &startz, &hx, &hy, &hz, &nExtrax, &nExtray, &nExtraz);
        DMStagGetGlobalSizes(dm_p, &Nx, &Ny, PETSC_NULLPTR);

        // Create DMs for all levels
        sys->mg_params->levels[finest_level]->dm_center= dm_p;
        DM *dm_coarse, *dm_fine;
        PetscInt NxC, NyC;
        for (PetscInt level=finest_level; level > 0; --level){
            dm_fine  = &sys->mg_params->levels[level]->dm_center;
            dm_coarse= &sys->mg_params->levels[level-1]->dm_center;
            PetscCall(DMCoarsen(*dm_fine, PETSC_COMM_WORLD, dm_coarse));
            PetscCall(DMStagSetUniformCoordinatesProduct(*dm_coarse, sys->lx, sys->Lx, sys->ly, sys->Ly, 0.0, 0.0));
            DMStagGetGlobalSizes(*dm_coarse, &NxC,&NyC, PETSC_NULLPTR);
            sys->mg_params->levels[level-1]->hx = (sys->Lx - sys->lx)/(NxC);
            sys->mg_params->levels[level-1]->hy = (sys->Ly - sys->ly)/(NyC);
        }


        //KSP Solve
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetType(ksp,KSPCG);
        PetscCall(KSPSetOperators(ksp, A_p,A_p));
        KSPSetDM(ksp, dm_p);
        KSPSetDMActive(ksp,PETSC_FALSE);
        KSPSetOptionsPrefix(ksp,"pre_mgcg_");
        KSPSetFromOptions(ksp);
        PC pc;
        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCMG);
        PetscCall(PCMGSetLevels(pc, sys->mg_params->n_levels, NULL));
        PetscCall(PCMGSetCycleType(pc, PC_MG_CYCLE_V));
        PetscCall(PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH));

        //mg settings
        for (PetscInt level=0; level< sys->mg_params->n_levels; ++level){
            KSP ksp_level;
            PC pc_level;

            // Smoothers
            PetscCall(PCMGGetSmoother(pc,level,&ksp_level));
            PetscCall(KSPGetPC(ksp_level,&pc_level));

            if (level>0) {
                PetscCall(KSPSetType(ksp_level, KSPRICHARDSON));
                PetscCall(PCSetType(pc_level, PCSOR));
                PetscCall(PCSORSetOmega(pc_level,1.85));
            }

            // Transfer Operators
            if(level>0){
                Mat restriction, interpolation;
                DM dm_level = sys->mg_params->levels[level]->dm_center;
                DM dm_coarser = sys->mg_params->levels[level-1] ->dm_center;
                PetscCall(DMCreateInterpolation(dm_coarser, dm_level, &interpolation, NULL));
                PetscCall(PCMGSetInterpolation(pc,level,interpolation));
                //dumpMatToFile(interpolation,"interp_"+ std::to_string(level));
                PetscCall(MatDestroy(&interpolation));
                PetscCall(DMCreateRestriction(dm_coarser, dm_level, &restriction));
                PetscCall(PCMGSetRestriction(pc, level, restriction));
                //dumpMatToFile(restriction,"restrict_"+ std::to_string(level));
                PetscCall(MatDestroy(&restriction));
            }
        }
        // Coarse Grid Solver
        KSP ksp_coarse;
        PC pc_coarse;

        //Solve
        KSPSetUp(ksp);
        PetscCall(DMCreateGlobalVector(dm_p, &sol_center));
        if (debug) DumpPressureData(dm_p, sol_center, "p_0");
        KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
        PetscCall(KSPSolve(ksp,rhs,sol_center));
        if (debug) DumpPressureData(dm_p, rhs, "p_rhs");
        DumpPressureData(dm_p, sol_center, "p_sol");
        if (debug){
            KSPConvergedReason reason;
            PetscCall(KSPGetConvergedReason(ksp, &reason));
            //PetscCall(PetscInfo(reason,"KSP Converged Reason: %d",reason));
            PetscCheck(reason >= 0, PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED, "Linear solve failed on Pressure");
        }
        PetscFunctionReturn(0);

    }
    PetscErrorCode InterpolateVelocityU(Vec* uc, Vec *uf, PetscInt nxc, PetscInt nyc, PetscInt nxf, PetscInt nyf){
        // This function was written before using DMStag -- probably won't be used now.
        // Add extra zero row above and below ? Add conditional for first and last row
        // No need, error on dirichlet is zero - no prolongation!
        // entries 0,1,2,... fill entries 0, 2, 4,...
        // So we go from j=0 to Nx
        // and fill entries 2*j of lower and upper
        // Then we use 1/2 * alternate entries to also fill center values
        // Each row has Nx+1 elements
        PetscInt row_size_c = nxc+1;
        PetscInt row_size_f = nxf+1;
        PetscInt num_rows_c = nyc;
        PetscInt row_ip1,lower_row_f,upper_row_f;

        PetscInt row_i_indices[row_size_c], row_ip1_indices[row_size_c];
        PetscInt lower_row_indices[row_size_c], upper_row_indices[row_size_c];
        PetscInt lower_bdry_indices[row_size_c], upper_bdry_indices[row_size_c];
        PetscInt lower_middle_indices[row_size_c - 1], upper_middle_indices[row_size_c - 1];
        PetscInt lower_bdry_middle_indices[row_size_c], upper_bdry_middle_indices[row_size_c];

        PetscScalar row_i_vals[row_size_c], row_ip1_vals[row_size_c];
        PetscScalar lower_row_vals[row_size_c], upper_row_vals[row_size_c];
        PetscScalar lower_bdry_vals[row_size_c], upper_bdry_vals[row_size_c];
        PetscScalar lower_middle_vals[row_size_c-1], upper_middle_vals[row_size_c-2];
        PetscScalar lower_bdry_middle_vals[row_size_c], upper_bdry_middle_vals[row_size_c];
        // Go over each Row for rows 0 to Ny-1,
        for (PetscInt row_i=0; row_i < num_rows_c-1; row_i++){
            // Grab the next Row
            row_ip1 = row_i + 1;
            // Corresponding Rows of uf to fill
            // When you have row i and i+1 for coarse, fill 2*(i+1) and 2*(i+1)-1 for fine
            upper_row_f = 2*(row_ip1);
            lower_row_f = upper_row_f-1;
            for (PetscInt j=0; j< row_size_c; j++){
                row_i_indices[j] = j + row_i*row_size_c;
                row_ip1_indices[j] = j + row_ip1 * row_size_c;
                lower_row_indices[j] = lower_row_f*row_size_f + j*2;
                upper_row_indices[j] = upper_row_f*row_size_f + j*2;
                if ( j != row_size_c-1){
                    lower_middle_indices[j]= lower_row_indices[j] + 1;
                    upper_middle_indices[j]= upper_row_indices[j] + 1;
                }
                if (row_i == 0){
                    lower_bdry_indices[j] = j*2;
                    if ( j != row_size_c-1) {
                        lower_bdry_middle_indices[j] = lower_bdry_indices[j] + 1;
                    }
                }
                if (row_ip1 == num_rows_c-1){
                    upper_bdry_indices[j] = upper_row_indices[j] + row_size_f;
                    if ( j != row_size_c-1) {
                        upper_bdry_middle_indices[j] = upper_bdry_indices[j] + 1;
                    }
                }
            }
            VecGetValues(*uc,row_size_c,row_i_indices,row_i_vals);
            VecGetValues(*uc,row_size_c,row_ip1_indices,row_ip1_vals);

            for (PetscInt j=0; j< row_size_c; j++){
                // In each row 3/4 * ith + 1/4* (i+1) fills lower of fine
                // In each row 1/4 * ith + 2/4* (i+1) fills upper of fine
                if(j!=0 && j!= row_size_c-1){
                    lower_row_vals[j] = 3.0/4.0 * row_i_vals[j] + 1.0/4.0*row_ip1_vals[j];
                    upper_row_vals[j] = 1.0/4.0 * row_i_vals[j] + 3.0/4.0*row_ip1_vals[j];
                }

                if(j != 0){
                    lower_middle_vals[j-1] = 1.0/2.0 * (lower_row_vals[j] + lower_row_vals[j-1]) ;
                    upper_middle_vals[j-1] = 1.0/2.0 * (upper_row_vals[j] + upper_row_vals[j-1]);
                }
                if (row_i == 0){
                    lower_bdry_vals[j] = 0.5 * lower_row_vals[j];
                    if(j != 0){
                        lower_bdry_middle_vals[j-1] = 1.0/2.0 * (lower_bdry_vals[j] + lower_bdry_vals[j-1]) ;
                    }
                }
                if (row_ip1 == num_rows_c-1){
                    upper_bdry_vals[j] = 0.5* upper_row_vals[j] ;
                    if(j != 0){
                        upper_bdry_middle_vals[j-1] = 1.0/2.0 * (upper_bdry_vals[j] + upper_bdry_vals[j-1]) ;
                    }
                }
            }
            VecSetValues(*uf, row_size_c,lower_row_indices,lower_row_vals,INSERT_VALUES);
            VecSetValues(*uf, row_size_c,upper_row_indices,upper_row_vals,INSERT_VALUES);
            VecSetValues(*uf, row_size_c-1,lower_middle_indices,lower_middle_vals,INSERT_VALUES);
            VecSetValues(*uf, row_size_c-1,upper_middle_indices,upper_middle_vals,INSERT_VALUES);
        }
        VecSetValues(*uf, row_size_c,lower_bdry_indices,lower_bdry_vals,INSERT_VALUES);
        VecSetValues(*uf, row_size_c,upper_bdry_indices,upper_bdry_vals,INSERT_VALUES);
        VecSetValues(*uf, row_size_c-1,lower_bdry_middle_indices,lower_bdry_vals,INSERT_VALUES);
        VecSetValues(*uf, row_size_c-1,upper_bdry_middle_indices,upper_bdry_vals,INSERT_VALUES);
        PetscFunctionReturn(0);
    }
    PetscErrorCode InterpolateVelocityV(Vec* vc, Vec *vf, PetscInt nxc, PetscInt nyc, PetscInt nxf, PetscInt nyf){
        // This function was written before using DMStag -- probably won't be used now.
        //TODO: Too many conditionals?
        PetscInt row_size_c = nxc;
        PetscInt row_size_f = nxf;
        PetscInt num_rows_c = nyc+1;
        PetscInt col_ip1,right_col_f,left_col_f;

        PetscInt col_i_indices[num_rows_c], col_ip1_indices[num_rows_c];
        PetscInt right_col_indices[num_rows_c], left_col_indices[num_rows_c];
        PetscInt right_middle_indices[num_rows_c - 1], left_middle_indices[num_rows_c - 1];

        PetscScalar col_i_vals[num_rows_c], col_ip1_vals[num_rows_c];
        PetscScalar right_col_vals[num_rows_c], left_col_vals[num_rows_c];
        PetscScalar right_middle_vals[num_rows_c-1], left_middle_vals[num_rows_c-2];

        for (PetscInt col_i=0; col_i < row_size_c-1; col_i++){
            // Grab the next col
            col_ip1 = col_i + 1;
            // Corresponding cols of vf to fill
            left_col_f = 2*(col_ip1);
            right_col_f = left_col_f-1;
            // j is going up rows
            for (PetscInt j=0; j< num_rows_c; j++){
                col_i_indices[j] = col_i + j*row_size_c;
                col_ip1_indices[j] = col_ip1 + j* row_size_c;
                left_col_indices[j] = left_col_f + j*2*row_size_f ;
                right_col_indices[j]= right_col_f + j*2*row_size_f ;
                if ( j != num_rows_c-1){
                    right_middle_indices[j]= right_col_indices[j] + row_size_f;
                    left_middle_indices[j]= left_col_indices[j] + row_size_f;
                }
            }
            VecGetValues(*vc,num_rows_c,col_i_indices,  col_i_vals);
            VecGetValues(*vc,num_rows_c,col_ip1_indices,col_ip1_vals);

            for (PetscInt j=0; j< num_rows_c; j++){
                // In each row 3/4 * ith + 1/4* (i+1) fills lower of fine
                // In each row 1/4 * ith + 2/4* (i+1) fills upper of fine
                if(j!=0 && j!=num_rows_c-1){
                    right_col_vals[j] = 3.0/4.0 * col_i_vals[j] + 1.0/4.0*col_ip1_vals[j];
                    left_col_vals[j] = 1.0/4.0 * col_i_vals[j] + 3.0/4.0*col_ip1_vals[j];
                }

                if(j != 0){
                    right_middle_vals[j-1] = 1.0/2.0 * (right_col_vals[j] + right_col_vals[j-1]) ;
                    left_middle_vals[j-1] = 1.0/2.0 * (left_col_vals[j] + left_col_vals[j-1]);
                }
            }
            VecSetValues(*vf, num_rows_c,right_col_indices,right_col_vals,INSERT_VALUES);
            VecSetValues(*vf, num_rows_c,left_col_indices,left_col_vals,INSERT_VALUES);
            VecSetValues(*vf, num_rows_c-1,right_middle_indices,right_middle_vals,INSERT_VALUES);
            VecSetValues(*vf, num_rows_c-1,left_middle_indices,left_middle_vals,INSERT_VALUES);
        }
        PetscFunctionReturn(0);
    }
}
