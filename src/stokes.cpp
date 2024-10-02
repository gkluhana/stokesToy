//
// Created by gkluhana on 16/07/24.
//
#include <stokes.h>
namespace ibImplicit{
    PetscErrorCode StokesSolver::updateSys(SystemParameters sys){
        PetscFunctionBeginUser;
        d_sys = sys;
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::linearInterpGhostCells(Vec x, Vec y){
        PetscFunctionBeginUser;
        Vec             stokesLocal;
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        PetscInt count;
        PetscFunctionBeginUser;
        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        DMBoundaryType xtype,ytype,ztype;
        PetscCall(DMStagGetBoundaryTypes(dm_stokes,&xtype,&ytype,&ztype));
        PetscCall(DMCreateGlobalVector(dm_stokes, &y));
        PetscCall(DMGetLocalVector(dm_stokes, &stokesLocal));
        PetscCall(DMGlobalToLocal(dm_stokes, x, INSERT_VALUES, stokesLocal));


        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);
                DMStagStencil from[2], ghost;
                PetscScalar valFrom[2], valGhost;
                if(right_boundary){
                    from[0].i = ex;
                    from[0].j = ey;
                    from[0].loc = DMSTAG_RIGHT;
                    from[0].c = 0;
                    from[1].i = ex;
                    from[1].j = ey;
                    from[1].loc = DMSTAG_LEFT;
                    from[1].c = 0;
                    PetscCall(DMStagVecGetValuesStencil(dm_stokes,stokesLocal,2,from,valFrom));
                    ghost.i = ex+1;
                    ghost.j = ey;
                    ghost.loc = DMSTAG_RIGHT;
                    ghost.c = 0;
                    valGhost = 2*valFrom[0] - valFrom[1];
                    PetscCall(DMStagVecSetValuesStencil(dm_stokes, y, 1, &ghost, &valGhost, INSERT_VALUES));
                }
                else if(left_boundary){

                }
                else if(top_boundary){

                }
                else if(bottom_boundary){

                }
            }
        }
        VecAssemblyBegin(y);
        VecAssemblyEnd(y);
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::initialize(){
        PetscFunctionBeginUser;
        const PetscInt dof0 = 0, dof1 = 1, dof2 = 1; /* 1 dof on each edge and element center */
        const PetscInt stencilWidth = 1;
        PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, d_sys->Nx, d_sys->Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm_stokes));
        PetscCall(DMSetFromOptions(dm_stokes));
        PetscCall(DMSetUp(dm_stokes));
        PetscCall(DMStagSetUniformCoordinatesProduct(dm_stokes, d_sys->lx, d_sys->Lx, d_sys->ly, d_sys->Ly, 0.0, 0.0));
        pinPressure = PETSC_FALSE; // If this is set to False we need to attach a constant nullspace with the DM
        d_sys->mg_params->faces_only=PETSC_FALSE;
        PetscCall(this->CreateStokesSystem(&A_s, &rhs_stokes));
        PetscCall(AttachStokesNullspace(dm_stokes, &A_s));
        //Setup data for preconditioner
        dt = d_sys->dt;
        rho = d_sys->rho;
        mu = d_sys->mu;

        PetscCall(DMStagCreateCompatibleDMStag(dm_stokes,0,1,0,0,&ctx_dm_uv));
        PetscCall(DMStagCreateCompatibleDMStag(dm_stokes,0,0,1,0,&ctx_dm_p));
        PetscCall(DMStagSetUniformCoordinatesProduct(ctx_dm_uv, d_sys->lx, d_sys->Lx, d_sys->ly , d_sys->Ly, 0.0, 0.0));
        PetscCall(DMStagSetUniformCoordinatesProduct(ctx_dm_p, d_sys->lx, d_sys->Lx, d_sys->ly , d_sys->Ly, 0.0, 0.0));

        // Set up the work vectors needed for pc
        PetscCall(DMCreateGlobalVector(dm_stokes, &y_eul_p));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &u_star));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &Btphi));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &Bty_p));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &x_uv));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &y_S1_uv));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &f_uv));
        PetscCall(DMCreateGlobalVector(ctx_dm_uv, &u_tilde));

        PetscCall(DMCreateGlobalVector(ctx_dm_p, &Bu_star));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &Bu_tilde));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &phi));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &x_p));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &y_p));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &f_p));
        PetscCall(DMCreateGlobalVector(ctx_dm_p, &y_S1_p));

        PetscCall(DMCreateGlobalVector(dm_stokes,&y_uv));

        // Get IS from a dummy fieldsplit KSP -- should find a more efficient alternative
        PC pc_dummy;
        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp_dummy));
        PetscCall(KSPSetOperators(ksp_dummy, A_s, A_s));
        PetscCall(KSPSetDM(ksp_dummy,dm_stokes));
        PetscCall(KSPSetDMActive(ksp_dummy,PETSC_FALSE));
        PetscCall(KSPGetPC(ksp_dummy,&pc_dummy));
        PetscCall(KSPSetFromOptions(ksp_dummy));
        PetscCall(PCSetType(pc_dummy,PCFIELDSPLIT));
        PetscCall(PCFieldSplitSetType(pc_dummy,PC_COMPOSITE_SCHUR));
        PetscCall(KSPSetUp(ksp_dummy));


        // Get Index Sets for velocity and pressure fields
        IS is_vel, is_p;
        PetscCall(PCFieldSplitGetISByIndex(pc_dummy,0,&is_vel));
        PetscCall(PCFieldSplitGetISByIndex(pc_dummy,1,&is_p));
        // Get the B and Bt operators
        PetscCall(MatCreateSubMatrix(A_s,is_p,is_vel,MAT_INITIAL_MATRIX,&B));
        PetscCall(MatCreateSubMatrix(A_s,is_vel,is_p,MAT_INITIAL_MATRIX,&Bt));
        if (debug) dumpMatToFile(B,"B_dm");
        if (debug) dumpMatToFile(Bt,"Bt_dm");

        // Generate the matrix operators for pressure and velocity subfields
        PetscCall(this->CreateVelocitySystem(&A_uv,nullptr));

        // Set up MG-CG KSP for Velocity
        // Create DMs for all levels -- only needed for velocity sub solves with MG...
        PetscInt finest_level = d_sys->mg_params->n_levels-1;
        d_sys->mg_params->levels[finest_level]->dm_faces= ctx_dm_uv;
        DM *dm_coarse, *dm_fine;
        PetscInt NxC, NyC;
        for (PetscInt level=finest_level; level > 0; --level){
            dm_fine  = &d_sys->mg_params->levels[level]->dm_faces;
            dm_coarse= &d_sys->mg_params->levels[level-1]->dm_faces;
            PetscCall(DMCoarsen(*dm_fine, PETSC_COMM_WORLD, dm_coarse));
            PetscCall(DMStagSetUniformCoordinatesProduct(*dm_coarse, d_sys->lx, d_sys->Lx, d_sys->ly, d_sys->Ly, 0.0, 0.0));
            PetscCall(DMStagGetGlobalSizes(*dm_coarse, &NxC,&NyC, PETSC_NULLPTR));
            d_sys->mg_params->levels[level-1]->hx = (d_sys->Lx - d_sys->lx)/(NxC);
            d_sys->mg_params->levels[level-1]->hy = (d_sys->Ly - d_sys->ly)/(NyC);
        }

        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp_uv));
        PetscCall(KSPSetType(ksp_uv, KSPCG));
        PetscCall(KSPSetOperators(ksp_uv, A_uv, A_uv));
        PetscCall(KSPSetDM(ksp_uv, ctx_dm_uv));
        PetscCall(KSPSetDMActive(ksp_uv, PETSC_FALSE));
        PetscCall(KSPSetOptionsPrefix(ksp_uv, "vel_mgcg_"));
        PetscCall(KSPSetFromOptions(ksp_uv));
        PetscCall(KSPSetNormType(ksp_uv, KSP_NORM_NATURAL));

        PC pc_uv;
        PetscCall(KSPGetPC(ksp_uv, &pc_uv));
        PetscCall(PCSetType(pc_uv,PCMG));
        PetscCall(PCMGSetLevels(pc_uv, d_sys->mg_params->n_levels, NULL));
        PetscCall(PCMGSetCycleType(pc_uv, PC_MG_CYCLE_V));
        PetscCall(PCMGSetGalerkin(pc_uv,PC_MG_GALERKIN_BOTH));

        //mg settings
        for (PetscInt level=0; level< d_sys->mg_params->n_levels; ++level){
            KSP ksp_level;
            PC pc_level;

            // Smoothers
            PetscCall(PCMGGetSmoother(pc_uv,level,&ksp_level));
            PetscCall(KSPGetPC(ksp_level,&pc_level));

            if (level>0) {
                PetscCall(KSPSetType(ksp_level, KSPRICHARDSON));
                PetscCall(PCSetType(pc_level, PCSOR));
                PetscCall(PCSORSetOmega(pc_level,2.0/3.0));
            }

            // Transfer Operators
            if(level>0){
                Mat restriction, interpolation;
                DM dm_level = d_sys->mg_params->levels[level]->dm_faces;
                DM dm_coarser = d_sys->mg_params->levels[level-1] ->dm_faces;
                PetscCall(DMCreateInterpolation(dm_coarser, dm_level, &interpolation, NULL));
                PetscCall(PCMGSetInterpolation(pc_uv,level,interpolation));
                //dumpMatToFile(&interpolation,"interp_"+ std::to_string(level));
                PetscCall(MatDestroy(&interpolation));
                PetscCall(DMCreateRestriction(dm_coarser, dm_level, &restriction));
                PetscCall(PCMGSetRestriction(pc_uv, level, restriction));
                //dumpMatToFile(&restriction,"restrict_"+ std::to_string(level));
                PetscCall(MatDestroy(&restriction));
                PetscCall(DMDestroy(&dm_coarser));
            }
        }
        PetscCall(KSPSetUp(ksp_uv));

        // Set up MG-CG for Pressure
        // Create DMs for all levels
        PetscCall(this->CreatePressureCenterPoissonSystem(ctx_dm_p, &A_p, NULL));
        PetscCall(this->CreatePressureCenterUnsteadySystem(ctx_dm_p, &A_c, NULL));
        PetscCall(AttachPressureNullspace(ctx_dm_p, &A_p));
        d_sys->mg_params->levels[finest_level]->dm_center = ctx_dm_p;
        for (PetscInt level=finest_level; level > 0; --level){
            dm_fine  = &d_sys->mg_params->levels[level]->dm_center;
            dm_coarse= &d_sys->mg_params->levels[level-1]->dm_center;
            PetscCall(DMCoarsen(*dm_fine, PETSC_COMM_WORLD, dm_coarse));
            PetscCall(DMStagSetUniformCoordinatesProduct(*dm_coarse, d_sys->lx, d_sys->Lx, d_sys->ly, d_sys->Ly, 0.0, 0.0));
            PetscCall(DMStagGetGlobalSizes(*dm_coarse, &NxC,&NyC, PETSC_NULLPTR));
            d_sys->mg_params->levels[level-1]->hx = (d_sys->Lx - d_sys->lx)/(NxC);
            d_sys->mg_params->levels[level-1]->hy = (d_sys->Ly - d_sys->ly)/(NyC);
        }
        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp_p));
        PetscCall(KSPSetType(ksp_p,KSPCG));
        PetscCall(KSPSetOperators(ksp_p, A_p, A_p));
        if (debug) dumpMatToFile(A_p,"A_p");

        PetscCall(KSPSetDM(ksp_p, ctx_dm_p));
        PetscCall(KSPSetDMActive(ksp_p,PETSC_FALSE));
        PetscCall(KSPSetOptionsPrefix(ksp_p, "pre_mgcg_"));
        PetscCall(KSPSetFromOptions(ksp_p));
        PC pc_p;
	PetscCall(KSPGetPC(ksp_p,&pc_p));
        PetscCall(PCSetType(pc_p,PCMG));
        PetscCall(PCMGSetLevels(pc_p, d_sys->mg_params->n_levels, NULL));
        PetscCall(PCMGSetCycleType(pc_p, PC_MG_CYCLE_V));
        PetscCall(PCMGSetGalerkin(pc_p,PC_MG_GALERKIN_BOTH));

        //mg settings
        for (PetscInt level=0; level< d_sys->mg_params->n_levels; ++level){
            KSP ksp_level;
            PC pc_level;

            // Smoothers
            PetscCall(PCMGGetSmoother(pc_p,level,&ksp_level));
            PetscCall(KSPGetPC(ksp_level,&pc_level));

            if (level>0) {
                PetscCall(KSPSetType(ksp_level, KSPRICHARDSON));
                PetscCall(PCSetType(pc_level, PCSOR));
                PetscCall(PCSORSetOmega(pc_level,1.85));
            }

            // Transfer Operators
            if(level>0){
                Mat restriction, interpolation;
                DM dm_level = d_sys->mg_params->levels[level]->dm_center;
                DM dm_coarser = d_sys->mg_params->levels[level-1] ->dm_center;
                PetscCall(DMCreateInterpolation(dm_coarser, dm_level, &interpolation, NULL));
                PetscCall(PCMGSetInterpolation(pc_p,level,interpolation));
                if(debug) dumpMatToFile(interpolation,"interp_"+ std::to_string(level));
                PetscCall(MatDestroy(&interpolation));
                PetscCall(DMCreateRestriction(dm_coarser, dm_level, &restriction));
                PetscCall(dumpMatToFile(restriction,"restrict_"+ std::to_string(level)));
                PetscCall(PCMGSetRestriction(pc_p, level, restriction));
                if(debug) dumpMatToFile(restriction,"restrict_"+ std::to_string(level));
                PetscCall(MatDestroy(&restriction));
                PetscCall(DMDestroy(&dm_coarser));
            }
        }
        // Create the S1 matrix shell
        Mat S1;
        PetscInt n_local;
        PetscCall(VecGetLocalSize(x_p, &n_local));
        PetscCall(MatCreateShell(PETSC_COMM_WORLD,n_local,n_local, PETSC_DETERMINE,PETSC_DETERMINE,this,&S1));
        PetscCall(MatShellSetOperation(S1, MATOP_MULT, reinterpret_cast<void (*)(void)>(S1_Apply_Wrapper)));

        //Create the S1 solver
        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp_S1));
        PetscCall(KSPSetOperators(ksp_S1, S1, S1));
        PetscCall(KSPSetOptionsPrefix(ksp_S1, "ib_S1_"));
        PetscCall(KSPSetReusePreconditioner(ksp_S1, PETSC_FALSE));
        PetscCall(KSPGetPC(ksp_S1, &S1_pc));
        PetscCall(PCSetType(S1_pc,PCNONE));
        PetscCall(KSPSetFromOptions(ksp_S1));
        PetscCall(KSPSetUp(ksp_S1));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::clean(){
        PetscFunctionBeginUser;
        PetscCall(MatDestroy(&A_s));
        PetscCall(MatDestroy(&A_p));
        PetscCall(MatDestroy(&A_uv));
        PetscCall(MatDestroy(&B));
        PetscCall(MatDestroy(&Bt));
        PetscCall(DMDestroy(&dm_stokes));
        PetscCall(DMDestroy(&ctx_dm_p));
        PetscCall(DMDestroy(&ctx_dm_uv));
        PetscCall(KSPDestroy(&ksp_uv));
        PetscCall(KSPDestroy(&ksp_p));
        PetscCall(KSPDestroy(&ksp_dummy));

        if(rhs_assembled){PetscCall(VecDestroy(&rhs_stokes));}
        PetscCall(VecDestroy(&u_star));
        PetscCall(VecDestroy(&Bu_star));
        PetscCall(VecDestroy(&Btphi));
        PetscCall(VecDestroy(&phi));
        PetscCall(VecDestroy(&x_uv));
        PetscCall(VecDestroy(&x_p));
        PetscCall(VecDestroy(&y_uv));
        PetscCall(MatNullSpaceDestroy(&matNullSpacePressure));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::solve(Vec x, Vec y){
        PetscFunctionBeginUser;

        // Create KSP
        KSP ksp_stokes;
        PC pc_stokes;
        PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp_stokes));
        PetscCall(KSPSetType(ksp_stokes,KSPFGMRES));
        PetscCall(KSPSetOperators(ksp_stokes, A_s, A_s));
        //PetscCall(KSPSetDM(ksp_stokes, dm_stokes));
        //PetscCall(KSPSetDMActive(ksp_stokes,PETSC_FALSE));
        std::string prefix;
        if (inner_solve){
            prefix = "stokes_inner_";
        }
        else{
            prefix = "stokes_";
        }
        PetscCall(KSPSetOptionsPrefix(ksp_stokes,prefix.c_str()));
        PetscCall(KSPSetFromOptions(ksp_stokes));

        // Create PC Shell
        PetscCall(KSPGetPC(ksp_stokes,&pc_stokes));
        PetscCall(PCSetType(pc_stokes,PCSHELL));

        // Set PC Shell Apply, Context, Setup,
        PetscCall(PCShellSetContext(pc_stokes,this));
        PetscCall(PCShellSetApply(pc_stokes,PCApply_Projection_Wrapper));
        KSPSetUp(ksp_stokes);

        KSPSolve(ksp_stokes,x,y);
        PetscInt its_stokes;
        PetscCall(KSPGetIterationNumber(ksp_stokes,&its_stokes));
        its += its_stokes;
        if (inner_solve){
            inner_its += its_stokes;
            total_inner_its += its_stokes;
            //if(debug) std::cout << "total_inner_its: " << total_inner_its << std::endl;
        }
        else{
            total_its += its_stokes;
            //if(debug) std::cout << "total_its: " << total_its << std::endl;
        }

        // Print Solution
        if(debug) DumpStokesData(dm_stokes, y, "stokes");
        if(debug) dumpVecToFile(y,"dm_stokes");

        // Destroy Objects
        PetscCall(KSPDestroy(&ksp_stokes));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::solve(){
        PetscFunctionBeginUser;
        Vec sol_stokes;
        DMCreateGlobalVector(dm_stokes, &sol_stokes);

        this->solve(rhs_stokes,sol_stokes);
        VecDestroy(&sol_stokes);
        PetscFunctionReturn(0);
    }

    PetscErrorCode StokesSolver::PCApply_Projection(Vec x, Vec y){
        PetscFunctionBeginUser;
        // Get required context and work vectors

        // Get the vectors for u,v,p
        PetscCall(DMStagMigrateVec(dm_stokes, x, ctx_dm_p, x_p));
        PetscCall(DMStagMigrateVec(dm_stokes, x, ctx_dm_uv, x_uv));


        // Solve A_vel * u_star = (fu,fv) with mg-cg
        KSPSolve(ksp_uv,x_uv,u_star);

        // Solve -L_p * phi =   rho/dt * (B*u_star - fp) , where B = -D dot
        // Ap * phi = rho/dt * (Bu_star - fp)

        MatMult(B, u_star, Bu_star);
        VecAXPY(Bu_star,-1.0,x_p);
        VecScale(Bu_star,rho/dt);

        // solve with mgcg to get phi = Ap_inv(rho/dt (Bu_star - fp))
        KSPSolve(ksp_p,Bu_star,phi);

        // u = u_star - dt / rho * grad (phi)
        MatMult(Bt,phi,Btphi);
        VecAXPY(u_star,-dt/rho,Btphi);

        // p = phi + dt/rho * mu * (rho/dt * (Bu* - x_p))
        // Bu_star is rho/dt * (Bu* - x_p)
        VecAXPY(phi,mu*dt/rho ,Bu_star);

        // Migrate the vectors to y
        PetscCall(DMStagMigrateVec(ctx_dm_p, phi, dm_stokes, y));
        PetscCall(DMStagMigrateVec(ctx_dm_uv, u_star, dm_stokes, y_uv));
        PetscCall(VecAXPY(y,1.0,y_uv));

        //Destroy temporary vectors
        PetscFunctionReturn(0);
    }
    PetscErrorCode PCApply_None(PC pc, Vec x, Vec y){
        PetscFunctionBeginUser;
        VecCopy(x,y);
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateUnSteadyStokesSystem(Mat *pA, Vec *pRhs)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt level = d_sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;

        PetscCall(DMCreateMatrix(dm_stokes, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        PetscCheck(!(d_sys->mg_params->faces_only && build_rhs), PetscObjectComm((PetscObject)dm_stokes), PETSC_ERR_SUP, "RHS for faces-only not supported");
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_stokes, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        rho = d_sys->rho;
        mu = d_sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        hx = d_sys->hx;
        hy = d_sys->hy;
        PetscInt truncated_sys = d_sys->truncated_sys;

        /* Loop over all local elements. Note that it may be more efficient in real
           applications to loop over each boundary separately */
        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);
                if (right_boundary) { /* Right Boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = RIGHT;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                }
                if (top_boundary) {
                    /* Top boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = UP;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                if (bottom_boundary) {
                    /* Bottom boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                } else {
                    /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y */
                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
                    PetscInt      nEntries;

                    row.i   = ex;
                    row.j   = ey;
                    row.loc = DOWN;
                    row.c   = 0;
                    if (ex == 0) {
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1){
                            valA[count] = 0.0;
                            // Bottom adjacent boundary
                            valRhs =  mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        }
                        else {
                            valA[count] = - mu * 1.0 / (hy * hy);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (top_boundary){
                            valA[count]    = 0.0;
                            // Top adjacent Boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hy * hy);
                        }
                        /* Missing left element */
                        // Left near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][iprev], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hy;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hy;
                            ++count;
                        }
                    }else if (right_boundary) {
                        /* Right boundary y velocity stencil */
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1){
                            valA[count]    = 0.0;
                            // Bottom adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu  * 1.0 / (hy * hy);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (top_boundary){
                            valA[count]    = 0.0;
                            // Top adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hy * hy);
                        }
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        /* Missing right element */
                        // Right near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][inext], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hy;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = +1.0 / hy;
                            ++count;
                        }
                    }else if(ey==1){ //interior bottom boundary adjacent points
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0;
                        valRhs     = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hy;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hy;
                            ++count;
                        }
                    }else if(top_boundary){ //interior top boundary adjacent points
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0 ;
                        ++count;
                        valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hy;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hy;
                            ++count;
                        }
                    } else {
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hy;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hy;
                            ++count;
                        }
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                if (ex == 0) {
                    /* Left velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = LEFT;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                } else {
                    /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = LEFT;
                    row.c   = 0;

                    if (ey == 0) {
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
                        ++count;
                        /* missing term from element below */
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex == 1){
                            valA[count]    = 0.0;
                            // Left adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex == N[0] -1 ){
                            valA[count]    = 0.0;
                            // Right adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                        }
                        else {
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        // Bottom near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex - 1;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hx;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hx;
                            ++count;
                        }
                    } else if (top_boundary) {
                        /* Top boundary x velocity stencil */
                        count =0;
                        row.i      = ex;
                        row.j      = ey;
                        row.loc    = LEFT;
                        row.c      = 0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        /* Missing element above term */
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex==1){
                            valA[count]    = 0.0;
                            // Left adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex==N[0]-1){
                            valA[count]    = 0.0;
                            // Right adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }

                        // Up near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][inext],top_boundary,truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex - 1;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hx;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hx;
                            ++count;
                        }
                    } else {
                        // interior excluding ey=0,ey=N[1]-1
                        /* Note how this is identical to the stencil for U_y, with "DOWN" replaced by "LEFT" and the pressure derivative in the other direction */
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    =  rho/dt + mu *(2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;

                        if (ex == 1){
                            valA[count]    = 0.0;
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else {
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex == N[0] - 1){
                            valA[count]    = 0.0;
                            valRhs =  mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],top_boundary,truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                        if(!d_sys->mg_params->faces_only){
                            col[count].i   = ex - 1;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = -1.0 / hx;
                            ++count;
                            col[count].i   = ex;
                            col[count].j   = ey;
                            col[count].loc = ELEMENT;
                            col[count].c   = 0;
                            valA[count]    = 1.0 / hx;
                            ++count;
                        }
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                /* P equation : u_x + v_y = g
                   Note that this includes an explicit zero on the diagonal. This is only needed for
                   direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
                if(!d_sys->mg_params->faces_only){
                    if (pinPressure && ex == 0 && ey == 0) { /* Pin the first pressure node, if requested */
                        DMStagStencil row;
                        PetscScalar   valA, valRhs;
                        row.i   = ex;
                        row.j   = ey;
                        row.loc = ELEMENT;
                        row.c   = 0;
                        valA    = 1.0;
                        PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                        valRhs = pRef(cArrX[ex][icenter], cArrY[ey][icenter],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                    } else {
                        DMStagStencil row, col[5];
                        PetscScalar   valA[5], valRhs;

                        row.i      = ex;
                        row.j      = ey;
                        row.loc    = ELEMENT;
                        row.c      = 0;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].loc = LEFT;
                        col[0].c   = 0;
                        if (ex==0){
                            valA[0]    = 0.0;
                            // Left Boundary
                            valRhs = - (1.0 / hx) * uxRef(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[0]    = 1.0 / hx;
                        }
                        col[1].i   = ex;
                        col[1].j   = ey;
                        col[1].loc = RIGHT;
                        col[1].c   = 0;
                        if (ex == N[0]-1){
                            valA[1]    = 0.0;
                            // Right Boundary
                            valRhs = (1.0 / hx) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[1]    = -1.0 / hx;
                        }
                        col[2].i   = ex;
                        col[2].j   = ey;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        if (ey == 0){
                            valA[2]    = 0.0;
                            // Bottom Boundary
                            valRhs = -(1.0 / hy) * uyRef(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[2]     = 1.0 / hy;
                        }
                        col[3].i   = ex;
                        col[3].j   = ey;
                        col[3].loc = UP;
                        col[3].c   = 0;
                        if (ey == N[1]-1){
                            valA[3]    = 0.0;
                            // Top Boundary
                            valRhs = (1.0 / hy) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[3]    = -1.0 / hy;
                        }
                        col[4]     = row;
                        valA[4]    = 0.0;
                        PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 5, col, valA, INSERT_VALUES));
                        valRhs = g(cArrX[ex][icenter], cArrY[ey][icenter],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                    }
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        if(build_rhs) rhs_assembled = PETSC_TRUE;
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateTruncatedUnSteadyStokesSystem_NEUMANN(Mat *pA)
    {
        //  some code is repeated here
        //  calls original version and changes the boundary stencils only
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt level = d_sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;

        this->CreateUnSteadyStokesSystem(pA,nullptr);
        A = *pA;

        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        rho = d_sys->rho;
        mu = d_sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        hx = d_sys->hx;
        hy = d_sys->hy;

        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y */
                DMStagStencil row, col[7];
                PetscScalar   valA[7];
                PetscInt      nEntries;

                if(top_boundary){ //top boundary points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = UP;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    if(left_boundary||right_boundary){
                        valA[count]    = rho/dt + mu * ( 1.5 / (hx * hx) + 0.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;
                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);


                    ++count;

                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
                if(right_boundary){ //interior top boundary adjacent points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = RIGHT;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    if(top_boundary||bottom_boundary){
                        valA[count]    = rho/dt + mu * ( 0.5 / (hx * hx) + 1.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;

                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                row.i   = ex;
                row.j   = ey;
                row.loc = DOWN;
                row.c   = 0;
                count=0;

                if(bottom_boundary){ //on bottom boundary points
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    if(left_boundary||right_boundary){
                        valA[count]    = rho/dt + mu * ( 0.5 / (hx * hx) + 1.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);

                    ++count;
                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                }
                else if (left_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    /* Missing left element */
                    // Left near boundary

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hy;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hy;
                        ++count;
                    }
                }
                else if (right_boundary) {
                    /* Right boundary y velocity stencil */
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu  * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;
                    /* Missing right element */
                    // Right near boundary
                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hy;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = +1.0 / hy;
                        ++count;
                    }
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }


                /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                row.i   = ex;
                row.j   = ey;
                row.loc = LEFT;
                row.c   = 0;
                count =0;

                if (left_boundary) {
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    if(top_boundary||bottom_boundary){
                        valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 0.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    /* missing term from element below */
                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;

                    }
                }
                else if (bottom_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;
                    /* missing term from element below */

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hx;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hx;
                        ++count;
                    }
                }
                else if (top_boundary) {
                    /* Top boundary x velocity stencil */
                    row.i      = ex;
                    row.j      = ey;
                    row.loc    = LEFT;
                    row.c      = 0;

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;
                    /* Missing element above term */
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hx;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hx;
                        ++count;
                    }
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateTruncatedUnSteadyStokesSystem_TRUNCATED(Mat *pA)
    {
        //  some code is repeated here
        //  calls Neumann version and changes the boundary stencils only
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt level = d_sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;

        this->CreateTruncatedUnSteadyStokesSystem_NEUMANN(pA);
        A = *pA;

        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        rho = d_sys->rho;
        mu = d_sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        hx = d_sys->hx;
        hy = d_sys->hy;



        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y */
                DMStagStencil row, col[7];
                PetscScalar   valA[7];
                PetscInt      nEntries;

                if(top_boundary){ //top boundary points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = UP;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);


                    ++count;

                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
                if(right_boundary){ //interior right boundary adjacent points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = RIGHT;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;

                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                row.i   = ex;
                row.j   = ey;
                row.loc = DOWN;
                row.c   = 0;
                count=0;

                if(bottom_boundary){ //on bottom boundary points
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);

                    ++count;
                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }
                }
                else if (left_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    /* Missing left element */
                    // Left near boundary

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hy;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hy;
                        ++count;
                    }
                }
                else if (right_boundary) {
                    /* Right boundary y velocity stencil */
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu  * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;
                    /* Missing right element */
                    // Right near boundary
                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hy;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = +1.0 / hy;
                        ++count;
                    }
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                row.i   = ex;
                row.j   = ey;
                row.loc = LEFT;
                row.c   = 0;
                count =0;
                if (left_boundary) {
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    /* missing term from element below */
                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                    }
                } else if (bottom_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    /* missing term from element below */

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hx;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hx;
                        ++count;
                    }
                } else if (top_boundary) {
                    /* Top boundary x velocity stencil */
                    row.i      = ex;
                    row.j      = ey;
                    row.loc    = LEFT;
                    row.c      = 0;

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;
                    /* Missing element above term */
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    if(!d_sys->mg_params->faces_only){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -1.0 / hx;
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = 1.0 / hx;
                        ++count;
                    }
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateTruncatedVelocityFacesUnsteadySystem_NEUMANN(DM dm_uv, Mat *pA)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        Mat             A;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt        count;

        PetscFunctionBeginUser;
        this->CreateVelocityFacesUnsteadySystem(dm_uv, pA, nullptr);
        A = *pA;

        PetscCall(DMStagGetCorners(dm_uv, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_uv, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        mu = d_sys->mu;
        rho = d_sys->rho;

        hx = d_sys->hx;
        hy = d_sys->hy;

        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y */
                DMStagStencil row, col[5];
                PetscScalar   valA[5];
                PetscInt      nEntries;

                if(top_boundary){ //top boundary points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = UP;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    if(left_boundary||right_boundary){
                        valA[count]    = rho/dt + mu * ( 1.5 / (hx * hx) + 0.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;
                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);


                    ++count;

                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
                if(right_boundary){ //interior top boundary adjacent points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = RIGHT;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    if(top_boundary||bottom_boundary){
                        valA[count]    = rho/dt + mu * ( 0.5 / (hx * hx) + 1.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;

                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                row.i   = ex;
                row.j   = ey;
                row.loc = DOWN;
                row.c   = 0;
                count=0;

                if(bottom_boundary){ //on bottom boundary points
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    if(left_boundary||right_boundary){
                        valA[count]    = rho/dt + mu * ( 0.5 / (hx * hx) + 1.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);

                    ++count;
                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                }
                else if (left_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 1.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    /* Missing left element */
                    // Left near boundary

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                }
                else if (right_boundary) {
                    /* Right boundary y velocity stencil */
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu  * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;
                    /* Missing right element */
                    // Right near boundary
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }


                /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                row.i   = ex;
                row.j   = ey;
                row.loc = LEFT;
                row.c   = 0;
                count =0;

                if (left_boundary) {
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    if(top_boundary||bottom_boundary){
                        valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 0.5 / (hy * hy));
                    }
                    else{
                        valA[count]    = rho/dt + mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                    }
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    /* missing term from element below */
                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 0.5 / (hy * hy);
                        ++count;

                    }
                }
                else if (bottom_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;
                    /* missing term from element below */

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                }
                else if (top_boundary) {
                    /* Top boundary x velocity stencil */
                    row.i      = ex;
                    row.j      = ey;
                    row.loc    = LEFT;
                    row.c      = 0;

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;
                    /* Missing element above term */
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
            }
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateTruncatedVelocityFacesUnsteadySystem_TRUNCATED(DM dm_uv, Mat *pA)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        Mat             A;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt        count;
        /* Create the Velocity System Matrix for velocities defined on the edges */

        PetscFunctionBeginUser;
        this->CreateTruncatedVelocityFacesUnsteadySystem_NEUMANN(dm_uv,pA);
        A = *pA;

        PetscCall(DMStagGetCorners(dm_uv, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_uv, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        rho = d_sys->rho;
        mu = d_sys->mu;
        hx = d_sys->hx;
        hy = d_sys->hy;



        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                /* Y-momentum equation : (u_xx + u_yy) - p_y = f^y */
                DMStagStencil row, col[5];
                PetscScalar   valA[5];
                PetscInt      nEntries;

                if(top_boundary){ //top boundary points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = UP;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = UP;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);


                    ++count;

                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = UP;
                        col[count].c   = 0;
                        valA[count]    = -mu * 0.5 / (hx * hx);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
                if(right_boundary){ //interior right boundary adjacent points
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = RIGHT;
                    row.c   = 0;

                    count=0;
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = RIGHT;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;

                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = RIGHT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                row.i   = ex;
                row.j   = ey;
                row.loc = DOWN;
                row.c   = 0;
                count=0;

                if(bottom_boundary){ //on bottom boundary points
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hy * hy);

                    ++count;
                    if(!left_boundary){
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }

                    if(!right_boundary){
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }
                }
                else if (left_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    /* Missing left element */
                    // Left near boundary

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                }
                else if (right_boundary) {
                    /* Right boundary y velocity stencil */
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu  * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = DOWN;
                    col[count].c   = 0;
                    valA[count]    = -mu * 1.0 / (hx * hx);
                    ++count;
                    /* Missing right element */
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }

                /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                row.i   = ex;
                row.j   = ey;
                row.loc = LEFT;
                row.c   = 0;
                count =0;
                if (left_boundary) {
                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    /* missing term from element below */
                    if(!top_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                    }
                    if(!bottom_boundary){
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                    }
                } else if (bottom_boundary) {

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;
                    /* missing term from element below */

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                } else if (top_boundary) {
                    /* Top boundary x velocity stencil */
                    row.i      = ex;
                    row.j      = ey;
                    row.loc    = LEFT;
                    row.c      = 0;

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;
                    /* Missing element above term */
                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = LEFT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                }
                if(count){
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                }
            }
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreatePressureCenterUnsteadySystem(DM dm_p, Mat *pA, Vec *pRhs) {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscFunctionBeginUser;

        build_rhs = (PetscBool)(pRhs!=NULL);
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_p, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        mu = d_sys->mu;
        dt = d_sys->dt;
        rho = d_sys->rho;
        hx = d_sys->hx;
        hy = d_sys->hy;
        // Construct -L
        CreatePressureCenterPoissonSystem(dm_p, pA, pRhs);
        A = *pA;
        // A = rho/dt  I - mu L
        PetscCall(MatScale(A,mu));
        PetscCall(MatShift(A,rho/dt));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreatePressureCenterPoissonSystem(DM dm_p, Mat *pA, Vec *pRhs)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy;
        PetscScalar     **cArrX, **cArrY;

        PetscFunctionBeginUser;
        /* Probably don't need level info, creating operators with galerkin approach*/
        PetscInt level = d_sys->mg_params->level_curr;
        PetscInt count=0;
        /* Assuming dm_p has dof on cell elements only */
        PetscCall(DMCreateMatrix(dm_p, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_p, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_p, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_p, &N[0], &N[1], NULL));

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_p, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, RIGHT, &inext));
        hx = d_sys->hx;
        hy = d_sys->hy;

        PetscInt truncated_sys = d_sys->truncated_sys;
        /* Loop over all local elements. Note that it may be more efficient in real
           applications to loop over each boundary separately */
        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                DMStagStencil row, col[5];
                PetscScalar   valA[5], valRhs;

                row.i   = ex;
                row.j   = ey;
                row.loc = ELEMENT;
                row.c   = 0;

                count =0;
                /* p poisson equation : -(p_xx + p_yy) = f^p */
                if (left_boundary) { //p point is near left boundary
                    if (bottom_boundary) {
                        // Bottom Left
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    =  (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        // Missing Bottom Entry

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hy * hy);
                        ++count;

                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hx * hx);
                        ++count;

                    } else if (top_boundary) {
                        // Top Left

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    =  (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hy * hy);
                        ++count;

                        // Missing Top Entry


                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hx * hx);
                        ++count;
                    } else {
                        // Left

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    =  (1.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hy * hy);
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hy * hy);
                        ++count;


                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hx * hx);
                        ++count;
                    }
                }
                else if(right_boundary) {
                    if (bottom_boundary) {
                        // Bottom Right

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    =  (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        // Missing Bottom Entry

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = -  1.0 / (hy * hy);
                        ++count;

                        /* Missing right element */

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hx * hx);
                        ++count;
                    } else if (top_boundary) {
                        // Top Right


                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hy * hy);
                        ++count;

                        // Missing Top Entry
                        // Fill in boundary contribution in rhs
                       // valRhs = ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);


                        /* Missing right element */
                        // Fill in boundary contribution in rhs
                        //valRhs += ( 2.0/  (hx * hx) ) * pRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hx * hx);
                        ++count;
                    } else {
                        // Right

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = (1.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hy * hy);
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hy * hy);
                        ++count;


                        /* Missing right element */
                        // Fill in boundary contribution in rhs
                        //valRhs = ( 2.0/  (hx * hx) ) * pRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - 1.0 / (hx * hx);
                        ++count;
                    }
                }
                else if(top_boundary) {
                    // Top excluding corners

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hy * hy);
                    ++count;

                    // Missing Top Entry
                    // Fill in boundary contribution in rhs
                    //valRhs = ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys), truncated_sys;
                    //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));


                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;
                }
                else if(bottom_boundary) {
                    // Bottom excluding corners

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hy * hy);
                    ++count;

                    // Missing Bottom Entry
                    // Fill in boundary contribution in rhs
                    //valRhs = ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
                    //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;
                }
                else {
                    // Interior Points

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hy * hy);
                    ++count;


                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hy * hy);
                    ++count;


                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - 1.0 / (hx * hx);
                    ++count;
                }
                PetscCall(DMStagMatSetValuesStencil(dm_p, A, 1, &row, count, col, valA, INSERT_VALUES));
                valRhs = g(cArrX[ex][icenter], cArrY[ey][icenter],ex,ey, truncated_sys);
                if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_p, &cArrX, &cArrY, NULL));
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode StokesSolver::CreateVelocityFacesUnsteadySystem(DM dm_uv, Mat *pA, Vec *pRhs)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;

        /* Create the Velocity System Matrix for velocities defined on the edges */

        PetscFunctionBeginUser;
        PetscInt level = d_sys->mg_params->level_curr;
        PetscInt count=0;
        PetscCall(DMCreateMatrix(dm_uv, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_uv, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_uv, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_uv, &N[0], &N[1], NULL));
        dt = d_sys->dt;
        mu = d_sys->mu;
        rho = d_sys->rho;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_uv, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_uv, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_uv, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_uv, RIGHT, &inext));
        hx = d_sys->hx;
        hy = d_sys->hy;
        PetscInt truncated_sys = d_sys->truncated_sys;
        /* Loop over all local elements. Note that it may be more efficient in real
           applications to loop over each boundary separately */
        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {

                const PetscBool left_boundary   = (PetscBool)(ex == 0);
                const PetscBool right_boundary  = (PetscBool)(ex == N[0] - 1);
                const PetscBool bottom_boundary = (PetscBool)(ey == 0);
                const PetscBool top_boundary    = (PetscBool)(ey == N[1] - 1);

                if (right_boundary) {
                    /* Right Boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = RIGHT;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                }
                if (top_boundary) {
                    /* Top boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = UP;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                if (bottom_boundary) {
                    /* Bottom boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                } else {
                    /* Y-momentum equation : -(u_xx + u_yy) = f^y */
                    DMStagStencil row, col[5];
                    PetscScalar   valA[5], valRhs;
                    PetscInt      nEntries;

                    row.i   = ex;
                    row.j   = ey;
                    row.loc = DOWN;
                    row.c   = 0;
                    if (left_boundary) { //y point is near left boundary
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1) { //y point below is on bottom boundary
                            valA[count] = 0.0;
                            // Bottom adjacent boundary
                            valRhs =  mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        }
                        else {
                            valA[count] = - mu * 1.0 / (hy * hy);
                            valRhs = 0;
                        }
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (top_boundary){
                            valA[count]    = 0.0;
                            // Top adjacent Boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hy * hy);
                        }

                        /* Missing left element */
                        // Left near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][iprev], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                    }else if (right_boundary) {
                        /* Right boundary y velocity stencil */
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;

                        if (!(ey==1)){
                            col[count].i   = ex;
                            col[count].j   = ey - 1;
                            col[count].loc = DOWN;
                            col[count].c   = 0;
                            valA[count]    = - mu  * 1.0 / (hy * hy);
                            valRhs = 0.0;
                            ++count;
                        }
                        else{
                            // Bottom adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        }

                        if (!top_boundary){
                            col[count].i   = ex;
                            col[count].j   = ey + 1;
                            col[count].loc = DOWN;
                            col[count].c   = 0;
                            valA[count]    = - mu * 1.0 / (hy * hy);
                            ++count;
                        }
                        else{
                            // Top adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        }

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        /* Missing right element */
                        // Right near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][inext], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;

                    }else if(ey==1){ //interior bottom boundary adjacent points
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0;
                        valRhs     = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }else if(top_boundary){ //interior top boundary adjacent points
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0 ;
                        ++count;
                        valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    } else {
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                if (left_boundary) {
                    /* Left velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.loc                = LEFT;
                    row.c                  = 0;
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                } else {
                    /* X-momentum equation : (u_xx + u_yy) - p_x = f^x */
                    DMStagStencil row, col[5];
                    PetscScalar   valA[5], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.loc = LEFT;
                    row.c   = 0;

                    if (bottom_boundary) { //x point is near bottom boundary
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
                        ++count;
                        /* missing term from element below */
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex == 1){
                            valA[count]    = 0.0;
                            // Left adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex == N[0] -1 ){
                            valA[count]    = 0.0;
                            // Right adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                        }
                        else {
                            valA[count]    = -mu*1.0 / (hx * hx);
                        }
                        // Bottom near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                    } else if (top_boundary) {
                        /* Top boundary x velocity stencil */
                        count =0;
                        row.i      = ex;
                        row.j      = ey;
                        row.loc    = LEFT;
                        row.c      = 0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = rho/dt + mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        /* Missing element above term */
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (ex==1){
                            valA[count]    = 0.0;
                            // Left adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                            valRhs = 0;
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (right_boundary){
                            valA[count]    = 0.0;
                            // Right adjacent boundary
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }

                        // Up near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][inext],top_boundary,truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                    } else {
                        // interior points excluding near boundary points
                        /* Note how this is identical to the stencil for U_y, with "DOWN" replaced by "LEFT" and the pressure derivative in the other direction */
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    =  rho/dt + mu *( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;
                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;

                        if (ex == 1){//left boundary adjacent
                            valA[count]    = 0.0;
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else {
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        if (right_boundary){//right boundary adjacent
                            valA[count]    = 0.0;
                            valRhs =  mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dm_uv, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_uv, rhs, 1, &row, &valRhs, ADD_VALUES));
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_uv, &cArrX, &cArrY, NULL));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }

    PetscErrorCode StokesSolver::CreateVelocitySystem(Mat *pA, Vec *pRhs){
        PetscFunctionBeginUser;
        if(d_sys->truncated_sys){
            switch (d_sys->truncated_BCtype){
                case DIRICHLET:
                    PetscCall(this->CreateVelocityFacesUnsteadySystem(ctx_dm_uv, pA,nullptr));
                    break;
                case NEUMANN:
                    PetscCall(this->CreateTruncatedVelocityFacesUnsteadySystem_NEUMANN(ctx_dm_uv,pA));
                    break;
                case TRUNCATED_BC:
                    PetscCall(this->CreateTruncatedVelocityFacesUnsteadySystem_TRUNCATED(ctx_dm_uv,pA));
                    break;
                default:
                    std::cerr << "Unknown boundary condition type" << std::endl;
                    break;
            }
            if (debug) dumpMatToFile(*pA, "A_uv_truncated");
        }
        else{
            PetscCall(this->CreateVelocityFacesUnsteadySystem(ctx_dm_uv, pA, nullptr));
            if (debug) dumpMatToFile(*pA, "A_uv");
        }
        PetscCall(MatAssemblyBegin(*pA, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(*pA, MAT_FINAL_ASSEMBLY));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::CreateStokesSystem(Mat *pA, Vec *pRhs){
        PetscFunctionBeginUser;
        if(d_sys->truncated_sys){
            switch (d_sys->truncated_BCtype){
                case DIRICHLET:
                    PetscCall(this->CreateUnSteadyStokesSystem(pA, nullptr));
                    break;
                case NEUMANN:
                    PetscCall(this->CreateTruncatedUnSteadyStokesSystem_NEUMANN(pA));
                    break;
                case TRUNCATED_BC:
                    PetscCall(this->CreateTruncatedUnSteadyStokesSystem_TRUNCATED(pA));
                    break;
                default:
                    std::cerr << "Unknown boundary condition type" << std::endl;
                    break;
            }
            if (debug) dumpMatToFile(*pA, "A_stokes_truncated");
        }
        else{
            PetscCall(this->CreateUnSteadyStokesSystem(pA, pRhs));
            if (debug) dumpMatToFile(*pA, "A_stokes");
            if (debug) dumpVecToFile(*pRhs, "rhs_stokes");
        }
        PetscCall(MatAssemblyBegin(*pA, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(*pA, MAT_FINAL_ASSEMBLY));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::UpperPreconditionerApply(Vec x_eul, Vec f_s, Vec y_eul){
        PetscFunctionBeginUser;
        // Apply Stokes part of the Upper Triangular Preconditioner
        // Get the vectors for u,v,p
        PetscCall(DMStagMigrateVec(dm_stokes, f_s, ctx_dm_p, f_p));
        PetscCall(DMStagMigrateVec(dm_stokes, f_s, ctx_dm_uv, f_uv));
        PetscCall(DMStagMigrateVec(dm_stokes, x_eul, ctx_dm_p, x_p));
        PetscCall(DMStagMigrateVec(dm_stokes, x_eul, ctx_dm_uv, x_uv));


        // u_tilde = A_inv * (f_u)
        PetscCall(KSPSolve(ksp_uv,f_uv,u_tilde));

        // yp = -S1_inv (x_p + Bu_tilde)
        PetscCall(MatMult(B,u_tilde,Bu_tilde));
        PetscCall(VecAXPY(x_p,1.0,Bu_tilde));
        switch(d_sys->s1_type){
            case LSC:
                PetscCall(S1_Solve_LSC(x_p, y_p));
                break;
            case EXACT:
                PetscCall(KSPSolve(ksp_S1, x_p, y_p));
                break;
            default:
                std::cerr << "Unknown S1 type" << std::endl;
                break;
        }
        PetscCall(VecScale(y_p,-1.0));

        // yu = A_inv (x_uv - Bt * yp) - u_tilde
        // yu = u_star - u_tilde
        PetscCall(MatMult(Bt, y_p, Bty_p));
        PetscCall(VecAXPY(x_uv, -1.0, Bty_p));
        PetscCall(KSPSolve(ksp_uv, x_uv, u_star));
        PetscCall(VecAXPY(u_star, -1.0, u_tilde));

        PetscCall(DMStagMigrateVec(ctx_dm_p, y_p, dm_stokes, y_eul_p));
        PetscCall(DMStagMigrateVec(ctx_dm_uv, u_star, dm_stokes, y_eul));
        PetscCall(VecAXPY(y_eul,1.0,y_eul_p));

        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::S1_Solve_LSC(Vec x, Vec y){
        PetscFunctionBeginUser;
        // LSC approximation of S1_inv;
        // S1_inv = (BBt)_inv * BAB^t * (BBt)_inv
        // S1_inv = Ac (BBt)^-1
        PetscCall(KSPSolve(ksp_p, x, y_S1_p));
        PetscCall(MatMult(A_c, y_S1_p, y));
        //PetscCall(MatMult(A_uv,y_S1_uv,u_star));
        //PetscCall(MatMult(B, u_star, y_S1_p));
        //PetscCall(KSPSolve(ksp_p, y_S1_p, y));

        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::S1_Apply(Vec x, Vec y){
        PetscFunctionBeginUser;
        // S1 = B * A_inv * B^t
        PetscCall(MatMult(Bt, x, y_S1_uv));
        PetscCall(VecZeroEntries(u_star));
        PetscReal mean;
        PetscCall(VecMean(y_S1_uv,&mean));
        PetscCall(VecShift(y_S1_uv,-mean));
        PetscCall(KSPSolve(ksp_uv, y_S1_uv, u_star));
        PetscCall(MatMult(B, u_star, y));

        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::solveVelocity(Vec x, Vec y) {
        PetscFunctionBeginUser;

        // [f_uv, f_p] = x_eul
        PetscCall(DMStagMigrateVec(dm_stokes, x, ctx_dm_uv, f_uv));
        PetscCall(DMStagMigrateVec(dm_stokes, x, ctx_dm_p , f_p));

        // u_star = A_inv * (f_uv)
        PetscCall(KSPSolve(ksp_uv,f_uv,u_star));
        // y_uv = [u_star, 0]
        PetscCall(DMStagMigrateVec(ctx_dm_uv,u_star, dm_stokes, y_uv));
        // y = [0, f_p]
        PetscCall(DMStagMigrateVec(ctx_dm_p,f_p, dm_stokes, y));
        // y = [u_star, f_p]
        PetscCall(VecAXPY(y,1.0,y_uv));

        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode StokesSolver::AttachPressureNullspace(DM dmPressure, Mat *A)
    {
        Vec          constantPressure, basis;
        PetscReal    nrm;
        PetscFunctionBeginUser;
        //PetscCall(DMGetGlobalVector(dmPressure, &constantPressure));
        //PetscCall(VecSet(constantPressure, 1.0));
        //PetscCall(VecNorm(constantPressure, NORM_2, &nrm));
        //PetscCall(VecScale(constantPressure, 1.0 / nrm));
        //PetscCall(DMCreateGlobalVector(dmPressure, &basis));
        //PetscCall(VecCopy(constantPressure, basis));
        PetscCall(MatNullSpaceCreate(PetscObjectComm((PetscObject)dmPressure), PETSC_TRUE, 0, nullptr, &matNullSpacePressure));
        PetscCall(MatSetNullSpace(*A, matNullSpacePressure));
        //PetscCall(VecDestroy(&basis));
        //PetscCall(DMRestoreGlobalVector(dmPressure, &constantPressure));
	// Test Nullspace
	PetscBool isNull;
	PetscCall(MatNullSpaceTest(matNullSpacePressure, *A, &isNull));
	std::cout << "isNull : " << isNull << std::endl;
	// View Nullspace
	PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	PetscViewerSetType(viewer, PETSCVIEWERASCII);
	PetscCall(MatNullSpaceView(matNullSpacePressure, viewer));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
}
