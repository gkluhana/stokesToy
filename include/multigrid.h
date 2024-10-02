//
// Created by gkluhana on 29/04/24.
//
#include <staggered.h>
#ifndef IBIMPLICIT_MULTIGRID_H
#define IBIMPLICIT_MULTIGRID_H
namespace ibImplicit{
    // This functionality was created before using Petsc DM+multigrid
    inline PetscErrorCode LevelCtxCreate(LevelCtx *p_level_ctx) {
        LevelCtx level_ctx;
        PetscFunctionBeginUser;
        PetscCall(PetscMalloc1(1, p_level_ctx));
        level_ctx = *p_level_ctx;
        level_ctx->dm_level= NULL;
        level_ctx->dm_faces= NULL;
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode MultigridOnPressurePoissonDM(SystemParameters sys);
    PetscErrorCode MultigridOnVelPoissonDM(SystemParameters sys);
    PetscErrorCode MultigridOnVelDM(SystemParameters sys);
    PetscErrorCode SolveDMStokes(SystemParameters sys);
    PetscErrorCode InterpolateVelocityU(Vec* uc, Vec *uf, PetscInt nxc, PetscInt nyc, PetscInt nxf, PetscInt nyf);
    PetscErrorCode InterpolateVelocityV(Vec* vc, Vec *vf, PetscInt nxc, PetscInt nyc, PetscInt nxf, PetscInt nyf);
}
#endif //IBIMPLICIT_MULTIGRID_H
