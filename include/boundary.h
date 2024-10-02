#ifndef IBIMPLICIT_BOUNDARY_H
#define IBIMPLICIT_BOUNDARY_H
#include <iostream>
#include <petscsys.h>
#include <petscmat.h>
#include <utilities.h>
namespace ibImplicit{


    PetscErrorCode getBdryIndxLex(BCIndices bc_idx, PetscInt Nx, PetscInt Ny );
    PetscErrorCode modifyBoundaryRows(Mat *pA,
                                             std::vector<int> bc_idx,
                                             BCType bc_type,
                                             int Nx,
                                             int nfluid,
                                             char loc
                                            );
    PetscErrorCode modifyNearBoundaryRows(Mat *pA,
                                          std::vector<int> bc_idx,
                                          BCType bc_type,
                                          PetscInt Nx,
                                          PetscReal fac_h,
                                          PetscReal fac_v,
                                          PetscReal shft,
                                          char loc);
    PetscErrorCode determineRowsToRemove(PetscInt* num_rows,
                                         BCStaggeredType bc_stag_type,
                                         PetscInt Nx,
                                         PetscInt Ny);
}
#endif // IBIMPLICIT_BOUNDARY_H