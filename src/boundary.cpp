//
// Created by gkluhana on 16/04/24.
//
#include "../include/boundary.h"
namespace ibImplicit{
    PetscErrorCode getBdryIndxLex(BCIndices bc_idx, PetscInt Nx, PetscInt Ny ){
        PetscInt u_total = Ny* (Nx+1);
        PetscInt v_total = Nx* (Ny+1);
        PetscInt f_total = u_total+v_total;
        // East/Right Boundary Starts at Nx, incremented by Nx+1 each row
        for (int i = Nx ; i <= u_total-1 ; i += Nx + 1) {
            bc_idx->u_e.push_back(i);
        }
        //West Boundary Starts at 0, incremented by Nx+1
        for (int i = 0; i <= (Ny - 1) * (Nx + 1) ; i += Nx + 1) {
            bc_idx->u_w.push_back(i);
        }
        //South Boundary v
        for (int i = 0; i <= Nx-1; ++i) {
            bc_idx->v_s.push_back(i + u_total);
        }
        // North Boundary v
        for (int i = Nx * Ny; i < Nx * (Ny + 1); ++i) {
            bc_idx->v_n.push_back(i + u_total);
        }
        // Near Boundary Indices
        for (int i = 1 ; i <= Nx-1 ; ++i ) {
            bc_idx->u_s_near.push_back(i);
        }

        for (int i = (Ny-1)*(Nx+1)+1; i <= (Nx+1)*Ny-2 ; ++i ) {
            bc_idx->u_n_near.push_back(i);
        }

        for (int i = 2*Nx -1 ; i <= Nx*Ny -1  ; i+=Nx) {
            bc_idx->v_e_near.push_back(i + u_total);
        }

        for (int i = Nx ; i <= Nx * (Ny - 1); i+=Nx) {
            bc_idx->v_w_near.push_back(i + u_total);
        }

        // Pressure Boundary Indices
        for (int i = 0; i < Nx ;i++){
            bc_idx->p_s_near.push_back(i + f_total);
            bc_idx->p_n_near.push_back( i + Nx*(Ny-1) + f_total);
        }
        for(int i = 0; i<Ny ; i++ ){
            bc_idx->p_w_near.push_back(i * Nx + f_total);
            bc_idx->p_e_near.push_back((i+1)*Nx - 1 +f_total);
        }
        return 0;
    }//getBdryIndxLex
    PetscErrorCode determineRowsToRemove(PetscInt* num_rows,
                                         BCStaggeredType bc_stag_type, PetscInt Nx, PetscInt Ny)
                                         {
        PetscFunctionBeginUser;

            if (bc_stag_type->u_e == DIRICHLET) *num_rows += Ny;
            if (bc_stag_type->u_w == DIRICHLET) *num_rows += Ny;
        PetscFunctionReturn(0);
    }
    PetscErrorCode modifyBoundaryRows(Mat *pA,
                                             std::vector<int> bc_idx,
                                             BCType bc_type,
                                             int Nx,
                                             int nfluid,
                                             char loc
                                             ) {
        PetscFunctionBeginUser;
        Mat A = *pA;
        switch (bc_type) {
            case DIRICHLET:
                if (loc == 'e' || loc == 'w') {
                    for (int idx: bc_idx) {
                        MatSetValue(A, idx, idx, 1.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx + 1, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx - 1, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx - Nx - 1, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx + Nx + 1, 0.0, INSERT_VALUES);
                        //Remove corresponding entries from adjacent nodes to preserve symmetry
                        MatSetValue(A, idx+1, idx, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx - 1, idx, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx - Nx - 1, idx, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx + Nx + 1, idx, 0.0, INSERT_VALUES);
                    }
                } else if (loc == 's') {
                    for (int idx: bc_idx) {
                        MatSetValue(A, idx, idx, 1.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx + 1, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx - 1, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx + Nx, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx, idx - Nx, 0.0, INSERT_VALUES);
                        //Remove corresponding entries from adjacent nodes to preserve symmetry
                        MatSetValue(A,  idx + 1,idx, 0.0, INSERT_VALUES);
                        MatSetValue(A,  idx - 1,idx, 0.0, INSERT_VALUES);
                        MatSetValue(A,  idx + Nx,idx, 0.0, INSERT_VALUES);
                        MatSetValue(A,  idx - Nx,idx, 0.0, INSERT_VALUES);
                    }
                } else if (loc == 'n') {
                    for (int idx: bc_idx) {
                        MatSetValue(A, idx, idx, 1.0, INSERT_VALUES);
                        if (idx < nfluid - 1) {
                            MatSetValue(A, idx, idx + 1, 0.0, INSERT_VALUES);
                            MatSetValue(A,  idx + 1, idx,0.0, INSERT_VALUES);
                        }
                        MatSetValue(A, idx, idx - 1, 0.0, INSERT_VALUES);
                        MatSetValue(A,  idx - 1, idx,0.0, INSERT_VALUES);
                        if (idx < nfluid - Nx - 2) {
                            MatSetValue(A, idx, idx + Nx + 1, 0.0, INSERT_VALUES);
                            MatSetValue(A, idx + Nx + 1,idx, 0.0, INSERT_VALUES);
                        }
                        MatSetValue(A, idx, idx - Nx, 0.0, INSERT_VALUES);
                        MatSetValue(A, idx - Nx,idx , 0.0, INSERT_VALUES);
                    }
                }
                break;
            case NEUMANN:
                //TODO: TREATMENT FOR NEUMANN BOUNDARY NODES
                break;
            default:
                std::cout << "Invalid BC Type, should throw error" << std::endl;
                break;
        }
        PetscFunctionReturn(0);
    }//modifyBoundaryRows
    PetscErrorCode modifyNearBoundaryRows(Mat *pA,
                                          std::vector<int> bc_idx,
                                          BCType bc_type,
                                          PetscInt Nx,
                                          PetscReal fac_h,
                                          PetscReal fac_v,
                                          PetscReal shft,
                                          char loc) {
        PetscFunctionBeginUser;
        Mat A = *pA;
        switch (bc_type){
            case DIRICHLET:
                if (loc == 's' ){
                    for (int idx: bc_idx){
                        MatSetValue(A, idx, idx, -(2*fac_h + 3*fac_v)+shft, INSERT_VALUES);
                    }
                }
                else if (loc == 'n'){
                    for (int idx: bc_idx){
                        MatSetValue(A, idx, idx, -(2*fac_h + 3*fac_v)+shft, INSERT_VALUES);
                    }
                }
                else if(loc == 'e') {
                    for (int idx: bc_idx) {
                        MatSetValue(A, idx, idx,  -(3*fac_h + 2*fac_v)+ shft, INSERT_VALUES);
                        MatSetValue(A, idx, idx + 1, 0.0, INSERT_VALUES);
                    }
                }
                else if (loc == 'w'){
                     for (int idx: bc_idx){
                         MatSetValue(A, idx, idx, -(3*fac_h + 2*fac_v)+shft, INSERT_VALUES);
                         MatSetValue(A, idx, idx-1, 0.0, INSERT_VALUES);
                     }

                }
                break;

            case NEUMANN:
                //TODO: TREATMENT FOR NEUMANN NEAR BOUNDARY NODES
                break;
            default:
                std::cout << "Invalid BC Type, should throw error" << std::endl;
                break;
        }
            PetscFunctionReturn(0);
    } //modifyNearBoundaryRows
}