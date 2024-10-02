//
// Created by gkluhana on 12/04/24.
//
#include "../include/staggered.h"

namespace ibImplicit {
    const double PI = 3.141592653589793;
    PetscScalar uxRef(PetscScalar x, PetscScalar y,PetscBool is_top, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;

        if (is_top){
            return 0.5*(1.0-std::cos(2*PI*x));
            //return 1 - x*x*x*x;
        }else
        {
            return 0;
        }
        //return 20*x*y*y*y;
    } /* no x-dependence  */
    PetscScalar uxRef(PetscScalar x, PetscScalar y,PetscInt truncated_sys)
    {
        if(truncated_sys) return 0;
        return 0;
        //return 20*x*y*y*y;
    } /* no x-dependence  */
    PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscInt truncated_sys)
    {
        //return 5.0 * (x*x*x*x - y*y*y*y);
        if (truncated_sys) return 0;
        return 0;
    } /* no y-dependence  */
    PetscScalar pRef(PetscScalar x, PetscScalar y, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;
        return 60*x*x*y  - 20*y*y*y;
        //return 0;
    } /* zero integral    */
    PetscScalar fx(PetscScalar x, PetscScalar y, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;
        return 0.0 * x * y;
        //    return 1 ;
    } /* no x-dependence  */
    PetscScalar fy(PetscScalar x, PetscScalar y, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;
        return 0* (x*x - y*y);
        //return 1 ;
    }
    PetscScalar g(PetscScalar x, PetscScalar y, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;
        return 0.0 ;
    } /* identically zero */
    PetscScalar g(PetscScalar x, PetscScalar y, PetscInt ex, PetscInt ey, PetscInt truncated_sys)
    {
        if (truncated_sys) return 0;
        if (ex==16 & ey==16){
            return -5;
        }
        else if(ex==48 & ey==48 ){
            return 5;
        }
        else{
            return 0.0 * x * y;
        }
    } /* identically zero */
    PetscErrorCode DumpStokesData(DM dm, Vec x, std::string var_name) {
        PetscInt ex, e, startx, starty, nx, ny;
        Vec vec_c;
        DM da_center;
        Vec vec_u,vec_v;
        DM da_u_left, da_v_down;
        PetscFunctionBeginUser;
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_ELEMENT, 0, &da_center, &vec_c));
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_LEFT, 0, &da_u_left, &vec_u));
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_DOWN, 0, &da_v_down, &vec_v));
        PetscCall(PetscObjectSetName((PetscObject)vec_c, ("p_"+var_name).c_str()));
        PetscCall(PetscObjectSetName((PetscObject)vec_u, ("u_"+var_name).c_str()));
        PetscCall(PetscObjectSetName((PetscObject)vec_v, ("v_"+var_name).c_str()));
        /* Dump edge-based field to .vtr file */
        {
            PetscViewer viewer;
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_u_left), ("./dump/u_"+var_name+".vtr").c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_u, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        /* Dump edge-based field to a second .vtr file */
        {
            PetscViewer viewer;
            std::string filename = "./dump/v_" + var_name + ".vtr";
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_v_down), filename.c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_v, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        /* Dump center-based field to .vtr file */
        {
            PetscViewer viewer;
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_center), ("./dump/p_"+ var_name+".vtr").c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_c, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        PetscCall(DMDestroy(&da_v_down));
        PetscCall(VecDestroy(&vec_v));
        PetscCall(DMDestroy(&da_u_left));
        PetscCall(VecDestroy(&vec_u));
        PetscCall(DMDestroy(&da_center));
        PetscCall(VecDestroy(&vec_c));
        PetscFunctionReturn(0);
    }
    PetscErrorCode DumpPressureData(DM dm, Vec x, std::string var_name) {
        PetscInt ex, e, startx, starty, nx, ny;
        Vec vec_c;
        DM da_center;
        PetscFunctionBeginUser;
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_ELEMENT, 0, &da_center, &vec_c));
        PetscCall(PetscObjectSetName((PetscObject)vec_c, var_name.c_str()));
        /* Dump edge-based field to .vtr file */
        {
            PetscViewer viewer;
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_center), ("./dump/"+var_name+".vtr").c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_c, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        PetscCall(DMDestroy(&da_center));
        PetscCall(VecDestroy(&vec_c));
        PetscFunctionReturn(0);
    }
    PetscErrorCode DumpVelocityData(DM dm, Vec x, std::string var_name) {
        PetscInt ex, e, startx, starty, nx, ny;
        Vec vec_u,vec_v;
        DM da_u_left, da_v_down;
        PetscFunctionBeginUser;
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_LEFT, 0, &da_u_left, &vec_u));
        PetscCall(DMStagVecSplitToDMDA(dm, x, DMSTAG_DOWN, 0, &da_v_down, &vec_v));
        PetscCall(PetscObjectSetName((PetscObject)vec_u, (var_name+"_left").c_str()));
        PetscCall(PetscObjectSetName((PetscObject)vec_v, (var_name+"_down").c_str()));
        /* Dump edge-based field to .vtr file */
        {
            PetscViewer viewer;
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_u_left), ("./dump/"+var_name+"_left.vtr").c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_u, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        /* Dump edge-based field to a second .vtr file */
        {
            PetscViewer viewer;
            PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_v_down), ("./dump/"+var_name+"_down.vtr").c_str(), FILE_MODE_WRITE, &viewer));
            PetscCall(VecView(vec_v, viewer));
            PetscCall(PetscViewerDestroy(&viewer));
        }
        PetscCall(DMDestroy(&da_v_down));
        PetscCall(VecDestroy(&vec_v));
        PetscCall(DMDestroy(&da_u_left));
        PetscCall(VecDestroy(&vec_u));
        PetscFunctionReturn(0);
    }

    PetscErrorCode StaggeredGrid::ApplyStokes2D(Vec* x, Vec* y){
        PetscFunctionBeginUser;
        PetscInt Nx = d_sys-> Nx;
        PetscInt Ny =   d_sys-> Ny;
        PetscReal dt =  d_sys->dt;
        PetscReal mu =  d_sys->mu;
        PetscReal rho = d_sys->rho;
        PetscInt nu =   d_sys->nu;
        PetscInt nv =   d_sys->nv;
        PetscReal hx =  d_sys->hx;
        PetscReal hy =  d_sys->hy;
        PetscReal fac_c = - (1/(hx*hx)+ 1/(hy*hy))* mu ;
        PetscReal fac_h = - 1.0 / (hx*hx) * mu ;
        PetscReal fac_v = - 1.0 / (hy*hy) * mu ;
        PetscReal shft = rho / dt ;
        BCStaggeredType bc_stag_type = d_sys->bcSpec->bc_stag_type;
        Vec u = *x ;
        Vec Au = *y ;
        //Apply on Interior Nodes
        PetscScalar u_ii;
        PetscScalar u_left;
        PetscScalar u_right;
        PetscScalar u_top;
        PetscScalar u_bottom;
        PetscInt ileft;
        PetscInt iright;
        PetscInt itop;
        PetscInt ibottom;
        PetscScalar Au_ii;
        // Apply to horizontal velocity u DOF
        for (PetscInt i=0; i < nu ; i ++){
            ileft = i-1;
            iright = i+1;
            itop = i+Nx+1;
            ibottom = i-Nx-1;
            VecGetValues(u, 1, &i, &u_ii);
            if (i % (Nx+1) != Nx) VecGetValues(u, 1, &iright, &u_right);
            if (i % (Nx+1) != 0)  VecGetValues(u, 1, &ileft, &u_left);
            if (i <= (Nx +1)*(Ny-1) ) VecGetValues(u, 1, &itop, &u_top);
            if (i > Nx) VecGetValues(u, 1, &ibottom, &u_bottom);

            // West Boundary
            if (i % (Nx+1) == 0) {
                // Treat based on Boundary Type
                if (bc_stag_type->u_w == DIRICHLET){ Au_ii =u_ii;}// don't do anything
                else if (bc_stag_type->u_w == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // East Boundary
            else if (i % (Nx+1) == Nx) {
                if (bc_stag_type->u_e == DIRICHLET){ Au_ii =u_ii;}// don't do anything
                else if (bc_stag_type->u_e == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // West Adjacent Nodes Symmetric
            else if (i % (Nx+1) ==1) {
                if (bc_stag_type->u_w == DIRICHLET) {
                    //Corner Nodes
                    if (i == 1) {
                        // South West Node
                        Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*u_right + fac_v*u_top + shft;
                    } else if (i == (Nx + 1) * (Ny - 1) + 1) {
                        // North West Node
                        Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*u_right + fac_v*u_bottom + shft;
                    } else {
                        Au_ii = -2.0*fac_c * u_ii + fac_h*u_right + fac_v*(u_top+u_bottom) + shft;
                    }
                }
            }
            // East Adjacent Nodes Symmetric
            else if ( ( i % (Nx+1) ) == (Nx-1)){
                if (bc_stag_type->u_e == DIRICHLET){
                    //Corner Nodes
                    if (i == Nx-1 ){
                        //South East Node
                        Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*u_left+ fac_v*u_top + shft;
                    }
                    else if (i== (Nx+1)*(Ny-1)+ Nx-1){
                        //North East Node
                        Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*u_left+ fac_v*u_bottom+ shft;
                    }
                    else{
                        Au_ii = -2.0*fac_c * u_ii + fac_h*u_left+ fac_v*(u_top+u_bottom) + shft;
                    }
                }
            }
            // South Boundary
            else if (i < Nx){
                if(bc_stag_type->u_s== DIRICHLET){
                    Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*(u_left+u_right)+ fac_v*u_top + shft;
                }
                else if (bc_stag_type->u_s == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // North Boundary
            else if (i > (Nx +1)*(Ny-1) ){
                if(bc_stag_type->u_n== DIRICHLET){
                    Au_ii = -(2.0*fac_h + 3*fac_v) * u_ii + fac_h*(u_left+u_right)+ fac_v*u_bottom+ shft;
                }
                else if (bc_stag_type->u_n == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            //Interior Node
            else {
                Au_ii = -2.0*fac_c * u_ii + fac_h*(u_left+ u_right)+ fac_v*(u_top+u_bottom) + shft;
            }
            VecSetValue(Au, i, Au_ii, INSERT_VALUES );
        }
        // Apply to vertical velocity u DOF
        PetscInt offset = nu;
        PetscInt ii=0;
        for (PetscInt i=0; i < nv ; i ++){
            ii = i + offset;
            ileft = i-1 +offset;
            iright = i+1 + offset;
            itop = i+Nx + offset;
            ibottom = i-Nx + offset;
            VecGetValues(u, 1, &ii, &u_ii);

            if (i % (Nx) != Nx-1 ) VecGetValues(u, 1, &iright, &u_right);
            if (i % (Nx) != 0) VecGetValues(u, 1, &ileft, &u_left);
            if (i >= Nx) VecGetValues(u, 1, &ibottom, &u_bottom);
            if (i < (Nx)*(Ny) ) VecGetValues(u, 1, &itop, &u_top);

            // South Boundary
            if (i < Nx){
                if(bc_stag_type->u_s== DIRICHLET){
                    Au_ii = u_ii;
                }
                else if (bc_stag_type->u_s == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // South Adjacent Nodes
            else if (i < 2*Nx){
                if (i==Nx){
                    //South East Node
                    Au_ii = -(3.0*fac_h + 2*fac_v)* u_ii + fac_h*u_right + fac_v*u_top + shft;
                }
                else if(i== 2*Nx-1){
                    //South West Node
                    Au_ii = -(3.0*fac_h + 2*fac_v)* u_ii + fac_h*u_left+ fac_v*u_top + shft;
                }
                else {
                    Au_ii = -2.0*fac_c * u_ii + fac_h*(u_right +u_left)+ fac_v*u_top+ shft;
                }
            }
            // North Boundary
            else if (i >= (Nx)*(Ny) ){
                if(bc_stag_type->u_n== DIRICHLET){
                    Au_ii = u_ii;
                }
                else if (bc_stag_type->u_n == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // North Adjacent Nodes
            else if (i >= Nx*(Ny-1)){
                if (i==Nx*(Ny-1)){
                    Au_ii =  -(3.0*fac_h + 2*fac_v)* u_ii + fac_h*u_right + fac_v*u_bottom+ shft;
                }
                else if(i== Nx*Ny -1){
                    Au_ii =  -(3.0*fac_h + 2*fac_v)* u_ii + fac_h*u_left+ fac_v*u_bottom+ shft;
                }
                else {
                    Au_ii = -2.0*fac_c * u_ii + fac_h*(u_left + u_right) + fac_v*u_bottom+ shft;
                }
            }
            // West Boundary
            else if (i % (Nx) == 0) {
                // Treat based on Boundary Type
                if (bc_stag_type->u_w == DIRICHLET){
                    Au_ii = -(3.0*fac_h + 2*fac_v) * u_ii + fac_h*u_right+ fac_v*(u_top + u_bottom) + shft;
                }//
                else if (bc_stag_type->u_w == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            // East Boundary
            else if (i % (Nx) == Nx-1) {
                if (bc_stag_type->u_e == DIRICHLET){
                    Au_ii = -(3.0*fac_h + 2*fac_v) * u_ii + fac_h*u_left+ fac_v*(u_top + u_bottom) + shft;
                }
                else if (bc_stag_type->u_e == NEUMANN) {
                    // TODO: Treat Neumann Boundary
                }
            }
            //Interior Node
            else {
                Au_ii = -2.0*fac_c * u_ii + fac_h*(u_left+ u_right)+ fac_v*(u_top+u_bottom) + shft;
            }

            VecSetValue(Au, i+offset, Au_ii, INSERT_VALUES );
        }
        VecAssemblyBegin(Au);
        VecAssemblyEnd(Au);

        PetscFunctionReturn(0);
    }//ApplyStokes2D
    PetscErrorCode StaggeredGrid::ApplyStaggeredDiv2D(Vec*x , Vec* y){
        PetscFunctionBeginUser;
        PetscInt Nx = d_sys-> Nx;
        PetscInt Ny = d_sys-> Ny;
        PetscInt nu = d_sys->nu;
        PetscReal hx = d_sys->hx;
        PetscReal hy = d_sys->hy;
        Vec u = *x ;
        Vec Bu = *y ;
        // Traverse over cells, Nx*Ny
        // For each cell, apply central differencing in both directions
        PetscScalar u_left, u_right,  v_top,  v_bottom;
        PetscInt ileft, iright, itop, ibottom;
        PetscInt u_offset, v_offset, cell ;
        PetscScalar Bu_ii;
        PetscInt j =0;
        for (PetscInt i =0; i<Nx*Ny; i++){
            u_offset = j*(Nx+1);
            v_offset = j*(Nx) + nu;
            cell = i % Nx;
            ileft = u_offset + cell ;
            iright = ileft + 1;
            ibottom= v_offset + cell ;
            itop= ibottom+ Nx;
            Bu_ii = 0;
            if (i % Nx != (Nx-1)) {
                VecGetValues(u, 1, &iright, &u_right);
                Bu_ii += 1/hx *(u_right);
            }
            if (i % Nx != 0)  {
                VecGetValues(u, 1, &ileft, &u_left);
                Bu_ii -= 1/hx *(u_left);
            }
            if (i < Nx*(Ny-1)) {
                VecGetValues(u, 1, &itop, &v_top);
                Bu_ii += 1/hy *(v_top);
            }
            if (i >= Nx) {
                VecGetValues(u, 1, &ibottom, &v_bottom);
                Bu_ii -= 1/hy * (v_bottom);
            }
            VecSetValue(Bu,i,Bu_ii,INSERT_VALUES);
            if (i % Nx == (Nx-1)) j++;
        }
        VecAssemblyBegin(Bu);
        VecAssemblyEnd(Bu);
        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::AssembleVelocityRHS(Vec*x , Vec* y){
        PetscFunctionBeginUser;
        PetscInt Nx =   d_sys-> Nx;
        PetscReal dt =  d_sys->dt;
        PetscReal mu =  d_sys->mu;
        PetscReal rho = d_sys->rho;
        PetscInt nu =  d_sys->nu;
        PetscReal hx =  d_sys->hx;
        PetscReal hy =  d_sys->hy;
        PetscReal fac_v =2.0/(hy*hy);
        PetscReal fac_h =2.0/(hx*hx);
        PetscReal shft = rho / dt ;
        BCIndices bc_idx= d_sys->bcSpec->bc_idx;
        Vec u = *x ;
        Vec rhs = *y ;
        VecCopy(u,rhs);
        VecScale(rhs, shft);
        PetscScalar bc_val;
        PetscScalar x_coord, y_coord;
        // Near Boundary Nodes
        for (auto idx: bc_idx->v_e_near){
            x_coord = d_sys->Lx;
            y_coord = d_sys->ly + ( (idx-nu)/(Nx+1) ) * hy;
            d_sys->bcSpec->bc_funcs->ge(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= fac_h;
            VecSetValue(rhs,idx,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->v_w_near){ x_coord =0;
            x_coord = d_sys->lx;
            y_coord = d_sys->ly + ( (idx -nu)/(Nx) ) * hy;
            d_sys->bcSpec->bc_funcs->gw(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= fac_h;
            VecSetValue(rhs,idx,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->u_n_near){
            x_coord = d_sys->lx + (idx % (Nx+1)) * hx;
            y_coord = d_sys->Ly;
            d_sys->bcSpec->bc_funcs->fn(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= fac_v;
            VecSetValue(rhs,idx,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->u_s_near){
            x_coord = d_sys->lx + (idx % (Nx+1)) * hx;
            y_coord = d_sys->ly ;
            d_sys->bcSpec->bc_funcs->fs(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *=fac_v;
            VecSetValue(rhs,idx,bc_val,ADD_VALUES);
        }
        fac_v = 1/(hy*hy);
        fac_h = 1/(hx*hx);
        // Boundary Nodes and Boundary Adjacent Nodes
        for (auto idx : bc_idx->u_e){
            x_coord = d_sys->Lx;
            y_coord = d_sys->ly + (idx/(Nx)-1) * hy + hy/2.0;
            d_sys->bcSpec->bc_funcs->fe(x_coord,y_coord,mu,hx,hy,&bc_val);
            VecSetValue(rhs,idx,bc_val,INSERT_VALUES);
            bc_val *= fac_h;
            VecSetValue(rhs,idx-1,bc_val,ADD_VALUES);
        }
        for (auto idx : bc_idx->u_w){
            x_coord = d_sys->lx;
            y_coord = d_sys->ly + (idx/(Nx+1)) * hy + hy/2.0;
            d_sys->bcSpec->bc_funcs->fw(x_coord,y_coord,mu,hx,hy,&bc_val);
            VecSetValue(rhs,idx,bc_val,INSERT_VALUES);
            bc_val *= fac_h;
            VecSetValue(rhs,idx+1,bc_val,ADD_VALUES);
        }
        for (auto idx : bc_idx->v_s){
            x_coord = d_sys->lx +((idx) % (Nx+1)) * hx + hx/2.0;
            y_coord = d_sys->ly;
            d_sys->bcSpec->bc_funcs->gs(x_coord,y_coord,mu,hx,hy,&bc_val);
            VecSetValue(rhs,idx,bc_val,INSERT_VALUES);
            bc_val *= fac_v;
            VecSetValue(rhs,idx+Nx,bc_val,ADD_VALUES);
        }
        for (auto idx : bc_idx->v_n){
            x_coord = d_sys->lx + ( (idx % Nx) ) * hx + hx/2.0;
            y_coord = d_sys->Ly;
            d_sys->bcSpec->bc_funcs->gn(x_coord,y_coord,mu,hx,hy,&bc_val);
            VecSetValue(rhs,idx,bc_val,INSERT_VALUES);
            bc_val *= fac_v;
            VecSetValue(rhs,idx-Nx,bc_val,ADD_VALUES);
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::AssemblePressureRHS( Vec* y){
        PetscFunctionBeginUser;
        PetscInt Nx =   d_sys-> Nx;
        PetscInt Ny =   d_sys-> Ny;
        PetscReal mu =  d_sys->mu;
        PetscInt nu =  d_sys->nu;
        PetscInt nf =  d_sys->nfluid;
        PetscReal hx =  d_sys->hx;
        PetscReal hy =  d_sys->hy;
        PetscReal fac_v = 1.0/(hy);
        PetscReal fac_h = 1.0/(hx);
        BCIndices bc_idx= d_sys->bcSpec->bc_idx;
        PetscScalar bc_val;
        PetscScalar x_coord, y_coord;
        Vec rhs = *y ;
        // Set to 0
        VecAssemblyBegin(rhs);
        VecAssemblyEnd(rhs);
        PetscScalar zero = 0.0;
        VecSet(rhs,zero);
        // Near Boundary Nodes
        for (auto idx: bc_idx->p_e_near){
            x_coord = d_sys->Lx;
            y_coord = d_sys->ly + ( (idx-nf)/(Nx) ) * hy + hy/2.0;
            d_sys->bcSpec->bc_funcs->fe(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= fac_h;
            VecSetValue(rhs,idx - nf,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->p_w_near){
            x_coord = d_sys->lx;
            y_coord = d_sys->ly + ( (idx-nf)/(Nx) ) * hy + hy/2.0;
            d_sys->bcSpec->bc_funcs->fw(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= -fac_h;
            VecSetValue(rhs,idx - nf,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->p_s_near){
            x_coord = d_sys->lx+ ( (idx-nf) ) * hy + hy/2.0;
            y_coord = d_sys->ly;
            d_sys->bcSpec->bc_funcs->gs(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= -fac_v;
            VecSetValue(rhs,idx - nf,bc_val,ADD_VALUES);
        }
        for (auto idx: bc_idx->p_n_near){
            x_coord = d_sys->lx+ ( (idx-nf- Nx*(Ny-1) ) ) * hy + hy/2.0;
            y_coord = d_sys->Ly;
            d_sys->bcSpec->bc_funcs->gn(x_coord,y_coord,mu,hx,hy,&bc_val);
            bc_val *= fac_v;
            VecSetValue(rhs,idx - nf,bc_val,ADD_VALUES);
        }


        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::ApplyCellCenteredGrad2D(Vec*x , Vec* y){
        PetscFunctionBeginUser;
        PetscInt Nx = d_sys-> Nx;
        PetscInt Ny = d_sys-> Ny;
        PetscInt nu = d_sys->nu;
        PetscReal hx = d_sys->hx;
        PetscReal hy = d_sys->hy;
        Vec p = *x ;
        Vec Btp = *y ;
        PetscScalar p_left, p_right,  p_top,  p_bottom;
        PetscInt ileft, iright, itop, ibottom;
        PetscInt u_offset, v_offset, u_edge , v_edge;
        PetscScalar Btp_ii ;
        BCSpec bcSpec = d_sys->bcSpec;
        BCStaggeredType bc_stag_type = bcSpec->bc_stag_type;
        PetscInt j =0;
        // Traverse over vertical edges
        // For each edge, apply horizontal central differencing
        for (PetscInt i =0; i<(Nx+1)*Ny; i++){
            u_offset = j*(Nx);
            u_edge = i % (Nx + 1);
            ileft = u_offset + u_edge - 1;
            iright = ileft + 1;
            Btp_ii = 0;
            if (i % (Nx+1) != (Nx)) {
                VecGetValues(p, 1, &iright, &p_right);
                if(i %(Nx+1)==0 && bc_stag_type->u_e == DIRICHLET){}
                else{
                    Btp_ii+= 1/hx *(p_right);
                }
            }
            if (i % (Nx+1) != 0)  {
                VecGetValues(p, 1, &ileft, &p_left);
                if(i %(Nx+1)==Nx && bc_stag_type->u_w == DIRICHLET){}
                else{
                    Btp_ii -= 1/hx *(p_left);
                }
            }
            VecSetValue(Btp,i,Btp_ii,INSERT_VALUES);
            if (i % (Nx+1) == (Nx)) j++;
        }
        // Traverse over horizontal edges
        // For each edge, apply vertical central differencing
        j=0;
        PetscInt offset = nu;
        for (PetscInt i =0; i<(Nx+1)*Ny; i++){
            v_offset = j*(Nx) ;
            v_edge = i % (Nx);
            itop= v_offset+v_edge;
            ibottom=  itop- Nx;
            Btp_ii = 0;
            if (i >= Nx) {
                VecGetValues(p, 1, &ibottom, &p_bottom);
                if(bc_stag_type->u_n == DIRICHLET && i >= Nx*(Ny-1)){
                }
                else {
                    Btp_ii-= 1/hy *(p_bottom);
                }
            }
            if (i < Nx*(Ny-1))  {
                VecGetValues(p, 1, &itop, &p_top);
                if(bc_stag_type->u_s == DIRICHLET && i < Nx){
                }
                else{
                    Btp_ii += 1/hy *(p_top);
                }
            }
            VecSetValue(Btp,i+offset,Btp_ii,INSERT_VALUES);
            if (i % (Nx) == (Nx-1)) j++;
        }
        VecAssemblyBegin(Btp);
        VecAssemblyEnd(Btp);
        PetscFunctionReturn(0);

    }
    PetscErrorCode StaggeredGrid::AssembleStokes2D() {
        PetscFunctionBeginUser;
        PetscInt Nx = d_sys-> Nx;
        PetscInt Ny = d_sys-> Ny;
        std::cout << "Assemble Stokes 2D called for Nx, Ny :"<< Nx << " , " << Ny << std::endl;
        PetscReal dt = d_sys->dt;
        PetscReal mu = d_sys->mu;
        PetscReal rho = d_sys->rho;
        PetscInt nfluid = d_sys->nfluid;

        PetscInt nu = d_sys->nu;
        PetscInt nv = d_sys->nv;
        PetscReal hx = d_sys->hx;
        PetscReal hy = d_sys->hy;
        PetscReal fac_c = - (1/(hx*hx)+ 1/(hy*hy))* mu ;
        PetscReal fac_h = - 1.0 / (hx*hx) * mu ;
        PetscReal fac_v = - 1.0 / (hy*hy) * mu ;
        PetscReal shft = rho / dt ;
        BCIndices bc_idx = d_sys->bcSpec->bc_idx;
        BCStaggeredType bc_stag_type = d_sys->bcSpec->bc_stag_type;
        MatCreate(PETSC_COMM_WORLD, &d_A);

        // Decide Av DOF
        MatSetSizes(d_A, PETSC_DECIDE, PETSC_DECIDE, nfluid, nfluid);
        MatSetUp(d_A);
        for (PetscInt i=0; i < nu ; i ++){
            MatSetValue(d_A, i, i, -2.0*fac_c + shft, INSERT_VALUES);
            if (i != 0) MatSetValue(d_A, i, i-1, 1.0*fac_h, INSERT_VALUES);
            if (i != (nu - 1) ) MatSetValue(d_A, i, i+1, 1.0*fac_h, INSERT_VALUES);
            if (i > Nx) MatSetValue(d_A, i, i-Nx-1, 1.0*fac_v, INSERT_VALUES);
            if (i < (nu- Nx - 1)) MatSetValue(d_A, i, i+Nx+1, 1.0*fac_v, INSERT_VALUES);
        }
        PetscInt offset = nu;
        for (PetscInt i=0; i < nv ; i ++){
            MatSetValue(d_A, i+offset, i+offset, -2.0*fac_c+ shft, INSERT_VALUES);
            if (i != 0) MatSetValue(d_A, i + offset, i-1 + offset, 1.0*fac_h, INSERT_VALUES);
            if (i != (nv - 1) ) MatSetValue(d_A, i + offset, i+1 + offset, 1.0*fac_h, INSERT_VALUES);
            if (i > Nx-1) MatSetValue(d_A, i + offset, i-Nx + offset, 1.0*fac_v, INSERT_VALUES);
            if (i < (nv- Nx)) MatSetValue(d_A, i + offset, i+Nx + offset, 1.0*fac_v, INSERT_VALUES);
        }
        // Dirichlet on Boundaries
        // Call Function to Fill in Boundaries
        std::vector<PetscInt> *rows_to_remove;
        //if (bc_stag_type->u_w == DIRICHLET) nAu -= Ny;
        modifyBoundaryRows(&d_A, bc_idx->u_e, bc_stag_type->u_e, Nx, nfluid, 'e');
        modifyBoundaryRows(&d_A, bc_idx->u_w, bc_stag_type->u_w, Nx, nfluid, 'w');
        modifyBoundaryRows(&d_A, bc_idx->v_s, bc_stag_type->u_s, Nx, nfluid, 's');
        modifyBoundaryRows(&d_A, bc_idx->v_n, bc_stag_type->u_n, Nx, nfluid, 'n');
        // Dirichlet on Near Boundaries
        modifyNearBoundaryRows(&d_A, bc_idx->u_s_near, bc_stag_type->u_s, Nx, fac_h, fac_v, shft , 's');
        modifyNearBoundaryRows(&d_A, bc_idx->u_n_near, bc_stag_type->u_n, Nx, fac_h, fac_v, shft , 'n');
        modifyNearBoundaryRows(&d_A, bc_idx->v_e_near, bc_stag_type->u_e, Nx, fac_h, fac_v, shft , 'e');
        modifyNearBoundaryRows(&d_A, bc_idx->v_w_near, bc_stag_type->u_w, Nx, fac_h, fac_v, shft , 'w');
        // Remove Dirichlet Rows
        // Assemble the matrix
        PetscCall(MatAssemblyBegin(d_A, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(d_A, MAT_FINAL_ASSEMBLY));
        PetscFunctionReturn(0);
    }
    PetscErrorCode GenerateFluidGridMatrices(SystemParameters sys){
        PetscFunctionBeginUser;
        StaggeredGrid grid(sys);
        Mat A = grid.getA();
        Mat Bt = grid.getB();
        Mat B = grid.getBt();
        dumpMatToFile(A, "A_MAT");
        std::cout << "\n Printed A Matrix in computeMatrix \n" << std::endl;
        dumpMatToFile(B, "B_MAT");
        std::cout << "\n Printed B Matrix in computeMatrix \n" << std::endl;
        dumpMatToFile(Bt, "Bt_MAT");
        std::cout << "\n Printed Bt Matrix in computeMatrix \n" << std::endl;
        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::AssembleStaggeredDiv2D(){
        PetscFunctionBeginUser;
        PetscInt Nx = d_sys->Nx;
        PetscInt Ny = d_sys->Ny;
        PetscInt np = d_sys->Nx*Ny;
        PetscInt nu = d_sys->nu;
        PetscInt nv = d_sys-> nv;
        PetscReal hx = d_sys-> hx;
        PetscReal hy = d_sys-> hy;
        MatCreate(PETSC_COMM_WORLD, &d_B);
        MatSetSizes(d_B, PETSC_DECIDE, PETSC_DECIDE, np, nu+nv);
        MatSetUp(d_B);
        PetscInt j = 0;
        for (PetscInt i=0; i < Nx*Ny ; i++) {
            if(i % Nx !=0) MatSetValue(d_B, i, i+j, -1.0 / hx , INSERT_VALUES);
            if(i % Nx != (Nx-1))MatSetValue(d_B, i, i+j+1, 1.0 / hx, INSERT_VALUES);
            if(i >= Nx) MatSetValue(d_B, i, i+nu, -1.0 / hy , INSERT_VALUES);
            if(i < Nx*(Ny-1))MatSetValue(d_B, i, i+nu+Nx, 1.0 / hy, INSERT_VALUES);
            if ( (i % Nx) == (Nx-1) ) j++;
        }
        MatAssemblyBegin(d_B,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(d_B,MAT_FINAL_ASSEMBLY);
        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::AssembleCellCenteredGrad2D(){
        PetscFunctionBeginUser;
        std::cout << "Assemble Cell Centered Grad Called" << std::endl;
        PetscInt Nx = d_sys->Nx;
        PetscInt Ny = d_sys->Ny;
        PetscInt nu = d_sys->nu;
        PetscInt nv = d_sys-> nv;

        PetscScalar zero = 0.0;
        AssembleStaggeredDiv2D();
        MatTranspose(d_B, MAT_INITIAL_MATRIX, &d_Bt);
        PetscInt row = 0;
        // TODO: Check if this is really needed?
        MatZeroRows(d_Bt,1,&row,zero,PETSC_NULLPTR,PETSC_NULLPTR );
        row = (Nx+1)*Ny -1 ;
        MatZeroRows(d_Bt,1,&row,zero,PETSC_NULLPTR,PETSC_NULLPTR );
        for (PetscInt j=1; j < Ny ; j++) {
            row =j*(Nx+1);
            MatZeroRows(d_Bt,1,&row,zero,PETSC_NULLPTR,PETSC_NULLPTR );
            row -= 1;
            MatZeroRows(d_Bt,1,&row,zero,PETSC_NULLPTR,PETSC_NULLPTR );
        }
        row = nu;
        PetscInt rowNorth = nu+nv-1;
        for (PetscInt i =0; i< Nx; i++){
            MatZeroRows(d_Bt,1,&row,zero, PETSC_NULLPTR,PETSC_NULLPTR);
            row++;
            MatZeroRows(d_Bt,1,&rowNorth,zero, PETSC_NULLPTR,PETSC_NULLPTR);
            rowNorth--;
        }
        MatScale(d_Bt,-1.0);
        PetscFunctionReturn(0);
    }
    PetscErrorCode CheckFluidSystem(StaggeredGrid *grid){
        PetscFunctionBeginUser;
        SystemParameters sys = grid->getSys();
        Mat A = grid->getA();
        Mat Bt = grid->getBt();
        Mat B = grid->getB();;
        dumpMatToFile(A, "A_MAT");
        dumpMatToFile(B, "B_MAT");
        dumpMatToFile(Bt, "Bt_MAT");
        Mat Z22;
        MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, sys->np, sys->np, 0.0,&Z22);
        MatSetUp(Z22);
        MatAssemblyBegin(Z22,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Z22,MAT_FINAL_ASSEMBLY);
        Mat Z23;
        MatCreate(PETSC_COMM_WORLD,&Z23);
        MatSetSizes(Z23, PETSC_DECIDE, PETSC_DECIDE, sys->np, 2*sys->ns );
        MatSetUp(Z23);
        MatZeroEntries(Z23);
        Mat Z32;
        MatCreate(PETSC_COMM_WORLD,&Z32);
        MatSetSizes(Z32, PETSC_DECIDE, PETSC_DECIDE, 2*sys->ns, sys->np );
        MatSetUp(Z32);
        MatZeroEntries(Z32);
        Mat I;
        MatCreateConstantDiagonal(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, 2*sys->ns, 2*sys->ns, 1.0,&I);
        MatSetUp(I);
        MatAssemblyBegin(I,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(I,MAT_FINAL_ASSEMBLY);

        Mat A_s;
        Mat array[2*2];
        array[0] = A;
        array[1] = Bt;

        array[2] = B;
        array[3] =Z22;
        //Assemble Block System
        MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL,array,&A_s);
        MatAssemblyBegin(A_s,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A_s,MAT_FINAL_ASSEMBLY);
        Mat A_s_Mat;
        MatConvert(A_s,MATSEQAIJ,MAT_INITIAL_MATRIX,&A_s_Mat);
        dumpMatToFile(A_s_Mat,"A_s");
        Vec rhs,x,rhsf,rhsp;

        MatCreateVecs(A_s_Mat,&rhs,&x);
        IS rows[2];
        MatNestGetISs(A_s,rows,NULL);
        VecGetSubVector(rhs,rows[0],&rhsf);
        VecGetSubVector(rhs,rows[1],&rhsp);
        grid->AssembleFluidRHS(&rhsf,&rhsp);
        VecRestoreSubVector(rhs,rows[0],&rhsf);
        VecRestoreSubVector(rhs,rows[1],&rhsp);
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp,A_s_Mat,A_s_Mat);
        KSPSetUp(ksp);
        KSPSolve(ksp,rhs,x);
        dumpVecToFile(x,"u_sol");
        KSPDestroy(&ksp);
        PetscFunctionReturn(0);

    }
    PetscErrorCode CheckMatrixFreeOps(StaggeredGrid *grid){
        PetscFunctionBeginUser;
        // Test Matrix Free Stokes 2D
        SystemParameters sys = grid->getSys();
        Vec u;
        Vec v;
        Vec p;
        Vec Bu;
        Vec Btp;
        VecCreate(PETSC_COMM_WORLD, &u);
        VecSetSizes(u, PETSC_DECIDE, (sys)->nfluid);
        VecSetType(u, VECSEQ);
        VecAssemblyBegin(u);
        VecAssemblyEnd(u);
        PetscScalar one = 1.0;
        VecSet(u,one);
        VecCreate(PETSC_COMM_WORLD, &v);
        VecSetSizes(v, PETSC_DECIDE, (sys)->nfluid);
        VecSetType(v, VECSEQ);
        grid->ApplyStokes2D(&u,&v);
        dumpVecToFile(v, "Au_VEC");

        VecCreate(PETSC_COMM_WORLD, &Bu);
        VecSetSizes(Bu, PETSC_DECIDE, (sys)->np);
        VecSetType(Bu, VECSEQ);
        grid->ApplyStaggeredDiv2D(&u ,  &Bu);
        dumpVecToFile(Bu, "Bu_VEC");

        VecCreate(PETSC_COMM_WORLD, &p);
        VecSetSizes(p, PETSC_DECIDE, (sys)->np);
        VecSetType(p, VECSEQ);
        VecAssemblyBegin(p);
        VecAssemblyEnd(p);
        VecSet(p,one);

        VecCreate(PETSC_COMM_WORLD, &Btp);
        VecSetSizes(Btp, PETSC_DECIDE, (sys)->nfluid);
        VecSetType(Btp, VECSEQ);
        grid->ApplyCellCenteredGrad2D(&p ,  &Btp);
        ibImplicit::dumpVecToFile(Btp, "Btp_VEC");
        VecDestroy(&u);
        VecDestroy(&v);
        VecDestroy(&p);
        VecDestroy(&Bu);
        VecDestroy(&Btp);
        PetscFunctionReturn(0);
    }
    PetscErrorCode StaggeredGrid::AssembleFluidRHS(Vec *prhsf, Vec* prhsp){
        PetscFunctionBeginUser;
       Vec un;
       Vec pn;
       Vec rhsf = *prhsf;
       Vec rhsp = *prhsp;
       VecCreate(PETSC_COMM_WORLD, &un);
       VecSetSizes(un, PETSC_DECIDE, (d_sys)->nfluid);
       VecSetType(un, VECSEQ);
       VecAssemblyBegin(un);
       VecAssemblyEnd(un);
       PetscScalar zero = 0.0;
       VecSet(un,zero);
       AssembleVelocityRHS(&un, &rhsf);
       ibImplicit::dumpVecToFile(rhsf, "rhsf");
       AssemblePressureRHS( &rhsp);
       ibImplicit::dumpVecToFile(rhsp, "rhsp");
       PetscFunctionReturn(0);
    }
    PetscErrorCode CreateSystem(DM dmSol, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        /* Here, we showcase two different methods for manipulating local vector entries.
           One can use DMStagStencil objects with DMStagVecSetValuesStencil(),
           making sure to call VecAssemble[Begin/End]() after all values are set.
           Alternately, one can use DMStagVecGetArray[Read]() and DMStagVecRestoreArray[Read]().
           The first approach is used to build the rhs, and the second is used to
           obtain coordinate values. Working with the array is almost certainly more efficient,
           but only allows setting local entries, requires understanding which "slot" to use,
           and doesn't correspond as precisely to the matrix assembly process using DMStagStencil objects */
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count;


        PetscFunctionBeginUser;
        PetscCall(DMCreateMatrix(dmSol, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        PetscCheck(!(sys->mg_params->faces_only && build_rhs), PetscObjectComm((PetscObject)dmSol), PETSC_ERR_SUP, "RHS for faces-only not supported");
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dmSol, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dmSol, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dmSol, &N[0], &N[1], NULL));
        rho = sys->rho;
        dt = sys->dt;
        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, RIGHT, &inext));
        PetscReal Lx = cArrX[0][2*N[0]];
        PetscReal Ly = cArrY[0][2*N[1]];
        hx = Lx / N[0];
        hy = Ly / N[1];
        PetscInt truncated_sys = sys->truncated_sys;

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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if (ey==N[1]-1){
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        }
                    }else if (ex == N[0] - 1) {
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
                        if (ey == N[1]-1){
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
                        ++count;
                        /* Missing right element */
                        // Right near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][inext], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));

                        if(!sys->mg_params->faces_only){
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
                        }
                    }else if(ey==1){
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
                        ++count;
                        valRhs     = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        }
                    }else if(ey== N[1] - 1 ){
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        }
                    }
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            valA[count]    = -1.0 / (hx * hx);
                        }
                        // Bottom near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][iprev],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;

                        if(!sys->mg_params->faces_only){
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
                        }
                    } else if (ey == N[1] - 1) {
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));

                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        }
                    } else {
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            valRhs =  mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        if(!sys->mg_params->faces_only){
                            ++count;
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
                        }
                    }
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                /* P equation : u_x + v_y = g
                   Note that this includes an explicit zero on the diagonal. This is only needed for
                   direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
                if(!sys->mg_params->faces_only){
                    if (pinPressure && ex == 0 && ey == 0) { /* Pin the first pressure node, if requested */
                        DMStagStencil row;
                        PetscScalar   valA, valRhs;
                        row.i   = ex;
                        row.j   = ey;
                        row.loc = ELEMENT;
                        row.c   = 0;
                        valA    = 1.0;
                        PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                        valRhs = pRef(cArrX[ex][icenter], cArrY[ey][icenter],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[3]    = -1.0 / hy;
                        }
                        col[4]     = row;
                        valA[4]    = 0.0;
                        PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 5, col, valA, INSERT_VALUES));
                        valRhs = g(cArrX[ex][icenter], cArrY[ey][icenter],truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                    }
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreateGradientCenterOperator(DM dm_stokes, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys){
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        /* Here, we showcase two different methods for manipulating local vector entries.
           One can use DMStagStencil objects with DMStagVecSetValuesStencil(),
           making sure to call VecAssemble[Begin/End]() after all values are set.
           Alternately, one can use DMStagVecGetArray[Read]() and DMStagVecRestoreArray[Read]().
           The first approach is used to build the rhs, and the second is used to
           obtain coordinate values. Working with the array is almost certainly more efficient,
           but only allows setting local entries, requires understanding which "slot" to use,
           and doesn't correspond as precisely to the matrix assembly process using DMStagStencil objects */
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;
        PetscCall(DMCreateMatrix(dm_stokes, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        PetscCheck(!(sys->mg_params->faces_only && build_rhs), PetscObjectComm((PetscObject)dm_stokes), PETSC_ERR_SUP, "RHS for faces-only not supported");
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_stokes, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        dt = sys->dt;
        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        PetscReal Lx = *cArrX[N[0]];
        PetscReal Ly = *cArrY[N[1]];
        PetscReal lx = *cArrX[0];
        PetscReal ly = *cArrY[0];
        hx = (Lx-lx)/N[0];
        hy = (Ly-ly)/N[1];
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreateUnSteadyStokesSystem(DM dm_stokes, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;

        PetscCall(DMCreateMatrix(dm_stokes, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        PetscCheck(!(sys->mg_params->faces_only && build_rhs), PetscObjectComm((PetscObject)dm_stokes), PETSC_ERR_SUP, "RHS for faces-only not supported");
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_stokes, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));
        dt = sys->dt;
        rho = sys->rho;
        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;
        PetscInt truncated_sys = sys->truncated_sys;

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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                        if(!sys->mg_params->faces_only){
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
                if(!sys->mg_params->faces_only){
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
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreateSteadyStokesSystem(DM dm_stokes, Mat *pA, Vec *pRhs, PetscBool pinPressure, SystemParameters sys)
    {
        PetscInt        N[2];
        PetscInt        ex, ey, startx, starty, nx, ny;
        PetscInt        iprev, icenter, inext;
        Mat             A;
        Vec             rhs;
        PetscBool       build_rhs;
        PetscReal       hx, hy, rho, dt, mu;
        PetscScalar     **cArrX, **cArrY;
        /* Here, we showcase two different methods for manipulating local vector entries.
           One can use DMStagStencil objects with DMStagVecSetValuesStencil(),
           making sure to call VecAssemble[Begin/End]() after all values are set.
           Alternately, one can use DMStagVecGetArray[Read]() and DMStagVecRestoreArray[Read]().
           The first approach is used to build the rhs, and the second is used to
           obtain coordinate values. Working with the array is almost certainly more efficient,
           but only allows setting local entries, requires understanding which "slot" to use,
           and doesn't correspond as precisely to the matrix assembly process using DMStagStencil objects */
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count;
        PetscFunctionBeginUser;
        PetscCall(DMCreateMatrix(dm_stokes, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        PetscCheck(!(sys->mg_params->faces_only && build_rhs), PetscObjectComm((PetscObject)dm_stokes), PETSC_ERR_SUP, "RHS for faces-only not supported");
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dm_stokes, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dm_stokes, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_stokes, &N[0], &N[1], NULL));

        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_stokes, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;
        PetscInt truncated_sys = sys->truncated_sys;
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
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
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
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
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
                        valA[count]    = mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1){
                            valA[count] = 0.0;
                            // Bottom adjacent boundary
                            valRhs =  mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
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
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hy * hy);
                        }
                        /* Missing left element */
                        // Left near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][iprev], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                        if(!sys->mg_params->faces_only){
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
                        valA[count]    = mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1){
                            valA[count]    = 0.0;
                            // Bottom adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
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
                        if (ey == N[1]-1){
                            valA[count]    = 0.0;
                            // Top adjacent boundary
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
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
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][inext], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));

                        ++count;
                        if(!sys->mg_params->faces_only){
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
                        valA[count]    = mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0;
                        valRhs     = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
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
                        if(!sys->mg_params->faces_only){
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
                        valA[count]    = mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                        valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
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
                        if(!sys->mg_params->faces_only){
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                        if(!sys->mg_params->faces_only){
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
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
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
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter], truncated_sys);
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        }
                        else {
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        // Bottom near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        if(!sys->mg_params->faces_only){
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
                    } else if (ey == N[1] - 1) {
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }

                        // Up near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][inext],top_boundary, truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                        if(!sys->mg_params->faces_only){
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
                        valA[count]    =  mu *(2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs =  mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                        if(!sys->mg_params->faces_only){
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
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                }

                /* P equation : u_x + v_y = g
                   Note that this includes an explicit zero on the diagonal. This is only needed for
                   direct solvers (not required if using an iterative solver and setting the constant-pressure nullspace) */
                if(!sys->mg_params->faces_only){
                    if (pinPressure && ex == 0 && ey == 0) { /* Pin the first pressure node, if requested */
                        DMStagStencil row;
                        PetscScalar   valA, valRhs;
                        row.i   = ex;
                        row.j   = ey;
                        row.loc = ELEMENT;
                        row.c   = 0;
                        valA    = 1.0;
                        PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                        valRhs = pRef(cArrX[ex][icenter], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = - (1.0 / hx) * uxRef(cArrX[ex][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = (1.0 / hx) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = -(1.0 / hy) * uyRef(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
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
                            valRhs = (1.0 / hy) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[3]    = -1.0 / hy;
                        }
                        col[4]     = row;
                        valA[4]    = 0.0;
                        PetscCall(DMStagMatSetValuesStencil(dm_stokes, A, 1, &row, 5, col, valA, INSERT_VALUES));
                        valRhs = g(cArrX[ex][icenter], cArrY[ey][icenter], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_stokes, rhs, 1, &row, &valRhs, ADD_VALUES));
                    }
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_stokes, &cArrX, &cArrY, NULL));
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreateVelocityFacesUnsteadySystem(DM dmSol, Mat *pA, Vec *pRhs, SystemParameters sys)
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
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count=0;
        PetscCall(DMCreateMatrix(dmSol, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dmSol, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dmSol, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dmSol, &N[0], &N[1], NULL));
        dt = sys->dt;
        mu = sys->mu;
        rho = sys->rho;

        PetscCall(DMStagGetProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;
        PetscInt truncated_sys = sys->truncated_sys;



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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][inext], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));

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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter],truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreateVelocityFacesPoissonSystem(DM dmSol, Mat *pA, Vec *pRhs, SystemParameters sys)
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
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count=0;
        PetscCall(DMCreateMatrix(dmSol, pA));
        A = *pA;
        build_rhs = (PetscBool)(pRhs!=NULL);
        if(build_rhs){
            PetscCall(DMCreateGlobalVector(dmSol, pRhs));
            rhs = *pRhs;
        } else{
            rhs = NULL;
        }
        PetscCall(DMStagGetCorners(dmSol, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dmSol, &N[0], &N[1], NULL));
        dt = sys->dt;
        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dmSol, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;
        PetscInt truncated_sys = sys->truncated_sys;


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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uyRef(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        valA[count]    = mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        if (ey==1) { //y point below is on bottom boundary
                            valA[count] = 0.0;
                            // Bottom adjacent boundary
                            valRhs =  mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
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
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hy * hy);
                        }

                        /* Missing left element */
                        // Left near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][iprev], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));

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
                        valA[count]    = mu * (3.0 / (hx * hx) + 2.0 / (hy * hy));
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
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
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
                            valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                        }

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = -mu * 1.0 / (hx * hx);
                        /* Missing right element */
                        // Right near boundary
                        valRhs += mu * ( 2.0/  (hx * hx) ) * uyRef(cArrX[ex][inext], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;

                    }else if(ey==1){ //interior bottom boundary adjacent points
                        count=0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;
                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = DOWN;
                        col[count].c   = 0;
                        valA[count]    = 0.0;
                        valRhs     = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey-1][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        valA[count]    = mu * ( 2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                        valRhs = mu * ( 1.0 / (hy * hy) ) * uyRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fy(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES));
                    valRhs = uxRef(cArrX[ex][iprev], cArrY[ey][icenter], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        }
                        else {
                            valA[count]    = -mu*1.0 / (hx * hx);
                        }
                        // Bottom near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][iprev], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                        valA[count]    = mu * (2.0 / (hx * hx) + 3.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }

                        // Up near boundary
                        valRhs +=  mu * (2.0/(hy * hy)) * uxRef(cArrX[ex][iprev], cArrY[ey][inext], truncated_sys);
                        if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        ++count;
                    } else {
                        // interior points excluding near boundary points
                        /* Note how this is identical to the stencil for U_y, with "DOWN" replaced by "LEFT" and the pressure derivative in the other direction */
                        count =0;
                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = LEFT;
                        col[count].c   = 0;
                        valA[count]    =  mu *( 2.0 / (hx * hx) + 2.0 / (hy * hy));
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
                            valRhs = mu * (1.0/(hx * hx)) * uxRef(cArrX[ex-1][iprev], cArrY[ey][icenter], truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
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
                            valRhs =  mu * (1.0/(hx * hx)) * uxRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                            if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                        }
                        else{
                            valA[count]    = - mu * 1.0 / (hx * hx);
                        }
                        ++count;
                    }
                    PetscCall(DMStagMatSetValuesStencil(dmSol, A, 1, &row, count, col, valA, INSERT_VALUES));
                    valRhs = fx(cArrX[ex][iprev], cArrY[ey][icenter], truncated_sys);
                    if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, ADD_VALUES));
                }
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dmSol, &cArrX, &cArrY, NULL));
        PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyBegin(rhs));
        PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
        if(build_rhs) PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreatePressureCenterExactSolution(DM dm_p, Vec *pRhs, SystemParameters sys)
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
        /* Probably don't need level info, creating operators with galerkin approach*/
        PetscInt level = sys->mg_params->level_curr;
        PetscInt count=0;
        /* Assuming dm_p is on cell elements only */
        PetscCall(DMStagGetCorners(dm_p, &startx, &starty, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
        PetscCall(DMStagGetGlobalSizes(dm_p, &N[0], &N[1], NULL));
        dt = sys->dt;
        mu = sys->mu;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_p, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;

        PetscInt truncated_sys = sys->truncated_sys;

        PetscCall(DMCreateGlobalVector(dm_p, pRhs));
        rhs = *pRhs;
        DMStagStencil row;
        PetscScalar   valRhs;
        /* Loop over all local elements. Note that it may be more efficient in real
           applications to loop over each boundary separately */
        for (ey = starty; ey < starty + ny; ++ey) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ex = startx; ex < startx + nx; ++ex) {
                row.i   = ex;
                row.j   = ey;
                row.loc = ELEMENT;
                row.c   = 0;
                valRhs = pRef(cArrX[ex][icenter], cArrY[ey][icenter], truncated_sys);
                PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));
            }
        }
        PetscCall(DMStagRestoreProductCoordinateArraysRead(dm_p, &cArrX, &cArrY, NULL));
        PetscCall(VecAssemblyBegin(rhs));
        PetscCall(VecAssemblyEnd(rhs));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CreatePressureCenterPoissonSystem(DM dm_p, Mat *pA, Vec *pRhs, SystemParameters sys)
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
        /* Probably don't need level info, creating operators with galerkin approach*/
        PetscInt level = sys->mg_params->level_curr;
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
        mu = 1.0;

        PetscCall(DMStagGetProductCoordinateArraysRead(dm_p, &cArrX, &cArrY, NULL));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, ELEMENT, &icenter));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, LEFT, &iprev));
        PetscCall(DMStagGetProductCoordinateLocationSlot(dm_p, RIGHT, &inext));
        hx = sys->hx;
        hy = sys->hy;

        PetscInt truncated_sys = sys->truncated_sys;
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
                        valA[count]    = mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        // Missing Bottom Entry

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;

                    } else if (top_boundary) {
                        // Top Left

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        // Missing Top Entry


                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                    } else {
                        // Left

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = mu * (1.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;


                        /* Missing left element */

                        col[count].i   = ex + 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
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
                        valA[count]    = mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        // Missing Bottom Entry

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        /* Missing right element */

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                    } else if (top_boundary) {
                        // Top Right


                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = mu * (1.0 / (hx * hx) + 1.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        // Missing Top Entry
                        // Fill in boundary contribution in rhs
                        valRhs = mu * ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys);


                        /* Missing right element */
                        // Fill in boundary contribution in rhs
                        valRhs += mu * ( 2.0/  (hx * hx) ) * pRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                    } else {
                        // Right

                        col[count].i   = ex;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = mu * (1.0 / (hx * hx) + 2.0 / (hy * hy));
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey - 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;

                        col[count].i   = ex;
                        col[count].j   = ey + 1;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hy * hy);
                        ++count;


                        /* Missing right element */
                        // Fill in boundary contribution in rhs
                        valRhs = mu * ( 2.0/  (hx * hx) ) * pRef(cArrX[ex][inext], cArrY[ey][icenter], truncated_sys);
                        //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                        col[count].i   = ex - 1;
                        col[count].j   = ey;
                        col[count].loc = ELEMENT;
                        col[count].c   = 0;
                        valA[count]    = - mu * 1.0 / (hx * hx);
                        ++count;
                    }
                }
                else if(top_boundary) {
                    // Top excluding corners

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    // Missing Top Entry
                    // Fill in boundary contribution in rhs
                    valRhs = mu * ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][inext], truncated_sys), truncated_sys;
                    //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));


                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                }
                else if(bottom_boundary) {
                    // Bottom excluding corners

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = mu * (2.0 / (hx * hx) + 1.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;

                    // Missing Bottom Entry
                    // Fill in boundary contribution in rhs
                    valRhs = mu * ( 2.0/  (hy * hy) ) * pRef(cArrX[ex][icenter], cArrY[ey][iprev], truncated_sys);
                    //if(build_rhs) PetscCall(DMStagVecSetValuesStencil(dm_p, rhs, 1, &row, &valRhs, ADD_VALUES));

                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;
                }
                else {
                    // Interior Points

                    col[count].i   = ex;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = mu * (2.0 / (hx * hx) + 2.0 / (hy * hy));
                    ++count;

                    col[count].i   = ex;
                    col[count].j   = ey + 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;


                    col[count].i   = ex;
                    col[count].j   = ey - 1;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hy * hy);
                    ++count;


                    col[count].i   = ex - 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
                    ++count;

                    col[count].i   = ex + 1;
                    col[count].j   = ey;
                    col[count].loc = ELEMENT;
                    col[count].c   = 0;
                    valA[count]    = - mu * 1.0 / (hx * hx);
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
    PetscErrorCode AttachStokesNullspace(DM dmSol, Mat *A)
    {
        DM           dmPressure;
        Vec          constantPressure, basis;
        PetscReal    nrm;
        MatNullSpace matNullSpace;
        PetscFunctionBeginUser;
        PetscCall(DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &dmPressure));
        PetscCall(DMGetGlobalVector(dmPressure, &constantPressure));
        PetscCall(VecSet(constantPressure, 1.0));
        PetscCall(VecNorm(constantPressure, NORM_2, &nrm));
        PetscCall(VecScale(constantPressure, 1.0 / nrm));
        PetscCall(DMCreateGlobalVector(dmSol, &basis));
        PetscCall(DMStagMigrateVec(dmPressure, constantPressure, dmSol, basis));
        PetscCall(MatNullSpaceCreate(PetscObjectComm((PetscObject)dmSol), PETSC_FALSE, 1, &basis, &matNullSpace));
        PetscCall(VecDestroy(&basis));
        PetscCall(MatSetNullSpace(*A, matNullSpace));
        PetscCall(MatNullSpaceDestroy(&matNullSpace));
        PetscCall(DMRestoreGlobalVector(dmPressure, &constantPressure));
        PetscCall(DMDestroy(&dmPressure));
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    PetscErrorCode CheckDMStag(SystemParameters sys){
        PetscFunctionBeginUser;
        /* This is a testing/experiments function that creates a DMStag object, checks DMCoarsen, creates the Steady State Stokes
         system and solves it with a direct solver to check the correctness of the code.  */
        DM dm0, dmc, dmcc;
        Mat A, Ac, Acc;
        Vec rhs, rhsc, rhscc;
        PetscBool pinPressure;
        KSP ksp;
        Vec sol;

        const PetscInt dof0 = 0, dof1 = 1, dof2 = 1; /* 1 dof on each edge and element center */
        const PetscInt stencilWidth = 1;
        PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, sys->Nx, sys->Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm0));
        PetscCall(DMSetFromOptions(dm0));
        PetscCall(DMSetUp(dm0));
        PetscCall(DMStagSetUniformCoordinatesProduct(dm0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));
        pinPressure = PETSC_FALSE; // If this is set to False we need to attach a constant nullspace with the DM
        PetscCall(CreateSteadyStokesSystem(dm0, &A, &rhs, pinPressure, sys));
        dumpMatToFile(A,"A_dmstag");
        dumpVecToFile(rhs,"rhsDM");
        PetscViewer     viewer;
        DMCoarsen(dm0, PETSC_COMM_WORLD, &dmc);
        DMCoarsen(dmc, PETSC_COMM_WORLD, &dmcc);
        PetscInt startx, starty, startz, hx, hy, hz, nExtrax, nExtray, nExtraz;
        PetscInt Nx, Ny;
        DMStagGetCorners(dm0, &startx, &starty, &startz, &hx, &hy,&hz, &nExtrax, &nExtray, &nExtraz);
        DMStagGetGlobalSizes(dm0, &Nx,&Ny, PETSC_NULLPTR);
        DMStagGetGlobalSizes(dmc, &Nx,&Ny, PETSC_NULLPTR);
        DMStagGetGlobalSizes(dmcc, &Nx,&Ny, PETSC_NULLPTR);
        PetscCall(CreateSteadyStokesSystem(dmc, &Ac, &rhsc, pinPressure, sys));
        dumpMatToFile(Ac,"A_dmstagc");
        dumpVecToFile(rhsc,"rhsDMc");
        PetscCall(CreateSteadyStokesSystem(dmcc, &Acc, &rhscc, pinPressure, sys));
        dumpMatToFile(Acc,"A_dmstagcc");
        dumpVecToFile(rhscc,"rhsDMcc");

        //PetscCall(AttachPressureNullspace(dm0, &A));
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetType(ksp,KSPMINRES);
        PC pc;
        KSPSetOperators(ksp,A,A);
        KSPSetDM(ksp,dm0);
        KSPSetDMActive(ksp,PETSC_FALSE);
        KSPSetUp(ksp);

        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCNONE);

        PetscCall(DMCreateGlobalVector(dm0, &sol));
        KSPSolve(ksp,rhs,sol);
        {
            KSPConvergedReason reason;
            PetscCall(KSPGetConvergedReason(ksp, &reason));
            //PetscCheck(reason >= 0, PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED, "Linear solve failed");
        }
        dumpVecToFile(sol,"dmsol");
        PetscFunctionReturn(0);
    }
}
