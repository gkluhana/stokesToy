//
// Created by gkluhana on 25/03/24.
//
#include "../include/utilities.h"

namespace ibImplicit{
    bool debug = false;
    PetscErrorCode getTruncatedMGLevels(PetscInt powerOfTwo, PetscInt *result) {
        PetscFunctionBegin;

        if (powerOfTwo <= 0 || (powerOfTwo & (powerOfTwo - 1)) != 0) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Input must be a power of 2.");
        }

        PetscInt n = 0;
        while (powerOfTwo >>= 1) {
            ++n;
        }

        *result = n-3;
        if(*result==0) *result=1;

        PetscFunctionReturn(0);
    }
    PetscInt RoundUpToNextPowerOf2(PetscInt n) {
        if (n <= 0) return 1;

        --n;

        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;

        return ++n;
    }
    PetscInt RoundDownToPowerOf2(PetscInt n) {
        if (n == 0) return 0;

        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;

        return n - (n >>1);
    }

    PetscErrorCode SystemParametersData::clean(){
        PetscFunctionBeginUser;
        for (int i=0;i<mg_params->n_levels;i++){
            PetscCall(DMDestroy(&mg_params->levels[i]->dm_level));
            PetscCall(DMDestroy(&mg_params->levels[i]->dm_level));
            PetscCall(PetscFree(mg_params->levels[i]));
        }

        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode dumpISToFile(ISLocalToGlobalMapping is, std::string str){
        PetscFunctionBeginUser;
        PetscCall(PetscInfo(is,"Printed IS: %s\n",str.c_str()));
        std::ostringstream is_filename;
        is_filename << "./MatrixDump/"<< str;
        PetscViewer matlab_viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, is_filename.str().c_str(),  &matlab_viewer);
        PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_ASCII_MATLAB);
        ISLocalToGlobalMappingView(is, NULL);
        PetscCall(PetscViewerDestroy(&matlab_viewer));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode dumpVecToFile(Vec v, std::string str){
        PetscFunctionBeginUser;
        PetscCall(PetscInfo(v,"Printed Vector: %s\n",str.c_str()));
        std::ostringstream v_filename;
        v_filename << "./MatrixDump/"<< str;
        PetscViewer matlab_viewer;
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, v_filename.str().c_str(), FILE_MODE_WRITE, &matlab_viewer));
        PetscCall(PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_BINARY_MATLAB));
        PetscCall(VecView(v, matlab_viewer));
        PetscCall(PetscViewerDestroy(&matlab_viewer));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode dumpMatToFile(Mat M, std::string str){
        PetscFunctionBeginUser;
        PetscCall(PetscInfo(M,"Printed Matrix : %s\n",str.c_str()));
        std::ostringstream v_filename;
        v_filename << "./MatrixDump/"<< str;
        PetscViewer matlab_viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, v_filename.str().c_str(), FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_BINARY_MATLAB);
        PetscCall(MatView(M, matlab_viewer));
        PetscCall(PetscViewerDestroy(&matlab_viewer));
        PetscFunctionReturn(0);
    }
    PetscErrorCode dumpMatNestToFile(Mat M, std::string str){
        PetscFunctionBeginUser;
        PetscCall(PetscInfo(M,"Printed Matrix : %s\n",str.c_str()));
        std::ostringstream v_filename;
        v_filename << "./MatrixDump/"<< str;
        PetscViewer ascii_viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, v_filename.str().c_str(),  &ascii_viewer);
        PetscViewerPushFormat(ascii_viewer, PETSC_VIEWER_ASCII_MATLAB);
        PetscCall(MatView(M, ascii_viewer));
        PetscCall(PetscViewerDestroy(&ascii_viewer));
        PetscFunctionReturn(0);
    }
}