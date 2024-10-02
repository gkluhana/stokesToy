//
// Created by gkluhana on 25/03/24.
//
#ifndef IBIMPLICIT_INITIALIZE_H
#define IBIMPLICIT_INITIALIZE_H
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmstag.h>

#include <staggered.h>
#include <lagrangian.h>
#include <multigrid.h>
#include <stokes.h>

namespace ibImplicit {


    /// Initialize the Problem parameters with Input File
    /// \param argc
    /// \param argv
    /// \param data
    PetscErrorCode initializeProblem(int argc, char **argv, std::unordered_map <std::string, std::string> *data);

    /// Parse the Input File and populate the problem parameters data
    /// \param filename
    /// \param data
    PetscErrorCode parseFile(const std::string &filename, std::unordered_map <std::string, std::string> *data);


    /// Convert string to a BCType
    /// \param[in] str
    /// \return BC Type
    BCType getBCTypeFromString(const char *str);

    /// Populate the fields of a SystemParameters with the data map
    /// \param[in,out] sys
    /// \param[in] data
    /// \return
    PetscErrorCode createSystemParameters(SystemParameters *sys, std::unordered_map<std::string,std::string>* data);

    PetscErrorCode populateBCStaggeredType(BCStaggeredType  bc_type, std::unordered_map<std::string,std::string>* data);
    PetscErrorCode populateBCFunctions(BCFunctions2D  bc_funcs, std::unordered_map<std::string,std::string>* data);
}
#endif