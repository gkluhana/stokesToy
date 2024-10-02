//
// Created by gkluhana on 25/03/24.
//
#include "../include/initialize.h"

namespace ibImplicit{
    const double PI = 3.14159265358979323846;
    void fe(double x, double y,double mu, double hx, double hy, double *f){
        x = 2*x-1;
        y = 2*y-1;
        *f = 20*x*y*y*y;
    }

    void fw(double x, double y,double mu, double hx, double hy, double* f) {
        x = 2*x-1;
        y = 2*y-1;
        *f = 20*x*y*y*y;
    }

    void fs(double x, double y,double mu, double hx, double hy, double* f) {
        x = 2*x-1;
        y = 2*y-1;
        *f = 20*x*y*y*y;
    }

    void fn(double x, double y,double mu, double hx, double hy, double* f) {
        x = 2*x-1;
        y = 2*y-1;
        *f = 20*x*y*y*y;
        //*f = mu*(0.5 - 0.5*std::cos(2*PI*x));
    }

    void ge(double x, double y,double mu, double hx, double hy, double *g){
        x = 2*x-1;
        y = 2*y-1;
        *g = 5*(x*x*x*x - y*y*y*y) ;
    }

    void gw(double x, double y,double mu, double hx, double hy, double* g) {
        x = 2*x-1;
        y = 2*y-1;
        *g = 5*(x*x*x*x - y*y*y*y) ;
    }

    void gs(double x, double y,double mu, double hx, double hy, double* g) {
        x = 2*x-1;
        y = 2*y-1;
        *g = 5*(x*x*x*x - y*y*y*y) ;
    }

    void gn(double x, double y,double mu, double hx, double hy, double* g) {
        x = 2*x-1;
        y = 2*y-1;
        *g = 5*(x*x*x*x - y*y*y*y) ;
    }
    SCType getSCTypeFromString(const std::string& str) {
        if (str == "EXACT") {
            return EXACT;
        } else if (str == "LSC") {
            return LSC;
        } else if (str == "TRUNCATED_SC") {
            return TRUNCATED_SC;
        } else {
            // Handle the error case: unknown string
            std::cerr << "Unknown SCType: " << str << std::endl;
            // Return a default value or handle the error appropriately
            return LSC; // Or any other default/error value
        }
    }
    BCType getBCTypeFromString(const std::string& str) {
        if (str == "DIRICHLET") {
            return DIRICHLET;
        } else if (str == "NEUMANN") {
            return NEUMANN;
        } else if (str == "TRUNCATED_BC") {
            return TRUNCATED_BC;
        } else {
            // Handle the error case: unknown string
            std::cerr << "Unknown BCType: " << str << std::endl;
            // Return a default value or handle the error appropriately
            return DIRICHLET; // Or any other default/error value
        }
    }
    PetscErrorCode populateBCStaggeredType(BCStaggeredType  bc_type, std::unordered_map<std::string,std::string>* data){
        PetscFunctionBeginUser;
        for (const auto& entry : *data) {
            // Velocity Boundaries
            if (entry.first == "bc_vel_type_west")
                bc_type->u_w = (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_vel_type_east")
                bc_type->u_e= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_vel_type_south")
                bc_type->u_s= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_vel_type_north")
                bc_type->u_n= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            // Pressure Boundaries
            else if (entry.first == "bc_p_type_north")
                bc_type->p_n= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_p_type_east")
                bc_type->p_e= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_p_type_west")
                bc_type->p_w= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
            else if (entry.first == "bc_p_type_south")
                bc_type->p_s= (entry.second == "dirichlet") ? DIRICHLET : NEUMANN;
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode populateBCFunctions(BCFunctions2D  bc_funcs, std::unordered_map<std::string,std::string>* data){
        PetscFunctionBeginUser;
       // This is hardcoded for now
        bc_funcs->fe = fe;
        bc_funcs->fw = fw;
        bc_funcs->fs = fs;
        bc_funcs->fn = fn;

        bc_funcs->ge = ge;
        bc_funcs->gw = gw;
        bc_funcs->gs = gs;
        bc_funcs->gn = gn;
        PetscFunctionReturn(0);
    }
    PetscErrorCode initializeProblem(int argc, char** argv, std::unordered_map<std::string, std::string>* data){
        PetscFunctionBeginUser;
        if (argc <= 1) {
            throw ibError("Invalid command line arguments. Usage: " + std::string(argv[0]) + " <filename>");
        }

        std::string filename = argv[1];

        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        if (rank == 0) {
            PetscCall(PetscInfo(NULL,"Reading Input from filename: %s \n",filename.c_str()));
        }

        PetscCall(parseFile(filename, data));

        // Display parsed data
        for (const auto& pair : *data) {
            //if (rank == 0) std::cout << pair.first.c_str() << ": " << pair.second.c_str() << std::endl;
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode createSystemParameters(SystemParameters *sys, std::unordered_map<std::string,std::string>* data){
        PetscFunctionBeginUser;
        PetscCall(PetscMalloc1(1, sys));
        //TODO : This function needs cleaning up
        SystemParameters SYS=*sys;
        // Physical Parameters
        SYS->lx = std::stod(data->at("lx"));
        SYS->ly = std::stod(data->at("ly"));
        SYS->Lx = std::stod(data->at("Lx"));
        SYS->Ly = std::stod(data->at("Ly"));
        SYS->mu = std::stod(data->at("mu"));
        SYS->rho = std::stod(data->at("rho"));

        // Discretization Parameters
        SYS->Nx = std::stoi(data->at("Nx"));
        SYS->Ny = std::stoi(data->at("Ny"));
        SYS->hx = ((*sys)->Lx - (*sys)->lx) / (*sys)->Nx;
        SYS->hy = ((*sys)->Ly - (*sys)->ly)/ (*sys)->Ny;
        SYS->dt = std::stod(data->at("dt"));
        SYS->nfluid= (*sys)->Ny * ((*sys)->Nx+1) + (*sys)->Nx * ((*sys)->Ny + 1);
        SYS->np= (*sys)->Nx * (*sys)->Ny;
        SYS->nu = ((*sys)->Nx+1) * (*sys)->Ny;
        SYS->nv = ((*sys)->Ny+1) * (*sys)->Nx;



        //Solver Parameters
        SYS->mg_params = new MGParametersData;
        SYS->mg_params->n_levels = std::stoi(data->at("mg_levels"));
        SYS->mg_params->level_curr= std::stoi(data->at("mg_levels"))-1;
        SYS->mg_params->RefRatio= std::stoi(data->at("mg_RefRatio"));
        SYS->mg_params->levels = new LevelCtx[SYS->mg_params->n_levels];

        SYS->pre_block_type = std::stoi(data->at("pre_block_type"));
        SYS->truncated_sys = 0;
        std::string *s1_type_string = &data->at("S1_type");
        SYS->s1_type = getSCTypeFromString(*s1_type_string);

        PetscCall(PetscMalloc1(SYS->mg_params->n_levels, &SYS->mg_params->levels));
        for (PetscInt i = 0; i <SYS->mg_params->n_levels; ++i) PetscCall(LevelCtxCreate(&SYS->mg_params->levels[i]));
        PetscFunctionReturn(0);
    }

    // Functions to parse input files
    PetscErrorCode parseFile(const std::string& filename, std::unordered_map<std::string, std::string>* data) {
        PetscFunctionBeginUser;
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open()) {
            throw ibError("Error opening file: " + filename);
        }
        std::unordered_map<std::string, std::string> d_data;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss >> std::ws, key, '=')) {
                if (std::getline(iss >> std::ws, value)) {
                    // Trim leading and trailing whitespace
                    key.erase(0, key.find_first_not_of(" \t\r\n"));
                    key.erase(key.find_last_not_of(" \t\r\n") + 1);
                    value.erase(0, value.find_first_not_of(" \t\r\n"));
                    value.erase(value.find_last_not_of(" \t\r\n") + 1);
                    d_data[key] = value;
                }
            }
        }

        *data = d_data;
        PetscFunctionReturn(0);
}
    PetscErrorCode readVertexFile(const std::string& filename, int& ns, Vec*X) {
        PetscFunctionBeginUser;
        Vec X0 = *X;
        std::string vertex_filename = filename;
        vertex_filename.append(".vertex");
        std::ifstream file(vertex_filename);
        if (!file.is_open()) {
            std::cerr << "Error: Failed to open file " << vertex_filename<< std::endl;
            return false;
        }

        // Read ns from the first line
        if (!(file >> ns)) {
            std::cerr << "Error: Failed to read ns from file " << filename << std::endl;
            return false;
        }
        PetscCall(VecSetSizes(X0,PETSC_DECIDE,2*ns));
        PetscCall(VecSetType(X0,VECSTANDARD));

        // Create an array to store x and y coordinates
        PetscScalar *coordinates;
        PetscMalloc1(2 * ns, &coordinates);

        // Read x and y coordinates from the file and store them in the array
        for (int i = 0; i < ns; ++i) {
            if (!(file >> coordinates[2*i] >> coordinates[2*i + 1])) {
                std::cerr << "Error: Failed to read coordinates from file " << filename << std::endl;
                PetscFree(coordinates);
                return false;
            }

            // Set the values of the Petsc Vec X
            PetscCall(VecSetValue(X0, 2*i,coordinates[2*i] , INSERT_VALUES));
            PetscCall(VecSetValue(X0, 2*i+1,coordinates[2*i+1] , INSERT_VALUES));
        }

        PetscCall(VecAssemblyBegin(X0));
        PetscCall(VecAssemblyEnd(X0));

        // Free the memory allocated for the array
        PetscFree(coordinates);

        // Close the file
        file.close();
        ibImplicit::dumpVecToFile(X0,"X0");

        return true;
    }
    PetscErrorCode readSpringFile(const std::string& filename, SpringParameters springData) {
        PetscFunctionBeginUser;
        // TODO: Adjust code so that vectors are Petsc Vecs

        std::string shell_spring_filename = filename;
        shell_spring_filename.append(".spring");
        std::ifstream file(shell_spring_filename);

        if (!file.is_open()) {
            std::cerr << "Error: Failed to open file " << shell_spring_filename << std::endl;
            PetscFunctionReturn(1);
        }

        // Read the number of edges (first line)
        file >> springData->num_edges;

        // Resize vectors to store data
        springData->master_idx.resize(springData->num_edges);
        springData->slave_idx.resize(springData->num_edges);
        springData->kappa.resize(springData->num_edges);
        springData->rest_length.resize(springData->num_edges);

        // Read the rest of the file
        for (int i = 0; i < springData->num_edges; ++i) {
            // Read the current line
            file >> springData->master_idx[i] >> springData->slave_idx[i] >> springData->kappa[i] >> springData->rest_length[i];
        }

        // Close the file
        file.close();

        PetscFunctionReturn(0);
    }
    void SystemParametersData::recomputeGridParameters(){
        //Call this function when changing any of the system parameters
        PetscReal wx = Lx-lx;
        PetscReal wy = Ly-ly;
        hx  = wx / Nx;
        hy = wy / Ny;
        nfluid= Ny * (Nx+1) + Nx * (Ny + 1);
        np= Nx * Ny;
        nu = (Nx+1) * Ny;
        nv = (Ny+1) * Nx;
        free(bcSpec->bc_idx);
        bcSpec->bc_idx = new BCIndicesData;
        getBdryIndxLex(bcSpec->bc_idx,Nx,Ny);
    }
    void SystemParametersData::recomputeGridParametersTruncated(PetscReal lx_f, PetscReal ly_f){
        //Call this function when creating the domain for the truncation preconditioner


        PetscReal wx = Lx-lx;
        PetscReal wy = Ly-ly;

        // Add margin on both sides of the structure
        lx -= 0.5 * wx ;
        Lx += 0.5 * wx ;
        ly -= 0.5 * wy ;
        Ly += 0.5 * wy ;

        // Align lx,ly with the closest cell of full discretization
        lx = std::floor( (lx - lx_f) / hx) * hx;
        ly = std::floor( (ly - ly_f) / hy) * hy;

        // Check the width
        wx = Lx-lx;
        wy = Ly-ly;
        PetscReal w_max = std::max(wx,wy);

        // Compute the max number of cells in each direction
        Nx =  w_max/hx;
        Ny =  w_max/hy;

        // Expand Nx,Ny to a power of 2
        Nx = RoundUpToNextPowerOf2(Nx);
        Ny = RoundUpToNextPowerOf2(Ny);

        // Expand the box
        Lx = lx + Nx*hx;
        Ly = ly + Ny*hy;

        // Compute multigrid levels
        getTruncatedMGLevels(Nx,&mg_params->n_levels);

        // Update remaining parameters
        nfluid= Ny * (Nx+1) + Nx * (Ny + 1);
        np= Nx * Ny;
        nu = (Nx+1) * Ny;
        nv = (Ny+1) * Nx;
        free(bcSpec->bc_idx);
        bcSpec->bc_idx = new BCIndicesData;
        getBdryIndxLex(bcSpec->bc_idx,Nx,Ny);
    }

}
