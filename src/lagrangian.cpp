//
// Created by gkluhana on 23/04/24.
//
#include "../include/lagrangian.h"
namespace ibImplicit{
    PetscErrorCode DumpLagrangianToFile(Vec X,  std::string var_name){
        PetscScalar *coords;
        PetscInt numPoints;
        PetscFunctionBeginUser;
        // Extract data from the vector
        PetscCall(VecGetArray(X, &coords));
        PetscCall(VecGetSize(X,&numPoints));
        numPoints /=2;
        // Open the file for writing
        std::string filename = "./dump/"+ var_name + ".vtu";
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            PetscPrintf(PETSC_COMM_WORLD, "Error opening file for writing.\n");
            PetscFinalize();
            return -1;
        }

        // Write the VTK header
        outFile << "<?xml version=\"1.0\"?>\n";
        outFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        outFile << "  <UnstructuredGrid>\n";
        outFile << "    <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"0\">\n";

        // Write the points
        outFile << "      <Points>\n";
        outFile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (PetscInt i = 0; i < numPoints; ++i) {
            outFile << coords[2*i] << " " << coords[2*i+1] << " 0.0\n";
        }
        outFile << "        </DataArray>\n";
        outFile << "      </Points>\n";

        // Write the cells (empty in this example)
        outFile << "      <Cells>\n";
        outFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";
        outFile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";
        outFile << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\"/>\n";
        outFile << "      </Cells>\n";

        // Close the VTK file
        outFile << "    </Piece>\n";
        outFile << "  </UnstructuredGrid>\n";
        outFile << "</VTKFile>\n";

        outFile.close();

        // Clean up
        PetscCall(VecRestoreArray(X, &coords));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode LagrangianGrid::AssembleForceDensityOp(){
        PetscMPIInt rank;

        PetscFunctionBeginUser;
        if(!K_assembled) {
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


            std::vector<int> master_idx = d_springData->master_idx;
            std::vector<int> slave_idx = d_springData->slave_idx;
            std::vector<float> kappa = d_springData->kappa;
            std::vector<float> rest_length = d_springData->rest_length;
            PetscScalar kappa_uni = -1;
            if (d_springData->uniform_spring_stiffness) {
                kappa_uni = d_springData->uniform_spring_stiffness;
            }
            PetscInt ne = d_springData->num_edges;
            PetscInt ns;
            PetscCall(VecGetSize(d_X, &ns));
            ns /= 2;
            PetscCall(MatCreate(PETSC_COMM_WORLD, &d_K));
            PetscCall(MatSetSizes(d_K, PETSC_DECIDE, PETSC_DECIDE, 2 * ns, 2 * ns));
            PetscCall(MatSetUp(d_K));
            // Ensure Mat K is initialized
            // Ensure Mat K size is (2*ns)^2
            // Ensure Vec X contains current structural configuration
            PetscInt m_idx, s_idx;
            PetscInt m_indices[2], s_indices[2];
            PetscScalar D[2], M[2], S[2], dF_dX[4], kappa_e;
            PetscReal R, T, dT_dR;


            for (int e = 0; e < ne; e++) {
                m_idx = 2 * master_idx[e];
                s_idx = 2 * slave_idx[e];
                m_indices[0] = m_idx;
                m_indices[1] = m_idx + 1;
                s_indices[0] = s_idx;
                s_indices[1] = s_idx + 1;

                PetscCall(VecGetValues(Xlocal, 2, m_indices, M));
                PetscCall(VecGetValues(Xlocal, 2, s_indices, S));

                D[0] = S[0] - M[0];
                D[1] = S[1] - M[1];
                R = std::sqrt(D[0] * D[0] + D[1] * D[1]);
                if (kappa_uni != -1) { kappa_e = kappa_uni; }
                else { kappa_e = kappa[e]; }
                T = kappa_e * (R - rest_length[e]);
                dT_dR = kappa_e;
                for (int j = 0; j < 2; j++) {
                    for (int i = 0; i < 2; i++) {
                        if (i == j) {
                            dF_dX[i + j * 2] = (T / R) +
                                               (dT_dR - T / R) * D[i] * D[j] / (R * R);
                        } else {
                            dF_dX[i + j * 2] = (dT_dR - T / R) * D[i] * D[j] / (R * R);
                        }
                        if (dF_dX[i + j * 2] < 1e-12) dF_dX[i + j * 2] = 0.0;
                    }
                }
                if (!rank) {
                    PetscCall(MatSetValues(d_K, 2, s_indices, 2, m_indices, dF_dX, ADD_VALUES));
                    PetscCall(MatSetValues(d_K, 2, m_indices, 2, s_indices, dF_dX, ADD_VALUES));
                }
                for (int i = 0; i < 4; i++) {
                    dF_dX[i] *= -1.0;
                }
                if (!rank) {
                    PetscCall(MatSetValues(d_K, 2, s_indices, 2, s_indices, dF_dX, ADD_VALUES));
                    PetscCall(MatSetValues(d_K, 2, m_indices, 2, m_indices, dF_dX, ADD_VALUES));
                }
            }
            PetscCall(MatAssemblyBegin(d_K, MAT_FINAL_ASSEMBLY));
            PetscCall(MatAssemblyEnd(d_K, MAT_FINAL_ASSEMBLY));
            K_assembled = PETSC_TRUE;
        }
        else{
            if(debug) std::cerr << "Lagrangian Force Density Operator Already Assembled" << std::endl;
        }
        PetscFunctionReturn(0);
    }
    PetscErrorCode LagrangianGrid::clean() {
        PetscFunctionBeginUser;
        PetscCall(MatDestroy(&d_K));
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode LagrangianGrid::initialize(){
	PetscMPIInt rank;
        PetscFunctionBeginUser;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//std::cout << "Scattering Lagrangian Vector, rank:" << rank << std::endl;
        VecScatter scatter_ctx;
        PetscCall(VecScatterCreateToAll(d_X, &scatter_ctx, &Xlocal));
        VecScatterBegin(scatter_ctx, d_X, Xlocal, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter_ctx, d_X, Xlocal, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter_ctx);
	//std::cout << "Completed Scattering " << std::endl;
        AssembleForceDensityOp();
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    PetscErrorCode LagrangianGrid::ApplyForceDensity(Vec *pX, Vec *pF){
        VecScatter scatter_ctx;
        Vec Xlocal;
        PetscFunctionBeginUser;
        PetscCall(VecScatterCreateToAll(d_X,&scatter_ctx,&Xlocal));
        std::vector<int> master_idx = d_springData->master_idx;
        std::vector<int> slave_idx= d_springData->slave_idx;
        std::vector<float> kappa= d_springData->kappa;
        std::vector<float> rest_length= d_springData->rest_length;
        PetscScalar kappa_uni = -1;
        if(d_springData->uniform_spring_stiffness){
            kappa_uni = d_springData->uniform_spring_stiffness;
        }
        PetscInt ne = d_springData->num_edges;
        PetscInt ns;
        PetscCall(VecGetSize(*pX,&ns));
        ns /=2;
        Vec X = *pX;
        PetscCall(VecCreate(PETSC_COMM_WORLD,pF));
        Vec F = *pF;
        PetscCall(VecSetSizes(F,PETSC_DECIDE,2*ns));
        PetscCall(VecSetType(F,VECSTANDARD));
        PetscCall(VecZeroEntries(F));
        // Ensure Vec X contains current structural configuration
        PetscInt m_idx, s_idx;
        PetscInt m_indices[2], s_indices[2];
        PetscScalar D[2], M[2],S[2], dF_dX[4], kappa_e;
        PetscReal R, T, dT_dR;
        PetscScalar Xs[2],Xm[2];
        for (int e =0; e<ne; e++) {
            m_idx = 2 * master_idx[e];
            s_idx = 2 * slave_idx[e];
            m_indices[0] = m_idx;
            m_indices[1] = m_idx + 1;
            s_indices[0] = s_idx;
            s_indices[1] = s_idx + 1;

            PetscCall(VecGetValues(Xlocal, 2, m_indices, M));
            PetscCall(VecGetValues(Xlocal, 2, s_indices, S));

            D[0] = S[0] - M[0];
            D[1] = S[1] - M[1];
            R = std::sqrt(D[0] * D[0] + D[1] * D[1]);
            if (kappa_uni != -1){kappa_e = kappa_uni;}
            else {kappa_e = kappa[e];}
            T = kappa_e * (R - rest_length[e]);
            dT_dR = kappa_e;
            for (int j = 0; j < 2; j++){
                for (int i = 0; i < 2; i++) {
                    if(i==j) {
                        dF_dX[i + j * 2] = (T / R) +
                                           (dT_dR - T / R) * D[i] * D[j] / (R * R);
                    }
                    else {
                        dF_dX[i + j * 2] =(dT_dR - T / R) * D[i] * D[j] / (R * R);
                    }
                    if(dF_dX[i+ j*2] < 1e-12) dF_dX[i+j*2] = 0.0;
                }
            }
            // Off-Diagonal Effects Zero, i.e. dimensions x,y decoupled
            PetscCall(VecSetValue(F,s_indices[0],dF_dX[0]*M[0],ADD_VALUES));
            PetscCall(VecSetValue(F,s_indices[1],dF_dX[3]*M[1],ADD_VALUES));
            PetscCall(VecSetValue(F,m_indices[0],dF_dX[0]*S[0],ADD_VALUES));
            PetscCall(VecSetValue(F,m_indices[1],dF_dX[3]*S[1],ADD_VALUES));
            for (int i =0; i<4 ; i++){
                dF_dX[i] *= -1.0;
            }
            // Diagonal
            PetscCall(VecSetValue(F,s_indices[0],dF_dX[0]*S[0],ADD_VALUES));
            PetscCall(VecSetValue(F,s_indices[1],dF_dX[3]*S[1],ADD_VALUES));
            PetscCall(VecSetValue(F,m_indices[0],dF_dX[0]*M[0],ADD_VALUES));
            PetscCall(VecSetValue(F,m_indices[1],dF_dX[3]*M[1],ADD_VALUES));

        }
        PetscCall(VecAssemblyBegin(F));
        PetscCall(VecAssemblyEnd(F));
        PetscFunctionReturn(0);
    }
    PetscErrorCode CheckForceDensityOp(Vec *pX, SpringParameters springData){
        LagrangianGrid lgrid(springData, pX);
        lgrid.AssembleForceDensityOp();
        Mat K = lgrid.getK();
        ibImplicit::dumpMatToFile(K,"K_MAT");
        PetscFunctionReturn(0);
    }
}
