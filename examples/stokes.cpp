static char help[] = "Solve a toy 2D problem on a staggered grid\n\n";
#include <initialize.h>
 int main(int argc, char **argv) {
         ibImplicit::debug = true;

         // Initialize System Parameters
         PetscCall(PetscInitialize(&argc, &argv, (char *)"PetscOptions.dat", help));
         std::unordered_map<std::string, std::string>* data = new std::unordered_map<std::string,std::string>();
         PetscCall(ibImplicit::initializeProblem(argc,argv,data));
         ibImplicit::SystemParameters sys = new ibImplicit::SystemParametersData ;
         ibImplicit::createSystemParameters(&sys, data);

         //Initialize Eulerian Grids and Data Structures
         ibImplicit::StokesSolver stokes_solver(sys);
         stokes_solver.initialize();


         //Solve the system
         stokes_solver.solve();

         //Call clean-up functions
         PetscCall(stokes_solver.clean());
         PetscCall(sys->clean());
         PetscFree(sys);

         PetscCall(PetscFinalize());
         return 0;
}




