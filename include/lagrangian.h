//
// Created by gkluhana on 22/04/24.
//

#ifndef IBIMPLICIT_LAGRANGIAN_H
#define IBIMPLICIT_LAGRANGIAN_H
#include <utilities.h>
#include <vector>
namespace ibImplicit{
    PetscErrorCode DumpLagrangianToFile(Vec X,  std::string var_name);
    class LagrangianGrid{
    public:
        LagrangianGrid(){
            //Default Constructor
        }
        LagrangianGrid(SystemParameters sys, Vec *pX){
            d_springData = sys->springData;
            d_X = *pX;
            common_init();
        }
        LagrangianGrid(SpringParameters springData, Vec *pX): d_springData(springData){
            d_X = *pX;
            common_init();
        }
        PetscErrorCode AssembleForceDensityOp();
        PetscErrorCode ApplyForceDensity(Vec *X, Vec *KX);
        PetscErrorCode initialize();
        PetscErrorCode clean();
        Mat getK(){
            return d_K;
        }
        Vec getX(){
            return d_X;
        }
        Vec getXlocal(){
            return Xlocal;
        }
    private:
        SpringParameters d_springData;
        Mat d_K;
        Vec d_X,Xlocal;
        PetscBool K_assembled;
        void common_init(){
            K_assembled = PETSC_FALSE;
        }
    };

    PetscErrorCode CheckForceDensityOp(Vec *pX, SpringParameters springData);
}
#endif //IBIMPLICIT_LAGRANGIAN_H
