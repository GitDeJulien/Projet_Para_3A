#ifndef _FUNCTION_CPP

#include <cmath>

#include "functions.h"


Function::Function()
{

};


double Function::InitialCondition(const Data* data, const double x, const double y) const 
{

    if (data->Get_key_InitialCondition() == 1) return 100*ExactSolution(data, x, y);
    else if (data->Get_key_InitialCondition() == 2) return cos(x);
    else if (data->Get_key_InitialCondition() == 3) return 0.0;
    else {
        std::cerr << "Error: The initial condition key " << data->Get_key_InitialCondition() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

}

double Function::SourceFunction(const Data* data, const double x, const double y, const double t) const 
{
    if (data->Get_key_SourceTerme() == 1) return 2*(x-(x*x)+y-(y*y));
    else if (data->Get_key_SourceTerme() == 2) return sin(x)+cos(y);
    else if (data->Get_key_SourceTerme() == 3) return exp(-(x-data->Get_Lx()/2.)*(x-data->Get_Lx()/2.))*exp(-(y-data->Get_Ly()/2.)*(y-data->Get_Ly()/2.))*cos(M_PI/2.*t);
    else {
        std::cerr << "Error: The source terme key " << data->Get_key_SourceTerme() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::ExactSolution(const Data* data, const double x, const double y) const 
{
    if (data->Get_key_SourceTerme() == 1 && data->Get_key_RightBoundCond() == 1 && data->Get_key_LeftBoundCond() == 1 && data->Get_key_DownBoundCond() == 1 && data->Get_key_UpBoundCond() == 1){
        return x * (1-x) * y*(1-y);
    }
    else if (data->Get_key_SourceTerme() == 2 && data->Get_key_RightBoundCond() == 2 && data->Get_key_LeftBoundCond() == 2 && data->Get_key_UpBoundCond() == 2 && data->Get_key_DownBoundCond() == 2){
        return sin(x) + cos(y);
    }
    else {
        std::cout << "Exacte solution havn't been determined analyticaly" << std::endl;
        return 0.0;
    }
}

double Function::BoundaryCondition_Left(const Data* data, const double x, const double y) const 
{
    if (data->Get_key_LeftBoundCond() == 1) return 0.0;
    else if (data->Get_key_LeftBoundCond() == 2) return cos(x) + sin(y);
    else if (data->Get_key_LeftBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Left boundary condition key " << data->Get_key_LeftBoundCond() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

}

double Function::BoundaryCondition_Right(const Data* data, const double x, const double y) const 
{
    if (data->Get_key_RightBoundCond() == 1) return 0.0;
    else if (data->Get_key_RightBoundCond() == 2) return cos(x) + sin(y);
    else if (data->Get_key_RightBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Right boundary condition key" << data->Get_key_RightBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::BoundaryCondition_Up(const Data* data, const double x, const double y) const 
{
    if (data->Get_key_UpBoundCond() == 1) return 0.0;
    else if (data->Get_key_UpBoundCond() == 2) return cos(x) + sin(y);
    else if (data->Get_key_UpBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Up boundary condition key" << data->Get_key_UpBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::BoundaryCondition_Down(const Data* data, const double x, const double y) const 
{
    if (data->Get_key_DownBoundCond() == 1) return 0.0;
    else if (data->Get_key_DownBoundCond() == 2) return cos(x) + sin(y);
    else if (data->Get_key_DownBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Down boundary condition key" << data->Get_key_DownBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}


#define _FUNCTION_CPP
#endif