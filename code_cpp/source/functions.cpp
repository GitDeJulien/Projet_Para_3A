#ifndef _FUNCTION_CPP

#include <cmath>

#include "functions.h"


Function::Function(Data* data): _data(data)
{

};


double Function::InitialCondition(const double x, const double y) const 
{

    if (_data->Get_key_InitialCondition() == 1) return 20*ExactSolution(x, y);
    else if (_data->Get_key_InitialCondition() == 2) return cos(x);
    else if (_data->Get_key_InitialCondition() == 3) return 0.0;
    else {
        std::cerr << "Error: The initial condition key " << _data->Get_key_InitialCondition() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

}

double Function::SourceFunction(const double x, const double y, const double t) const 
{
    if (_data->Get_key_SourceTerme() == 1) return 2*(x-(x*x)+y-(y*y));
    else if (_data->Get_key_SourceTerme() == 2) return sin(x)+cos(y);
    else if (_data->Get_key_SourceTerme() == 3) return exp(-(x-_data->Get_Lx()/2.)*(x-_data->Get_Lx()/2.))*exp(-(y-_data->Get_Ly()/2.)*(y-_data->Get_Ly()/2.))*cos(M_PI/2.*t);
    else {
        std::cerr << "Error: The source terme key " << _data->Get_key_SourceTerme() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::ExactSolution(const double x, const double y) const 
{
    if (_data->Get_key_SourceTerme() == 1 && _data->Get_key_RightBoundCond() == 1 && _data->Get_key_LeftBoundCond() == 1 && _data->Get_key_DownBoundCond() == 1 && _data->Get_key_UpBoundCond() == 1){
        return x*(1-x) * y*(1-y);
    }
    else if (_data->Get_key_SourceTerme() == 2 && _data->Get_key_RightBoundCond() == 2 && _data->Get_key_LeftBoundCond() == 2 && _data->Get_key_UpBoundCond() == 2 && _data->Get_key_DownBoundCond() == 2){
        return sin(x) + cos(y);
    }
    else {
        std::cout << "Exacte solution havn't been determined analyticaly" << std::endl;
        return 0.0;
    }
}

double Function::BoundaryCondition_Left(const double x, const double y) const 
{
    if (_data->Get_key_LeftBoundCond() == 1) return 0.0;
    else if (_data->Get_key_LeftBoundCond() == 2) return cos(x) + sin(y);
    else if (_data->Get_key_LeftBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Left boundary condition key " << _data->Get_key_LeftBoundCond() << " is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }

}

double Function::BoundaryCondition_Right(const double x, const double y) const 
{
    if (_data->Get_key_RightBoundCond() == 1) return 0.0;
    else if (_data->Get_key_RightBoundCond() == 2) return cos(x) + sin(y);
    else if (_data->Get_key_RightBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Right boundary condition key" << _data->Get_key_RightBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::BoundaryCondition_Up(const double x, const double y) const 
{
    if (_data->Get_key_UpBoundCond() == 1) return 0.0;
    else if (_data->Get_key_UpBoundCond() == 2) return cos(x) + sin(y);
    else if (_data->Get_key_UpBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Up boundary condition key" << _data->Get_key_UpBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}

double Function::BoundaryCondition_Down(const double x, const double y) const 
{
    if (_data->Get_key_DownBoundCond() == 1) return 0.0;
    else if (_data->Get_key_DownBoundCond() == 2) return cos(x) + sin(y);
    else if (_data->Get_key_DownBoundCond() == 3) return 1.0;
    else {
        std::cerr << "Error: The Down boundary condition key" << _data->Get_key_DownBoundCond() << "is not referenced" << std::endl;
        exit (EXIT_FAILURE);
    }
}


#define _FUNCTION_CPP
#endif