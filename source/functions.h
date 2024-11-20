#ifndef _FUNCTION_H

#include <fstream>
#include <iostream>
#include <cmath>

#include "data.h"

class Function {

    // private:


    //     double _Lx;
    //     double _Ly;
    //     int _key_InitialCondition;
    //     int _key_SourceTerme;
    //     int _key_LeftRightBoundCond;
    //     int _key_UpDownBoundCond;



    public:
        //Constructor
        Function();

        // Initial condition
        double InitialCondition(const Data* data, const double x, const double y) const;

        // Source terme
        double SourceFunction(const Data* data, const double x, const double y, const double t) const;

        // Exacte solution if it known (uusefull for validation)
        double ExactSolution(const Data* data, const double x, const double y) const;

        //Boundary condition
        double BoundaryCondition_Right(const Data* data, const double x, const double y) const;

        double BoundaryCondition_Left(const Data* data, const double x, const double y) const;

        double BoundaryCondition_Up(const Data* data, const double x, const double y) const;

        double BoundaryCondition_Down(const Data* data, const double x, const double y) const;


        


};

#define _FUNCTION_H
#endif