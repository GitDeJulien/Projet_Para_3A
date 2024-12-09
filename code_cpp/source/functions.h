#ifndef _FUNCTION_H

#include <fstream>
#include <iostream>
#include <cmath>

#include "data.h"

class Function {

    private:

    Data* _data;

    //     double _Lx;
    //     double _Ly;
    //     int _key_InitialCondition;
    //     int _key_SourceTerme;
    //     int _key_LeftRightBoundCond;
    //     int _key_UpDownBoundCond;



    public:
        //Constructor
        Function(Data* data);

        // Initial condition
        double InitialCondition(const double x, const double y) const;

        // Source terme
        double SourceFunction(const double x, const double y, const double t) const;

        // Exacte solution if it known (uusefull for validation)
        double ExactSolution(const double x, const double y) const;

        //Boundary condition
        double BoundaryCondition_Right(const double x, const double y) const;

        double BoundaryCondition_Left(const double x, const double y) const;

        double BoundaryCondition_Up(const double x, const double y) const;

        double BoundaryCondition_Down(const double x, const double y) const;


        


};

#define _FUNCTION_H
#endif