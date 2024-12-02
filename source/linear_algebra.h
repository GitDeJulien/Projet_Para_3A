#ifndef _LINEAR_ALGEBRA_H
#define _LINEAR_ALGEBRA_H

#include <string>
#include <iostream>
#include <vector>

#include "data.h"
#include "functions.h"
#include "space_scheme.h"
#include "matrix.h"

class LinearAlgebra {

    private: 

    Data* _data;
    SpaceScheme* _ssch;

    static double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            result += a[i] * b[i];
        return result;
    }

    static double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }

    public:

    LinearAlgebra(SpaceScheme* ssch, Data* data);

    //Resolve the linear systeme AX=b with LU decomposition and backward, forward methods
    std::vector<double> LU(const Matrix& A, const std::vector<double> b);


    //Resolve the linear systeme directly with the Laplacian matrix-vector product
    std::vector<double> Lap_BiCGStab(const std::vector<double> b, int maxIterations, double tol);

};

#endif