#ifndef _SPACE_SCHEME_H
#define _SPACE_SCHEME_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>

#include "data.h"
#include "functions.h"
#include "matrix.h"

class SpaceScheme {

    private:

    Data* _data;
    Function* _fct;

    int _Nx, _Ny, _N_pts;
    double _hx, _hy;
    double _D, _dt;

    public:

        //Constructor
        SpaceScheme(Data* data, Function* fct);

        //Transforme matrix index (i,j) in vector index l
        int index_MatToVect(const int i, const int j);

        //Transforme matrix index (i,j) in vector index l
        std::pair<int, int> index_VectToMat(const int l);

        //Initialize the solution vector at t=0.0
        std::vector<double> Initialize();

        //Building the general Matrix M depending on the scheme
        Matrix BuildMatrix(Data* data);

        //Mat Vect product with the Laplacian matrix
        std::vector<double> Lap_MatVectProduct(std::vector<double> U);

        //Building the Source terme S depending on the border conditions
        std::vector<double> SourceTerme(const double t, std::vector<double> U);

};





#endif