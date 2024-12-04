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

    int _Nx, _Ny, _N_pts;
    double _hx, _hy;
    double _D, _dt;

    public:

        //Constructor
        SpaceScheme(Data* data);

        //Transforme matrix index (i,j) in vector index l
        int index_MatToVect(Data* data, const int i, const int j);

        //Transforme matrix index (i,j) in vector index l
        std::pair<int, int> index_VectToMat(Data* data, const int l);

        //Initialize the solution vector at t=0.0
        std::vector<double> Initialize(Data* data, Function* function);

        //Building the general Matrix M depending on the scheme
        // Matrix BuildMatrix(Data* data);

        //Mat Vect product with the Laplacian matrix
        std::vector<double> Lap_MatVectProduct(Data* data, std::vector<double> U);

        //Building the Source terme S depending on the border conditions
        std::vector<double> SourceTerme(Data* data, Function* function, const double t);

};





#endif