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

    double* _sol0; 
    double* _Matrix;

    public:

        //Constructor
        SpaceScheme();

        //Transforme matrix index (i,j) in vector index l
        int index_MatToVect(Data* data, const int i, const int j);

        //Transforme matrix index (i,j) in vector index l
        std::pair<int, int> index_VectToMat(Data* data, const int l);

        //Initialize the solution vector at t=0.0
        std::vector<double> Initialize(Data* data, Function* function);

        //Building the general Matrix M depending on the scheme
        Matrix BuildMatrix(Data* data);

        //Building the Source terme S depending on the border conditions
        std::vector<double> SourceTerme(Data* data, Function* function, const double t);

};





#endif