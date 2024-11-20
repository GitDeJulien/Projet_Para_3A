#ifndef _TIME_SCHEME_H
#define _TIME_SCHEME_H

#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <fstream>
#include <errno.h>
#include <cstring>

#include "data.h"
#include "space_scheme.h"
#include "linear_algebra.h"
#include "matrix.h"


class TimeScheme {

    private:

        Data* _data;
        LinearAlgebra* _lin;
        Function* _fct;
        SpaceScheme* _ssch;

        int _key_TimeScheme;
        double _dt;

    public:

        TimeScheme(Data* data, LinearAlgebra* lin, Function* fct, SpaceScheme* ssch);

        std::vector<double> EulerExplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn);

        std::vector<double> EulerImplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bnp1);

        std::vector<double> CranckNicholson(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn, const std::vector<double> bnp1);

        std::vector<double> Advance(const Matrix& A, const std::vector<double> Un, const double tn);

        void SaveSol(const std::vector<double>& sol, const std::string& path, int n);

};



#endif