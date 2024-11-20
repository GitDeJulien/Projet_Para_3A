#ifndef _TIME_SCHEME_CPP
#define _TIME_SCHEME_CPP


#include "time_scheme.h"

TimeScheme::TimeScheme(Data* data, LinearAlgebra* lin, Function* fct, SpaceScheme* ssch) :_data(data), _lin(lin), _fct(fct), _ssch(ssch)
{
    _dt = data->Get_dt();
    _key_TimeScheme = data->Get_TimeScheme();

}

std::vector<double> TimeScheme::EulerExplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn)
{
    if(A.cols != static_cast<int>(Un.size()) || Un.size() != bn.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);

    Unp1 = ((A.ScalarMultiply(-this->_dt)).AddMatrix(I)).MatrixVectorProduct(Un);

    for(int i=0; i<N; ++i) {
        Unp1[i] += this->_dt*bn[i];
    }

    return Unp1;

}

std::vector<double> TimeScheme::EulerImplicite(const Matrix& A, const std::vector<double> Un , const std::vector<double> bnp1)
{

    if(A.cols != static_cast<int>(Un.size()) || Un.size() != bnp1.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);
    Matrix A_star = (A.ScalarMultiply(this->_dt)).AddMatrix(I);
    std::vector<double> U_star;
    U_star.resize(N);

    for(int i=0; i<N; ++i) {
        U_star[i] = Un[i] + this->_dt*bnp1[i];
    }

    //Unp1 = _lin->LU(A_star, U_star);
    Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);

    return Unp1;
}

std::vector<double> TimeScheme::CranckNicholson(const Matrix& A, const std::vector<double> Un , const std::vector<double> bn, const std::vector<double> bnp1){

    if(A.cols != static_cast<int>(Un.size()) || Un.size() != bn.size() || Un.size() != bnp1.size()) {
        std::cerr << "Error: dimensions don't agree in Euler Explicit function!" << std::endl;
        exit(EXIT_FAILURE);
    }    

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix
    Matrix I = Matrix::Identity(N);
    Matrix A_star = (A.ScalarMultiply(this->_dt/2.)).AddMatrix(I);
    std::vector<double> U_star;
    U_star.resize(N);

    U_star = ((A.ScalarMultiply(-this->_dt/2.)).AddMatrix(I)).MatrixVectorProduct(Un);

    for(int i=0; i<N; ++i) {
        U_star[i] += this->_dt/2.*(bnp1[i]+bn[i]);
    }

    Unp1 = _lin->BiCGStab(A_star, U_star, 10000, 1e-6);

    return Unp1;
}

std::vector<double> TimeScheme::Advance(const Matrix& A, const std::vector<double> Un, const double tn) 
{
    int N = Un.size();
    std::vector<double> Unp1;
    std::vector<double> Sn;
    std::vector<double> Snp1;
    Unp1.resize(N);
    double tnp1 = tn + _data->Get_dt();

    switch (_data->Get_TimeScheme())
        {
        case 1:
            Sn.resize(N);
            Sn = _ssch->SourceTerme(_data, _fct, tn);
            Unp1 = this->EulerExplicite(A, Un, Sn);
            break;

        case 2:
            Snp1.resize(N);
            Snp1 = _ssch->SourceTerme(_data, _fct, tnp1);
            Unp1 = this->EulerImplicite(A, Un, Snp1);
            break;
            
        case 3:
            Sn.resize(N);
            Snp1.resize(N);
            Sn = _ssch->SourceTerme(_data, _fct, tn);
            Snp1 = _ssch->SourceTerme(_data, _fct, tnp1);
            Unp1 = this->CranckNicholson(A, Un, Sn, Snp1);
            break;
        
        default:
            std::cout << "The time scheme key " << _data->Get_TimeScheme() << " is not valid !" << std::endl;
            exit(EXIT_FAILURE);
    }

    return Unp1;

}


void TimeScheme::SaveSol(const std::vector<double>& sol, const std::string& path, int n) 
{

    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        // Try to create the directory
        if (mkdir(path.c_str(), 0777) == -1) {
            std::cout << "Directory already exist." << std::endl;
        } else {
            std::cout << "Directory " << path << " created successfully." << std::endl;
        }
    }
    std::string n_file = path + "/sol_" + std::to_string(n) + ".vtk";
    std::ofstream solution;
    solution.open(n_file);
    if (solution.is_open()){
        solution << "# vtk DataFile Version 3.0" << std::endl;
        // solution << "#key time scheme: " << _data->Get_TimeScheme() << std::endl;
        // solution << "#key space scheme: " << _data->Get_SpaceScheme() << std::endl;
        // solution << "#key left/right boundary condition: " << _data->Get_key_LeftRightBoundCond() << std::endl;
        // solution << "#key up/down boundary condition: " << _data->Get_key_UpDownBoundCond() << std::endl;
        // solution << "#key initial condition: " << _data->Get_key_InitialCondition() << std::endl;
        // solution << "#key source terme: " << _data->Get_key_SourceTerme() << std::endl;
        solution << "sol" << std::endl;
        solution << "ASCII" << std::endl;
        solution << "DATASET STRUCTURED_POINTS" << std::endl;
        solution << "DIMENSIONS " << _data->Get_Nx() << " " << _data->Get_Ny() << " " << 1 << std::endl;
        solution << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        solution << "SPACING " << _data->Get_hx() << " " << _data->Get_hy() << " " << 1 << std::endl;;
        solution << "POINT_DATA " << _data->Get_Nx()*_data->Get_Ny() << std::endl;
        solution << "SCALARS sol float" << std::endl;
        solution << "LOOKUP_TABLE default" << std::endl;
        for(int j=0; j<_data->Get_Ny(); ++j)
        {
            for(int i=0; i<_data->Get_Nx(); ++i)
            {
                solution << sol[i+j*_data->Get_Nx()] << " ";
            }
            solution << std::endl;
        }
        solution.close();
    }
    else {
        std::cout << "Error opening results file!" << std::endl;
    }
}




#endif