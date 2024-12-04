#ifndef _TIME_SCHEME_CPP
#define _TIME_SCHEME_CPP


#include "time_scheme.h"

TimeScheme::TimeScheme(Data* data, LinearAlgebra* lin, Function* fct, SpaceScheme* ssch) :_data(data), _lin(lin), _fct(fct), _ssch(ssch)
{
    _dt = data->Get_dt();
    _key_TimeScheme = data->Get_TimeScheme();

}

std::vector<double> TimeScheme::EulerExplicite(const std::vector<double> Un , const std::vector<double> bn)
{
    if (Un.size() != bn.size()) {
        std::cout << "Error: Un and bn have differed shape!";
    }

    int N = Un.size();

    std::vector<double> Unp1, H;
    Unp1.resize(N);
    H.resize(N);

    H = _ssch->Lap_MatVectProduct(_data, Un);

    for(int l=0; l<N; ++l) {
        Unp1[l] = Un[l] - this->_dt*H[l] + this->_dt*bn[l];
    }

    return Unp1;

}

std::vector<double> TimeScheme::EulerImplicite(const std::vector<double> Un , const std::vector<double> bnp1)
{
    if (Un.size() != bnp1.size()) {
        std::cout << "Error: Un and bn have differed shape!";
    }

    int N = Un.size();

    std::vector<double> Unp1, H;
    Unp1.resize(N);

    std::vector<double> U_star;
    U_star.resize(N);

    for(int l=0; l<N; ++l) {
        U_star[l] = Un[l] + this->_dt*bnp1[l];
    }

    Unp1 = _lin->Lap_BiCGStab(U_star, 10000, 1e-6);

    return Unp1;
}

std::vector<double> TimeScheme::CranckNicholson(const std::vector<double> Un , const std::vector<double> bn, const std::vector<double> bnp1){

    if (Un.size() != bn.size()) {
        std::cout << "Error: Un and bn have differed shape!";
    }

    int N = Un.size();

    std::vector<double> Unp1;
    Unp1.resize(N);

    //Identity matrix

    std::vector<double> U_star;
    U_star.resize(N);

    Unp1 = _ssch->Lap_MatVectProduct(_data, Un); //il manque un 1/2 à ajouter dans le produit mat-vect

    for(int l=0; l<N; ++l) {
        U_star[l] += this->_dt/2.*(bnp1[l]+bn[l]);
    }

    Unp1 = _lin->Lap_BiCGStab(U_star, 10000, 1e-6); //de même ici

    return Unp1;
}

std::vector<double> TimeScheme::Advance(const std::vector<double> Un, const double tn) 
{
    int N = Un.size();
    std::vector<double> Unp1;
    std::vector<double> Sn;
    std::vector<double> Snp1;
    Unp1.resize(N);
    double tnp1 = tn + this->_dt;

    switch (_data->Get_TimeScheme())
        {
        case 1:
            Sn.resize(N);
            Sn = _ssch->SourceTerme(_data, _fct, tn);
            Unp1 = this->EulerExplicite(Un, Sn);
            break;

        case 2:
            Snp1.resize(N);
            Snp1 = _ssch->SourceTerme(_data, _fct, tnp1);
            Unp1 = this->EulerImplicite(Un, Snp1);
            break;
            
        case 3:
            Sn.resize(N);
            Snp1.resize(N);
            Sn = _ssch->SourceTerme(_data, _fct, tn);
            Snp1 = _ssch->SourceTerme(_data, _fct, tnp1);
            Unp1 = this->CranckNicholson(Un, Sn, Snp1);
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
        solution << "sol" << std::endl;
        solution << "ASCII" << std::endl;
        solution << "DATASET STRUCTURED_POINTS" << std::endl;
        solution << "DIMENSIONS " << _data->Get_Nx() << " " << _data->Get_Ny() << " " << 1 << std::endl;
        solution << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
        solution << "SPACING " << _data->Get_hx() << " " << _data->Get_hy() << " " << 1 << std::endl;;
        solution << "POINT_DATA " << _data->Get_N_pts() << std::endl;
        solution << "SCALARS sol float" << std::endl;
        solution << "LOOKUP_TABLE default" << std::endl;
        for(int j=1; j<=_data->Get_Ny(); ++j)
        {
            for(int i=1; i<=_data->Get_Nx(); ++i)
            {
                int l = _ssch->index_MatToVect(_data, i, j);
                solution << sol[l] << " ";
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