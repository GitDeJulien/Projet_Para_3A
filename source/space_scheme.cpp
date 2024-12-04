#ifndef _SPACE_SCHEME_CPP
#define _SPACE_SCHEME_CPP

#include <iostream>
#include <utility>

#include "space_scheme.h"

SpaceScheme::SpaceScheme(Data* data): _data(data) {
    this->_Nx = data->Get_Nx();
    this->_Ny = data->Get_Ny();
    this->_N_pts = data->Get_N_pts();

    this->_D = data->Get_diffusion_coeff();
    this->_dt = data->Get_dt();
    this->_hx = data->Get_hx();
    this->_hy = data->Get_hy();
}

int SpaceScheme::index_MatToVect(Data* data, const int i, const int j)
{
    int Nx(0);
    int l(0);

    Nx = data->Get_Nx();

    l = (j-1)*Nx + (i-1);

    return(l);
}

std::pair<int, int> SpaceScheme::index_VectToMat(Data* data, const int l) {

    int Ny(0);
    int i(0);
    int j(0);

    Ny = data->Get_Ny();

    i = l/Ny+1;
    j = l%Ny+1;

    return {i,j};
}

std::vector<double> SpaceScheme::Initialize(Data* data, Function* function)
{
    std::vector<double> U0;
    double x(0.0), y(0.0);

    int N = _N_pts;
    U0.resize(N);

    for(int i=1; i<=_Nx; ++i){
        for(int j=1; j<=_Ny; ++j){
            int l = this->index_MatToVect(data, i, j);

            x = i*data->Get_hx();
            y = j*data->Get_hy();

            U0[l] = function->InitialCondition(data, x, y);
        }
    }


    return U0;
}

Matrix SpaceScheme::BuildMatrix(Data* data) 
{
    std::pair<int, int> indMat;

    int Nx(0);
    int Ny(0);

    double alpha(0.0);
    double beta(0.0);
    double gamma(0.0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = Nx*Ny;

    alpha = data->Get_diffusion_coeff()*(2/(data->Get_hx()*data->Get_hx()) + 2/(data->Get_hy()*data->Get_hy()));
    beta = - data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
    gamma = - data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());

    Matrix matrix(N,N);

    for (int I=0; I<N; ++I) {
        for (int J=0; J<N; ++J){
            matrix(I,J) = 0.0;
        }
    }

    if(data->Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        for(int i=1; i<=Nx; ++i){
            for(int j=1; j<=Ny; ++j){
                int l = this->index_MatToVect(data, i, j);

                matrix(l,l) = alpha;
                if (i>1) {
                    matrix(l,l-1) = beta;
                }
                if (i<Nx){
                    matrix(l,l+1) = beta;
                }
                if (j>1) {
                    matrix(l,l-Nx) = gamma;
                }
                if (j<Ny) {
                    //matrix(l,l+Nx) = gamma;
                    matrix(l,l+Nx) = gamma;
                }
            }
        }
    }
    else {
        std::cerr << "The Space Scheme key" << data->Get_SpaceScheme() << "is not define. Please change it in the 'input/data.dat' file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return matrix;     
}

std::vector<double> SpaceScheme::Lap_MatVectProduct(Data* data, std::vector<double> U) {

    /*
    This function compute direclty the matrix-vector product
    --------------- U_res = (I+/-dt*A)*U ---------------------

    And return the resulting vector U_res
    expl_impl == 1 explicite time scheme
    expl_impl == 2 implicit time scheme
    
    - This is usefull to not stock the full matrix A
    - It work with Euler Implicite/Explicite and Crank-Nicholson scheme
    */

    
    double alpha = _D*(2/(_hx*_hx) + 2/(_hy*_hy));
    double beta = - _D/(_hx*_hx);
    double gamma = - _D/(_hy*_hy);


    std::vector<double> U_res;
    U_res.resize(this->_N_pts);
    int l(0);
    if(data->Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        for(int i=1; i<=this->_Nx; ++i){
            for(int j=1; j<=this->_Ny; ++j){
                l = this->index_MatToVect(data, i, j);

                if(data->Get_key_Schwarz_Bounds() == 1){

                    U_res[l] = alpha*U[l];

                    if (i>1) U_res[l] += beta*U[l-1];

                    if (i<_Nx) U_res[l] += beta*U[l+1];

                    if (j>1) U_res[l] += gamma*U[l-_Nx];

                    if (j<_Ny) U_res[l] += gamma*U[l+_Nx];
                }
            }
        }
    }
    else {
        std::cerr << "The Space Scheme key" << data->Get_SpaceScheme() << "is not define. Please change it in the 'input/data.dat' file" << std::endl;
        exit(EXIT_FAILURE);
    }
    return U_res;  
}


std::vector<double> SpaceScheme::SourceTerme(Data* data, Function* function, const double t){

    std::vector<double> S;
    double x(0.0), y(0.0);

    int N = this->_N_pts;
    S.resize(N);

    double beta = - _D/(_hx*_hx);
    double gamma = - _D/(_hy*_hy);

    for(int i=1; i<=this->_Nx; ++i){
        for(int j=1; j<=this->_Ny; ++j){
            int l = index_MatToVect(data, i, j);

            x = i*this->_hx;
            y = j*this->_hy;
            
            S[l] = function->SourceFunction(data, x, y, t);

            if (i==1) { 
                //Derichlet non homogène
                S[l] += -beta * function->BoundaryCondition_Left(data, x -_hx, y);
                // std::cout << S[l] << std::endl;
            }
            if (i==_Nx) { 
                //Derichlet non homogène
                S[l] += -beta * function->BoundaryCondition_Right(data, x + _hx, y);
                // std::cout << S[l] << std::endl;
            }
            if (j==1) {
                //Derichlet non homogène
                S[l] += -gamma * function->BoundaryCondition_Down(data, x, y -_hy);
                // std::cout << S[l] << std::endl;
            }
            if (j==_Ny) {
                //Derichlet non homogène
                S[l] += -gamma * function->BoundaryCondition_Up(data, x, y + _hy);
                // std::cout << S[l] << std::endl;
            }
        }
    }

    // TODO!
    // if there is other border condition than homogenus
    // need to modify the source terme in appropriate spotes...
    // make so if, else if, else case do so...

    return S;

}


#endif
