#ifndef _SPACE_SCHEME_CPP
#define _SPACE_SCHEME_CPP

#include <iostream>
#include <utility>

#include "space_scheme.h"

SpaceScheme::SpaceScheme(){}

int SpaceScheme::index_MatToVect(Data* data, const int i, const int j)
{
    int Nx(0);
    int l(0);

    Nx = data->Get_Nx();

    l = j*(Nx+1) + i;

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
    int Nx(0), Ny(0);
    int l(0);
    double x(0.0), y(0.0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = (Nx+1)*(Ny+1);
    U0.resize(N);

    for(int i=0; i<=Nx; ++i){
        for(int j=0; j<=Ny; ++j){
            l = index_MatToVect(data, i, j);

            x = i*data->Get_hx();
            y = j*data->Get_hy();

            U0[l] = function->InitialCondition(data, x, y);
        }
    }


    return U0;
}

// Matrix SpaceScheme::BuildMatrix(Data* data) 
// {
//     std::pair<int, int> indMat;

//     int Nx(0);
//     int Ny(0);

//     double alpha(0.0);
//     double beta(0.0);
//     double gamma(0.0);

//     Nx = data->Get_Nx();
//     Ny = data->Get_Ny();

//     int N = Nx*Ny;

//     alpha = data->Get_diffusion_coeff()*(2/(data->Get_hx()*data->Get_hx()) + 2/(data->Get_hy()*data->Get_hy()));
//     beta = - data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
//     gamma = - data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());

//     Matrix matrix(N,N);

//     for (int I=0; I<N; ++I) {
//         for (int J=0; J<N; ++J){
//             matrix(I,J) = 0.0;
//         }
//     }

//     if(data->Get_SpaceScheme() == 1){ //Laplacian centered discretisation
//         for(int i=1; i<=Nx; ++i){
//             for(int j=1; j<=Ny; ++j){
//                 int l = this->index_MatToVect(data, i, j);

//                 matrix(l,l) = alpha;
//                 if (i>1) {
//                     matrix(l,l-1) = beta;
//                 }
//                 if (i<Nx){
//                     matrix(l,l+1) = beta;
//                 }
//                 if (j>1) {
//                     matrix(l,l-Nx) = gamma;
//                 }
//                 if (j<Ny) {
//                     //matrix(l,l+Nx) = gamma;
//                     matrix(l,l+Nx) = gamma;
//                 }
//             }
//         }
//     }
//     else {
//         std::cerr << "The Space Scheme key" << data->Get_SpaceScheme() << "is not define. Please change it in the 'input/data.dat' file" << std::endl;
//         exit(EXIT_FAILURE);
//     }
//     return matrix;     
// }

std::vector<double> SpaceScheme::Lap_MatVectProduct(Data* data, std::vector<double> U, int expl_impl)
{

    /*
    This function compute direclty the matrix-vector product
    --------------- U_res = (I+/-dt*A)*U ---------------------

    And return the resulting vector U_res
    expl_impl == 1 explicite time scheme
    expl_impl == 2 implicit time scheme
    
    - This is usefull to not stock the full matrix A
    - It work with Euler Implicite/Explicite and Crank-Nicholson scheme
    */

    int Nx(0);
    int Ny(0);
    double dt(0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();
    dt = data->Get_dt();

    double alpha(0.0);
    double beta(0.0);
    double gamma(0.0);

    if (expl_impl == 1) {
        alpha = -dt*data->Get_diffusion_coeff()*(2/(data->Get_hx()*data->Get_hx()) + 2/(data->Get_hy()*data->Get_hy()));
        beta = dt*data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
        gamma = dt*data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());
    }

    else if (expl_impl == 2) {
        alpha = dt*data->Get_diffusion_coeff()*(2/(data->Get_hx()*data->Get_hx()) + 2/(data->Get_hy()*data->Get_hy()));
        beta = -dt*data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx());
        gamma = -dt*data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy());
    }

    else {
        std::cout << "No other value for expl_impl than 1 or 2" << std::endl;
    }

    std::vector<double> U_res;
    U_res.resize((Nx+1)*(Ny+1));

    if(data->Get_SpaceScheme() == 1){ //Laplacian centered discretisation
        for (int i=0; i<=Nx; i++) {
            U_res[i] = U[i];
            U_res[i+Ny*(Nx+1)] = U[i+Ny*(Nx+1)];
        }
        for(int j=1; j<Ny; j++){ //probleme d'indice
            U_res[j*(Nx+1)] = U[j*(Nx+1)];
            U_res[Nx+j*(Nx+1)] = U[Nx+j*(Nx+1)];

            for(int ij=1+j*(Nx+1); ij<Nx+j*(Nx+1); ij++){
                //int l = this->index_MatToVect(data, i, j);
                //int l = (j-1) * Nx + (i-1);
                //Different case for condition Derichlet or Robin
                if(data->Get_key_Schwarz_Bounds() == 1){ //Dirichlet homogene

                    U_res[ij] = 1.0 + alpha * U[ij] + beta*U[ij-1] + beta*U[ij+1] + gamma*U[ij-(Nx+1)] + gamma*U[ij+Nx+1];
                    



                    // if (ij>1+j*(Nx+1)) U_res[ij] += beta*U[ij-1];

                    // if (ij<Nx-1+j*(Nx+1)) U_res[ij] += beta*U[ij+1];

                    // if (j>1) U_res[ij] += gamma*U[ij-(Nx+1)];

                    // if (j<Ny-1) U_res[ij] += gamma*U[ij+Nx+1];
                }
            }
        }
    }

    return U_res;

}


std::vector<double> SpaceScheme::SourceTerme(Data* data, Function* function, const double t){

    std::vector<double> S;
    int Nx(0), Ny(0);

    Nx = data->Get_Nx();
    Ny = data->Get_Ny();

    int N = (Nx+1)*(Ny+1);
    S.resize(N);

    // for(int i=1; i<=Nx; ++i){
    //     for(int j=1; j<=Ny; ++j){
    //         l = index_MatToVect(data, i, j);

    //         x = i*data->Get_hx();
    //         y = j*data->Get_hy();

    //         S[l] = function->SourceFunction(data, x, y, t);

    //         if (i==1) { 
    //             //Derichlet non homogène
    //             S[l] += data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx()) * function->BoundaryCondition_Left(data, x-data->Get_hx(), y);
    //         }
    //         if (i==Nx) { 
    //             //Derichlet non homogène
    //             S[l] += data->Get_diffusion_coeff()/(data->Get_hx()*data->Get_hx()) * function->BoundaryCondition_Right(data, x+data->Get_hx(), y);
    //         }
    //         if (j==1) {
    //             //Derichlet non homogène
    //             S[l] += data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy()) * function->BoundaryCondition_Down(data, x, y-data->Get_hy());
    //         }
    //         if (j==Ny) {
    //             //Derichlet non homogène
    //             S[l] += data->Get_diffusion_coeff()/(data->Get_hy()*data->Get_hy()) * function->BoundaryCondition_Up(data, x, y+ data->Get_hy());
    //         }
    //     }
    // }

    for (int i=0; i<=Nx; i++) {
        double x =i*data->Get_hx();
        S[i] = function->BoundaryCondition_Down(data, x, 0.0);
        S[i+Ny*(Nx+1)] = function->BoundaryCondition_Up(data, x, data->Get_Ly());
    }
    for(int j=1; j<Ny; j++){ //probleme d'indice
        double y = j*data->Get_hy();
        S[j*(Nx+1)] = function->BoundaryCondition_Left(data, 0.0, y);
        S[Nx+j*(Nx+1)] = function->BoundaryCondition_Left(data, data->Get_Lx(), y);

        for(int ij=1+j*(Nx+1); ij<Nx+j*(Nx+1); ij++){
            //Different case for condition Derichlet or Robin
            double x =(ij-j*(Nx-1))*data->Get_hx();

            if(data->Get_key_Schwarz_Bounds() == 1){ //Dirichlet homogene

                S[ij] = function->SourceFunction(data, x, y, t);



                // if (ij>1+j*(Nx+1)) U_res[ij] += beta*U[ij-1];

                // if (ij<Nx-1+j*(Nx+1)) U_res[ij] += beta*U[ij+1];

                // if (j>1) U_res[ij] += gamma*U[ij-(Nx+1)];

                // if (j<Ny-1) U_res[ij] += gamma*U[ij+Nx+1];
            }
        }
    }

    return S;

}


#endif
