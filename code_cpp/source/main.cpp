#include <string>
#include <vector>
#include<time.h>

#include "data.h"
#include "functions.h"
#include "matrix.h"
#include "linear_algebra.h"
#include "space_scheme.h"
#include "time_scheme.h"



int main(int argc, char** argv) {

    if (argc < 2)
    {
        std::cout << "Please, enter the name of your data file." << std::endl;
        exit (EXIT_FAILURE);
    }

    //std::string filename = "./input/data.dat";
    std::string filename = argv[1];

    //Pointer to Data class
    Data* data = new Data(filename);

    //Pointer to Function class
    Function* function = new Function(data);

    //Pointer to Space Scheme class
    SpaceScheme* ssch = new SpaceScheme(data, function);

    //Pointer to Linear Algebra class
    LinearAlgebra* lin = new LinearAlgebra(data, ssch);

    //Pointer to Time Scheme class
    TimeScheme* tsch = new TimeScheme(data, lin, function, ssch);

    //Display all the parameters and conditions used for computation
    data->display_parameters();


    int Nx = data->Get_Nx();
    int Ny = data->Get_Ny();
    int N_pts = data->Get_N_pts();
    double hx = data->Get_hx();
    double hy = data->Get_hy();

    std::vector<double> Un(N_pts);
    std::vector<double> Unp1(N_pts);
    std::vector<double> Sn(N_pts);

    std::cout << "N_pts: " << N_pts << std::endl;

    //Initial solution
    //Un = ssch->Initialize();
    //Unp1 = ssch->Initialize(data, function);

    tsch->SaveSol(Un, data->Get_outputPath(), 0);

    double tn = data->Get_t0();
    double nb_iteration = data->Get_niter();
    double dt = data->Get_dt();

    std::vector<double> U_exact(N_pts);
    double error(0.0);

    for(int l=0; l<N_pts; ++l){

        int j = l/Nx+1;
        int i = l%Nx+1;
        double x = i*hx;
        double y = j*hy;

        Un[l] = function->ExactSolution(x,y);
        U_exact[l] = function->ExactSolution(x, y);
    }

    tsch->SaveSol(U_exact, "output/Exact1", 0);

    std::cout << "Starting time loop" << std::endl;
    tn += dt;
    for(int iter = 1; iter<nb_iteration; ++iter) {

        Sn = ssch->SourceTerme(tn+dt, Un);
        Unp1 = lin->Lap_BiCGStab(Sn, 10000, 1e-6);

        //Advance of a time step with the chosen time scheme
        //Unp1 = tsch->Advance(Un, tn);

        //Download result in vtk files
        tsch->SaveSol(Unp1, data->Get_outputPath(), iter);

        //Error computation
        for (int l=0; l<N_pts; ++l){
            error += pow(fabs(Unp1[l] - U_exact[l]),2);
        }

        //Update
        Un = Unp1;
        tn += dt;
        std::cout << "tn: " << tn << ", " << "error: " << 1./N_pts*sqrt(error) << std::endl;
        error = 0.0;
    }

    //Pointer deletion
    delete data, delete function, delete lin, delete ssch, delete tsch;
    return 0;
}


//############################//
//########## TEST ############//
//############################//

/* TEST PROD MAT-VECT
    Matrix matrix = ssch->BuildMatrix(data);
    std::cout << matrix.rows << matrix.cols << std::endl;
    std::vector<double> U_test(N_pts, 1.0);
    std::vector<double> U_test2(N_pts, 1.0);
    U_test = ssch->Lap_MatVectProduct(data, U_test);
    U_test2 = matrix.MatrixVectorProduct(U_test2);

    std::cout << "U_test - U_test2" << std::endl;
    for (int l=0; l<N_pts; ++l) {
        std::cout << U_test[l] - U_test2[l] << std::endl;
    }
*/

/*
TEST MATRIX CLASS

Matrix A(2, 3); // 2x3 matrix
A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;

std::vector<double> x = {1, 2, 3};
std::vector<double> result_vector = A.MatrixVectorProduct(x);

std::cout << "Matrix-Vector Product:" << std::endl;
for (double val : result_vector) {
    std::cout << val << " ";
}
std::cout << std::endl;

Matrix B(3, 2); // 3x2 matrix
B(0, 0) = 7; B(0, 1) = 8;
B(1, 0) = 9; B(1, 1) = 10;
B(2, 0) = 11; B(2, 1) = 12;

Matrix C = A.MatrixMatrixProduct(B);

std::cout << "Matrix-Matrix Product:" << std::endl;
for (int i = 0; i < C.rows; ++i) {
    for (int j = 0; j < C.cols; ++j) {
        std::cout << C(i, j) << " ";
    }
    std::cout << std::endl;
}

Matrix A(2, 2); // 2x2 matrix
A(0, 0) = 1; A(0, 1) = 2;
A(1, 0) = 3; A(1, 1) = 4;

Matrix B(2, 2); // Another 2x2 matrix
B(0, 0) = 5; B(0, 1) = 6;
B(1, 0) = 7; B(1, 1) = 8;

Matrix C = A.Add(B);

std::cout << "Matrix Sum:" << std::endl;
for (int i = 0; i < C.rows; ++i) {
    for (int j = 0; j < C.cols; ++j) {
        std::cout << C(i, j) << " ";
    }
    std::cout << std::endl;
}

*/