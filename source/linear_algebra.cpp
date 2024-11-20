#ifndef _LINEAR_ALGEBRA_CPP
#define _LINEAR_ALGEBRA_CPP

#include "linear_algebra.h"


LinearAlgebra::LinearAlgebra(){};

std::vector<double> LinearAlgebra::LU(const Matrix& A, const std::vector<double> b){

    int N = b.size();

    Matrix L(N, N);
    Matrix U(N, N);

    U = A;

    std::vector<double> X;
    X.resize(N);

    if (A.cols != N || A.cols != A.rows) {
        std::cerr << "Error: dimensions don't agree in LU function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        L(i,i) = 1.;
    }

    //building L and U via gauss elimination
    for (int k = 0; k < N; ++k) {
        for (int i = k+1; i < N; i++) {
            double multiplier = U(i,k)/U(k,k);
            L(i,k) = multiplier;
            U(i,k) = 0.;
            for (int j = k+1; j < N; j++) {
                U(i,j) -= multiplier*U(k,j);
            }
        }
    }

    //backward solving Ux=y
    for (int i = 1; i <= N; ++i) {
        double s = 0.;
        for (int j = 1; j < i; ++j){
            s += U(N-i,N-j)*X[N-j];
        }
        X[N-i] = (b[N-i] - s)/U(N-i,N-i);
    }

    //forward solving Ly = b
    for (int i = 0; i < N; ++i) {
        double s = 0.;
        for (int j = 0; j < i; ++j)
            s += L(i,j)*X[j];
        X[i] = (b[i] - s)/L(i,i);
    }

    return X;

}

std::vector<double> LinearAlgebra::BiCGStab(const Matrix& A, const std::vector<double> b, int maxIterations = 1000, double tol = 1e-6){
    
    int N = b.size();

    std::vector<double> X(N, 0.0);             // Initial guess (zero vector)
    std::vector<double> r = b;                 // Residual vector r = b - A * X (initially, b)
    std::vector<double> r0 = r;                // Copy of the initial residual
    std::vector<double> p(N, 0.0), v(N, 0.0), s(N, 0.0), t(N, 0.0), h(N,0.0);
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_prev = 0.0, beta = 0.0;
    
    double normB = norm(b);
    if (normB < tol) return X; // If b is small enough, return X = 0
    
    p = r;           
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        rho_prev = rho;
        rho = dot(r0, r);
        beta = (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < N; ++i){
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        if (std::abs(rho_prev) < tol) break;       // Breakdown check
        
        v = A.MatrixVectorProduct(p);
        alpha = rho / dot(r0, v);
        
        for (int i = 0; i < N; ++i){
            h[i] = X[i] + alpha * p[i];
            s[i] = r[i] - alpha * v[i];
        }
        
        if (norm(h) < tol) {
            for (int i = 0; i < N; ++i)
                X[i] += alpha * p[i];
            break;
        }
        
        t = A.MatrixVectorProduct(s);
        omega = dot(t, s) / dot(t, t);
        
        for (int i = 0; i < N; ++i)
            X[i] = h[i] + omega * s[i];
        
        for (int i = 0; i < N; ++i)
            r[i] = s[i] - omega * t[i];
        
        //normB = norm(b);
        if (norm(r) < tol*normB) break;
        
    }

    return X;
}


#endif