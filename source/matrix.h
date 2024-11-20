#ifndef _MATRIX_H
#define _MATRIX_H

#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>

class Matrix {

    public:

        std::vector<double> data;
        int rows, cols;

        // Constructor prototype
        Matrix(int r, int c);

        //Matrix(int r, int c) : rows(r), cols(c), data(r * c, 0) {};

        double& operator()(int i, int j) { return data[i * cols + j]; };
        const double& operator()(int i, int j) const { return data[i * cols + j]; };

        // Method to multiply a matrix with a vector
        std::vector<double> MatrixVectorProduct(const std::vector<double>& x) const;

        // Method to multiply two matrices
        Matrix MatrixMatrixProduct(const Matrix& B) const;

        // Method to add two matrices
        Matrix AddMatrix(const Matrix& B) const;

        // Method to multiply matrix by a scalar value
        Matrix ScalarMultiply(double lambda) const;

        // Static method to generate an identity matrix
        static Matrix Identity(int size);


};


#endif