#ifndef _MATRIX_CPP
#define _MATRIX_CPP

#include "matrix.h"


Matrix::Matrix(int r, int c) : data(r * c, 0), rows(r), cols(c) {}

std::vector<double> Matrix::MatrixVectorProduct(const std::vector<double>& x) const {
    if (cols != static_cast<int>(x.size())) {
        throw std::invalid_argument("Matrix columns must match vector size.");
    }

    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += (*this)(i, j) * x[j];
        }
    }
    return result;
}

Matrix Matrix::MatrixMatrixProduct(const Matrix& B) const {
    if (cols != B.rows) {
        throw std::invalid_argument("Number of columns of A must match number of rows of B.");
    }

    Matrix C(rows, B.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < B.cols; ++j) {
            C(i, j) = 0;
            for (int k = 0; k < cols; ++k) {
                C(i, j) += (*this)(i, k) * B(k, j);
            }
        }
    }
    return C;
}

Matrix Matrix::AddMatrix(const Matrix& B) const {
    if (rows != B.rows || cols != B.cols) {
        throw std::invalid_argument("Matrices must have the same dimensions to be added.");
    }

    Matrix C(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            C(i, j) = (*this)(i, j) + B(i, j);
        }
    }

    return C;
}

Matrix Matrix::ScalarMultiply(double lambda) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) * lambda;
        }
    }
    return result;
}

Matrix Matrix::Identity(int size) {
    Matrix identity(size, size);
    for (int i = 0; i < size; ++i) {
        identity(i, i) = 1.0;
    }
    return identity;
}


#endif