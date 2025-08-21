#include "matrix.hpp"
#include <random>
#include <thread>


Matrix& Matrix::operator*(const Matrix& rhs){
    if (this->cols == rhs.rows) {
        Matrix* result = new Matrix(this->rows, rhs.cols);
        for(int i = 0; i < this->rows; i++) {
            for (int k = 0; k < rhs.cols; k++) {
                for(int j = 0; j < rhs.rows; j++) {
                    result->matrix[i][k] += this->matrix[i][j] * rhs.matrix[j][k];
                }
            }
        }
        return result;
    } else {
        std::cerr << "Error: Inner dimensions do not match.\n";
        return NULL;
    }
}

Matrix::~Matrix(){
    for (int i = 0; i < rows; i++){
        delete [] matrix[i];
    }

    delete [] matrix;
}

void Matrix::setDataRand() const {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            this->matrix[i][j] = rand()%9 + 1; //Random integer between 1-10
        }
    }
}

void Matrix::setDataFile(const char* file) const {
    std::ifstream myFile(file);
    while(!myFile.eof()){
        int i = 0;
        int j = 0;
        myFile >> matrix[i][j++];
        if (j == cols) {
            i++;
            j = 0;
        }
    }

}

double Matrix::determinant(const Matrix& A){
    if (A.rows == 2){
        return (A.matrix[0][0] * A.matrix[1][1]) - (A.matrix[1][0] * A.matrix[0][1]);
    } else {
        int i = 0;
        double det = 0;
        for(int j = 0; j < A.cols; j++){
            det += pow(-1, i+j+2) * A.matrix[i][j] * determinant(subMatrix(A,i,j));
        }
        return det;
    }
}

Matrix Matrix::inverse(const Matrix &A){
    //Need to fix for integers and floats or treat solely as floats
    Matrix inverse(A);
    Matrix identity('I',A.rows);
    if (determinant(A) != 0) {
        for (int i = 0; i < A.rows; i++){
            for(int j = 0; j < A.rows; j++){
                if(inverse.matrix[i][0] != 0){
                    inverse.matrix[i][j] = inverse.matrix[i][j]/inverse.matrix[i][0]; //only works for floats
                }
            }
        }
        //A -> identity using Gauss-Jordan Elimination
        //Do same row operations to identity as A
    }
    return inverse;
}

void Matrix::matrixToI(const Matrix& A){
    Matrix tmp(A.rows);
    Matrix identity('I',A.rows);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < rows; j++){
            tmp.matrix[i][j] = identity.matrix[i][j] - A.matrix[i][j];
        }
    }
    tmp.print();
}

Matrix Matrix::cofactor(const Matrix& A){
    Matrix cofactorM(A.rows);
    for(int i = 0; i < A.rows; i++){
        for(int j = 0; j < A.rows; j++){
            cofactorM.matrix[i][j] = pow(-1,i+j) * determinant(subMatrix(A,i,j));
        }
    }

    return cofactorM;
}

Matrix Matrix::transpose(const Matrix& A){
    Matrix transposeM(A.rows);
    for(int i = 0; i < A.rows; i++){
        for(int j = 0; j < A.rows; j++){
            transposeM.matrix[i][j] = A.matrix[j][i];
        }
    }

    return transposeM;
}

const Matrix Matrix::subMatrix(const Matrix& A, int i, int j){
    //Only call if n x n and n >= 3;
    Matrix sM(A.rows - 1);
    for(int m = 0; m < A.rows; m++){
        for(int n = 0; n < A.rows; n++){
            if (m != i && n != j){
                int row = m,col = n;
                if (n > j){
                    col = n - 1;
                }
                if (m > i){
                    row = m - 1;
                }
                sM.matrix[row][col] = A.matrix[m][n];
            }
        }
    }

    return sM;
}

void Matrix::subMatrix(Matrix& sM, const Matrix& A, const int i, const int j){
    // Same as above but useful if subMatrix matrix should be usable and not const
    for(int m = 0; m < A.rows; m++){
        for(int n = 0; n < A.rows; n++){
            if (m != i && n != j){
                int row = m,col = n;
                if (n > j){
                    col = n - 1;
                }
                if (m > i){
                    row = m - 1;
                }
                sM.matrix[row][col] = A.matrix[m][n];
            }
        }
    }
}

void Matrix::print() {
    using std::cout;
    cout << "[\n";
    for (int i = 0; i < rows; i++) {
        cout << "\t";
        for (int j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "]\n";
}

/*
 *TODO:
    IMPLEMENT:
        Implement Gauss-Jordan method for inverse function
        Define row operations
 */
