#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <concepts>
#include <fstream>
#include <optional>
#include <random>
#include <thread>

template <std::integral T>
class Matrix{
    T** matrix;
	size_t rows;
	size_t cols;

    Matrix subMatrix(const Matrix& A, const int i, const int j)
    {
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
    };

    void subMatrix(Matrix& cM, const Matrix& A, const int i, const int j)
    {
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
                    cM.matrix[row][col] = A.matrix[m][n];
                }
            }
        }
    };


public:
    // ReSharper disable once CppNonExplicitConvertingConstructor
    Matrix(const size_t n):rows(n), cols(n) {
        srand(clock());
        //Construction of a square matrix of size n
        this->matrix = new T*[rows];
        for (size_t i = 0; i < rows; i++){
            this->matrix[i] = new T[cols];
        }
    };

    Matrix(const size_t m, const size_t n): rows(m), cols(n){
        //Construction of a matrix of size m(rows) x n(columns)
        srand(clock());
        this->matrix = new T*[rows];
        for(size_t i = 0; i < rows; i++) {
            this->matrix[i] = new T[cols];
        }
    };

    Matrix(const char I, const size_t n): rows(n), cols(n){
        //Creates nxn Identity matrix
        if (I == 'I'){
            srand(clock());
            //Construction of a square matrix of size n
            this->matrix = new T*[rows];
            for (size_t i = 0; i < rows; i++){
                this->matrix[i] = new T[cols];
            }

            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < n; j++) {
                    if (i == j){
                        matrix[i][j] = 1;
                    } else {
                        matrix[i][j] = 0;
                    }
                }
            }
        }
    };

    Matrix(const Matrix& copy): rows(copy.rows), cols(copy.cols){
        srand(clock());
        this->matrix = new T*[this->rows];
        for(size_t i = 0; i < this->rows; i++) {
            this->matrix[i] = new T[this->cols];
            for(size_t j = 0; j < this->cols; j++){
                this->matrix[i][j] = copy.matrix[i][j];
            }
        }
    };

    ~Matrix() {
        {
            for (int i = 0; i < rows; i++){
                delete [] matrix[i];
            }

            delete [] matrix;
        }
    };

    Matrix& operator=(const Matrix& copy) {
        srand(clock());
        if(this != &copy){
            this->rows = copy.rows;
            this->cols = copy.cols;
            this->matrix = new T*[this->rows];
            for(size_t i = 0; i < this->rows; i++) {
                this->matrix[i] = new T[this->cols];
                for(size_t j = 0; j < this->cols; j++){
                    this->matrix[i][j] = copy.matrix[i][j];
                }
            }
        }
        return *this;
    }

    constexpr bool operator==(const Matrix& rhs) {
        if(this->rows == rhs.rows && this->cols == rhs.cols) {
            for(size_t i = 0; i < this->rows; i++) {
                for(size_t j = 0; j < this->cols; j++) {
                    if(this->matrix[i][j] != rhs.matrix[i][j]) return false;
                }
            }
            return true;
        }
        return false;
    }

    constexpr bool operator!=(const Matrix& rhs) {
        return !(*this == rhs);
    }

    constexpr std::optional<Matrix> operator+(const Matrix& rhs){
        if(this->rows == rhs.rows && this->cols == rhs.cols){
            Matrix result(rows, cols);
            for(int i = 0; i < this->rows; i++) {
                for (int k = 0; k < this->cols; k++) {
                    result.matrix[i][k] = this->matrix[i][k] + rhs.matrix[i][k];
                }
            }
            return result;
        }
        return std::nullopt;
    };

    constexpr std::optional<Matrix> operator-(const Matrix& rhs){
        if(this->rows == rhs.rows && this->cols == rhs.cols){
            Matrix result(rows, cols);
            for(int i = 0; i < this->rows; i++) {
                for (int k = 0; k < this->cols; k++) {
                    result.matrix[i][k] = this->matrix[i][k] - rhs.matrix[i][k];
                }
            }
            return result;
        }
        return std::nullopt;
    };

    constexpr std::optional<Matrix> operator*(const Matrix& rhs){
        if (this->cols == rhs.rows) {
            Matrix result{this->rows, rhs.cols};
            for(size_t i = 0; i < this->rows; i++) {
                for (size_t k = 0; k < rhs.cols; k++) {
                    for(size_t j = 0; j < rhs.rows; j++) {
                        result.matrix[i][k] += this->matrix[i][j] * rhs.matrix[j][k];
                    }
                }
            }
            return result;
        }
        return std::nullopt;
    };

    constexpr std::optional<Matrix> operator/(const Matrix& divisor) {
        /*
         *"X = B / A"
         *AX = B
         *A^-1 * AX = A^-1 * B
         *IX = A^-1 * B
         *X = A^-1 * B
         */
        return inverse(divisor) * this;
    };

    void setDataRand() const{
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                this->matrix[i][j] = rand()%9 + 1; //Random integer between 1-10
            }
        }
    };

    void setDataFile(const char* file) const{
        std::ifstream myFile(file);
        while(!myFile.eof()){
            size_t i = 0;
            size_t j = 0;
            myFile >> matrix[i][j++];
            if (j == cols) {
                i++;
                j = 0;
            }
        }
    };

    constexpr double determinant(const Matrix& A){
        if (A.rows == 2){
            return (A.matrix[0][0] * A.matrix[1][1]) - (A.matrix[1][0] * A.matrix[0][1]);
        }
        size_t i = 0;
        double det = 0;
        for(size_t j = 0; j < A.cols; j++){
            det += pow(-1, i+j+2) * A.matrix[i][j] * determinant(subMatrix(A,i,j));
        }
        return det;
    };

    constexpr auto inverse(const Matrix& A){
        Matrix inverse(A);
        Matrix identity('I',A.rows);
        if (determinant(A) != 0) {
            for (size_t i = 0; i < A.rows; i++){
                for(size_t j = 0; j < A.rows; j++){
                    if(inverse.matrix[i][0] != 0){
                        inverse.matrix[i][j] = inverse.matrix[i][j]/inverse.matrix[i][0]; //only works for floats
                    }
                }
            }
            //TODO:
            //A -> identity using Gauss-Jordan Elimination
            //Do same row operations to identity as A
        }
        return &inverse;
    };

    constexpr auto cofactor(const Matrix& A){
        Matrix cofactorM(A.rows);
        for(size_t i = 0; i < A.rows; i++){
            for(int j = 0; j < A.rows; j++){
                cofactorM.matrix[i][j] = pow(-1,i+j) * determinant(subMatrix(A,i,j));
            }
        }

        return &cofactorM;
    };

    constexpr auto transpose(const Matrix& cofactorM) {
        Matrix transposeM(cofactorM.rows);
        for(int i = 0; i < cofactorM.rows; i++){
            for(int j = 0; j < cofactorM.rows; j++){
                transposeM.matrix[i][j] = cofactorM.matrix[j][i];
            }
        }
        return &transposeM;
    };

    constexpr auto adjoint(const Matrix& A) {
        return transpose(cofactor(A));
    }

    void matrixToI(const Matrix& A); //Make private

    constexpr bool isSquare() const {
        return rows == cols;
    }
    void print();
};

#endif
