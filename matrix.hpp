#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <concepts>
#include <fstream>
#include <optional>
#include <random>
#include <thread>
namespace Neo {
    struct Dims {
        size_t rows;
        size_t cols;
        bool operator==(const Dims& other) const {if (rows == other.rows && cols == other.cols) return true; return false;}
        [[nodiscard]] bool canMultiply(const Dims& other) const {if (cols == other.rows) return true;  return false;}
    };

    template <std::integral T>
    class Matrix{
        T** matrix;
        Dims dims{};

        Matrix subMatrix(const Matrix &A, size_t i, size_t j);
        static void subMatrix(Matrix& cM, const Matrix& A, size_t i, size_t j);


    public:
        explicit Matrix(size_t n);
        Matrix(size_t m, size_t n);
        Matrix(T** data, size_t m, size_t n);
        Matrix(std::initializer_list<std::initializer_list<T>> data, size_t m, size_t n);
        explicit Matrix(const std::string& fileName = "matrix.csv"); // First line contains dimensions, remaining file contents are data points //
        Matrix(const std::string& I, size_t n);
        Matrix(const Matrix& copy);
        ~Matrix();

        [[nodiscard]] const Dims& getDims() const {return dims;}
        Matrix& operator=(Matrix &&other) noexcept {
            //If other is destroyed this object well lose its values?
            if (this == &other)
                return *this;
            matrix = other.matrix;
            dims = other.dims;
            for (size_t i = 0; i < dims.rows; i++) {
                for (size_t j = 0; j < dims.cols; j++) {
                    this->matrix[i][j] = std::move(other.matrix[i][j]);
                }
            }
            return *this;
        }

        Matrix& operator=(const Matrix& copy) {
            srand(clock());
            if(this != &copy){
                this->dims = copy.dims;
                this->matrix = new T*[this->dims.rows];
                for(size_t i = 0; i < this->dims.rows; i++) {
                    this->matrix[i] = new T[this->dims.cols];
                    for(size_t j = 0; j < this->dims.cols; j++){
                        this->matrix[i][j] = copy.matrix[i][j];
                    }
                }
            }
            return *this;
        }

        constexpr bool operator==(const Matrix<T>& rhs) {
            if(dims == rhs.dims) {
                for(size_t i = 0; i < this->dims.rows; i++) {
                    for(size_t j = 0; j < this->dims.cols; j++) {
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
            if(dims == rhs.dims) {
                Matrix result(dims.rows, dims.cols);
                for(int i = 0; i < this->dims.rows; i++) {
                    for (int k = 0; k < this->dims.cols; k++) {
                        result.matrix[i][k] = this->matrix[i][k] + rhs.matrix[i][k];
                    }
                }
                return result;
            }
            return std::nullopt;
        }

        constexpr std::optional<Matrix> operator-(const Matrix& rhs){
            if(dims == rhs.dims) {
                Matrix result(dims.rows, dims.cols);
                for(int i = 0; i < this->dims.rows; i++) {
                    for (int k = 0; k < this->dims.cols; k++) {
                        result.matrix[i][k] = this->matrix[i][k] - rhs.matrix[i][k];
                    }
                }
                return result;
            }
            return std::nullopt;
        }

        constexpr std::optional<Matrix> operator*(const Matrix& rhs){
            if (dims.canMultiply(rhs.getDims())) {
                Matrix result{this->dims.rows, rhs.dims.cols};
                for(size_t i = 0; i < this->dims.rows; i++) {
                    for (size_t k = 0; k < rhs.dims.cols; k++) {
                        for(size_t j = 0; j < rhs.dims.rows; j++) {
                            result.matrix[i][k] += this->matrix[i][j] * rhs.matrix[j][k];
                            // std::fputs((std::to_string(this->matrix[i][j] * rhs.matrix[j][k]) + "+").data(), stdout);
                        }
                        // std::fputs(("=" + std::to_string(result.matrix[i][k]) + " ").data(), stdout);
                    }
                    // std::puts("");
                }
                return result;
            }
            return std::nullopt;
        }

        constexpr std::optional<Matrix> operator/(const Matrix& divisor) {
            /*
             *"X = B / A"
             *AX = B
             *A^-1 * AX = A^-1 * B
             *IX = A^-1 * B
             *X = A^-1 * B
             */
            Matrix inverseA = *inverse(divisor);
            if (inverseA.dims.cols == this->dims.cols) return inverseA * (*this);
            return std::nullopt;
        }

        constexpr T* operator[](const size_t m) {
            return this->matrix[m];
        }
        //Only Compatible with std=c++23 or greater
        // constexpr T operator[](const size_t m, const size_t n) {
        //     return this->matrix[m][n];
        // }

        void setDataRand() const{
            for (size_t i = 0; i < dims.rows; i++) {
                for (size_t j = 0; j < dims.cols; j++) {
                    this->matrix[i][j] = rand()%9 + 1; //Random integer between 1-10
                }
            }
        }

        // void setDataFile(const char* file) const{
        //     //Rudimentary file io
        //     std::ifstream myFile(file);
        //     while(!myFile.eof()){
        //         size_t i = 0;
        //         size_t j = 0;
        //         myFile >> matrix[i][j++];
        //         if (j == cols) {
        //             i++;
        //             j = 0;
        //         }
        //     }
        // }

        constexpr double determinant(const Matrix& A){
            //Need to develop error handling here
            if (A.dims.rows == 2){
                return (A.matrix[0][0] * A.matrix[1][1]) - (A.matrix[1][0] * A.matrix[0][1]);
            }
            size_t i = 0;
            double det = 0;
            //Consider refactoring as sufficiently large matrices could cause a stack overflow with recursion
            for(size_t j = 0; j < A.dims.cols; j++){
                //parallize by splitting matrix into 4 submatrices and solve determinant for each quadrant in sep thread
                det += pow(-1, i+j+2) * A.matrix[i][j] * determinant(subMatrix({},A,i, j));
            }
            return det;
        }

        constexpr auto inverse(const Matrix& A){
            //Not currently Functional
            Matrix inverse(A);
            Matrix identity('I',A.dims.rows);
            if (determinant(A) != 0) {
                for (size_t i = 0; i < A.dims.rows; i++){
                    for(size_t j = 0; j < A.dims.rows; j++){
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
        }

        constexpr auto cofactor(const Matrix& A){
            Matrix cofactorM(A.dims.rows);
            for(size_t i = 0; i < A.dims.rows; i++){
                for(int j = 0; j < A.dims.rows; j++){
                    cofactorM.matrix[i][j] = pow(-1,i+j) * determinant(subMatrix({},A,i, j));
                }
            }

            return &cofactorM;
        }

        constexpr auto transpose(const Matrix& cofactorM) {
            Matrix transposeM(cofactorM.dims.rows);
            for(int i = 0; i < cofactorM.dims.rows; i++){
                for(int j = 0; j < cofactorM.dims.rows; j++){
                    transposeM.matrix[i][j] = cofactorM.matrix[j][i];
                }
            }
            return &transposeM;
        }

        constexpr auto adjoint(const Matrix& A) {
            return transpose(cofactor(A));
        }

        // void matrixToI(const Matrix& A);

        [[nodiscard]] constexpr bool isSquare() const {
            return dims.rows == dims.cols;
        }
        void print();
        void serialize(const std::string& fileName = "matrix.csv") {
            std::fstream newFile;
            newFile.open(fileName, std::ios::out | std::ios::app);

            if (!newFile.is_open()) {
                std::puts(std::string_view("Error opening file " + fileName).data());
            } else {
                newFile << dims.rows << ',' << dims.cols << std::endl;
                for(size_t i = 0; i < dims.rows; i++) {
                    for(size_t j = 0; j < dims.cols; j++) {
                        newFile << this->matrix[i][j];
                        if(j != dims.cols-1) newFile << ',';
                    }
                    if (i != dims.rows-1) newFile << std::endl;
                }
                newFile.close();
            }
        }
    };

    //Big Five
    template<std::integral T>
    Matrix<T>::Matrix(const size_t n):dims(n, n) {
        srand(clock());
        //Construction of a square matrix of size n
        this->matrix = new T*[n];
        for (size_t i = 0; i < n; i++){
            this->matrix[i] = new T[n];
        }
    };

    template<std::integral T>
    Matrix<T>::Matrix(const size_t m, const size_t n): dims(m,n) {
        //Construction of a matrix of size m(rows) x n(columns)
        srand(clock());
        this->matrix = new T *[m];
        for (size_t i = 0; i < m; i++) {
            this->matrix[i] = new T[n];
            for(size_t j = 0; j < n; j++) {
                this->matrix[i][j] = 0;
            }
        }
    };

    template<std::integral T>
    Matrix<T>::Matrix(T** data, const size_t m, const size_t n): dims(m,n) {
        srand(clock());
        this->matrix = new T*[m];
        for(size_t i = 0; i < m; i++) {
            this->matrix[i] = new T[n];
            for(size_t j = 0; j < n; j++) {
                this->matrix[i][j] = data[i][j];
            }
        }
    }

    template<std::integral T>
    Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> data,const size_t m, const size_t n): dims(m,n) {
        srand(clock());
        this->matrix = new T*[m];
        auto it = data.begin();
        for(size_t i = 0; i < m; i++) {
            this->matrix[i] = new T[n];
            for(size_t j = 0; j < n; j++) {
                this->matrix[i][j] = it->begin()[j];
            }
            ++it;
        }
    }

    template<std::integral T>
    Matrix<T>::Matrix(const std::string& I, const size_t n): dims(n,n){
        //Creates nxn Identity matrix
        if (I == "I"){
            srand(clock());
            this->matrix = new T*[n];
            for(size_t i = 0; i < n; i++) {
                this->matrix[i] = new T[n];
                for(size_t j = 0; j < n; j++) {
                    if (i == j){
                        matrix[i][j] = 1;
                    } else {
                        matrix[i][j] = 0;
                    }
                }
            }
        }
    }

    template<std::integral T>
    Matrix<T>::Matrix(const Matrix& copy): dims(copy.dims) {
        srand(clock());
        this->matrix = new T*[this->dims.rows];
        for(size_t i = 0; i < this->dims.rows; i++) {
            this->matrix[i] = new T[this->dims.cols];
            for(size_t j = 0; j < this->dims.cols; j++){
                this->matrix[i][j] = copy.matrix[i][j];
            }
        }
    }

    template<std::integral T>
    Matrix<T>::~Matrix()
    {
        for (int i = 0; i < dims.rows; i++){
            delete [] matrix[i];
        }

        delete [] matrix;
    }

    template<std::integral T>
    Matrix<T>::Matrix(const std::string& fileName ): dims(0,0)
    {
        std::fstream newFile;
        newFile.open(fileName, std::ios::in);
        if (!newFile.is_open()) {
            std::puts(std::string_view("Error opening file " + fileName).data());
        } else {
            newFile >> dims.rows >> dims.cols;
            matrix = new T*[dims.rows];
            for(size_t i = 0; i < dims.rows; i++) {
                matrix[i] = new T[dims.cols];
                for(size_t j = 0; j < dims.cols; j++) {
                    newFile >> this->matrix[i][j];
                }
            }
            newFile.close();
        }
    }

    template<std::integral T>
    Matrix<T> Matrix<T>::subMatrix(const Matrix &A, const size_t i, const size_t j)
    {
        //Only call if n x n and n >= 3;
        Matrix sM(A.dims.rows - 1);
        for(int m = 0; m < A.dims.rows; m++){
            for(int n = 0; n < A.dims.rows; n++){
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


    template<std::integral T>
    void Matrix<T>::subMatrix(Matrix &cM, const Matrix &A, const size_t i, const size_t j)
    {
        // Same as above but useful if subMatrix matrix should be usable and not const
        for(int m = 0; m < A.dims.rows; m++){
            for(int n = 0; n < A.dims.rows; n++){
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
    }
    /*
    template<std::integral T>
    void Matrix<T>::matrixToI(const Matrix &A)
    {
        //Incomplete
        Matrix tmp(A.rows);
        Matrix identity('I',A.rows);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < rows; j++){
                tmp.matrix[i][j] = identity.matrix[i][j] - A.matrix[i][j];
            }
        }
        tmp.print();
    }
    */
    template<std::integral T>
    void Matrix<T>::print()
    {
        std::puts("[");
        for (int i = 0; i < dims.rows; i++) {
            std::fputs("\t", stdout);
            for (int j = 0; j < dims.cols; j++) {
                std::fputs(std::to_string(matrix[i][j]).data(), stdout);
                std::fputs(" ",stdout);
            }
            std::puts("");
        }
        std::puts("]");
    }
}
#endif
