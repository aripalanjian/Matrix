# Neo - Header Only Linear Algebra Library
- Description: A header only linear algebra library
- Contributor(s): Ari Palanjian
- [GitHub](https://github.com/aripalanjian/Matrix.git)
## About
Began as an exercise to develop mathematical equations while studying linear algebra at university to reinforce
the covered material. This project is in transition from a compiled version to a header only library and thus at this time 
is not guaranteed to work. As I am able, I will update this library to more modern c++20 standards and to add multithreading.
## Getting Started
- Build from source
  - Download repository
  - In your CMake project structure in the directory where the library is stored add 
  ```add_subdirectory(Neo)``` to your CMakeLists.txt file.
  - Reload or Sync your integrated development environment.
## Interface
### Constructors
 - Matrix(const size_t n):                 Construction of a square matrix of size n
 - Matrix(const size_t m, const size_t n): Construction of a matrix of size m(rows) x n(columns)
 - Matrix(const char I, const size_t n):   Creates Identity matrix of size n
 - Matrix(const Matrix& copy):             Copies a given matrix into a new object

### Overloaded Operators
- operator= : Working on implementation 
- operator==: Returns true if every element in lhs matrix is the same as the rhs
- operator!=: Return the opposite of operator==
- operator\+: Applies matrix addition to two matrices
- operator\-: Applies matrix subtraction to two matrices
- operator\*: Applies matrix multiplication to two matrices
- operator\/: Applies matrix multiplication to two matrices to perform division in the form $X = A^{-1} * B$

### Operations
- determinant
- inverse
- cofactor
- transpose
- adjoint
- isSquare
- print

## Future
- Update codebase fully to -std=c++20
- Add multithreaded operation
