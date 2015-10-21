#ifndef __Matrix_CRS__
#define __Matrix_CRS__

#include "../include/Matrix_base.hpp"
#include "../include/Matrix_COO.hpp"
#include "../include/Matrix_CSC.hpp"

class Matrix_CRS : public Matrix_base {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> col;
    std::vector<int> irow;
    std::vector<int> idiag;
    int nnz;

    /****Métodos*******/
    Matrix_CRS(void): Matrix_base() {
        data.resize(0);
        col.resize(0);
        irow.resize(0);
    }
    Matrix_CRS(int);
    ~Matrix_CRS() {
    }
    void initialize();
    void print();
    void convert(const Matrix_COO&);
    void convertCCStoCRS(const Matrix_CSC&);
    std::vector<double> operator * (const std::vector<double>&);
    Matrix_CRS const &operator=(Matrix_CRS const&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
    void idiagCalculate();
    Matrix_CRS diag();
    Matrix_CRS ILU();
    Matrix_CRS MILU();
    Matrix_CRS ICHOL();
    void JacobiSolve(std::vector<double>& , const std::vector<double>&);
    void LUSolve(std::vector<double>& , const std::vector<double>&);
};

#endif
