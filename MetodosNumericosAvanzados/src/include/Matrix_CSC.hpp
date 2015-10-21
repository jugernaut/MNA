#ifndef __Matrix_CSC__
#define __Matrix_CSC__

#include "../include/Matrix_base.hpp"
#include "../include/Matrix_COO.hpp"

class Matrix_CSC : public Matrix_base {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> row;
    std::vector<int> icol;
    int nnz;

    /****Métodos*******/
    Matrix_CSC(void): Matrix_base() {
        data.resize(0);
        row.resize(0);
        icol.resize(0);
    }
    Matrix_CSC(int);
    ~Matrix_CSC() {
    }
    void convert(const Matrix_COO &);
    void initialize();
    void print();
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
};

#endif
