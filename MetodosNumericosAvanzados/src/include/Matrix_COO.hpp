#ifndef __Matrix_COO__
#define __Matrix_COO__

#include "../include/Matrix_base.hpp"
#include "../include/Matrix_dense.hpp"

class Matrix_COO : public Matrix_base {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> row;
    std::vector<int> col;
    int nnz;

    /****Métodos*******/
    Matrix_COO(void): Matrix_base() {
        data.resize(0);
        row.resize(0);
        col.resize(0);
        nnz=0;
    }
    Matrix_COO(int);
    ~Matrix_COO() {
    }
    void convert(const Matrix_dense &);
    void initialize();
    void print();
    void insert(int, int, double);
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
};

#endif
