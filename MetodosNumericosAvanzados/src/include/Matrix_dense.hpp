#ifndef __Matrix_dense__
#define __Matrix_dense__

#include "../include/Matrix_base.hpp"
#include <iostream>

class Matrix_dense : public Matrix_base {
public:
    /****Atributos*******/
    std::vector<std::vector<double> > A;

    /****Métodos*******/
    Matrix_dense(void): Matrix_base() {
        A.resize( Row , std::vector<double>( Col) );
    }
    ~Matrix_dense() {
    }
    void inicializa(int, int, double);
    void print();
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
    Matrix_dense const &operator=(Matrix_dense const&);
    Matrix_dense diag();
    void forwardBackward(std::vector<double>&, const std::vector<double>&);
};

#endif

