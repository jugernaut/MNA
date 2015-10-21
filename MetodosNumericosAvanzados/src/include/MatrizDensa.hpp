#ifndef __MatrizDensa__
#define __MatrizDensa__

#include "../include/MatrizBase.hpp"
#include <iostream>

class MatrizDensa : public MatrizBase {
public:
    /****Atributos*******/
    std::vector<std::vector<double> > A;

    /****Métodos*******/
    MatrizDensa(void): MatrizBase() {
        A.resize( Row , std::vector<double>( Col) );
    }
    ~MatrizDensa() {
    }
    void inicializa(int, int, double);
    void print();
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
    MatrizDensa const &operator=(MatrizDensa const&);
    MatrizDensa diag();
    void forwardBackward(std::vector<double>&, const std::vector<double>&);
};

#endif

