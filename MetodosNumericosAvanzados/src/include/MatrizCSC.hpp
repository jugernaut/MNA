#ifndef __MatrizCSC__
#define __MatrizCSC__

#include "../include/MatrizBase.hpp"
#include "../include/MatrizCOO.hpp"

class MatrizCSC : public MatrizBase {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> row;
    std::vector<int> icol;
    int nnz;

    /****Métodos*******/
    MatrizCSC(void): MatrizBase() {
        data.resize(0);
        row.resize(0);
        icol.resize(0);
    }
    MatrizCSC(int);
    ~MatrizCSC() {
    }
    void convert(const MatrizCOO &);
    void initialize();
    void print();
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
};

#endif
