#ifndef __MatrizCRS__
#define __MatrizCRS__

#include "../include/MatrizBase.hpp"
#include "../include/MatrizCOO.hpp"
#include "../include/MatrizCSC.hpp"

class MatrizCRS : public MatrizBase {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> col;
    std::vector<int> irow;
    std::vector<int> idiag;
    int nnz;

    /****Métodos*******/
    MatrizCRS(void): MatrizBase() {
        data.resize(0);
        col.resize(0);
        irow.resize(0);
    }
    MatrizCRS(int);
    ~MatrizCRS() {
    }
    void initialize();
    void print();
    void convert(const MatrizCOO&);
    void convertCCStoCRS(const MatrizCSC&);
    std::vector<double> operator * (const std::vector<double>&);
    MatrizCRS const &operator=(MatrizCRS const&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
    void idiagCalculate();
    MatrizCRS diag();
    MatrizCRS ILU();
    MatrizCRS MILU();
    MatrizCRS ICHOL();
    void JacobiSolve(std::vector<double>& , const std::vector<double>&);
    void LUSolve(std::vector<double>& , const std::vector<double>&);
};

#endif
