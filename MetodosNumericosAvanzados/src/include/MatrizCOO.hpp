#ifndef __MatrizCOO__
#define __MatrizCOO__

#include "../include/MatrizBase.hpp"
#include "../include/MatrizDensa.hpp"

class MatrizCOO : public MatrizBase {
public:
    /****Atributos*******/
    std::vector<double> data;
    std::vector<int> row;
    std::vector<int> col;
    int nnz;

    /****Métodos*******/
    MatrizCOO(void): MatrizBase() {
        data.resize(0);
        row.resize(0);
        col.resize(0);
        nnz=0;
    }
    MatrizCOO(int);
    ~MatrizCOO() {
    }
    void convert(const MatrizDensa &);
    void initialize();
    void print();
    void insert(int, int, double);
    std::vector<double> operator * (const std::vector<double>&);
    void JacobiIter(std::vector<double>&, const std::vector<double>&);
};

#endif
