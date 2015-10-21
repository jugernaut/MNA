#ifndef __MatrizCDS__
#define __MatrizCDS__

#include "../include/MatrizBase.hpp"

class MatrizCDS : public MatrizBase {
public:
    /****Atributos*******/
    double **A;
    int p, q, w,dimension,n;
    std::vector<int> identificacionbandas;

    /****Métodos*******/
    MatrizCDS(void): MatrizBase() {
        dimension=0;
        Col=0;
        Row=0;
        identificacionbandas.resize(0);
    }
    MatrizCDS(int ,int , int , std::vector<int>&);
    MatrizCDS(const MatrizCDS&);
    ~MatrizCDS();
    void inicializa();
    void print();
    void escribir(int,int,double);
    std::vector<double> operator * (const std::vector<double>&);
    MatrizCDS& operator =(MatrizCDS& );
    double obtener(int i,int j);
    void JacobiIter(std::vector<double>& x, const std::vector<double>& b);
};

#endif
