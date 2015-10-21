#include "../include/MatrizDensa.hpp"
#include <iostream>
#include <stdio.h>

void MatrizDensa::inicializa(int Rowglones, int columnas, double valor=0) {
    Row=Rowglones;
    Col=columnas;
    A.resize( Row , std::vector<double>( Col  , valor) );
}

void MatrizDensa::print() {
    for(int i=0; i<Row; i++) {
        std::cout<<std::endl;
        for(int j=0; j<Col; j++)
            printf("%.3f ",A[i][j]);
    }
    std::cout<<std::endl;
}

std::vector<double> MatrizDensa::operator * (const std::vector<double>& v) {
    if (Row!=v.size()) {
        std::cout<<"Error en las dimensiones"<<std::endl;
        exit (1);
    }

    std::vector<double> y;
    y.resize(Row);
    for (int i = 0; i < Row; ++i) {
        y[i]=0;
        for (int j = 0; j < Col; ++j) {
            y[i]+=A[i][j]*v[j];
        }
    }
    return y;
}

void MatrizDensa::JacobiIter(std::vector<double>& x, const std::vector<double>& b) {
    std::cout<<"JacobiIter densa"<<std::endl;
}

MatrizDensa const& MatrizDensa::operator=(MatrizDensa const &rhs) {
    if (this != &rhs) {
        A=rhs.A;
        Col=rhs.Col;
        Row=rhs.Row;
    }
    return *this;
}

MatrizDensa MatrizDensa::diag() {
    MatrizDensa D;
    D.inicializa(Row,Col,0.0);
    for (int i = 0; i < Row; ++i) {
        D.A[i][i]=A[i][i];
    }
    return D;
}

void MatrizDensa::forwardBackward(std::vector<double>&x, const std::vector<double>&b) {
    std::vector<double> y(Row);
    for (int i = 0; i < Row; ++i) {
        y[i]=b[i];
        for (int j = 0; j<i; ++j) {
            y[i]-=A[i][j]*y[j];
        }
        //y[i]/=A[i][i];
    }
    std::cout<<"forward y"<<std::endl;
    printVector(y);

    for (int i = Row-1; i >=0; --i) {
        x[i]=y[i];
        for (int j = i+1; j < Row; ++j) {
            x[i]-=A[i][j]*x[j];
        }
        x[i]/=A[i][i];
    }
    std::cout<<"backward x"<<std::endl;
    printVector(x);

}
