#include "../include/MatrizCDS.hpp"
#include <algorithm>
#include <iostream>

MatrizCDS::MatrizCDS(int NuevoN, int numBandasSup, int numBandasInf, std::vector<int> &indices): MatrizBase() {

    p = numBandasSup; //numero de diagonales sobre la diagonal principal distintas de cero
    q = numBandasInf; //numero de diagonales bajo la diagonal principal distintas de cero

    identificacionbandas=indices;

    w = p + q + 1; //w=indices.size();

    Col=Row=NuevoN; //pensando que es una Matrix cuadrada de NuevoN x NuevoN

    n=NuevoN;

    dimension = n*w;
    A = new double*[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[w];
}

MatrizCDS::MatrizCDS(const MatrizCDS &m) {
    p=m.p;
    q=m.q;
    w=m.w;

    Col=m.Col;
    Row=m.Row;
    n=m.n;

    dimension=m.dimension;
    identificacionbandas = m.identificacionbandas;

    A = new double*[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[w];

    for(int i = 0; i < n; i++)
        for(int j = 0; j < w; j++)
            A[i][j] =  m.A[i][j];
}

MatrizCDS& MatrizCDS::operator =(MatrizCDS& m) {
    p=m.p;
    q=m.q;
    w=m.w;

    Col=m.Col;
    Row=m.Row;
    n=m.n;

    dimension=m.dimension;
    identificacionbandas = m.identificacionbandas;

    A = new double*[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[w];
    for(int i = 0; i < n; i++)
        for(int j = 0; j < w; j++)
            A[i][j] =  m.A[i][j];
}

double MatrizCDS::obtener(int i,int j) {
    int ind;
    ind=i-j;
    std::vector<int>::iterator it;
    int k;
    it = find (identificacionbandas.begin(), identificacionbandas.end(),ind);
    if (it!=identificacionbandas.end()) {
        k = std::distance( identificacionbandas.begin(), it );
        return A[i][k];
    } else
        return 0.0;
}

void MatrizCDS::escribir(int i, int j, double valor) {
    int ind;
    ind=j-i;
    std::vector<int>::iterator it;
    int k;
    it = find (identificacionbandas.begin(), identificacionbandas.end(),ind);
    if (it!=identificacionbandas.end()) {
        k = std::distance( identificacionbandas.begin(), it );
        A[i][k]=valor;
    } else
        std::cout<<"No se puede escribir ahí"<<std::endl;
}

MatrizCDS::~MatrizCDS () {
    for( int i = 0; i < n; i++ ) {
        delete [] A[i];
    }
    delete [] A;
    //A = NULL;
}

void MatrizCDS::print() {
    std::cout<<std::endl<<"vector de indices: [ ";
    for(int i = 0; i < w; i++)
        std::cout << identificacionbandas[i] << ",";
    std::cout<<" ]"<<std::endl;

    for(int i = 0; i < n ; i++) {
        for(int j = 0; j < w; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout<<std::endl;
    }
}

std::vector<double> MatrizCDS::operator * (const std::vector<double>& V) {
    std::vector<double> y;
    if(n!=V.size()) {
        std::cout<<"Error de dimensiones"<<std::endl;
        exit(0);
    } else {
        y.resize(n,0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < w; ++j) {
                if(! (i+identificacionbandas[j]<0 && i+identificacionbandas[j]>n-1) ) {
                    y[i]+=A[i][j]*V[i+identificacionbandas[j]];
                }
            }
        }
    }

    return y;
}

void MatrizCDS::JacobiIter(std::vector<double>& x, const std::vector<double>& b) {
    std::cout<<"JacobiIter cds"<<std::endl;
}

