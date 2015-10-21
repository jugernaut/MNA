#include "../include/Matrix_CSC.hpp"
#include <iostream>

Matrix_CSC::Matrix_CSC(int n) {
    nnz=0;
    Col=Row=n;
}

void Matrix_CSC::initialize() {
}

void Matrix_CSC::print() {
    std::cout<<std::endl<<"Data: [";
    for (std::vector<double>::iterator it = data.begin() ; it != data.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"row: [";
    for (std::vector<int>::iterator it = row.begin() ; it != row.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"icol: [";
    for (std::vector<int>::iterator it = icol.begin() ; it != icol.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;
}

std::vector<double> Matrix_CSC::operator * (const std::vector<double>& V) {
    if (Row!=V.size()) {
        std::cout<<"Error en las dimensiones"<<std::endl;
        exit (1);
    }
    std::vector<double> y;
    y.resize(Row,0.0);
    for (int j = 0; j < Row; ++j) {
        for (int l = icol[j]; l <= icol[j+1]-1; ++l) {
            y[row[l]]+=data[l]*V[j];
        }
    }
    return y;
}

void Matrix_CSC::convert(const Matrix_COO & coo) {
    nnz=coo.nnz;
    Col=Row=coo.Col;
    data.resize(nnz);
    row.resize(nnz);
    icol.resize(Col+1,0);
    for (int l = 0; l < nnz; ++l) {
        icol[coo.col[l]]++;
    }
    int sum=0;
    int tmp;
    for (int i = 0; i < Col+1; ++i) {
        tmp=icol[i];
        icol[i]=sum;
        sum+=tmp;
    }
    icol[Col]=nnz;
    std::vector<int> icolAux;
    icolAux=icol;
    int i,j,d;
    double val;
    for (int l = 0; l < nnz; ++l) {
        i=coo.row[l];
        j=coo.col[l];
        val=coo.data[l];
        d=icolAux[j];
        data[d]=val;
        row[d]=i;
        icolAux[j]++;
    }
}

void Matrix_CSC::JacobiIter(std::vector<double>& x, const std::vector<double>& b) {
    std::cout<<"JacobiIter csc"<<std::endl;
}
