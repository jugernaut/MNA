#include "../include/Matrix_COO.hpp"
#include <iostream>


Matrix_COO::Matrix_COO(int n) {
    nnz=0;
    Col=Row=n;
}

void Matrix_COO::initialize() {
    //parece que no es necesaria
}

void Matrix_COO::print() {
    std::cout<<std::endl<<"Data: [";
    for (std::vector<double>::iterator it = data.begin() ; it != data.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"Row: [";
    for (std::vector<int>::iterator it = row.begin() ; it != row.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"Col: [";
    for (std::vector<int>::iterator it = col.begin() ; it != col.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;
}

void Matrix_COO::insert(int i, int j, double val) {
    data.push_back(val);
    col.push_back(j);
    row.push_back(i);
    nnz++;
}

std::vector<double> Matrix_COO::operator * (const std::vector<double>& V) {
    if (Row!=V.size()) {
        std::cout<<"Error en las dimensiones"<<std::endl;
        exit (1);
    }
    std::vector<double> y;
    y.resize(Row,0.0);
    for (int l = 0; l < nnz; ++l) {
        y[row[l]]+=data[l]*V[col[l]];
    }
    return y;
}

void Matrix_COO::convert(const Matrix_dense& MD) { //opcional o de mejor entendimiento que es un constructor "convertidor"
    for (int i = 0; i < Row; ++i) {
        for (int j = 0; j < Col; ++j) {
            if(MD.A[i][j]!=0) {
                data.push_back(MD.A[i][j]);
                col.push_back(j);
                row.push_back(i);
                nnz++;
            }
        }
    }
}

void Matrix_COO::JacobiIter(std::vector<double>& x, const std::vector<double>& b) {
    std::cout<<"JacobiIter coo"<<std::endl;
}
