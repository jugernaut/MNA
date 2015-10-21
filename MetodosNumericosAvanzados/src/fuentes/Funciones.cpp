
#include "../include/Funciones.hpp"

void printVector(const std::vector<double> &r ) {
    for (int i = 0; i < r.size(); ++i) {
        std::cout<<r[i]<<", " ;
    }
    std::cout<<std::endl;
}

void printVector(const std::vector<int> &r ) {
    for (int i = 0; i < r.size(); ++i) {
        std::cout<<r[i]<<", " ;
    }
    std::cout<<std::endl;
}

double dotProduct(const std::vector<double> &a, const std::vector<double> &b ) {
    double aux=0;
    for (int i = 0; i < a.size(); ++i) {
        aux+=a[i]*b[i];
    }
    return aux;
}

void sub(const std::vector<double> &a,const std::vector<double> &b,std::vector<double> &c) { //c=a-b
    int n=std::min(a.size(),b.size());
    c.resize(n);
    for (int i = 0; i < n; ++i) {
        c[i]=a[i]-b[i];
    }
}

void sum(const std::vector<double> &a,const std::vector<double> &b,std::vector<double> &c) { //c=a+b
    int n=std::min(a.size(),b.size());
    c.resize(n);
    for (int i = 0; i < n; ++i) {
        c[i]=a[i]+b[i];
    }
}

void scalarTimesVector(double a, std::vector<double>& v,std::vector<double>& v2) {
    v2.resize(v.size());
    for (int i = 0; i < v.size(); ++i) {
        v2[i]=a*v[i];
    }
}

double evalError(std::vector<double> &x, const std::vector<double>& b) {
    double norm;
    sub(x,b,x);
    norm=dotProduct(x,x);
    return norm;
}
