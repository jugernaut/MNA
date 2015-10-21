
#ifndef __funcionesExternas__
#define __funcionesExternas__


#include <vector>
#include <iostream>
//#include "Matriz_base.hpp"

void printVector(const std::vector<double> &);
void printVector(const std::vector<int> &);
double dotProduct(const std::vector<double> &, const std::vector<double> &);
void sub(const std::vector<double> &,const std::vector<double> &,std::vector<double> &);
void sum(const std::vector<double> &,const std::vector<double> &,std::vector<double> &);
void scalarTimesVector(double, std::vector<double>&,std::vector<double>&);
double evalError(std::vector<double> &, const std::vector<double>&);


#endif
