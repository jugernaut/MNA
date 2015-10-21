#ifndef __Jacobi__
#define __Jacobi__

#include "../include/SolverIterativo.hpp"

class Jacobi : public SolverIterativo {
public:
    Jacobi(void):SolverIterativo() {
    }
    ~Jacobi() {
    }

    template <class T>
    void solve(T &A, std::vector<double> &x, std::vector<double> &b) {
        //aproximación inicial en x
        double norm;
        std::vector<double> aux;
        timer.tic();
        for(int k =1; k<=mmaxIts; ++k)
            //for(int k =0;k<1;++k)
        {
            mits++;
            A.JacobiIter(x,b);
            /*****************comprobar tolerancia**************/
            aux=A*x;
            norm=evalError(aux,b);
            if(norm<=mtol*mtol)
                break;
        }
        //mits++;
        merror = sqrt(norm);
        timer.toc();
        etime = timer.etime();
    }

    void report(std::string);
};

#endif
