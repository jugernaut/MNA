#ifndef __SolverIterativo__
#define __SolverIterativo__

#include "../include/Solver.hpp"
#include "../include/Timer.hpp"

class SolverIterativo : public Solver {
public:
    int mits;
    int mmaxIts;
    double mtol;
    double etime;
    double merror;
    bool precond=false;
    std::string namep;
    SolverIterativo(void) {
        mmaxIts = 2000;
        mtol = 1e-6;
        mits = 0;
        merror = 0;
        etime = 0;
    }
    ~SolverIterativo() {
    }
    void maxIts(int it) {
        mmaxIts = it;
    }
    void tol(double t) {
        mtol=t;
    }
};

#endif

