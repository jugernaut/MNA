#ifndef __solverIter__
#define __solverIter__

#include "../include/solver.hpp"
#include "../include/Timer.hpp"

class solverIter : public solver {
public:
    int mits;
    int mmaxIts;
    double mtol;
    double etime;
    double merror;
    bool precond=false;
    std::string namep;
    solverIter(void) {
        mmaxIts = 2000;
        mtol = 1e-6;
        mits = 0;
        merror = 0;
        etime = 0;
    }
    ~solverIter() {
    }
    void maxIts(int it) {
        mmaxIts = it;
    }
    void tol(double t) {
        mtol=t;
    }
};

#endif

