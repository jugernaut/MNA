#ifndef __SolverDirecto__
#define __SolverDirecto__

#include "../include/Solver.hpp"


class SolverDirecto : public Solver {
public:
    bool ModifyMatrix;
    SolverDirecto(void) {
    }
    ~SolverDirecto() {
    }
};

#endif
