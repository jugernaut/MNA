#ifndef __SolverDirect__
#define __SolverDirect__

#include "../include/Solver.hpp"


class SolverDirect : public Solver {
public:
    bool ModifyMatrix;
    SolverDirect(void) {
    }
    ~SolverDirect() {
    }
};

#endif
