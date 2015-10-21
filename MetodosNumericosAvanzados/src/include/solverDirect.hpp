#ifndef __solverDirect__
#define __solverDirect__

#include "../include/solver.hpp"


class solverDirect : public solver {
public:
    bool ModifyMatrix;
    solverDirect(void) {
    }
    ~solverDirect() {
    }
};

#endif
