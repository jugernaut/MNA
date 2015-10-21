#ifndef __JacobiPreconditioner__
#define __JacobiPreconditioner__

#include "../include/Preconditioner.hpp"

template <class T>
class JacobiPreconditioner : public Preconditioner {
public:
    T mat;
    JacobiPreconditioner(void) {
    }
    ~JacobiPreconditioner() {
    }
    void calculate(T &A) {
        mat = A.diag();
    }
    void solve(std::vector<double> &z, const std::vector<double> &r) {
        mat.JacobiSolve(z,r);
    }
    std::string name() {
        return "Jacobi Preconditioner";
    }
};

#endif
