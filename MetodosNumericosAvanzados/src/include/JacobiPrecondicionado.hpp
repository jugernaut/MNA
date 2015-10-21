#ifndef __JacobiPrecondicionado__
#define __JacobiPrecondicionado__

#include "../include/Precondicionador.hpp"

template <class T>
class JacobiPrecondicionado : public Precondicionador {
public:
    T mat;
    JacobiPrecondicionado(void) {
    }
    ~JacobiPrecondicionado() {
    }
    void calculate(T &A) {
        mat = A.diag();
    }
    void solve(std::vector<double> &z, const std::vector<double> &r) {
        mat.JacobiSolve(z,r);
    }
    std::string name() {
        return "Jacobi Precondicionador";
    }
};

#endif
