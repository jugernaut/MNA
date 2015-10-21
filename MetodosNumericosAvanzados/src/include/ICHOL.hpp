#ifndef __ICHOL__
#define __ICHOL__

#include "../include/Precondicionador.hpp"

template <class T>
class ICHOL : public Precondicionador {
public:
    T mat;
    ICHOL(void) {
    }
    ~ICHOL() {
    }
    void calculate(T &A) {
        mat = A.ICHOL();
    }
    void solve(std::vector<double> &z, const std::vector<double> &r) {
        mat.CholSolve(z,r);
    }
    std::string name() {
        return "ICHOL Precondicionador";
    }
};

#endif
