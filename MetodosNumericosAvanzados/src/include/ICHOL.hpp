#ifndef __ICHOL__
#define __ICHOL__

#include "../include/Preconditioner.hpp"

template <class T>
class ICHOL : public Preconditioner {
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
        return "ICHOL Preconditioner";
    }
};

#endif
