#ifndef __ILU__
#define __ILU__

#include "../include/Preconditioner.hpp"

template <class T>
class ILU : public Preconditioner {
public:
    T mat;
    ILU(void) {
    }
    ~ILU() {
    }
    void calculate(T &A) {
        mat = A.ILU();
    }
    void solve(std::vector<double> &z, const std::vector<double> &r) {
        mat.LUSolve(z,r);
    }
    std::string name() {
        return "ILU Preconditioner";
    }
};

#endif
