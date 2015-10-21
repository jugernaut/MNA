#ifndef __ILU__
#define __ILU__

#include "../include/Precondicionador.hpp"

template <class T>
class ILU : public Precondicionador {
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
        return "ILU Precondicionador";
    }
};

#endif
