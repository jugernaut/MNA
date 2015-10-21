#ifndef __MILU__
#define __MILU__

#include "../include/Preconditioner.hpp"

template <class T>
class MILU : public Preconditioner {
public:
    T mat;
    MILU(void) {
    }
    ~MILU() {
    }
    void calculate(T &A) {
        mat = A.MILU();
    }
    void solve(std::vector<double> &z, const std::vector<double> &r) {
        mat.LUSolve(z,r);
    }
    std::string name() {
        return "MILU Preconditioner";
    }
};

#endif
