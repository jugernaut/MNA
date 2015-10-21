#ifndef __Matrix_base__
#define __Matrix_base__

//#include <iostream>
#include <vector>
#include "../include/Funciones.hpp"

class Matrix_base {
public:
    int Col;
    int Row;

    Matrix_base(void) {
        Col=0;
        Row=0;
    }
};

#endif
