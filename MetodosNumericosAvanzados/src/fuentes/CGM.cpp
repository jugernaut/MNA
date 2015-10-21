#include "../include/CGM.hpp"
#include <iostream>

void CGM::report(std::string ruta) {
    std::cout<<"CGM report: "<<std::endl;
    if(precond)
        std::cout<<"Precondition with: "<< namep <<std::endl;
    std::cout<<"Iterations: "<<mits<<std::endl;
    std::cout<<"Error: "<< merror <<std::endl;
    std::cout<< "Time elapse: "<< etime << "msec" << std::endl;

    /**************************Escribir a a archivo*******/
    std::ofstream myfile;
    myfile.open(ruta+".txt",std::ios::out);

    if (myfile.is_open()) {
        myfile<<"CGM report"<<std::endl;
        if(precond)
            myfile<<"Precondition with: "<< namep <<std::endl;
        myfile<<"Iterations: "<<mits<<std::endl;
        myfile<<"Error: "<< merror <<std::endl;
        myfile<< "Time elapse: "<< etime <<" msec" << std::endl;
    }
    myfile.close();
}
