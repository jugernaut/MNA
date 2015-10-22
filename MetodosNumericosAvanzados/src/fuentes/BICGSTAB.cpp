
#include "../include/BICGSTAB.hpp"
#include <iostream>


void BICGSTAB::report(std::string ruta) {
    std::cout<<"BICGSTAB report: "<<std::endl;
    if(precond)
        std::cout<<"Precondicionado con: "<< namep <<std::endl;
    std::cout<<"Iteraciones: "<<mits<<std::endl;
    std::cout<<"Error: "<< merror <<std::endl;
    std::cout<< "Tiempo transcurrido: "<< etime <<" msec" << std::endl;

    /**************************Escribir a a archivo*******/
    std::ofstream myfile;
    myfile.open(ruta+".txt",std::ios::out);

    if (myfile.is_open()) {
        myfile<<"BICGSTAB report"<<std::endl;
        if(precond)
            myfile<<"Precondicionado con: "<< namep <<std::endl;
        myfile<<"Iteraciones: "<<mits<<std::endl;
        myfile<<"Error: "<< merror <<std::endl;
        myfile<< "Tiempo transcurrido: "<< etime <<" msec" << std::endl;
    }
    myfile.close();
}
