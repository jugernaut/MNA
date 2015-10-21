#include "../include/Matrix_dense.hpp"
#include "../include/Matrix_COO.hpp"
#include "../include/Matrix_CRS.hpp"
#include "../include/Matrix_CSC.hpp"
#include "../include/MatrizCDS.hpp"
#include "../include/CGM.hpp"
#include "../include/Jacobi.hpp"
#include "../include/BICGSTAB.hpp"
#include "../include/JacobiPreconditioner.hpp"
#include "../include/ILU.hpp"
#include "../include/MILU.hpp"
#include <iostream>
#include <string>


using namespace std;


int main(int argc, char const *argv[]) {
    Matrix_dense dense;
    int n=4;
    dense.inicializa(n,n,0.0);
    std::vector<double> v,v1;
    std::vector<double> r;

    v= {1,2,3,4};
    dense.A= {{10,-1,-2,0},{-1,11,-1,3},{-2,-1,10,-1},{0,3,-1,8}};

    cout<<"Dense matrix: "<<endl;
    dense.print();

    cout<<endl<<"Dense matrix x vector: "<<endl;
    r=dense*v;
    printVector(r);

    cout<<endl<<"COO matrix: "<<endl;
    Matrix_COO coo(n);
    coo.convert(dense);
    coo.print();
    std::vector<double> cooXv;
    cooXv=coo*v;
    cout<<"COO x v: "<<endl;
    printVector(cooXv);


    Matrix_CSC csc(n);
    csc.convert(coo);
    cout<<endl<<"CSC Matrix: "<<endl;
    csc.print();
    std::vector<double> cscXv;
    cscXv=csc*v;
    cout<<endl<<"csc x v: "<<endl;
    printVector(cscXv);

    Matrix_CRS crs;
    //crs.convertCCStoCRS(csc);
    crs.convert(coo);
    cout<<endl<<"CRS matrix: "<<endl;
    crs.print();
    printVector(crs.idiag);
    std::vector<double> crsXv;
    crsXv=crs*v;
    cout<<endl<<"crs x v: "<<endl;
    printVector(crsXv);

    std::vector<double> x,b;
    b.resize(n,0.0);
    x.resize(n,0.0);

    b[0]=6.0;
    b[1]=25.0;
    b[2]=-11;
    b[3]=15.0;

    ////////////////////////////////////////////////////SOLVERS/////////////////////////////////////////
    cout<<endl<<"----------------------SOLVERS----------------------------"<<endl;

    cout<<"b: ";
    printVector(b);

    Jacobi jac;
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    jac.solve(crs,x,b);
    jac.report("prueba");
    printVector(x);
    cout<<endl<<endl;

    CGM cgm;
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    cgm.solve(crs,x,b);
    cgm.report("pruebaCGM");
    printVector(x);
    cout<<endl<<endl;


    BICGSTAB big;
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    x.resize(n,0.0);
    big.solve(crs,x,b);
    big.report("pruebaBIG");
    printVector(x);
    //big.reset(1);
    cout<<endl<<endl;

    //////////////////////////////////////////////////precondicionadores////////////////////////////
    cout<<endl<<"----------------------Preconditioners  CGM Jacobi----------------------------"<<endl;
    JacobiPreconditioner<Matrix_CRS> jacp2;
    jacp2.calculate(crs);
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    //CGM cgm2;
    cgm.solve(crs,x,b,jacp2);
    cgm.report("pruebaCGMPReJac");
    printVector(x);
    cout<<endl<<endl;

    cout<<endl<<"----------------------Preconditioners  ilu CGM----------------------------"<<endl;

    ILU<Matrix_CRS> ilu;
    ilu.calculate(crs);
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    CGM cgm3;
    cgm3.solve(crs,x,b,ilu);
    cgm3.report("pruebaCGMPReILU");
    printVector(x);
    cout<<endl<<endl;

    cout<<endl<<"----------------------Preconditioners  milu CGM----------------------------"<<endl;
    MILU<Matrix_CRS> milu;
    milu.calculate(crs);
    milu.mat.print();
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    CGM cgm1;
    cgm1.solve(crs,x,b,milu);
    cgm1.report("pruebaCGMPReMILU");
    printVector(x);
    cout<<endl<<endl;

    cout<<endl<<"----------------------Preconditioners  BICGSTAB Jacobi----------------------------"<<endl;
    JacobiPreconditioner<Matrix_CRS> jacpBIG;
    jacpBIG.calculate(crs);
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    BICGSTAB big2;
    big2.solve(crs,x,b,jacpBIG);
    big2.report("pruebaBICGPReJac");
    printVector(x);
    cout<<endl<<endl;

    cout<<endl<<"----------------------Preconditioners  BICGSTAB ILU----------------------------"<<endl;

    ILU<Matrix_CRS> iluBIG;
    iluBIG.calculate(crs);
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    BICGSTAB big3;
    big3.solve(crs,x,b,iluBIG);
    big3.report("pruebaCGMPReILU");
    printVector(x);
    cout<<endl<<endl;

    cout<<endl<<"----------------------Preconditioners  BICGSTAB MILU----------------------------"<<endl;


    MILU<Matrix_CRS> miluBIG;
    miluBIG.calculate(crs);
    x[0]=.6;
    x[1]=2.2727;
    x[2]=-1.1;
    x[3]=1.875;
    BICGSTAB big4;
    big4.solve(crs,x,b,miluBIG);
    big4.report("pruebaCGMPReMILU");
    printVector(x);
    cout<<endl<<endl;

    return 0;
}
