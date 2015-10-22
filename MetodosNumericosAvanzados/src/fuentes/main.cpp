#include "../include/MatrizDensa.hpp"
#include "../include/MatrizCOO.hpp"
#include "../include/MatrizCRS.hpp"
#include "../include/MatrizCSC.hpp"
#include "../include/MatrizCDS.hpp"
#include "../include/CGM.hpp"
#include "../include/Jacobi.hpp"
#include "../include/BICGSTAB.hpp"
#include "../include/JacobiPrecondicionado.hpp"
#include "../include/ILU.hpp"
#include "../include/MILU.hpp"
#include "../include/ICHOL.hpp"
#include <iostream>
#include <string>


using namespace std;

   /*{
    Timer timer;                   //mide tiempo de ejecucion

   //llena matriz COO
   cout<< endl << "Tamanio de problema " << n << "x" <<n<<endl<<endl;
   int l = 0;
      COO Acoo(n);                 //matriz temporal en formato de coordenadas
      timer.tic();                 //comienza a medir tiempo
      for(int j=1;j<ny;++j)
      {
         for(int i=1;i<nx;++i)
         {
             // prototipo de funcion de insertar void insert(int i,int j,doble val);
             if(j>1)
                Acoo.insert(l, l-(nx-1),-1.);
             if(i>1)
                Acoo.insert(l, l-1,-1.);
             Acoo.insert(l, l,4.);
             if(i<nx-1)
                Acoo.insert(l, l+1,-1.);
             if(j<ny-1)
                Acoo.insert(l, l+(ny-1),-1.);
             ++l;
         }
      }
      timer.toc();                 //termina de medir tiempo
      std::cout << "Tiempo de llenado de matriz   COO: " << timer.etime() << " ms" << std::endl;

      timer.tic();
      A.convert(Acoo);             //convierte matriz COO a formato CSR
      timer.toc();
      std::cout << "Tiempo de conversion de COO a CSR " <<timer.etime() << " ms" << std::endl;
   }//destruye Acoo*/



int main(int argc, char const *argv[]) {
    MatrizDensa dense;
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
    MatrizCOO coo(n);
    coo.convert(dense);
    coo.print();
    std::vector<double> cooXv;
    cooXv=coo*v;
    cout<<"COO x v: "<<endl;
    printVector(cooXv);


    MatrizCSC csc(n);
    csc.convert(coo);
    cout<<endl<<"CSC Matrix: "<<endl;
    csc.print();
    std::vector<double> cscXv;
    cscXv=csc*v;
    cout<<endl<<"csc x v: "<<endl;
    printVector(cscXv);

    MatrizCRS crs;
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

    ////////////////////////////////////////////////////SolverS/////////////////////////////////////////
    cout<<endl<<"----------------------SolverS----------------------------"<<endl;

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
    cout<<endl<<"----------------------Precondicionadors  CGM Jacobi----------------------------"<<endl;
    JacobiPrecondicionado<MatrizCRS> jacp2;
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

    cout<<endl<<"----------------------Precondicionadors  ilu CGM----------------------------"<<endl;

    ILU<MatrizCRS> ilu;
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

    cout<<endl<<"----------------------Precondicionadors  milu CGM----------------------------"<<endl;
    MILU<MatrizCRS> milu;
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

    cout<<endl<<"----------------------Precondicionadors  CGM ICHOL----------------------------"<<endl;
	ICHOL<MatrizCRS> ichol;
	ichol.calculate(crs);
	x[0]=.6;
	x[1]=2.2727;
	x[2]=-1.1;
	x[3]=1.875;
	CGM cgm4;
	cgm4.solve(crs,x,b,ichol);
	cgm4.report("pruebaCGMPReJac");
	printVector(x);
	cout<<endl<<endl;

    cout<<endl<<"----------------------Precondicionadors  BICGSTAB Jacobi----------------------------"<<endl;
    JacobiPrecondicionado<MatrizCRS> jacpBIG;
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

    cout<<endl<<"----------------------Precondicionadors  BICGSTAB ILU----------------------------"<<endl;

    ILU<MatrizCRS> iluBIG;
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

    cout<<endl<<"----------------------Precondicionadors  BICGSTAB MILU----------------------------"<<endl;


    MILU<MatrizCRS> miluBIG;
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

    cout<<endl<<"----------------------Precondicionadors  BICGSTAB ICHOL----------------------------"<<endl;
	ICHOL<MatrizCRS> icholBIG;
	icholBIG.calculate(crs);
	x[0]=.6;
	x[1]=2.2727;
	x[2]=-1.1;
	x[3]=1.875;
	BICGSTAB big5;
	big5.solve(crs,x,b,icholBIG);
	big5.report("pruebaCGMPReJac");
	printVector(x);
	cout<<endl<<endl;

    return 0;
}
