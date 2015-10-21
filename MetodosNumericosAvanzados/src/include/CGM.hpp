#ifndef __CGM__
#define __CGM__

#include "../include/solverIter.hpp"

class CGM : public solverIter {
public:
    CGM(void):solverIter() {
    }
    ~CGM() {
    }

    template <class T>
    void solve(T &A, std::vector<double> &x, std::vector<double> &b) {
        std::cout<<"CGM solve"<<std::endl;
        timer.tic();
        double tol2 = mtol*mtol;
        int n = x.size();
        mits=0;

        /*Vector r(n);
        Vector p(n);
        Vector Ap(n);*/

        std::vector<double> r;
        r.resize(n);
        std::vector<double> p;
        p.resize(n);
        std::vector<double> Ap;
        Ap.resize(n);

        double alpha;
        double beta;
        double rr0,rr;
        int k;

        //r = b-A*x;
        sub(b,A*x,r);

        //rr0 = r*r;
        rr0=dotProduct(r,r);

        p=r;

        for(k =0; k< mmaxIts; ++k) {
            Ap = A*p;
            //alpha =  (rr0) / (Ap*p);
            alpha=(rr0)/dotProduct(Ap,p);

            //x = x + alpha*p;
            scalarTimesVector(alpha,p,p);
            sum(x,p,x);

            //r = r - alpha*Ap;
            scalarTimesVector(alpha,Ap,Ap);
            sub(r,Ap,r);

            //rr = r*r;
            rr=dotProduct(r,r);
            merror = rr;
            if(merror < tol2)
                break;
            beta = rr / rr0;
            //p = r + beta*p;
            scalarTimesVector(beta,p,p);
            sum(r,p,p);
            rr0 = rr;
        }

        mits = k+1;
        merror = sqrt(merror);
        timer.toc();
        etime = timer.etime();
        precond = false;
    }

    template<class T, class U>
    void solve(T &A, std::vector<double> &x, std::vector<double> &b,U &M)  //T: tipo de matriz de entrada,
    //U: tipo de precondicinador de entrada
    {
        mits=0;
        timer.tic();                 //comienza a medir tiempo de ejecucion
        double tol2 = mtol*mtol;
        int n = x.size();

        std::vector<double> r(n),z(n),p(n),Ap(n);


        double alpha;
        double beta;
        double rr0,rr;
        int k;

        //r = b-A*x;                    //primer residual operacion vector = vector - matriz*vector
        sub(b,A*x,r);
        //printVector(r);
        M.solve(z,r);                         //resuleve el sistema Mz=r
        //printVector(z);
        p=z;                          //copia de vector

        //rr0 = z*r;                    //producto interno de dos vectores
        rr0=dotProduct(z,r);

        //for(k =0;k<0;++k)
        for(k =0; k<mmaxIts; ++k) {
            //std::cout<<"iteracion: "<<k<<std::endl;
            Ap = A*p;                  //producto matriz vector

            //alpha =  (rr0) / (Ap*p);   //producto interno en el divisor
            alpha=(rr0)/dotProduct(Ap,p);

            //x = x + alpha*p;           //operacion vector = vector + escalar*vector (saxpy en BLAS)
            scalarTimesVector(alpha,p,p);
            sum(x,p,x);

            //r = r - alpha*Ap;
            scalarTimesVector(alpha,Ap,Ap);
            sub(r,Ap,r);

            M.solve(z,r);              //resuleve el sistema Mz=r

            //rr = z*r;                  //producto interno
            rr=dotProduct(z,r);

            //merror = r*r;              //producto interno
            merror=dotProduct(r,r);

            if(merror < tol2)
                break;
            beta = rr / rr0;
            //p = z + beta*p;            //operacion vector = vector + escalar*vector
            scalarTimesVector(beta,p,p);
            sum(z,p,p);

            rr0 = rr;
        }
        mits = k+1;
        merror = sqrt(merror);
        timer.toc();                 //termina de medir tiempo
        etime = timer.etime();       //guarda en etime el tiempo de ejecucion
        precond = true;
        namep = M.name();            //guarda el nombre del precondicionador
    }

    void report(std::string);
};

#endif
