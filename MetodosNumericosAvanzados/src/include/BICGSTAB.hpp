#ifndef __BiCGSTAB__
#define __BiCGSTAB__
#include "../include/SolverIterativo.hpp"

class BICGSTAB : public SolverIterativo {
public:
    BICGSTAB(void): SolverIterativo() {
    }
    ~BICGSTAB() {
    }
    template <class T>
    void solve(T &A, std::vector<double> &x, std::vector<double> &b) {
        timer.tic();
        int k, n = x.size();
        double beta= 1 , rho_2 = 1, rho_1= 1, alpha = 1, omega = 1;
        double tol2 = mtol * mtol;
        std::vector<double> p(n), s(n), t(n), v(n), r(n), rtilde(n), tem(n), tem2(n);

        //Comienza el algoritmo
        //r = b-A*x;
        sub(b, A*x, r);
        rtilde = r;
        p = r;

        //ciclo principal
        for (k = 1; k <= mmaxIts; k++) {
            rho_1 = dotProduct(rtilde, r);

            beta = (rho_1/rho_2) * (alpha/omega);
            //p = r + beta * (p - omega * v);
            scalarTimesVector(omega, v, tem);  //omega*v
            sub(p, tem, tem);				   //p-omega*v
            scalarTimesVector(beta, tem, tem); //beta*(p-omega*v)
            sum(tem, r, p);					   //p=r+beta*(p-omega*v)

            //phat = M.solve(p); //precondicionador
            v = A*p;
            alpha = rho_1 / dotProduct(rtilde, v);
            //s = r - alpha * v;
            scalarTimesVector(alpha, v, tem); //alpha * v
            sub(r, tem, s);					  //s=r-alpha*v;

            //shat = M.solve(s); //precondicionador
            t = A*s;
            omega = dotProduct(t,s) / dotProduct(t,t);
            //x += alpha * p + omega * s;
            scalarTimesVector(alpha, p, tem);  //alpha * p
            scalarTimesVector(omega, s, tem2); //omega * s
            sum(tem, tem2, tem2);              //alpha * p + omega * s
            sum(x, tem2, x);                   //x = x + alpha * p + omega * s

            merror = dotProduct(s, s);
            if(merror < tol2) {
                break;
            }

            //r = s - omega * t;
            scalarTimesVector(omega, t, tem); //omega * t
            sub(s, tem, r);					  //r = s - omega * t

            rho_2 = rho_1;
        }
        mits = k + 1;
        merror = sqrt(merror);
        timer.toc();
        etime = timer.etime();
        precond = false;
    }

    template<class T, class U>//T: tipo de matriz de entrada, //U: tipo de precondicinador de entrada
    void solve(T &A, std::vector<double> &x, std::vector<double> &b,U &M) { //T: tipo de matriz de entrada,
        timer.tic();
        int k, n = x.size();
        double beta= 1 , rho_2 = 1, rho_1= 1, alpha = 1, omega = 1;
        double tol2 = mtol * mtol;
        std::vector<double> p(n), s(n), t(n), v(n), r(n), rtilde(n), tem(n), tem2(n), y(n), o(n);

        //Comienza el algoritmo
        //r = b-A*x;
        sub(b, A*x, r);
        rtilde = r;
        p = r;

        //ciclo principal
        for (k = 1; k <= mmaxIts; k++) {
            rho_1 = dotProduct(rtilde, r);

            beta = (rho_1/rho_2) * (alpha/omega);
            //p = r + beta * (p - omega * v);
            scalarTimesVector(omega, v, tem);  //omega*v
            sub(p, tem, tem);				   //p-omega*v
            scalarTimesVector(beta, tem, tem); //beta*(p-omega*v)
            sum(tem, r, p);					   //p=r+beta*(p-omega*v)

            //precondicionador
            M.solve(y,p);                         //resuleve el sistema My=p
            v = A*y;
            alpha = rho_1 / dotProduct(rtilde, v);
            //s = r - alpha * v;
            scalarTimesVector(alpha, v, tem); //alpha * v
            sub(r, tem, s);					  //s=r-alpha*v;

            //precondicionador
            M.solve(o,s);                         //resuleve el sistema Mo=s
            t = A*o;
            omega = dotProduct(t,s) / dotProduct(t,t);
            //x += alpha * p + omega * s;
            scalarTimesVector(alpha, y, tem);  //alpha * p
            scalarTimesVector(omega, o, tem2); //omega * s
            sum(tem, tem2, tem2);              //alpha * p + omega * s
            sum(x, tem2, x);                   //x = x + alpha * p + omega * s

            merror = dotProduct(s, s);
            if(merror < tol2) {
                break;
            }

            //r = s - omega * t;
            scalarTimesVector(omega, t, tem); //omega * t
            sub(s, tem, r);					  //r = s - omega * t

            rho_2 = rho_1;
        }
        mits = k + 1;
        merror = sqrt(merror);
        timer.toc();
        etime = timer.etime();
        precond = false;
    }

    void report(std::string);
};
#endif
