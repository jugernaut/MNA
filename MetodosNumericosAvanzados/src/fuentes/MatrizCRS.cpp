#include "../include/MatrizCRS.hpp"
#include <iostream>
#include <cmath>

MatrizCRS::MatrizCRS(int n) {
    nnz=0;
    Col=Row=n;
}

void MatrizCRS::initialize() {

}
void MatrizCRS::print() {
    std::cout<<std::endl<<"Data: [";
    for (std::vector<double>::iterator it = data.begin() ; it != data.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"Col: [";
    for (std::vector<int>::iterator it = col.begin() ; it != col.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

    std::cout<<"irow: [";
    for (std::vector<int>::iterator it = irow.begin() ; it != irow.end(); ++it)
        std::cout <<""<< *it<<", ";
    std::cout<<" ]"<<std::endl;

}
void MatrizCRS::convert(const MatrizCOO& coo) {
    MatrizCSC aux(coo.Row);
    aux.convert(coo);
    convertCCStoCRS(aux);
}

void MatrizCRS::convertCCStoCRS(const MatrizCSC& csc) {
    nnz=csc.nnz;
    Col=Row=csc.Col;
    data.resize(nnz);
    col.resize(nnz);
    irow.resize(Col+1,0);
    for (int l = 0; l < nnz; ++l) {
        irow[csc.row[l]]++;
    }
    int sum=0;
    int tmp;
    for (int i = 0; i < Col+1; ++i) {
        tmp=irow[i];
        irow[i]=sum;
        sum+=tmp;
    }
    irow[Col]=nnz;

    std::vector<int> irowAux;
    irowAux=irow;
    int i,j,d,k;
    k=0;
    double val;
    for (int l = 0; l < nnz; ++l) {
        i=csc.row[l];
        if (l >= csc.icol[k] && l < csc.icol[k+1]) {
            j=k;
            if(l==csc.icol[k+1]-1)
                k++;
        }
        val=csc.data[l];
        d=irowAux[i];
        data[d]=val;
        col[d]=j;
        irowAux[i]++;
    }
    idiagCalculate();
}

std::vector<double> MatrizCRS::operator * (const std::vector<double>& V) {
    if (Row!=V.size()) {
        std::cout<<"Error en las dimensiones"<<std::endl;
        exit (1);
    }
    std::vector<double> y;
    y.resize(Row,0.0);
    for (int i = 0; i < Row; ++i) {
        for (int l = irow[i]; l <= irow[i+1]-1; ++l) {
            y[i]+=data[l]*V[col[l]];
        }
    }
    return y;
}

MatrizCRS const& MatrizCRS::operator=(MatrizCRS const &rhs) {
    if (this != &rhs) { // && Row==rhs.Row && Col==rhs.Col)
        Col=rhs.Col;
        Row=rhs.Row;
        data=rhs.data;
        col=rhs.col;
        irow=rhs.irow;
        idiag=rhs.idiag;
        nnz=rhs.nnz;
    }
    return *this;
}

void MatrizCRS::JacobiIter(std::vector<double>& x0, const std::vector<double>& b) {
    std::vector<double> x;
    x=b;
    int j;
    double diag;
    for (int i = 0; i < Row; ++i) {
        for (int l = irow[i]; l <= irow[i+1]-1; ++l) {
            j=col[l];
            if (i!=j) {
                x[i]-=data[l]*x0[j];
            } else {
                diag=data[l];
            }
        }
        x[i]/=diag;
    }
    x0=x;
}

void MatrizCRS::idiagCalculate() {
    int j;
    for (int i = 0; i < Row; ++i) {
        for (int l = irow[i]; l <= irow[i+1]-1; ++l) {
            j=col[l];
            if(j==i) {
                idiag.push_back(l);
            }
        }
    }
}

MatrizCRS MatrizCRS::diag() {
    MatrizCRS D;
    for (int i = 0; i < Row; ++i) {
        D.col.push_back(i);
        D.irow.push_back(i);
        D.data.push_back(data[idiag[i]]);
    }
    D.irow.push_back(Row+1);
    D.Row=D.Col=Row;
    D.idiagCalculate();
    return D;
}

MatrizCRS MatrizCRS::ILU() {
    //aqui se copia A a LU
    MatrizCRS LU;
    LU = *this;
    int k,lli,llk,ji,jk;
    double Lik;
    for (int i=0; i < Row; ++i) {
        for (int l = LU.irow[i]; l < LU.idiag[i]; l++ ) {
            k=LU.col[l];
            Lik=LU.data[l]/LU.data[idiag[k]];
            lli=LU.irow[i];
            llk=LU.irow[k];
            ji=LU.col[lli];
            jk=LU.col[llk];
            while (lli<LU.irow[i+1] && llk<LU.irow[k+1]) {
                if(ji==jk) {
                    if(ji>k && jk > k) {
                        LU.data[lli]=LU.data[lli]-Lik*LU.data[llk];
                    }
                    lli++;
                    llk++;
                    ji=LU.col[lli];
                    jk=LU.col[llk];
                } else if (ji<jk) {
                    lli++;
                    ji=LU.col[lli];
                } else {
                    llk++;
                    jk=LU.col[llk];
                }
            }
            LU.data[l]=Lik;
        }
    }
    return LU;
}

MatrizCRS MatrizCRS::MILU() {
    //aqui se copia A a LU
    MatrizCRS LU;
    LU = *this;
    int k,lli,llk,ji,jk;
    double Lik;
    for (int i=0; i < Row; ++i) {
        for (int l = LU.irow[i]; l < LU.idiag[i]; l++ ) {
            k=LU.col[l];
            Lik=LU.data[l]/LU.data[idiag[k]];
            lli=LU.irow[i];
            llk=LU.irow[k];
            ji=LU.col[lli];
            jk=LU.col[llk];
            while (lli<LU.irow[i+1] && llk<LU.irow[k+1]) {
                if(ji==jk) {
                    if(ji>k && jk > k) {
                        LU.data[lli]=LU.data[lli]-Lik*LU.data[llk];
                        LU.data[idiag[i]]=LU.data[idiag[i]]-Lik*LU.data[llk];
                    }
                    lli++;
                    llk++;
                    ji=LU.col[lli];
                    jk=LU.col[llk];
                } else if (ji<jk) {
                    lli++;
                    ji=LU.col[lli];
                } else {
                    llk++;
                    jk=LU.col[llk];
                }
            }
            LU.data[l]=Lik;
        }
    }
    return LU;
}

void MatrizCRS::JacobiSolve(std::vector<double> & z, const std::vector<double>& r) {
    for (int i = 0; i < Row; ++i) {
        //std::cout<<r[i]<<"/"<<data[i]<<std::endl;
        z[i]=r[i]/data[idiag[i]];

    }
}

void MatrizCRS::LUSolve(std::vector<double> & z, const std::vector<double>& r) {
    std::vector<double> y;
    y.resize(Row);

    for (int i = 0; i < Row; ++i) {
        y[i]=r[i];
        for (int l = irow[i]; l < idiag[i]; ++l) {
            y[i]=y[i]-data[l]*y[col[l]];
        }
    }
    for (int i = Row-1; i >= 0; --i) {
        z[i]=y[i];
        for (int l = idiag[i]+1; l <= irow[i+1]-1 ; ++l) {
            z[i]-=data[l]*z[col[l]];
        }
        z[i]/=data[idiag[i]];
    }
}

int MatrizCRS::transposedElementIndexSearchInCRS(int i, int j){
    if (i==j) return idiag[i];
   int lside, rside;
   int pivote;
   lside = irow[j];
   rside = idiag[j];
   while (lside != rside){
        pivote = (lside+rside)/2;
        if(i==col[pivote])
            return pivote;
        else if (i < col[pivote])
            rside = pivote;
        else
            lside = pivote;
   }
   return -1;
}

MatrizCRS MatrizCRS::ICHOL(){
    MatrizCRS L;
    L= *this;
    int j,lli,llj,ji,jj;
    double dot;
    L.data[idiag[0]]=sqrt(data[idiag[0]]);
    for (int i=0; i < Row; ++i){
        for (int l = L.irow[i]; l < L.idiag[i];l++ ){
            j=L.col[l];
            lli=L.irow[i];
            llj=L.irow[j];
            ji=L.col[lli];
            jj=L.col[llj];
            L.data[l]=data[l];
            while (lli<L.irow[i+1] && llj<L.irow[j+1]){
                if(ji==jj){
                    if(ji<j && jj < j){
                        L.data[l]=L.data[l]-L.data[lli]*L.data[llj];//hazard cambiar l por lli
                    }
                    lli++;
                    llj++;
                    ji=L.col[lli];
                    jj=L.col[llj];
                }
                else if (ji<jj){
                    lli++;
                    ji=L.col[lli];
                }
                else{
                    llj++;
                    jj=L.col[llj];
                }
            }
            L.data[l]=L.data[l]/L.data[idiag[j]];
        }
        dot=0;
        for (int l = L.irow[i]; l < L.idiag[i];l++ )
        {
            dot=dot+L.data[l]*L.data[l];
        }
        //std::cout<<data[idiag[i]]-dot<<"\n";
        L.data[idiag[i]]=sqrt(data[idiag[i]]-dot);
    }
    return L;
}

void MatrizCRS::CHOLSolve(std::vector<double> & z, const std::vector<double>& r){
    std::vector<double> y;
    y.resize(Row);
    //Forward substitution
    for (int i = 0; i < Row; ++i){
        y[i]=r[i];
        for (int l = irow[i]; l < idiag[i]; ++l){
            y[i]-=data[l]*y[col[l]];
        }
        y[i]/=1.0;
    }

    int l;
    for (int i = Row-1; i >= 0; i--){
        z[i]=y[i];
        for (int j = i+1; j < Row ; j++){
            l=transposedElementIndexSearchInCRS(i,j);
            if (l != -1)
                z[i]-=data[l]*z[j];
        }
        z[i]/=data[idiag[i]];
    }
}

