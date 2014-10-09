#ifndef NEWTONGMRES_H
#define NEWTONGMRES_H

#include "math.h"
#include <iostream>
using namespace std;
#include "GMRES.h"

//determine whether a number is too small, correspoding matrix is singular
#ifndef singular
#define singular
const double SINGULAR=pow(0.1,10);
#endif

const int NGMRES_ITR_MAX=10;
const int NEWTON_ITR_MAX=80;

template<typename T, typename E>//class T need the function T(const E & x0,const E & x, const double & tol, E & y): compute y=\nabla a|_{x0} x , nabla a |_{x0} is appximated by center derivative
void GMRES(const T & a, const E & x0, const E & b, const int & m, const double & tol, E & x) //solve  \nabla a|_{x0} x =b
{
    x=0.;
    E res(b);
    const double r0=norm(b);
    double err=norm(res)/r0;
    E * v= new E [m+1];

    double *c=new double [m];
    double *s=new double [m];
    double temp[2];

    int k;
    double **H=new double * [m+1];
    double *y=new double [m];
    double * beta=new double [m+1];
    for(int i=0;i<m+1;i++)
    {
        H[i]=new double [m];
        beta[i]=0;
        for(int j=0;j<m;j++)
        {
            H[i][j]=0.;
        }
    }
    int itr=0;
    while(err>tol)
    {
        k=0;
        itr++;
        beta[0]=norm(res);
        for(int i=1;i<m+1;i++)
        {
            beta[i]=0;
        }
        v[0]=res/beta[0];

        for(int i=0;i<m;i++)
        {
            k++;
            a(x0,v[i],tol,v[i+1]);
            for(int j=0;j<=i;j++)
            {
                H[j][i]=v[j]*v[i+1];
                v[i+1]-=H[j][i]*v[j];
            }
            H[i+1][i]=norm(v[i+1]);
            v[i+1]/=H[i+1][i];

            for(int j=0;j<i;j++)
            {
                temp[0]=c[j]*H[j][i]+s[j]*H[j+1][i];
                temp[1]=-s[j]*H[j][i]+c[j]*H[j+1][i];
                H[j][i]=temp[0];
                H[j+1][i]=temp[1];
            }
            c[i]=H[i][i]/sqrt(H[i][i]*H[i][i]+H[i+1][i]*H[i+1][i]);
            s[i]=H[i+1][i]/sqrt(H[i][i]*H[i][i]+H[i+1][i]*H[i+1][i]);
            temp[0]=c[i]*H[i][i]+s[i]*H[i+1][i];
            temp[1]=-s[i]*H[i][i]+c[i]*H[i+1][i];
            H[i][i]=temp[0];
            H[i+1][i]=temp[1];
            temp[0]=c[i]*beta[i]+s[i]*beta[i+1];
            temp[1]=-s[i]*beta[i]+c[i]*beta[i+1];
            beta[i]=temp[0];
            beta[i+1]=temp[1];
            if(fabs(beta[i+1])<SINGULAR)
            {
                break;
            }
        }
        for(int i=k-1;i>=0;i--)
        {
            for(int j=k-1;j>i;j--)
            {
                beta[i]-=H[i][j]*y[j];
            }
            y[i]=beta[i]/H[i][i];
            x+=v[i]*y[i];
        }
        a(x0,x,tol,res);
        res=b-res;
        err=norm(res)/r0;
        cout<<"        GMRES Itr number :"<<itr<<" error: "<<err<<endl;
        if(itr>NGMRES_ITR_MAX)
        {
            cout<<"        GMRES Iteration might not converge! The relative error norm is "<<err<<endl;
            break;
        }
    }
    delete [] v;
    delete [] y;
    delete [] beta;
    delete [] c;
    delete [] s;
    for(int i=0;i<m+1;i++)
    {
        delete [] H[i];
    }
    delete [] H;
}

template<typename T, typename E>
void NewtonGmres(T & a, const E & b, const int m, const double tol, E & x)
{
    x=0;
    E error;
    E dx(0.);
    while(norm(dx)<tol)
    {
        x=(error+x)/2.;
        a(x,error);
        error-=b;
        a(x,error,tol,dx);
    }
    int itr=0;
    while(norm(error)>tol)
    {
        itr++;
        GMRES( a,x,error,m,0.01,dx);
        x-=dx;
        a(x,error);
        error-=b;
        cout<<"        Newton itr number: "<< itr << " error :"<<norm(error)<<endl;
        if(itr>NEWTON_ITR_MAX)
        {
            cout<<"        Newton Iteration might not converge, the relative error norm is "<<norm(error)/norm(b)<<endl;
            break;
        }
    }
}

//solve f(x)=0 with error tolerance
template<typename E, typename T>
class NEWTONGMRES
{
public:
    NEWTONGMRES(T & f, E & x, const E & r,const double & t, const int & M):tol(t),rhs(r),m(M),func(f),app(x)
    {
        error=0;
        dx=0;
    }
    void operator()(const E &, E &);
    void run();
    void ExactLineSearch();
private:
    const double & tol;
    const E & rhs;
    const int & m;
    T & func;
    E & app;
    E appleft;
    E appright;
    E error;
    E errorleft;
    E errorright;
    E dx;
};

template<typename E, typename T>
void NEWTONGMRES<E,T>::operator()(const E & x, E & dx)
{
    E temp;
    E temp1;
    E temp2;
    double length=0.5;
    double error=1.;
    while(error>0.0001)
    {

        temp=app+x*0.25*tol*length;
        func(temp,temp1);
        temp=app-x*0.25*tol*length;
        func(temp,temp2);
        dx=(temp1-temp2)/(0.5*tol*length);

        temp=app+x*0.5*tol*length;
        func(temp,temp1);
        temp=app-x*0.5*tol*length;
        func(temp,temp2);
        temp=(temp1-temp2)/(tol*length);

        error=norm(temp-dx)/norm(dx);
        length/=2.;
    }
//    while(norm(dx)<0.001*norm(x))
//    {
//        length*=2.;
//
//        temp=app+x*0.25*tol*length;
//        func(temp,temp1);
//        temp=app-x*0.25*tol*length;
//        func(temp,temp2);
//        dx=(temp1-temp2)/(0.5*tol*length);
//
//    }
}

template<typename E, typename T>
void NEWTONGMRES<E,T>::ExactLineSearch()
{
    //cout<<"begin line search"<<endl;
    double length=1.;
    appleft=app-dx*length;
    appright=app;
    app-=dx*length*0.5;
    func(app,error);
    error-=rhs;
    func(appleft,errorleft);
    errorleft-=rhs;
    func(appright,errorright);
    errorright-=rhs;
    while(length>0.0000001)
    {
        //cout<<length<<"  "<<norm(errorleft)<<"   "<<norm(error)<<"   "<<norm(errorright)<<endl;
        if(norm(errorleft)<norm(error)&&norm(errorleft)<norm(errorright))
        {
            length*=0.5;
            app=appleft;
            appright=app+length*dx;
            appleft=app-length*dx;
            error=errorleft;
            func(appleft,errorleft);
            func(appright,errorright);
            errorleft-=rhs;
            errorright-=rhs;
        }
        else if(norm(errorright)<norm(error)&&norm(errorright)<norm(errorleft))
        {
            length*=0.5;
            app=appright;
            error=errorright;
            appleft=app-length*dx;
            appright=app+length*dx;
            func(appleft,errorleft);
            func(appright,errorright);
            errorleft-=rhs;
            errorright-=rhs;
        }
        else //if(norm(error)<norm(errorleft)&&norm(error)<norm(errorright))
        {
            length*=0.5;
            appleft=app-length*dx;
            appright=app+length*dx;
            func(appleft,errorleft);
            func(appright,errorright);
            errorleft-=rhs;
            errorright-=rhs;
        }
    }
}

template<typename E, typename T>
void NEWTONGMRES<E,T>::run()
{
    cout<<"  Begin Newton Iteration"<<endl;
    func(app,error);
    (*this)(error,dx);
    while(norm(dx)<SINGULAR)
    {
        app=(error+app)/2.;
        func(app,error);
        error-=rhs;
        (*this)(error,dx);
    }
    int itr=0;
    cout<<"    Newton itr number: "<< itr << " error :"<<normmax(error)<<endl;
    while(normmax(error)>tol)
    {
        itr++;
        GMRES( *this,error,m,0.001,dx);
        ExactLineSearch();
        cout<<"    Newton itr number: "<< itr << " error :"<<normmax(error)<<endl;
        if(itr>NEWTON_ITR_MAX)
        {
            cout<<"    Newton Iteration might not converge, the relative error norm is "<<norm(error)<<endl;
            break;
        }
    }
    cout<<"  End Newton Iteration"<<endl;
}


#endif
