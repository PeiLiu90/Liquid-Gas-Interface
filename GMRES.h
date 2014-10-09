#ifndef GMRES_H
#define GMRES_H

//GMRES method with restart
//GMRES_RESTART=m

#include "math.h"
#include <iostream>
using namespace std;

//determine whether a number is too small, correspoding matrix is singular
#ifndef singular
#define singular
const double SINGULAR=pow(0.1,10);
#endif

const int GMRES_ITR_MAX=10;

template<typename T, typename E>//class T need the function T(const E &,E &)
void GMRES(T & a, const E & b, const int m, const double tol, E & x) //solve a*x=b, where a is linear operator and b is an element in the range of a,
{
    x=0.;
    E res(b);
    const double r0=norm(b);
    double err=norm(res)/r0;
    double errold;
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
        errold=err;
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
            a(v[i],v[i+1]);
            for(int j=0;j<=i;j++)
            {
                H[j][i]=v[j]*v[i+1];
                v[i+1]-=H[j][i]*v[j];
            }
            H[i+1][i]=norm(v[i+1]);
            if(fabs(H[i+1][i])<SINGULAR)
            {
                x=v[i+1];
                cout<<"SINGULAR, NOT CONVERGE "<<endl;
                break;
            }
            else
            {
                k++;
                v[i+1]/=H[i+1][i];
            }
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
            //cout<<fabs(beta[i+1])<<"  "<<H[i][i]<<endl;
            if(fabs(beta[i+1])<tol*r0)
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
        a(x,res);
        res=b-res;
        err=norm(res)/r0;
        cout<<"        GMRES Itr number :"<<itr<<" error: "<<err<<endl;
        if(itr>GMRES_ITR_MAX||err>0.99*errold)
        {
            cout<<"        GMRES Iteration might not converge, the relative error norm is "<<err<<endl;
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

template<typename T, typename E>//class T need the function operator()(const E &,E &), P(E &)
void GMRES_P(T & a, const E & b, const int m, const double tol, E & x) //solve a*x=b, where a is linear operator and b is an element in the range of a,
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
        itr++;
        k=0;
        a.P(res);
        beta[0]=norm(res);
        for(int i=1;i<m+1;i++)
        {
            beta[i]=0;
        }
        v[0]=res/beta[0];

        for(int i=0;i<m;i++)
        {
            k++;
            a(v[i],v[i+1]);
            a.P(v[i+1]);
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
        a(x,res);
        res=b-res;
        err=norm(res)/r0;
        if(itr>GMRES_ITR_MAX)
        {
            cout<<"GMRES Iteration might not converge, the relative error norm is "<<err<<endl;
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
#endif // GMRES_H
