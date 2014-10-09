#ifndef LINEAREQUATION_H
#define LINEAREQUATION_H

#include "math.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

/*
Function: Solve Linear Equation Ax=b using Jacobi Iteration

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void Jacobi(double **a, const double *b, const int n, const double tol, double *x)
{
    double err=0.;
    double *res=new double [n];
    //#pragma omp parallel for
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
    }
    int ite_num=0;
    //FILE *eps=fopen("Jacobi.plt","w+");
    while(err>tol)
    {
        ite_num++;
        for(int i=0;i<n;i++)
        {
            x[i]+=res[i]/a[i][i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        //fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<err<<endl;
    }
    delete [] res;
}

/*
Function: Solve Linear Equation Ax=b using Gauss Iteration

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void GaussSeidel(double **a, const double *b, const int n, const double tol, double *x)
{
    double err=0.;
    double *res=new double [n];
    double *xold=new double [n];
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
    }
    int ite_num=0;
    //FILE *eps=fopen("GaussSeidel.plt","w+");
    while(err>tol)
    {
        ite_num++;
        for(int i=0;i<n;i++)
        {
            xold[i]=x[i];
            x[i]=xold[i]*a[i][i]+res[i];
            for(int j=0;j<i;j++)
            {
                x[i]+=(xold[j]-x[j])*a[i][j];
            }
            x[i]/=a[i][i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        //fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] xold;
}

/*
Function: Solve Linear Equation Ax=b using SOR

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' , Relaxation factor 'omega', solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void SOR(double **a, const double *b, const int n, const double tol, const double omega, double *x)
{
    double err=0.;
    double *res=new double [n];
    double *xold=new double [n];
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
    }
    int ite_num=0;
    FILE *eps=fopen("SOR.plt","w+");
    while(err>tol)
    {
        ite_num++;
        for(int i=0;i<n;i++)
        {
            xold[i]=x[i];
            x[i]=xold[i]*a[i][i]+res[i]*omega;
            for(int j=0;j<i;j++)
            {
                x[i]+=(xold[j]-x[j])*a[i][j]*omega;
            }
            x[i]/=a[i][i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] xold;
}

/*
Function: Solve Linear Equation Ax=b using Steepest Decent

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void SteepestDecent(double **a, const double *b, const int n, const double tol, double *x)
{
    double err=0.;
    double length;
    double U=0;
    double D=0;
    double *res=new double [n];
    double *ARes=new double [n];
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
        U+=res[i]*res[i];
    }
    for(int i=0;i<n;i++)
    {
        ARes[i]=0;
        for(int j=0;j<n;j++)
        {
            ARes[i]+=a[i][j]*res[j];
        }
        D+=ARes[i]*res[i];
    }
    length=U/D;
    int ite_num=0;
    FILE *eps=fopen("SteepestDecent.plt","w+");
    while(err>tol)
    {
        ite_num++;
        for(int i=0;i<n;i++)
        {
            x[i]+=length*res[i];
        }
        err=0;
        U=0;
        D=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
            U+=res[i]*res[i];
        }
        //#pragma omp parallel for
        for(int i=0;i<n;i++)
        {
            ARes[i]=0;
            for(int j=0;j<n;j++)
            {
                ARes[i]+=a[i][j]*res[j];
            }
            D+=ARes[i]*res[i];
        }
        length=U/D;
        fprintf(eps,"%d   %.15f\n",ite_num, err);
        //cout<<ite_num<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] ARes;
}

/*
Function: Solve Linear Equation Ax=b using Conjugate Gradient

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void CG(double **a, const double *b, const int n, const double tol, double *x)
{
    double err=0.;
    double length;
    double glength;
    double U=0;
    double D=0;
    double *res=new double [n];
    double *gra=new double [n];
    double *AGra=new double [n];
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        gra[i]=res[i];
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
        U+=res[i]*res[i];
    }
    for(int i=0;i<n;i++)
    {
        AGra[i]=0;
        for(int j=0;j<n;j++)
        {
            AGra[i]+=a[i][j]*gra[j];
        }
        D+=AGra[i]*gra[i];
    }
    length=U/D;
    int ite_num=0;
    FILE *eps=fopen("CG.plt","w+");
    while(err>tol)
    {
        ite_num++;
        for(int i=0;i<n;i++)
        {
            x[i]+=length*gra[i];
        }
        err=0;
        D=0;
        for(int i=0;i<n;i++)
        {
            res[i]-=length*AGra[i];
        }
        for(int i=0;i<n;i++)
        {
            D+=res[i]*res[i];
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        glength=D/U;
        for(int i=0;i<n;i++)
        {
            gra[i]=glength*gra[i]+res[i];
        }
        U=D;
        D=0;
        for(int i=0;i<n;i++)
        {
            AGra[i]=0;
            for(int j=0;j<n;j++)
            {
                AGra[i]+=a[i][j]*gra[j];
            }
            D+=AGra[i]*gra[i];
        }
        length=U/D;
        fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] gra;
    delete [] AGra;
}

/*
Function: Solve Linear Equation Ax=b using LU factorization

Imput: Matrix 'A', problem size 'n', solution 'x' which is initialized to be the RHS

Output: solution 'x', Note that 'a' will be changed in this subprogram
*/

void Gauss(double **a,double *x,int n)
{
    int pv;
    double t;
    for(int i=0;i<n;i++)
    {
        pv=i;
        for(int j=pv+1;j<n;j++)
        {
            if(fabs(a[j][i])>fabs(a[pv][i]))
            {
                pv=j;
            }
        }
        if(pv!=i)
        {
            for(int j=i;j<n;j++)
            {
                t=a[i][j];
                a[i][j]=a[pv][j];
                a[pv][j]=t;
            }
            t=x[i];
            x[i]=x[pv];
            x[pv]=t;
        }
        for(int j=i+1;j<n;j++)
        {
            t=a[j][i]/a[i][i];
            for(int k=i;k<n;k++)
            {
                a[j][k]-=t*a[i][k];
            }
            x[j]-=t*x[i];
        }
    }
    for(int i=n-1;i>-1;i--)
    {
        x[i]/=a[i][i];
        for(int j=i-1;j>-1;j--)
        {
            x[j]-=a[j][i]*x[i];
        }
    }
}

void ChebyshevJacobi(double **a, const double *b, const int n, const double tol, const double rho, double *x)
{
    double err=0.;
    double *res=new double [n];
    double *xold=new double [n];
    double *xnew=new double [n];
    double t=(1.-sqrt(1.-rho*rho))/rho;
    double mu0=1;
    double mu1=1./rho;
    double mu2;
    //#pragma omp parallel for
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
        xold[i]=x[i];
        x[i]+=res[i]/a[i][i];
    }
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
    }
    int ite_num=0;
    //FILE *eps=fopen("Jacobi.plt","w+");
    while(err>tol)
    {
        ite_num++;
        mu2=(2./rho)*mu1-mu0;
        mu0=mu1;
        mu1=mu2;
        for(int i=0;i<n;i++)
        {
            xnew[i]=x[i]+res[i]/a[i][i];
            xnew[i]*=2.*mu0/(rho*mu1);
            xnew[i]+=(1.-2.*mu0/(rho*mu1))*xold[i];
            xold[i]=x[i];
            x[i]=xnew[i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        mu0*=t;
        mu1*=t;
        mu2*=t;
        //fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<mu1<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] xold;
    delete [] xnew;
}

/*
Function: Solve Linear Equation Ax=b using Gauss Iteration

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void ChebyshevGaussSeidel(double **a, const double *b, const int n, const double tol,const double rho, double *x)
{
    double err=0.;
    double *res=new double [n];
    double *xold=new double [n];
    double *xnew=new double [n];
    double t=(1.-sqrt(1.-rho*rho))/rho;
    double mu0=1;
    double mu1=1./rho;
    double mu2;
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
    }
    for(int i=0;i<n;i++)
    {
        xold[i]=x[i];
        x[i]=xold[i]*a[i][i]+res[i];
        for(int j=0;j<i;j++)
        {
            x[i]+=(xold[j]-x[j])*a[i][j];
        }
        x[i]/=a[i][i];
    }
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
    }
    int ite_num=0;
    //FILE *eps=fopen("GaussSeidel.plt","w+");
    while(err>tol)
    {
        ite_num++;
        mu2=(2./rho)*mu1-mu0;
        mu0=mu1;
        mu1=mu2;
        for(int i=0;i<n;i++)
        {
            xnew[i]=x[i]*a[i][i]+res[i];
            for(int j=0;j<i;j++)
            {
                xnew[i]+=(x[j]-xnew[j])*a[i][j];
            }
            xnew[i]/=a[i][i];
        }
        for(int i=0;i<n;i++)
        {
            xnew[i]*=2.*mu0/(rho*mu1);
            xnew[i]+=(1.-2.*mu0/(rho*mu1))*xold[i];
            xold[i]=x[i];
            x[i]=xnew[i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        mu0*=t;
        mu1*=t;
        mu2*=t;
        //fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<mu2<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] xold;
}

/*
Function: Solve Linear Equation Ax=b using Gauss Iteration

Imput: Matrix 'A', RHS 'b', problem size 'n', tolerance 'tol' ,solution 'x' which should be initialized to be an initial guess of the iteration

Output: solution 'x'
*/
void ChebyshevSOR(double **a, const double *b, const int n, const double tol,const double omega, const double rho, double *x)
{
    double err=0.;
    double *res=new double [n];
    double *xold=new double [n];
    double *xnew=new double [n];
    double t=(1.-sqrt(1.-rho*rho))/rho;
    double mu0=1;
    double mu1=1./rho;
    double mu2;
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
        if(err<fabs(res[i]))
        {
            err=fabs(res[i]);
        }
    }
    for(int i=0;i<n;i++)
    {
        xold[i]=x[i];
        x[i]=xold[i]*a[i][i]+res[i]*omega;
        for(int j=0;j<i;j++)
        {
            x[i]+=(xold[j]-x[j])*a[i][j]*omega;
        }
        x[i]/=a[i][i];
    }
    for(int i=0; i<n; i++)
    {
        res[i]=b[i];
        for(int j=0;j<n;j++)
        {
            res[i]-=a[i][j]*x[j];
        }
    }
    int ite_num=0;
    //FILE *eps=fopen("GaussSeidel.plt","w+");
    while(err>tol)
    {
        ite_num++;
        mu2=(2./rho)*mu1-mu0;
        mu0=mu1;
        mu1=mu2;
        for(int i=0;i<n;i++)
        {
            xnew[i]=x[i]*a[i][i]+res[i]*omega;
            for(int j=0;j<i;j++)
            {
                xnew[i]+=(x[j]-xnew[j])*a[i][j]*omega;
            }
            xnew[i]/=a[i][i];
        }
        for(int i=0;i<n;i++)
        {
            xnew[i]*=2.*mu0/(rho*mu1);
            xnew[i]+=(1.-2.*mu0/(rho*mu1))*xold[i];
            xold[i]=x[i];
            x[i]=xnew[i];
        }
        err=0;
        //#pragma omp parallel for
        for(int i=0; i<n; i++)
        {
            res[i]=b[i];
            for(int j=0;j<n;j++)
            {
                res[i]-=a[i][j]*x[j];
            }
            if(err<fabs(res[i]))
            {
                err=fabs(res[i]);
            }
        }
        mu0*=t;
        mu1*=t;
        mu2*=t;
        //fprintf(eps,"%d   %.15f\n",ite_num, err);
        cout<<ite_num<<"  "<<mu2<<"  "<<err<<endl;
    }
    delete [] res;
    delete [] xold;
}
#endif // LINEAREQUATION_H
