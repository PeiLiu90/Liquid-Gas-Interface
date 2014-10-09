#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "math.h"

#ifndef M_PI
#define M_PI 3.1415927
#endif

#include <iostream>
using namespace std;

double legendre(const int n, const double x)
{
	if(n==0)
	{
		return 1.;
	}
	else if(n==1)
	{
		return x;
	}
	else
	{
	    double p1=x;
	    double p0=1.;
	    double p2;
	    for(int i=2;i<=n;i++)
	    {
	        p2=((2*i-1)*x*p1-(i-1)*p0)/i;
	        p0=p1;
	        p1=p2;
	    }
	    return p1;
		/*double y1=legendre(n-1,x);
		double y2=legendre(n-2,x);
		return (double(2*n-1)*x*y1-double(n-1)*y2)/double(n);*/
	}
}


double legendre_dri(const int n, const double x)
{
	return (x*legendre(n,x)-legendre(n-1,x))/(x*x-1)*double(n);
}


//x gauss location, w gauss weight

void gauleg(const double xmin, const double xmax, const int n, double * x, double * w)
{
	const double EPS=1.0e-10;
	double m=double(n+1)/2.; //The roots are symmetric in the interval
	double xmean=(xmin+xmax)/2.;
	double xlength=(xmax-xmin)/2.;
	double z;//roots in (-1,1)
	double z1=0;
	double f,df;
	for (int i=0;i<m;i++)
	{
		z=cos(M_PI*(i+0.75)/(n+0.5));//initial guess
		while (fabs(z-z1) >EPS)
		{
            f=legendre(n,z);
            df=legendre_dri(n,z);
			z1=z;
			z=z1-f/df;
		}
        f=legendre(n,z);
        df=legendre_dri(n,z);
		x[i]=xmean-xlength*z;
		x[n-1-i]=xmean+xlength*z;
		w[i]=2.0*xlength/((1.0-z*z)*df*df);
		w[n-1-i]=w[i];
	}
}

//x gauss location, w gauss weight, (-1,1)

void gauleg(const int n, double * x, double * w)
{
	const double EPS=1.0e-10;
	double m=double(n+1)/2.; //The roots are symmetric in the interval
	double z;//roots in (-1,1)
	double z1=0;
	double f,df;
	for (int i=0;i<m;i++)
	{
		z=cos(M_PI*(i+0.75)/(n+0.5));//initial guess
		while (fabs(z-z1) >EPS)
		{
            f=legendre(n,z);
            df=legendre_dri(n,z);
			z1=z;
			z=z1-f/df;
		}
        f=legendre(n,z);
        df=legendre_dri(n,z);
		x[i]=-z;
		x[n-1-i]=z;
		w[i]=2.0/((1.0-z*z)*df*df);
		w[n-1-i]=w[i];
	}
}

double GaussQuadrature(const int & N, const double & xmin, const double & xmax, double (*f)(const double &))
{
    double x[N];
    double w[N];
    gauleg(xmin,xmax, N, x, w);
    double intergal=0.;
    for(int i=0;i<N;i++)
    {
        intergal+=w[i]*f(x[i]);
    }
    return intergal;
}
#endif
