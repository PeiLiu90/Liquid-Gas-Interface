//***************************
//This File use Traplize Rule For Numerical Integral
//**************************



#ifndef FMT_H
#define FMT_H

#include "vector2d.h"
#include "vector3d.h"
#include "FFT.h"
#include "parameter.h"
#include "LinearEquation.h"
#include "Gauss.h"
#include "NewtonGmres.h"
//#include "mathtool.h"

#include <iostream>
using namespace std;

template<int N, int M>
class FMT
{
public:
    FMT(const double);
    FMT(const double, double * ,double * );

    void set(double * ,double * );

    void ComputeC();//compute dcf using analytical expression
    void ComputeCK();
    void ComputeCK2();
    void ComputeHK();
    void ComputeH();
    void ComputeB(const vector3d<double,N,N,M> &);
    void ComputeHB();
    void ComputeHB_Picard();
    double * _a;//diameter
    double * _rho;//homogeneous density
    double _n[4];//weighted density
    double * _r;//real space
    double * _w;//fourier space

    vector3d<double,N,N,M> _C;
    vector3d<double,N,N,M> _CK;
    vector3d<double,N,N,M> _HK;
    vector3d<double,N,N,M> _H;
    vector3d<double,N,N,M> _CH;
    vector3d<double,N,N,M> _B;
    vector3d<double,N,N,M> _ECP;
    void operator()(const vector2d<double,N,M> & , vector2d<double,N,M> & );
//private:
    void ComputeCH();
    void ComputeWD();
    double Density(const int , const int );
    double radius(const int );
    void ComputeDPhi();
    void ComputeLocalChemicalPotential();
    double DPhi(const int , const int );

    //used for computing TCF. Inhomogeneous system
    int FixHS;//fixed particle
    vector2d<double,N,M> _density;
    //vector2d<double,N,M> _density_old;
    vector3d<double,N,6,M> _wd;//weighted density for each species
    vector2d<double,6,M> _WD;
    vector3d<double,N,6,M> _lcp;
    vector2d<double,N,M> _LocalChemicalPotential;
    vector2d<double,6,M> _dPhi;
    //***************************************************
    //functions to compute DCF****************************
    double V(const int & , const int & ,const int &);
    double S(const int & , const int & ,const int &);
    double R(const int & , const int & ,const int &);
    double R2(const int & , const int & ,const int &);
    double T(const int & , const int & ,const int &);
    //****************************************************

    const double _L;
    const double _h;//spatial discretize
    const double _dk;
};



template<int N, int M>
FMT<N,M>::FMT(const double L, double * r, double * rho)
:_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _a=new double [N];
    _rho=new double [N];
    _n[0]=_n[1]=_n[2]=_n[3]=0.;
    for(int i=0;i<N;i++)
    {
        _a[i]=floor(r[i]/(2.*_h)+0.5)*_h*2.;
        //_a[i]=r[i];
        _rho[i]=rho[i];
        _n[0]+=_rho[i];
        _n[1]+=_rho[i]*_a[i]/2.;
        _n[2]+=_rho[i]*M_PI*_a[i]*_a[i];
        _n[3]+=_rho[i]*M_PI*_a[i]*_a[i]*_a[i]/6.;
    }
    cout<<"Packing Fraction: "<<_n[3]<<endl;
    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
}

template<int N, int M>
FMT<N,M>::FMT(const double L):_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _a=new double [N];
    _rho=new double [N];
    _n[0]=_n[1]=_n[2]=_n[3]=0.;

    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
}

template<int N, int M>
void FMT<N,M>::set(double * r, double * rho)
{
    for(int i=0;i<N;i++)
    {
        _a[i]=floor(r[i]/_h+0.5)*_h;
        //_a[i]=r[i];
        _rho[i]=rho[i];
        _n[0]+=_rho[i];
        _n[1]+=_rho[i]*_a[i]/2.;
        _n[2]+=_rho[i]*M_PI*_a[i]*_a[i];
        _n[3]+=_rho[i]*M_PI*_a[i]*_a[i]*_a[i]/6.;
    }
    cout<<"Packing Fraction: "<<_n[3]<<endl;
}

template<int N, int M>
void FMT<N,M>::ComputeCK()
{
    double w1[6];
    double w2[6];
    double xi[5];
    xi[0]=1./(1.-_n[3]);
    xi[1]=_n[2]/pow(1.-_n[3],2);
    xi[2]=_n[1]/pow(1.-_n[3],2)-_n[2]*_n[2]/(12.*M_PI*_n[3]*_n[3])*(2.*log(1.-_n[3])/_n[3]+(2.-5.*_n[3]+_n[3]*_n[3])/pow(1.-_n[3],3));
    xi[3]=_n[0]/pow(1.-_n[3],2)+2.*_n[1]*_n[2]/pow(1.-_n[3],3)+
                                    pow(_n[2],3)/(36.*M_PI)*(_n[3]*(-5.*pow(_n[3],3)+26.*pow(_n[3],2)-21.*_n[3]+6.)+
                                    6.*pow(1.-_n[3],4)*log(1.-_n[3]))/pow((1.-_n[3])*_n[3],4);
    xi[4]=(_n[3]+(1-_n[3])*(1-_n[3])*log(1-_n[3]))*_n[2]/(6.*M_PI*pow((1-_n[3])*_n[3],2));

    double temp;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<M;k++)
            {
                temp=_w[k]*_a[i]/2.;
                w1[0]=sin(temp)/temp;
                w1[1]=sin(temp)/_w[k];
                w1[2]=2.*M_PI*_a[i]*sin(temp)/_w[k];
                w1[3]=4.*M_PI*(sin(temp)-temp*cos(temp))/(_w[k]*_w[k]*_w[k]);
                w1[5]=-w1[3]*_w[k];
                w1[4]=-w1[3]*_w[k]/(2.*M_PI*_a[i]);

                temp=_w[k]*_a[j]/2.;
                w2[0]=sin(temp)/temp;
                w2[1]=sin(temp)/_w[k];
                w2[2]=2.*M_PI*_a[j]*sin(temp)/_w[k];
                w2[3]=4.*M_PI*(sin(temp)-temp*cos(temp))/(_w[k]*_w[k]*_w[k]);
                w2[5]=-w2[3]*_w[k];
                w2[4]=-w2[3]*_w[k]/(2.*M_PI*_a[j]);

                _CK(i,j,k)=-(w1[0]*w2[3]+w1[3]*w2[0]+w1[1]*w2[2]+w1[2]*w2[1]-w1[4]*w2[5]-w1[5]*w2[4])*xi[0];
                _CK(i,j,k)-=(w1[1]*w2[3]+w1[3]*w2[1])*xi[1];
                _CK(i,j,k)-=(w1[2]*w2[2]-w1[5]*w2[5])*xi[4];
                _CK(i,j,k)-=(w1[2]*w2[3]+w1[3]*w2[2])*xi[2];
                _CK(i,j,k)-=w1[3]*w2[3]*xi[3];
                //_CK(i,j,k)=-w1[0]*w2[3]-w1[3]*w1[0]-w1[2]*w2[1]-w1[1]*w2[2]+w1[4]*w2[5]+w1[5]*w2[4];//dilute limite
                }
        }
    }
}

//template<int N,int M>
//void FMT<N,M>::ComputeCK()
//{
//    double * fsin=new double [1+M];
//    double * fcos=new double [1+M];
//    for(int i=0;i<M+1;i++)
//    {
//        fsin[i]=sin(double(i)*M_PI/double(M+1));
//        fcos[i]=cos(double(i)*M_PI/double(M+1));
//    }
//
//    int size=(M+1)*2;
//    double *re=new double [size];
//    double *im=new double [size];
//
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            re[0]=0.;
//            im[0]=0.;
//            re[M+1]=0.;
//            im[M+1]=0.;
//            for(int k=1;k<M+1;k++)
//            {
//                re[k]=_C(i,j,k-1)*_r[k-1]*_h*4.*M_PI;
//                im[k]=0.;
//                re[size-k]=-re[k];
//                im[size-k]=0.;
//            }
//            FFT(re,im,fsin,fcos,size);
//            for(int k=0;k<M;k++)
//            {
//                _CK(i,j,k)=im[k+1]/(2.*_w[k]);
//            }
//        }
//    }
//    delete [] fsin;
//    delete [] fcos;
//    delete [] re;
//    delete [] im;
//}

template<int N, int M>
double FMT<N,M>::V(const int & i, const int & j, const int & k)
{
    if(_r[k]*2.<fabs(_a[i]-_a[j]))
    {
        return M_PI*pow(min(_a[i],_a[j]),3)/6.;
    }
    else if(_r[k]*2.<(_a[i]+_a[j]))
    {
        return -M_PI/(64.*_r[k])*(_a[i]*_a[i]-_a[j]*_a[j])*(_a[i]*_a[i]-_a[j]*_a[j])+
                (pow(_a[i],3)+pow(_a[j],3))*M_PI/12.-M_PI*_r[k]*(_a[i]*_a[i]+_a[j]*_a[j])/8.+M_PI*pow(_r[k],3)/12.;
    }
    else
    {
        return 0.;
    }
}

template<int N, int M>
double FMT<N,M>::S(const int & i, const int & j, const int & k)
{
    if(_r[k]*2.<fabs(_a[i]-_a[j]))
    {
        return M_PI*pow(min(_a[i],_a[j]),2);
    }
    else if(_r[k]*2.<(_a[i]+_a[j]))
    {
        return -M_PI/(8.*_r[k])*(_a[i]+_a[j])*(_a[i]-_a[j])*(_a[i]-_a[j])+
                (_a[i]*_a[i]+_a[j]*_a[j])*M_PI/2.-M_PI*_r[k]*(_a[i]+_a[j])/2.;
    }
    else
    {
        return 0.;
    }
}

template<int N, int M>
double FMT<N,M>::R(const int & i, const int & j, const int & k)
{
    if(_r[k]*2.<fabs(_a[i]-_a[j]))
    {
        return min(_a[i],_a[j])/2.;
    }
    else if(_r[k]*2.<(_a[i]+_a[j]))
    {
        return -(_a[i]-_a[j])*(_a[i]-_a[j])/(16.*_r[k])+(_a[i]+_a[j])/4.-_r[k]/4.;
    }
    else
    {
        return 0.;
    }
}

template<int N, int M>
double FMT<N,M>::R2(const int & i, const int & j, const int & k)
{
    if(_r[k]*2.<fabs(_a[i]-_a[j]))
    {
        return 0.;
    }
    else if(_r[k]*2.<(_a[i]+_a[j]))
    {
        return -M_PI*(_a[i]-_a[j])*(_a[i]-_a[j])/(4.*_r[k])+_r[k]*M_PI;
    }
    else
    {
        return 0.;
    }
}

template<int N, int M>
double FMT<N,M>::T(const int & i, const int & j, const int & k)
{
    if(_r[k]*2.<(_a[i]+_a[j]))
    {
        return 1.;
    }
    else
    {
        return 0.;
    }
}

template<int N, int M>
void FMT<N,M>::ComputeC()
{
    double xi[5];
    xi[0]=1./(1.-_n[3]);
    xi[1]=_n[2]/pow(1.-_n[3],2);
    xi[2]=_n[1]/pow(1.-_n[3],2)-_n[2]*_n[2]/(12.*M_PI*_n[3]*_n[3])*(2.*log(1.-_n[3])/_n[3]+(2.-5.*_n[3]+_n[3]*_n[3])/pow(1.-_n[3],3));
    xi[3]=_n[0]/pow(1.-_n[3],2)+2.*_n[1]*_n[2]/pow(1.-_n[3],3)+
                                    pow(_n[2],3)/(36.*M_PI)*(_n[3]*(-5.*pow(_n[3],3)+26.*pow(_n[3],2)-21.*_n[3]+6.)+
                                    6.*pow(1.-_n[3],4)*log(1.-_n[3]))/pow((1.-_n[3])*_n[3],4);
    xi[4]=(_n[3]+(1-_n[3])*(1-_n[3])*log(1-_n[3]))*_n[2]/(6.*M_PI*pow((1-_n[3])*_n[3],2));

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<M;k++)
            {
                _C(i,j,k)=-xi[3]*V(i,j,k)-xi[2]*S(i,j,k)-xi[1]*R(i,j,k)-(xi[4]-xi[1]/(4.*M_PI))*R2(i,j,k)-xi[0]*T(i,j,k);
            }
        }
    }
}

//template<int N, int M>
//void FMT<N,M>::ComputeC2()
//{
//    double * fsin=new double [1+M];
//    double * fcos=new double [1+M];
//    for(int i=0;i<M+1;i++)
//    {
//        fsin[i]=sin(double(i)*M_PI/double(M+1));
//        fcos[i]=cos(double(i)*M_PI/double(M+1));
//    }
//
//    int size=(M+1)*2;
//    double *re=new double [size];
//    double *im=new double [size];
//
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            re[0]=0.;
//            im[0]=0.;
//            re[M+1]=0.;
//            im[M+1]=0.;
//            for(int k=1;k<M+1;k++)
//            {
//                re[k]=_CK(i,j,k-1)*_w[k-1]*dk/(2.*M_PI*M_PI);
//                im[k]=0.;
//                re[size-k]=-re[k];
//                im[size-k]=0.;
//            }
//            FFT(re,im,fsin,fcos,size);
//            for(int k=0;k<M;k++)
//            {
//                _C(i,j,k)=im[k+1]/(2.*_r[k]);
//            }
//        }
//    }
//    delete [] fsin;
//    delete [] fcos;
//    delete [] re;
//    delete [] im;
//}

template<int N, int M>
void FMT<N,M>::ComputeHK()
{
    double **a=new double * [N*N];
    double *b=new double [N*N];
    for(int i=0;i<N*N;i++)
    {
        a[i]=new double [N*N];
    }
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N*N;j++)
        {
            for(int k=0;k<N*N;k++)
            {
                if(j==k)
                {
                    a[j][k]=1;
                }
                else
                {
                    a[j][k]=0;
                }
            }
        }
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<N;k++)
            {
                b[N*j+k]=_CK(j,k,i);
                for(int p=0;p<N;p++)
                {
                    a[N*j+k][N*p+k]-=_rho[p]*_CK(j,p,i);
                }
            }
        }
        Gauss(a,b,N*N);
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<N;k++)
            {
                _HK(j,k,i)=b[N*j+k];
            }
        }
    }
    for(int i=0;i<N*N;i++)
    {
        delete [] a[i];
    }
    delete [] b;
    delete [] a;
}

template<int N,int M>
void FMT<N,M>::ComputeH()
{
    double * fsin=new double [1+M];
    double * fcos=new double [1+M];
    for(int i=0;i<M+1;i++)
    {
        fsin[i]=sin(double(i)*M_PI/double(M+1));
        fcos[i]=cos(double(i)*M_PI/double(M+1));
    }

    int size=(M+1)*2;
    double *re=new double [size];
    double *im=new double [size];

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            re[0]=0.;
            im[0]=0.;
            re[M+1]=0.;
            im[M+1]=0.;
            for(int k=1;k<M+1;k++)
            {
                re[k]=_HK(i,j,k-1)*_w[k-1]*_dk/(2.*M_PI*M_PI);
                im[k]=0.;
                re[size-k]=-re[k];
                im[size-k]=0.;
            }
            FFT(re,im,fsin,fcos,size);
            for(int k=0;k<M;k++)
            {
                _H(i,j,k)=im[k+1]/(2.*_r[k]);
            }
        }
    }
    delete [] fsin;
    delete [] fcos;
    delete [] re;
    delete [] im;
}

//template<int N, int M>
//void FMT<N,M>::ComputeH2()
//{
//    double err;
//    for(int k=0;k<N;k++)//fixed particle
//    {
//        for(int i=0;i<N;i++)
//        {
//            for(int j=0;j<M;j++)
//            {
//                if(_r[j]<(_a[k]+_a[i])*0.5)
//                {
//                    _density(i,j)=0.;
//                }
//                else
//                {
//                    _density(i,j)=_rho[i];
//                }
//            }
//        }
//        for(int l=0;l<3000;l++)//iteration
//        {
//            _density_old=_density;
//            ComputeWD();
//            ComputeDPhi();
//            ComputeLocalChemicalPotential();
//            err=0.;
//            for(int i=0;i<N;i++)
//            {
//                for(int j=0;j<M;j++)
//                {
//                    if(_r[j]<(_a[k]+_a[i])*0.5)
//                    {
//                        _density(i,j)=0.;
//                    }
//                    else
//                    {
//                        _density(i,j)=pow(M_E,-_LocalChemicalPotential(i,j));
//                    }
//                }
//                for(int j=0;j<M;j++)
//                {
//                    _density(i,j)=_density(i,j)*_rho[i]/_density(i,M-1);
//                    err+=fabs(_density(i,j)-_density_old(i,j));
//                }
//            }
//            err/=double(N*M);
//            _density=_density*0.01+_density_old*0.99;
//        }
//        cout<<err<<endl;
//        for(int i=0;i<N;i++)
//        {
//            for(int j=0;j<M;j++)
//            {
//                _H(k,i,j)=_density(i,j)/_rho[i]-1.;
//            }
//        }
//    }
//}

template<int N, int M>
void FMT<N,M>::operator()(const vector2d<double,N,M> & x, vector2d<double,N,M> & ax)
{
    _density=x;
    ComputeWD();
    ComputeDPhi();
    ComputeLocalChemicalPotential();
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            if(2.*_r[j]<_a[FixHS]+_a[i])
            {
                _density(i,j)=0.;
            }
            else
            {
                _density(i,j)=pow(M_E,-_LocalChemicalPotential(i,j));
            }
        }
        for(int j=0;j<M;j++)
        {
            _density(i,j)=_density(i,j)*_rho[i]/_density(i,M-1);
        }
    }
    ax=_density-x;
}

template<int N, int M>
void FMT<N,M>::ComputeWD()
{
    _WD=0.;
    int xmin;
    int xmax;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            _wd(i,2,j)=0.;
            _wd(i,3,j)=0.;
            _wd(i,5,j)=0.;
            xmin=floor((_r[j]-_a[i]*0.5)/_h+0.5)-1;
            xmax=floor((_r[j]+_a[i]*0.5)/_h+0.5)-1;

            for(int k=xmin+1;k<xmax;k++)
            {
//                _wd(i,2,j)+=(radius(k)*Density(i,k)+radius(k+1)*Density(i,k+1))*h*0.5;
//                _wd(i,3,j)+=(radius(k)*Density(i,k)*(_a[i]*_a[i]/4.-(_r[j]-radius(k))*(_r[j]-radius(k)))
//                             +radius(k+1)*Density(i,k+1)*(_a[i]*_a[i]/4.-(_r[j]-radius(k+1))*(_r[j]-radius(k+1))))*h*0.5;
//                _wd(i,5,j)+=(radius(k)*Density(i,k)*(_r[j]*_r[j]+_a[i]*_a[i]/4.-radius(k)*radius(k))
//                            +radius(k+1)*Density(i,k+1)*(_r[j]*_r[j]+_a[i]*_a[i]/4.-radius(k+1)*radius(k+1)))*h*0.5;
                _wd(i,2,j)+=radius(k)*Density(i,k)*_h;
                _wd(i,3,j)+=radius(k)*Density(i,k)*(_a[i]*_a[i]/4.-(_r[j]-radius(k))*(_r[j]-radius(k)))*_h;
                _wd(i,5,j)+=(radius(k)*Density(i,k)*(_r[j]*_r[j]+_a[i]*_a[i]/4.-radius(k)*radius(k)))*_h;
            }
            _wd(i,2,j)+=radius(xmin)*Density(i,xmin)*_h*0.5;
            _wd(i,3,j)+=radius(xmin)*Density(i,xmin)*(_a[i]*_a[i]/4.-(_r[j]-radius(xmin))*(_r[j]-radius(xmin)))*_h*0.5;
            _wd(i,5,j)+=(radius(xmin)*Density(i,xmin)*(_r[j]*_r[j]+_a[i]*_a[i]/4.-radius(xmin)*radius(xmin)))*_h*0.5;
            _wd(i,2,j)+=radius(xmax)*Density(i,xmax)*_h*0.5;
            _wd(i,3,j)+=radius(xmax)*Density(i,xmax)*(_a[i]*_a[i]/4.-(_r[j]-radius(xmax))*(_r[j]-radius(xmax)))*_h*0.5;
            _wd(i,5,j)+=(radius(xmax)*Density(i,xmax)*(_r[j]*_r[j]+_a[i]*_a[i]/4.-radius(xmax)*radius(xmax)))*_h*0.5;
            _wd(i,2,j)*=(M_PI*_a[i])/_r[j];
            _wd(i,3,j)*=(M_PI)/_r[j];
            _wd(i,5,j)*=(M_PI)/(_r[j]*_r[j]);
            _wd(i,0,j)=_wd(i,2,j)/(M_PI*_a[i]*_a[i]);
            _wd(i,1,j)=_wd(i,0,j)*_a[i]/2.;
            _wd(i,4,j)=_wd(i,5,j)/(2.*M_PI*_a[i]);
            _WD(0,j)+=_wd(i,0,j);
            _WD(1,j)+=_wd(i,1,j);
            _WD(2,j)+=_wd(i,2,j);
            _WD(3,j)+=_wd(i,3,j);
            _WD(4,j)+=_wd(i,4,j);
            _WD(5,j)+=_wd(i,5,j);
        }
    }
}

template<int N, int M>
void FMT<N,M>::ComputeDPhi()
{
    for(int i=0;i<M;i++)
    {
        if(fabs(_WD(3,i))<pow(0.1,10))
        {
            _dPhi(0,i)=0;
            _dPhi(1,i)=0;
            _dPhi(2,i)=0;
            _dPhi(3,i)=0;
            _dPhi(4,i)=0;
            _dPhi(5,i)=0;
        }
        else
        {
            _dPhi(0,i)=-log(1.-_WD(3,i));
            _dPhi(1,i)=_WD(2,i)/(1.-_WD(3,i));
            _dPhi(2,i)=_WD(1,i)/(1.-_WD(3,i))+(_WD(3,i)+(1.-_WD(3,i))*(1.-_WD(3,i))*log(1.-_WD(3,i)))/
                        (12.*M_PI*_WD(3,i)*_WD(3,i)*(1.-_WD(3,i))*(1.-_WD(3,i)))*(_WD(2,i)*_WD(2,i)-_WD(5,i)*_WD(5,i));
            _dPhi(3,i)=_WD(0,i)/(1.-_WD(3,i))+(_WD(1,i)*_WD(2,i)-_WD(4,i)*_WD(5,i))/(1.-_WD(3,i))/(1.-_WD(3,i))-
                        (_WD(3,i)*(_WD(3,i)*_WD(3,i)-5.*_WD(3,i)+2.)+2.*(1.-_WD(3,i))*(1.-_WD(3,i))*(1.-_WD(3,i))*log(1.-_WD(3,i)))
                        /(36.*M_PI*(1.-_WD(3,i))*(1.-_WD(3,i))*(1.-_WD(3,i))*_WD(3,i)*_WD(3,i)*_WD(3,i))*
                        (_WD(2,i)*_WD(2,i)*_WD(2,i)-3.*_WD(2,i)*_WD(5,i)*_WD(5,i));
            _dPhi(4,i)=-_WD(5,i)/(1.-_WD(3,i));
            _dPhi(5,i)=-_WD(4,i)/(1.-_WD(3,i))-(_WD(3,i)+(1.-_WD(3,i))*(1.-_WD(3,i))*log(1.-_WD(3,i)))/
                        (6.*M_PI*_WD(3,i)*_WD(3,i)*(1.-_WD(3,i))*(1.-_WD(3,i)))*_WD(2,i)*_WD(5,i);
        }
    }
}

template<int N, int M>
void FMT<N,M>::ComputeLocalChemicalPotential()
{
    int xmin;
    int xmax;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            xmin=floor((_r[j]-_a[i]*0.5)/_h+0.5)-1;
            xmax=floor((_r[j]+_a[i]*0.5)/_h+0.5)-1;
            _lcp(i,0,j)=0.;
            _lcp(i,1,j)=0.;
            _lcp(i,2,j)=0.;
            _lcp(i,3,j)=0.;
            _lcp(i,4,j)=0.;
            _lcp(i,5,j)=0.;
            for(int k=xmin;k<xmax;k++)
            {
                _lcp(i,0,j)+=(radius(k)*DPhi(0,k)+radius(k+1)*DPhi(0,k+1))*_h*0.5;
                _lcp(i,1,j)+=(radius(k)*DPhi(1,k)+radius(k+1)*DPhi(1,k+1))*_h*0.5;
                _lcp(i,2,j)+=(radius(k)*DPhi(2,k)+radius(k+1)*DPhi(2,k+1))*_h*0.5;
                _lcp(i,3,j)+=(radius(k)*DPhi(3,k)*(_a[i]*_a[i]/4.-(_r[j]-radius(k))*(_r[j]-radius(k)))
                              +radius(k+1)*DPhi(3,k+1)*(_a[i]*_a[i]/4.-(_r[j]-radius(k+1))*(_r[j]-radius(k+1))))*_h*0.5;
                _lcp(i,4,j)+=(DPhi(4,k)*(-_r[j]*_r[j]+_a[i]*_a[i]/4.+radius(k)*radius(k))
                              +DPhi(4,k+1)*(-_r[j]*_r[j]+_a[i]*_a[i]/4.+radius(k+1)*radius(k+1)))*_h*0.5;
                _lcp(i,5,j)+=(DPhi(5,k)*(-_r[j]*_r[j]+_a[i]*_a[i]/4.+radius(k)*radius(k))
                              +DPhi(5,k+1)*(-_r[j]*_r[j]+_a[i]*_a[i]/4.+radius(k+1)*radius(k+1)))*_h*0.5;
            }
            _lcp(i,2,j)*=(M_PI*_a[i])/_r[j];
            _lcp(i,3,j)*=(M_PI)/_r[j];
            _lcp(i,5,j)*=(M_PI)/_r[j];
            _lcp(i,0,j)/=(_a[i]*_r[j]);
            _lcp(i,1,j)/=(2.*_r[j]);
            _lcp(i,4,j)/=(2.*_a[i]*_r[j]);
            _LocalChemicalPotential(i,j)=_lcp(i,0,j)+_lcp(i,1,j)+_lcp(i,2,j)+_lcp(i,3,j)+_lcp(i,4,j)+_lcp(i,5,j);
            //_LocalChemicalPotential(i,j)=_lcp(i,2,j)*M_PI*sig*2.*h/_r[j]+_lcp(i,1,j)*h/(2.*_r[j])+_lcp(i,0,j)*h/(sig*2.*_r[j])
              //                          +_lcp(i,3,j)*M_PI*h/_r[j]+_lcp(i,5,j)*(M_PI*h)/(_r[j])+_lcp(i,4,j)*h/(4.*_r[j]*sig);
        }
    }
}

template<int N, int M>
double FMT<N,M>::radius(const int i)
{
    return double(i+1)*_h;
}

template<int N, int M>
double FMT<N,M>::Density(const int i , const int j)
{
    if (j>=M)
    {
        return _density(i,M-1);
    }
    else if (j==-1)
    {
        return _density(i,0)*2.-_density(i,1);
    }
    else if (j<-1)
    {
        return _density(i,fabs(j+2));
    }
    else
    {
        return _density(i,j);
    }
}

template<int N, int M>
double FMT<N,M>::DPhi(const int i , const int j)
{
    if (j>=M)
    {
        return _dPhi(i,M-1);
    }
    else if (j==-1)
    {
        return _dPhi(i,0)*2.-_dPhi(i,1);
    }
    else if (j<-1)
    {
        return _dPhi(i,fabs(j+2));
    }
    else
    {
        return _dPhi(i,j);
    }
}

template<int N, int M>
void FMT<N,M>::ComputeHB_Picard()
{
    vector2d<double,N,M> _density_old;
    double err;
    for(FixHS=0;FixHS<N;FixHS++)
    {
        //initialize*****************************
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                if(2.*_r[j]<_a[FixHS]+_a[i])
                {
                    _density(i,j)=0.;
                }
                else
                {
                    _density(i,j)=_rho[i]*(1.+_H(FixHS,i,j));
                }
            }
        }
        err=1.;
        //******************************************
        while(err>pow(0.1,8))
        {
            _density_old=_density;
            ComputeWD();
            ComputeDPhi();
            ComputeLocalChemicalPotential();
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<M;j++)
                {
                    if(2.*_r[j]<_a[FixHS]+_a[i])
                    {
                        _density(i,j)=0.;
                    }
                    else
                    {
                        _density(i,j)=pow(M_E,-_LocalChemicalPotential(i,j));
                    }
                }
                for(int j=0;j<M;j++)
                {
                    _density(i,j)=_density(i,j)*_rho[i]/_density(i,M-1);
                }
            }
            _density=_density*0.01+_density_old*0.99;
            err=norm(_density-_density_old);
	    cout<<err<<endl;
        }
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _H(FixHS,i,j)=_density(i,j)/_rho[i]-1.;
            }
        }
    }
}

template<int N, int M>
void FMT<N,M>::ComputeCH()
{
    double * fsin=new double [1+M];
    double * fcos=new double [1+M];
    for(int i=0;i<M+1;i++)
    {
        fsin[i]=sin(double(i)*M_PI/double(M+1));
        fcos[i]=cos(double(i)*M_PI/double(M+1));
    }

    int size=(M+1)*2;
    double *re=new double [size];
    double *im=new double [size];

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            re[0]=0.;
            im[0]=0.;
            re[M+1]=0.;
            im[M+1]=0.;
            for(int k=1;k<M+1;k++)
            {
                re[k]=_H(i,j,k-1)*_r[k-1]*_h*4.*M_PI;
                im[k]=0.;
                re[size-k]=-re[k];
                im[size-k]=0.;
            }
            FFT(re,im,fsin,fcos,size);
            for(int k=0;k<M;k++)
            {
                _HK(i,j,k)=im[k+1]/(2.*_w[k]);
            }
        }
    }
//    double **a=new double * [N*N];
//    double *b=new double [N*N];
//    for(int i=0;i<N*N;i++)
//    {
//        a[i]=new double [N*N];
//    }
//    for(int i=0;i<M;i++)
//    {
//        for(int j=0;j<N*N;j++)
//        {
//            for(int k=0;k<N*N;k++)
//            {
//                if(j==k)
//                {
//                    a[j][k]=1;
//                }
//                else
//                {
//                    a[j][k]=0;
//                }
//            }
//        }
//        for(int j=0;j<N;j++)
//        {
//            for(int k=0;k<N;k++)
//            {
//                b[N*j+k]=_HK(j,k,i);
//                for(int p=0;p<N;p++)
//                {
//                    a[N*j+k][N*p+k]+=_rho[p]*_HK(p,k,i);
//                }
//            }
//        }
//        Gauss(a,b,N*N);
//        for(int j=0;j<N;j++)
//        {
//            for(int k=0;k<N;k++)
//            {
//                _CK(j,k,i)=b[N*j+k];
//            }
//        }
//    }
//    for(int i=0;i<N*N;i++)
//    {
//        delete [] a[i];
//    }
//    delete [] b;
//    delete [] a;

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<M;k++)
            {
                _CH(i,j,k)=0.;
                for(int l=0;l<N;l++)
                {
                    _CH(i,j,k)+=_CK(i,l,k)*_HK(l,j,k)*_rho[l];
                }
            }
        }
    }
    //for(int i=0;i<M+1;i++)
    //{
        //fsin[i]=-sin(double(i)*M_PI/double(M+1));
        //fcos[i]=cos(double(i)*M_PI/double(M+1));
    //}
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            re[0]=0.;
            im[0]=0.;
            re[M+1]=0.;
            im[M+1]=0.;
            for(int k=1;k<M+1;k++)
            {
                re[k]=_CH(i,j,k-1)*_w[k-1]*_dk/(2.*M_PI*M_PI);
                im[k]=0.;
                re[size-k]=-re[k];
                im[size-k]=0.;
            }
            FFT(re,im,fsin,fcos,size);
            for(int k=0;k<M;k++)
            {
                _CH(i,j,k)=im[k+1]/(2.*_r[k]);
            }
        }
    }

    delete [] fsin;
    delete [] fcos;
    delete [] re;
    delete [] im;
}

template<int N, int M>
void FMT<N,M>::ComputeHB()
{
//    ComputeC();
//    ComputeCK();
//    ComputeHK();
//    ComputeH();
    vector2d<double,N,M> x;
    vector2d<double,N,M> rhs;
    rhs=0.;
    const double tol=pow(0.1,9);
    NEWTONGMRES<vector2d<double,N,M>, FMT<N,M> > solve(*this,x,rhs,tol,M);
    for(FixHS=0;FixHS<N;FixHS++)
    {
        cout<<"Fixed Particle "<<FixHS<<endl;
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                if(_r[j]*2.<_a[FixHS]+_a[i])
                {
                    x(i,j)=0.;
                }
                else
                {
                    x(i,j)=_rho[i]*(1.+_H(FixHS,i,j));
                }
            }
        }
        solve.run();
        (*this)(x,rhs);
        cout<<"final error "<<norm(rhs)<<endl;
        _density=x;
        ComputeWD();
        ComputeDPhi();
        ComputeLocalChemicalPotential();
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _H(FixHS,i,j)=_density(i,j)/_rho[i]-1.;
                _B(FixHS,i,j)=_LocalChemicalPotential(i,j)-_LocalChemicalPotential(i,M-1);
            }
        }
    }
    ComputeCH();
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<M;k++)
            {
                _B(i,j,k)+=_CH(i,j,k);
                //_B(i,j,k)+=_H(i,j,k)-_C(i,j,k);
            }
        }
    }
}

template<int N, int M>
void FMT<N,M>::ComputeB(const vector3d<double,N,N,M> & H)
{
    _H=H;
    for(FixHS=0;FixHS<N;FixHS++)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _density(i,j)=_rho[i]*(1.+_H(FixHS,i,j));
                //_density(i,j)=_rho[i];
            }
        }
        ComputeWD();
        ComputeDPhi();
        ComputeLocalChemicalPotential();
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _ECP(FixHS,i,j)=_LocalChemicalPotential(i,j)-_LocalChemicalPotential(i,M-1);
                _B(FixHS,i,j)=_LocalChemicalPotential(i,j)-_LocalChemicalPotential(i,M-1);
            }
        }
    }
    ComputeCH();
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<M;k++)
            {
                _B(i,j,k)+=_CH(i,j,k);
            }
        }
    }
}
#endif
