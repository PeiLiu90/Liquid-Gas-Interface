#ifndef LIQUID_VAPOR_H
#define LIQUID_VAPOR_H

#include "stdlib.h"
#include "vector2d.h"
#include "RHNC.h"
#include "load.h"

template<int M>
class Liquid_vapor
{
public:
    Liquid_vapor(const vector2d<double,1,1> &, const vector2d<double,1,1> &, const double ,double *);

    void UpdateRho();
    void SolveLMBW();
    void ComputeRHS();

    //value of Homogeneous DCF in Cylindrical Coordinate
    double clow(const int , const int );
    double clow2(const int , const int );
    double chigh(const int , const int );
    double chigh2(const int , const int );
    //Compute Integrated DCF of Interface System
    void ComputeDCF();
    //Density Distribution
    double _rho_b[4];//Four Bulk Homogeneous System

    double * _rho;//Interface System
    double * _rho_new;

    double * _drho_dz;
    double * _drho_dz_new;
    //DCF of Homogeneous System, Used to Approximate DCF of Inhomogeneous System
    //NOTE COORDINATEs of _Cb ARE NOT THE SAME WITH _C
    vector3d<double,1,1,M-1>  _Cb_low;
    vector3d<double,1,1,M-1>  _Cb_low2;
    vector3d<double,1,1,M-1>  _Cb_high;
    vector3d<double,1,1,M-1>  _Cb_high2;

    //Integrated correlation function of Homogeneous System
    double * _C_low2;
    double * _C_low;
    double * _C_high;
    double * _C_high2;
    //Integrated Correlation Function of Interface System
    vector2d<double,M,M> _C;//integrated parallel to the interface

    //coordinate
    const double _L;
    const double _h;//spatial discretize
    const double _dk;
    double * _z;
    double * _r;

    //temporary
    double tempr;
    double tempc;
    int index;
};

template<int M>
Liquid_vapor<M>::Liquid_vapor(const vector2d<double,1,1> & sig, const vector2d<double,1,1> & eps, const double L, double * rho)
:_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _z=new double [M];
    _r=new double [M];
    for(int i=0;i<M;i++)
    {
        _z[i]=double(i)*_h;
        _r[i]=_z[i];
    }

    RHNC<1,M-1> LJ(sig,eps,L*sqrt(2.));

    double * rhobulk=new double [1];
//Four Homogeneours System for One Component, Liquid-Vapor Interface
    rhobulk[0]=rho[0];
    LJ.set(rhobulk);
    LJ.ComputeCorrelationFunction();
    _Cb_high2=LJ._C;//Highest

    rhobulk[0]=rho[1];
    LJ.set(rhobulk);
    LJ.ComputeCorrelationFunction();
    _Cb_high=LJ._C;//high

    rhobulk[0]=rho[2];
    LJ.set(rhobulk);
    LJ.ComputeCorrelationFunction();
    _Cb_low=LJ._C;//low

    rhobulk[0]=rho[3];
    LJ.set(rhobulk);
    LJ.ComputeCorrelationFunction();
    _Cb_low2=LJ._C;//lowest
//******************************************************************
    _C_low2=new double [M];
    _C_low=new double [M];
    _C_high=new double [M];
    _C_high2=new double [M];
//Compute Integrated DCF of Homogeneous System**********************
    for(int i=0;i<M;i++)
    {
        _C_low[i]=0.;
        _C_low2[i]=0.;
        _C_high[i]=0.;
        _C_high2[i]=0.;
        for(int j=1;j<M;j++)
        {
            _C_low[i]+=clow(i,j)*_r[j]*_h;
            _C_low2[i]+=clow2(i,j)*_r[j]*_h;
            _C_high[i]+=chigh(i,j)*_r[j]*_h;
            _C_high2[i]+=chigh2(i,j)*_r[j]*_h;
        }
        _C_low[i]*=2.*M_PI;
        _C_low2[i]*=2.*M_PI;
        _C_high[i]*=2.*M_PI;
        _C_high2[i]*=2.*M_PI;
    }
    ofstream dcf("dcf-int.txt");
    for(int i=0;i<M;i++)
    {
        dcf<<_r[i]<<" "<<_C_low2[i]<<"  "<<_C_low[i]<<"  "<<_C_high[i]<<"   "<<_C_high2[i]<<endl;
    }
    dcf.close();
//******************************************************************
    _rho_b[0]=rho[0];
    _rho_b[1]=rho[1];
    _rho_b[2]=rho[2];
    _rho_b[3]=rho[3];
//Initial Guess
    zdf data;
    _rho=new double [M];
    _rho_new=new double [M];
    for(int i=0;i<M;i++)
    {
        //_rho[i]=0.9;
        //_rho[M-i]=0.1;
        _rho[i]=data(_z[i]-_L/2.);
    }

    _drho_dz=new double [M];
    _drho_dz_new=new double [M];
}

template<int M>
double Liquid_vapor<M>::chigh(const int i, const int j)
{
    tempr=sqrt(_z[i]*_z[i]+_r[j]*_r[j]);
    index=floor(tempr/(_h*sqrt(2.)))-1;
    if(index==-1)
    {
        tempc=(_Cb_high(0,0,1)*(tempr-sqrt(2)*_h)+_Cb_high(0,0,0)*(2*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    else if (index>=M-1)
    {
        tempc=0.;
    }
    else
    {
        tempc=(_Cb_high(0,0,index+1)*(tempr-double(index+1)*sqrt(2)*_h)+
               _Cb_high(0,0,index)*(double(index+2)*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    return tempc;
}

template<int M>
double Liquid_vapor<M>::chigh2(const int i, const int j)
{
    tempr=sqrt(_z[i]*_z[i]+_r[j]*_r[j]);
    index=floor(tempr/(_h*sqrt(2.)))-1;
    if(index==-1)
    {
        tempc=(_Cb_high2(0,0,1)*(tempr-sqrt(2)*_h)+_Cb_high2(0,0,0)*(2*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    else if (index>=M-1)
    {
        tempc=0.;
    }
    else
    {
        tempc=(_Cb_high2(0,0,index+1)*(tempr-double(index+1)*sqrt(2)*_h)+
               _Cb_high2(0,0,index)*(double(index+2)*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    return tempc;
}

template<int M>
double Liquid_vapor<M>::clow(const int i, const int j)
{
    tempr=sqrt(_z[i]*_z[i]+_r[j]*_r[j]);
    index=floor(tempr/(_h*sqrt(2.)))-1;
    if(index==-1)
    {
        tempc=(_Cb_low(0,0,1)*(tempr-sqrt(2)*_h)+_Cb_low(0,0,0)*(2*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    else if (index>=M-1)
    {
        tempc=0.;
    }
    else
    {
        tempc=(_Cb_low(0,0,index+1)*(tempr-double(index+1)*sqrt(2)*_h)+
               _Cb_low(0,0,index)*(double(index+2)*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    return tempc;
}

template<int M>
double Liquid_vapor<M>::clow2(const int i, const int j)
{
    tempr=sqrt(_z[i]*_z[i]+_r[j]*_r[j]);
    index=floor(tempr/(_h*sqrt(2.)))-1;
    if(index==-1)
    {
        tempc=(_Cb_low2(0,0,1)*(tempr-sqrt(2)*_h)+_Cb_low2(0,0,0)*(2*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    else if (index>=M-1)
    {
        tempc=0.;
    }
    else
    {
        tempc=(_Cb_low2(0,0,index+1)*(tempr-double(index+1)*sqrt(2)*_h)+
               _Cb_low2(0,0,index)*(double(index+2)*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    return tempc;
}

template<int M>
void Liquid_vapor<M>::ComputeDCF()
{
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<M;j++)
        {
            _C(i,j)=((_rho[i]-_rho_b[1])*(_rho[i]-_rho_b[2])*(_rho[i]-_rho_b[3])
                        /(_rho_b[0]-_rho_b[1])/(_rho_b[0]-_rho_b[2])/(_rho_b[0]-_rho_b[3])*_C_high2[abs(i-j)]
                   +(_rho[i]-_rho_b[0])*(_rho[i]-_rho_b[2])*(_rho[i]-_rho_b[3])
                        /(_rho_b[1]-_rho_b[0])/(_rho_b[1]-_rho_b[2])/(_rho_b[1]-_rho_b[3])*_C_high[abs(i-j)]
                   +(_rho[i]-_rho_b[0])*(_rho[i]-_rho_b[1])*(_rho[i]-_rho_b[3])
                        /(_rho_b[2]-_rho_b[0])/(_rho_b[2]-_rho_b[1])/(_rho_b[2]-_rho_b[3])*_C_low[abs(i-j)]
                   +(_rho[i]-_rho_b[0])*(_rho[i]-_rho_b[1])*(_rho[i]-_rho_b[2])
                        /(_rho_b[3]-_rho_b[0])/(_rho_b[3]-_rho_b[1])/(_rho_b[3]-_rho_b[2])*_C_low2[abs(i-j)]
                   +(_rho[j]-_rho_b[1])*(_rho[j]-_rho_b[2])*(_rho[j]-_rho_b[3])
                        /(_rho_b[0]-_rho_b[1])/(_rho_b[0]-_rho_b[2])/(_rho_b[0]-_rho_b[3])*_C_high2[abs(i-j)]
                   +(_rho[j]-_rho_b[0])*(_rho[j]-_rho_b[2])*(_rho[j]-_rho_b[3])
                        /(_rho_b[1]-_rho_b[0])/(_rho_b[1]-_rho_b[2])/(_rho_b[1]-_rho_b[3])*_C_high[abs(i-j)]
                   +(_rho[j]-_rho_b[0])*(_rho[j]-_rho_b[1])*(_rho[j]-_rho_b[3])
                        /(_rho_b[2]-_rho_b[0])/(_rho_b[2]-_rho_b[1])/(_rho_b[2]-_rho_b[3])*_C_low[abs(i-j)]
                   +(_rho[j]-_rho_b[0])*(_rho[j]-_rho_b[1])*(_rho[j]-_rho_b[2])
                        /(_rho_b[3]-_rho_b[0])/(_rho_b[3]-_rho_b[1])/(_rho_b[3]-_rho_b[2])*_C_low2[abs(i-j)])*0.5;
        }
    }
}

template<int M>
void Liquid_vapor<M>::ComputeRHS()
{
    _drho_dz[0]=0.;
    for(int i=1;i<M-1;i++)
    {
        _drho_dz[i]=(_rho[i+1]-_rho[i-1])/(2.*_h);
    }
    _drho_dz[M-1]=0.;
    for(int i=0;i<M;i++)
    {
        _drho_dz_new[i]=0.;
        for(int j=0;j<M;j++)
        {
            _drho_dz_new[i]+=_C(i,j)*_drho_dz[j]*_h;
        }
    }
}

template<int M>
void Liquid_vapor<M>::UpdateRho()
{
    double average;
    double err;
    _rho_new[0]=log(_rho[0]);
    for(int i=1;i<M;i++)
    {
        _rho_new[i]=_rho_new[i-1]+_drho_dz_new[i-1]*_h;
    }
    average=0.;
    err=0.;
    for(int i=0;i<M;i++)
    {
        _rho_new[i]=pow(M_E,_rho_new[i]);
        err+=fabs(_rho[i]-_rho_new[i]);
        _rho[i]=_rho[i]*0.9+_rho_new[i]*0.1;
        if(_rho[i]>0.9999)
        {
            _rho[i]=0.9999;
        }
        else if(_rho[i]<0.00001)
        {
            _rho[i]=0.00001;
        }
        average+=_rho[i];
    }
    average/=double(M);
    err/=double(M);
    for(int i=0;i<M;i++)
    {
        _rho[i]*=(0.35/average);
    }
    cout<<err<<endl;
}


template<int M>
void Liquid_vapor<M>::SolveLMBW()
{
    for(int i=0;i<200;i++)
    {
        ComputeDCF();
        ComputeRHS();
        UpdateRho();
    }
}
#endif
