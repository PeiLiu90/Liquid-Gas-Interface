#ifndef LIQUID_VAPOR_H
#define LIQUID_VAPOR_H

#include "stdlib.h"
#include "vector2d.h"
#include "RHNC.h"
#include "load.h"

template<int N, int M>
class Liquid_liquid
{
public:
    Liquid_liquid(const vector2d<double,N,N> &, const vector2d<double,N,N> &, const vector2d<double,N,N> & ,const double ,double *);

    void UpdateRho();
    void SolveLMBW();
    void ComputeRHS();

    //value of Homogeneous DCF in Cylindrical Coordinate
    double chomo(const int , const int , const int , const int , const int );//firt one int: system; middle two int : component; last two int : coordinate
    //Compute Integrated DCF of Interface System
    void ComputeDCF();
    //Density Distribution
    double _rho_bulk[N+1];//Bulk Homogeneous System
    vector2d<double,N,M> _rho;

    //DCF of homogeneous system
    vector3d<double,N,N,M-1> _C_homo[N+1];
    //Integrated DCF of Homogeneous System
    vector3d<double,N,N,M> _C_int[N+1];
    //Integrated DCF of Interface System
    vector2d<vector2d<double,M,M>,N,N> _C;
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

template<int N, int M>
Liquid_liquid<N,M>::Liquid_liquid(const vector2d<double,N,N> & sig, const vector2d<double, N, N> & eps, const vector2d<double, N, N> & alp, const double L, double * rho)
:_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _z=new double [M];
    _r=new double [M];
    for(int i=0;i<M;i++)
    {
        _z[i]=double(i)*_h;
        _r[i]=_z[i];
    }

    RHNC<N,M-1> LJ(sig,eps,alp,L*sqrt(2.));

    for(int i=0;i<N+1;i++)
    {
        _rho_bulk[i]=rho[i];
    }

    double * rhobulk=new double [N];
    for(int i=0;i<N;i++)
    {
        rhobulk[i]=rho[0];
    }

    LJ.set(rhobulk);
    LJ.ComputeCorrelationFunction();
    _C_homo[0]=LJ._C;
    for(int i=0;i<N;i++)
    {
        rhobulk[i]=rho[i+1];
        LJ.set(rhobulk);
        LJ.ComputeCorrelationFunction();
        _C_homo[i+1]=LJ._C;
        rhobulk[i]=rho[0];
    }
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<N;k++)
            {
                for(int p=0;p<M;p++)
                {
                    _C_int[i](j,k,p)=0.;
                    for(int q=0;q<M;q++)
                    {
                        _C_int[i](j,k,p)+=chomo(i,j,k,p,q)*_r[q]*_h;
                    }
                }
            }
        }
    }
    ofstream dcf("dcf_homo.txt");
    for(int i=0;i<M-1;i++)
    {
        dcf<<LJ._r[i]<<"  ";
        for(int j=0;j<N+1;j++)
        {
            for(int k=0;k<N;k++)
            {
                for(int l=0;l<N;l++)
                {
                    dcf<<_C_int[j](k,l,i)<<"  ";
                }
            }
        }
        dcf<<endl;
    }
    dcf.close();
}

template<int N, int M>
double Liquid_liquid<N,M>::chomo(const int i, const int j ,const int k, const int p , const int q)
{
    tempr=sqrt(_z[p]*_z[p]+_r[q]*_r[q]);
    index=floor(tempr/(_h*sqrt(2.)))-1;
    if(index==-1)
    {
        tempc=(_C_homo[i](j,k,1)*(tempr-sqrt(2)*_h)+_C_homo[i](j,k,0)*(2*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    else if (index>=M-1)
    {
        tempc=0.;
    }
    else
    {
        tempc=(_C_homo[i](j,k,index+1)*(tempr-double(index+1)*sqrt(2)*_h)+
               _C_homo[i](j,k,index)*(double(index+2)*sqrt(2)*_h-tempr))/(sqrt(2)*_h);
    }
    return tempc;
}

template<int N, int M>
void Liquid_liquid<N,M>::ComputeDCF()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int p=0;p<M;p++)
            {
                for(int q=0;q<M;q++)
                {
                    _C(i,j)(p,q)=0.;
                    for(int k=0;k<N+1;k++)
                    {
                        _C(i,j)(p,q)+=_C_int[k](i,j,abs(p-q));
                    }
                }
            }
        }
    }
}
#endif
