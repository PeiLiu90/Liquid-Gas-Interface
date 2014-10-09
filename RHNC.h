#ifndef RHNC_H
#define RHNC_H

#include "FMT.h"

template<int N, int M>
class RHNC
{
public:
    RHNC(const vector2d<double,N,N> & , const vector2d<double,N,N> &, const double );
    RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const double L, const double *rho);
    RHNC(const vector2d<double,N,N> & , const vector2d<double,N,N> &, const vector2d<double , N, N> &, const double );
    RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const vector2d<double , N, N> & , const double L, const double *rho);

    void set(double* );

    void ComputeCK();
    void ComputeHK();
    void ComputeH();

    void ComputeCorrelationFunction();

    void operator()(vector3d<double,N,N,M> & ,  vector3d<double,N,N,M> & );
    void ComputeLado();

    FMT<N,M> HardSphere;//we need the bridge function
    double Lado;

    vector3d<double, N, N, M> _C;
    vector3d<double, N, N, M> _CK;
    vector3d<double, N, N, M> _HK;
    vector3d<double, N, N, M> _H;
//coordinate
    const double _L;
    const double _h;
    const double _dk;

    double * _r;
    double * _w;
//
    double Energy(const int  ,const int , const double);

    double * _rho;
    vector2d<double, N, N> _epsilon;//potential depth of pairwise interaction
    vector2d<double, N, N> _sigma;//effective diameter , they could be not additive?
    vector2d<double, N, N> _alpha;

    double u_int;//internal energy
    double xi;//compressibility
    bool R;//switch RHNC and HNC
};

template<int N, int M>
RHNC<N,M>::RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const double L, const double *rho)
:HardSphere(L),_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _epsilon=eps;
    _sigma=sig;
    _alpha=1.;
    _rho=new double [N];
    for(int i=0;i<N;i++)
    {
        _rho[i]=rho[i];
    }
    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
    double *BH=new double [N];
    for(int i=0;i<N;i++)
    {
        //BH
//        BH[i]=pow(2./(1.+pow(1.+(1./_epsilon(i,i)-0.05536/_epsilon(i,i)/_epsilon(i,i)
//                         +0.0007278/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))/1.1287,0.5)),1./6.)*_sigma(i,i);
        //BH[i]=_sigma(i,i);
        //WCA
        BH[i]=pow(2./(1.+0.8165*
                      pow((1./_epsilon(i,i)-0.03367/_epsilon(i,i)/_epsilon(i,i)
                           +0.0003935/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))
                            /(1.-0.09835*_rho[i]+0.04937*_rho[i]*_rho[i]
                            -0.1415*_rho[i]*_rho[i]*_rho[i]),0.5)),1./6.)*_sigma(i,i);
    }
    HardSphere.set(BH,_rho);
    HardSphere.ComputeC();
    HardSphere.ComputeCK();
    delete [] BH;
//    cout<<"***************Begin Computing Reference Hard Sphere System****************"<<endl;
//    cout<<"Diameters are ";
//    for(int i=0;i<N;i++)
//    {
//        cout<<HardSphere._a[i]<<"  ";
//    }
//    cout<<endl;
//    HardSphere.ComputeHB();
//    ofstream fout("HardSphereBridgeFunction.txt");
//    for(int k=0;k<M;k++)
//    {
//        fout<<_r[k];
//        for(int i=0;i<N;i++)
//        {
//            for(int j=0;j<N;j++)
//            {
//                fout<<"  "<<HardSphere._B(i,j,k);
//            }
//        }
//        fout<<endl;
//    }
//    fout.close();
//    cout<<"***************End Computing Reference Hard Sphere System******************"<<endl;
}

template<int N, int M>
RHNC<N,M>::RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const vector2d<double,N,N> & alp, const double L, const double *rho)
:HardSphere(L),_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _epsilon=eps;
    _sigma=sig;
    _alpha=alp;
    _rho=new double [N];
    for(int i=0;i<N;i++)
    {
        _rho[i]=rho[i];
    }
    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
    double *BH=new double [N];
    for(int i=0;i<N;i++)
    {
        //BH
//        BH[i]=pow(2./(1.+pow(1.+(1./_epsilon(i,i)-0.05536/_epsilon(i,i)/_epsilon(i,i)
//                         +0.0007278/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))/1.1287,0.5)),1./6.)*_sigma(i,i);
        //BH[i]=_sigma(i,i);
        //WCA
        BH[i]=pow(2./(1.+0.8165*
                      pow((1./_epsilon(i,i)-0.03367/_epsilon(i,i)/_epsilon(i,i)
                           +0.0003935/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))
                            /(1.-0.09835*_rho[i]+0.04937*_rho[i]*_rho[i]
                            -0.1415*_rho[i]*_rho[i]*_rho[i]),0.5)),1./6.)*_sigma(i,i);
    }
    HardSphere.set(BH,_rho);
    HardSphere.ComputeC();
    HardSphere.ComputeCK();
    delete [] BH;
//    cout<<"***************Begin Computing Reference Hard Sphere System****************"<<endl;
//    cout<<"Diameters are ";
//    for(int i=0;i<N;i++)
//    {
//        cout<<HardSphere._a[i]<<"  ";
//    }
//    cout<<endl;
//    HardSphere.ComputeHB();
//    ofstream fout("HardSphereBridgeFunction.txt");
//    for(int k=0;k<M;k++)
//    {
//        fout<<_r[k];
//        for(int i=0;i<N;i++)
//        {
//            for(int j=0;j<N;j++)
//            {
//                fout<<"  "<<HardSphere._B(i,j,k);
//            }
//        }
//        fout<<endl;
//    }
//    fout.close();
//    cout<<"***************End Computing Reference Hard Sphere System******************"<<endl;
}

template<int N, int M>
RHNC<N,M>::RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const double L)
:HardSphere(L),_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _epsilon=eps;
    _sigma=sig;
    _alpha=1.;
    _rho=new double [N];
    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
}

template<int N, int M>
RHNC<N,M>::RHNC(const vector2d<double,N,N> & sig, const vector2d<double,N,N> & eps, const vector2d<double,N,N> & alp, const double L)
:HardSphere(L),_L(L),_h(_L/double(M+1)),_dk(M_PI/_L)
{
    _epsilon=eps;
    _sigma=sig;
    _alpha=alp;
    _rho=new double [N];
    _r=new double [M];
    _w=new double [M];
    for(int i=0;i<M;i++)
    {
        _r[i]=double(i+1)*_h;
        _w[i]=double(i+1)*_dk;
    }
}

template<int N, int M>
void RHNC<N,M>::set(double * rho)
{
    for(int i=0;i<N;i++)
    {
        _rho[i]=rho[i];
    }
    double *BH=new double [N];
    for(int i=0;i<N;i++)
    {
        //BH
//        BH[i]=pow(2./(1.+pow(1.+(1./_epsilon(i,i)-0.05536/_epsilon(i,i)/_epsilon(i,i)
//                         +0.0007278/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))/1.1287,0.5)),1./6.)*_sigma(i,i);
        //BH[i]=_sigma(i,i);
        //WCA
        BH[i]=pow(2./(1.+0.8165*
                      pow((1./_epsilon(i,i)-0.03367/_epsilon(i,i)/_epsilon(i,i)
                           +0.0003935/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i)/_epsilon(i,i))
                            /(1.-0.09835*_rho[i]+0.04937*_rho[i]*_rho[i]
                            -0.1415*_rho[i]*_rho[i]*_rho[i]),0.5)),1./6.)*_sigma(i,i);
    }
    HardSphere.set(BH,_rho);
    HardSphere.ComputeC();
    HardSphere.ComputeCK();
}

template<int N,int M>
void RHNC<N,M>::ComputeCK()
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
                re[k]=_C(i,j,k-1)*_r[k-1]*_h*4.*M_PI;
                im[k]=0.;
                re[size-k]=-re[k];
                im[size-k]=0.;
            }
            FFT(re,im,fsin,fcos,size);
            for(int k=0;k<M;k++)
            {
                _CK(i,j,k)=im[k+1]/(2.*_w[k]);
            }
        }
    }
    delete [] fsin;
    delete [] fcos;
    delete [] re;
    delete [] im;
}

template<int N, int M>
void RHNC<N,M>::ComputeHK()
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
                    for(int q=0;q<N;q++)
                    {
                        if(k==q)
                        {
                            a[N*j+k][N*p+q]-=_rho[p]*_CK(j,p,i);
                        }
                    }
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
void RHNC<N,M>::ComputeH()
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

template<int N, int M>
double RHNC<N,M>::Energy(const int i, const int j, const double r)
{
    return 4.*_epsilon(i,j)*( pow(_sigma(i,j)/r,12.)- _alpha(i,j)*pow(_sigma(i,j)/r,6));
}

template<int N, int M>
void RHNC<N,M>::operator()(vector3d<double,N,N,M> & x, vector3d<double,N,N,M> & ax)
{
    _C=x;
    ComputeCK();
    ComputeHK();
    ComputeH();
    if(R)
    {
        HardSphere.ComputeB(_H);
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                for(int k=0;k<M;k++)
                {
                    ax(i,j,k)=_H(i,j,k)+1.-pow(M_E,-Energy(i,j,_r[k])+_H(i,j,k)-_C(i,j,k)-HardSphere._B(i,j,k));
                }
            }
        }
    }
    else
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                for(int k=0;k<M;k++)
                {
                    ax(i,j,k)=_H(i,j,k)+1.-pow(M_E,-Energy(i,j,_r[k])+_H(i,j,k)-_C(i,j,k));
                }
            }
        }
    }
}

template<int N, int M>
void RHNC<N,M>::ComputeCorrelationFunction()
{
    vector3d<double,N,N,M> rhs;
    rhs=0.;
    vector3d<double,N,N,M> C;
    C=HardSphere._C;
    NEWTONGMRES<vector3d<double,N,N,M>,RHNC<N,M> > solve(*this,C,rhs,tol,M);
    cout<<"*************SOLVING HNC+OZ*****************"<<endl;
    R=false;
    solve.run();
    cout<<"**************HNC FINISHED******************"<<endl;
    cout<<"**************SOLVING RHNC+OZ***************"<<endl;
    R=true;
    solve.run();
    cout<<"***************END*************************"<<endl;
//    xi=1.;
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            for(int k=0;k<M;k++)
//            {
//                xi+=_rho[j]*_H(i,j,k)*4.*M_PI*_h*_r[k]*_r[k];
//            }
//        }
//    }
//    cout<<"compressibility is "<<xi<<endl;
//    u_int=0.;
//    for(int i=0;i<N;i++)
//    {
//        for(int j=0;j<N;j++)
//        {
//            for(int k=0;k<M;k++)
//            {
//                u_int+=_rho[j]*(_H(i,j,k)+1.)*Energy(i,j,_r[k])*4.*M_PI*_h*_r[k]*_r[k];
//            }
//        }
//    }
//    cout<<"Internal Energy is "<<u_int<<endl;
}

template<int N, int M>
void RHNC<N,M>::ComputeLado()
{

}
#endif
