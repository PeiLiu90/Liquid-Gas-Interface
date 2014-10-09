#ifndef FFT_H
#define FFT_H

#include "math.h"

void FFT(double * re , double * im ,double * fsin, double * fcos, const int & size)
{
    int i,j,k,p,q;
    int N=size;
    int m=1;
    double temp[2];
// CooleyTukey FFT algorithm
    while(N!=1)
    {
        k=N>>1;
        for(i=0;i<size;i+=N)
        {
            for(j=0,p=i,q=i+k;j<size>>1;j+=m,p++,q++)
            {
                temp[0]=re[p]-re[q];
                temp[1]=im[p]-im[q];
                re[p]+=re[q];
                im[p]+=im[q];
                re[q]=temp[0]*fcos[j]-temp[1]*fsin[j];
                im[q]=temp[0]*fsin[j]+temp[1]*fcos[j];
            }
        }
        m<<=1;
        N=k;
    }
//bit reversal, reordering
    j=0;
    m= size>>1;
    for (i = 0;i <size-1; i++)
    {
        if (i < j)
        {
            temp[0]=re[i];
            temp[1]=im[i];
            im[i] = im[j];
            re[i] = re[j];
            re[j] = temp[0];
            im[j] = temp[1];
        }
        k = m;
        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
}

void FST(double * inout , const int & _size)
{
    double * fsin=new double [_size];
    double * fcos=new double [_size];
    for(int i=0;i<_size;i++)
    {
        fsin[i]=sin(double(i)*M_PI/_size);
        fcos[i]=cos(double(i)*M_PI/_size);
    }
    int size=_size<<1;
    double *re=new double [size];
    double *im=new double [size];
    re[0]=0.;
    im[0]=0.;
    re[_size]=0.;
    im[_size]=0.;
    for(int i=1;i<_size;i++)
    {
        re[i]=inout[i-1];
        im[i]=0.;
        re[size-i]=-re[i];
        im[size-i]=0.;
    }
    FFT(re,im,fsin,fcos,size);
    for(int i=0;i<_size-1;i++)
    {
        inout[i]=im[i+1]/2;
    }
    delete [] fsin;
    delete [] fcos;
    delete [] re;
    delete [] im;
}

//Use FFT to Compute Function Derivatives
void FFT_dir(double * inout, const int & _size, const double & L)
{
    double * fsin=new double [_size];
    double * fcos=new double [_size];
    for(int i=0;i<_size;i++)
    {
        fsin[i]=sin(double(i)*2.*M_PI/_size);
        fcos[i]=cos(double(i)*2.*M_PI/_size);
    }
    double *re=new double [_size];
    double *im=new double [_size];
    for(int i=0;i<_size;i++)
    {
        re[i]=inout[i];
        im[i]=0.;
    }
    FFT(re,im,fsin,fcos,_size);
    for(int i=0;i<_size;i++)
    {
        fsin[i]=-fsin[i];
    }
    for(int i=1;i<_size/2;i++)
    {
        re[i]*=-double(i);
        re[_size-i]*=double(i);
        im[i]*=double(i);
        im[_size-i]*=-double(i);
    }
    re[0]=0.;
    re[_size/2]=0.;
    im[0]=0.;
    im[_size/2]=0.;

    FFT(im,re,fsin,fcos,_size);
    for(int i=0;i<_size;i++)
    {
        inout[i]=im[i]*2.*M_PI/(double(_size)*L);
    }
    delete [] fsin;
    delete [] fcos;
    delete [] re;
    delete [] im;
}
#endif // FFT_H
