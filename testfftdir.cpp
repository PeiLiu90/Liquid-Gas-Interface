#include "FFT.h"
#include <fstream>
#include <iostream>
using namespace std;

int main()
{
    cout<<"Hello World"<<endl;
    int N=1024;
    double L=20.;
    double h=L/double(N);
    double * x =new double [N];
    double * y =new double [N];
    double * z =new double [N];
    for(int i=0;i<N;i++)
    {
        x[i]=double(i)*h-L/2.;
    }
    y[0]=0.;
    y[N/2]=0.;
    for(int i=1;i<N/2;i++)
    {
        y[i]=0.;
        y[N-i]=1.;
    }
    FFT_dir(y,N,L);
    ofstream dout("dir.txt");
    for(int i=0;i<N;i++)
    {
        dout<<x[i]<<"  "<<y[i]<<endl;
    }
    return 0;
}
