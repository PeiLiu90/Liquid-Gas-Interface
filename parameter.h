#ifndef PARAMETER_H
#define PARAMETER_H

#include "math.h"
#include <fstream>
#include <sstream>
const int N1=2;
const int M1=1024;//USE FST NEED M TO BE POWER OF 2 （MINUS 1）

double rho[2]={0.9,0.0};//bulk density, liquid phase & vapor phase

const double l=10.24;//computational domain
//const double h=L/double(M1+1);//spatial discretize
//const double dk=M_PI/L;//fourier space

const int N_gauss=100;//FMT2 do not need gauss integral

double tol=pow(0.1,9);

//READ Lenard Jones Parameter sigma and epsilon FROM FILE************
template<int N>
void Load(vector2d<double,N,N> & sig, vector2d<double,N,N> & eps)
{
    fstream fin("input.dat");
    char * buffer=new char [20];
    fin.getline(buffer,20);
    while(buffer[0]!='#')
    {
        fin.getline(buffer,20);
    }
    for(int i=0;i<N;i++)
    {
        fin.getline(buffer,20);
        stringstream oss( buffer);
        for(int j=0;j<N;j++)
        {
            oss>>sig(i,j);
        }
    }
    fin.getline(buffer,20);
    while(buffer[0]!='#')
    {
        fin.getline(buffer,20);
    }
    for(int i=0;i<N;i++)
    {
        fin.getline(buffer,20);
        stringstream oss( buffer);
        for(int j=0;j<N;j++)
        {
            oss>>eps(i,j);
        }
    }
    fin.close();
}

//READ Lenard Jones Parameter sigma and epsilon FROM FILE************
template<int N>
void Load(vector2d<double,N,N> & sig, vector2d<double,N,N> & eps, vector2d<double,N,N> & alp)
{
    fstream fin("input.dat");
    char * buffer=new char [20];
    fin.getline(buffer,20);
    while(buffer[0]!='#')
    {
        fin.getline(buffer,20);
    }
    for(int i=0;i<N;i++)
    {
        fin.getline(buffer,20);
        stringstream oss( buffer);
        for(int j=0;j<N;j++)
        {
            oss>>sig(i,j);
        }
    }
    fin.getline(buffer,20);
    while(buffer[0]!='#')
    {
        fin.getline(buffer,20);
    }
    for(int i=0;i<N;i++)
    {
        fin.getline(buffer,20);
        stringstream oss( buffer);
        for(int j=0;j<N;j++)
        {
            oss>>eps(i,j);
        }
    }
    fin.getline(buffer,20);
    while(buffer[0]!='#')
    {
        fin.getline(buffer,20);
    }
    for(int i=0;i<N;i++)
    {
        fin.getline(buffer,20);
        stringstream oss( buffer);
        for(int j=0;j<N;j++)
        {
            oss>>alp(i,j);
        }
    }
    fin.close();
}
//**************************************************************
#endif // PARAMETER_H
