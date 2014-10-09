#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

class zdf
{
public:
    zdf();
    double operator()(double );

    double *x;
    double *y;
    int index;
};

zdf::zdf()
{
    x=new double [45];
    y=new double [45];
    fstream fin("data.txt");
    char * buffer=new char [50];
    for(int i=0;i<45;i++)
    {
        fin.getline(buffer,50);
        stringstream oss( buffer);
        oss>>x[i];
        oss>>y[i];
    }
    fin.close();
}

double zdf::operator()( double xx)
{
    index=0;
    while(x[index]<xx)
    {
        index++;
    }
    if(index<=0)
    {
        return (y[0]*(x[1]-xx)+y[1]*(xx-x[0]))/(x[1]-x[0]);
    }
    else if (index>=43)
    {
        return (y[42]*(x[43]-xx)+y[43]*(xx-x[42]))/(x[43]-x[42]);
    }
    else
    {
        return (y[index-1]*(x[index]-xx)+y[index]*(xx-x[index-1]))/(x[index]-x[index-1]);
    }
}
