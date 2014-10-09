#ifndef vector3d_H
#define vector3d_H

#include "math.h"
#include <iostream>
using namespace std;

template <typename E, typename T, int Nx,int Ny,int Nz>
class vector3dExpression
{
public:
    const T operator()(const int & i, const int & j, const int k) const{return static_cast<const E&>(*this)(i,j,k);}
    operator E&() { return static_cast<E&>(*this);}
    operator E const&() const { return static_cast<const E&>(*this);}
};

template<typename T ,int Nx,int Ny,int Nz>
class vector3d: public vector3dExpression<vector3d<T,Nx,Ny,Nz>, T, Nx,Ny,Nz>
{
public:
    vector3d();
    vector3d(const vector3d &);
    ~vector3d();
    template<typename E> vector3d(const vector3dExpression<E,T,Nx,Ny,Nz> & m)
    {
        _data=new T ** [Nx];
        for(int i=0;i<Nx;i++)
        {
            _data[i]=new T * [Ny];
            for(int j=0;j<Ny;j++)
            {
                _data[i][j]=new T [Nz];
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]=m(i,j,k);
                }
            }
        }
    }

    T & operator()(const int & ,const int & , const int & );
    const T & operator()(const int & , const int & , const int & ) const;
    const vector3d<T,Nx,Ny,Nz> & operator=(const T &a)
    {
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<Ny;j++)
            {
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]=a;
                }
            }
        }
        return *this;
    }
    const vector3d<T,Nx,Ny,Nz> & operator=(const vector3d<T,Nx,Ny,Nz> & m)
    {
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<Ny;j++)
            {
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]=m(i,j,k);
                }
            }
        }
        return *this;
    }
    template<typename E> const vector3d<T,Nx,Ny,Nz> & operator=(const vector3dExpression<E,T,Nx,Ny,Nz> & m)
    {
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<Ny;j++)
            {
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]=m(i,j,k);
                }
            }
        }
        return *this;
    }
    template<typename E> const vector3d<T,Nx,Ny,Nz> & operator+=(const vector3dExpression<E,T,Nx,Ny,Nz> & m)
    {
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<Ny;j++)
            {
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]+=m(i,j,k);
                }
            }
        }
        return *this;
    }
    template<typename E> const vector3d<T,Nx,Ny,Nz> & operator-=(const vector3dExpression<E,T,Nx,Ny,Nz> & m)
    {
        for(int i=0;i<Nx;i++)
        {
            for(int j=0;j<Ny;j++)
            {
                for(int k=0;k<Nz;k++)
                {
                    _data[i][j][k]-=m(i,j,k);
                }
            }
        }
        return *this;
    }
    const vector3d<T,Nx,Ny,Nz> & operator/=(const T & );
    const vector3d<T,Nx,Ny,Nz> & operator*=(const T & );
private:
    T *** _data;
};

template<typename T, int Nx,int Ny,int Nz>
vector3d<T,Nx,Ny,Nz>::vector3d()
{
    _data=new T ** [Nx];
    for(int i=0;i<Nx;i++)
    {
        _data[i]=new T * [Ny];
        for(int j=0;j<Ny;j++)
        {
            _data[i][j]=new T [Nz];
            for(int k=0;k<Nz;k++)
            {
                _data[i][j][k]=0;
            }
        }
    }
}

template<typename T, int Nx,int Ny,int Nz>
vector3d<T,Nx,Ny,Nz>::vector3d(const vector3d<T,Nx,Ny,Nz> & m)
{
    _data=new T ** [Nx];
    for(int i=0;i<Nx;i++)
    {
        _data[i]=new T * [Ny];
        for(int j=0;j<Ny;j++)
        {
            _data[i][j]=new T [Nz];
            for(int k=0;k<Nz;k++)
            {
                _data[i][j][k]=m(i,j,k);
            }
        }
    }
}

template<typename T, int Nx,int Ny,int Nz>
vector3d<T,Nx,Ny,Nz>::~vector3d()
{
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            delete [] _data[i][j];
        }
        delete [] _data[i];
    }
    delete [] _data;
}

template<typename T, int Nx,int Ny,int Nz>
T & vector3d<T,Nx,Ny,Nz>::operator()(const int & i, const int & j, const int & k)
{
    return _data[i][j][k];
}

template<typename T, int Nx,int Ny,int Nz>
const T & vector3d<T,Nx,Ny,Nz>::operator()(const int & i, const int & j, const int & k) const
{
    return _data[i][j][k];
}

template<typename T, int Nx,int Ny,int Nz>
const vector3d<T,Nx,Ny,Nz> & vector3d<T,Nx,Ny,Nz>::operator*=(const T & a)
{
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            for(int k=0;k<Nz;k++)
            {
                _data[i][j][k]*=a;
            }
        }
    }
    return *this;
}

template<typename T, int Nx,int Ny,int Nz>
const vector3d<T,Nx,Ny,Nz> & vector3d<T,Nx,Ny,Nz>::operator/=(const T & a)
{
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            for(int k=0;k<Nz;k++)
            {
                _data[i][j][k]/=a;
            }
        }
    }
    return *this;
}

template <typename E1, typename E2, typename T, int Nx,int Ny,int Nz>
class vector3dplus : public vector3dExpression<vector3dplus<E1, E2, T, Nx,Ny,Nz>, T, Nx,Ny,Nz>
{
    public:
        vector3dplus(vector3dExpression<E1, T, Nx,Ny,Nz > const& u, vector3dExpression<E2, T, Nx,Ny,Nz> const& v) : _u(u), _v(v){};
        const T operator()(const int & i, const int & j, const int & k) const { return _u(i,j,k) + _v(i,j,k); }
    private:
        E1 const& _u;
        E2 const& _v;
};

template<typename E1, typename E2, typename T, int Nx,int Ny,int Nz>
const vector3dplus<E1,E2,T,Nx,Ny,Nz> operator+(const vector3dExpression<E1,T,Nx,Ny,Nz> & u , const vector3dExpression<E2,T,Nx,Ny,Nz> & v)
{
    return vector3dplus<E1,E2,T,Nx,Ny,Nz>(u,v);
}

template <typename E1, typename E2, typename T, int Nx,int Ny,int Nz>
class vector3ddifference : public vector3dExpression<vector3ddifference<E1, E2, T, Nx,Ny,Nz>, T, Nx,Ny,Nz>
{
    public:
        vector3ddifference(vector3dExpression<E1, T, Nx,Ny,Nz > const& u, vector3dExpression<E2, T, Nx,Ny,Nz> const& v) : _u(u), _v(v){};
        const T operator()(const int & i, const int & j,const int & k) const { return _u(i,j,k) - _v(i,j,k); }
    private:
        E1 const& _u;
        E2 const& _v;
};

template<typename E1, typename E2, typename T, int Nx,int Ny,int Nz>
const vector3ddifference<E1,E2,T,Nx,Ny,Nz> operator-(const vector3dExpression<E1,T,Nx,Ny,Nz> & u , const vector3dExpression<E2,T,Nx,Ny,Nz> & v)
{
    return vector3ddifference<E1,E2,T,Nx,Ny,Nz>(u,v);
}

template<typename E, typename T, int Nx,int Ny,int Nz>
class vector3dScalarProduct: public vector3dExpression<vector3dScalarProduct<E,T,Nx,Ny,Nz>,T,Nx,Ny,Nz>
{
public:
    vector3dScalarProduct(vector3dExpression<E,T,Nx,Ny,Nz> const & u, const T & a):_u(u),_a(a){};
    const T operator()(const int & i, const int & j,const int & k) const { return _u(i,j,k) *_a; }
private:
    E const & _u;
    T const & _a;
};

template<typename E, typename T, int Nx,int Ny,int Nz>
const vector3dScalarProduct<E,T,Nx,Ny,Nz> operator*(const vector3dExpression<E,T,Nx,Ny,Nz> & u, const T & a)
{
    return vector3dScalarProduct<E,T,Nx,Ny,Nz>(u,a);
}

template<typename E, typename T, int Nx,int Ny,int Nz>
const vector3dScalarProduct<E,T,Nx,Ny,Nz> operator*(const T & a, const vector3dExpression<E,T,Nx,Ny,Nz> & u)
{
    return vector3dScalarProduct<E,T,Nx,Ny,Nz>(u,a);
}

template<typename E, typename T, int Nx,int Ny,int Nz>
class vector3dScalarDivision: public vector3dExpression<vector3dScalarDivision<E,T,Nx,Ny,Nz>,T,Nx,Ny,Nz>
{
public:
    vector3dScalarDivision(vector3dExpression<E,T,Nx,Ny,Nz> const & u, const T & a):_u(u),_a(a){};
    const T operator()(const int & i, const int & j,const int & k) const { return _u(i,j,k) / _a; }
private:
    E const & _u;
    T const & _a;
};

template<typename E, typename T, int Nx,int Ny,int Nz>
const vector3dScalarDivision<E,T,Nx,Ny,Nz> operator/(const vector3dExpression<E,T,Nx,Ny,Nz> & u, const T & a)
{
    return vector3dScalarDivision<E,T,Nx,Ny,Nz>(u,a);
}

template<typename E1, typename E2, typename T, int Nx,int Ny,int Nz>
const T operator*(const vector3dExpression<E1,T,Nx,Ny,Nz> & u, const vector3dExpression<E2,T,Nx,Ny,Nz> & v)
{
    T temp=0;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            for(int k=0;k<Nz;k++)
            {
                temp+=u(i,j,k)*v(i,j,k);
            }
        }
    }
    return temp/(Nx*Ny*Nz);
}

template<typename E, typename T, int Nx,int Ny,int Nz>
const T norm(const vector3dExpression<E,T,Nx,Ny,Nz> & u)
{
    T temp=0;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            for(int k=0;k<Nz;k++)
            {
                temp+=u(i,j,k)*u(i,j,k);
            }
        }
    }
    return sqrt(temp/(Nx*Ny*Nz));
}

template<typename E, typename T, int Nx,int Ny,int Nz>
const T normmax(const vector3dExpression<E,T,Nx,Ny,Nz> & u)
{
    T temp=0;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0;j<Ny;j++)
        {
            for(int k=0;k<Nz;k++)
            {
                temp=((fabs(u(i,j,k))>temp)?fabs(u(i,j,k)):temp);
            }
        }
    }
    return temp;
}

#endif // vector3d_H
