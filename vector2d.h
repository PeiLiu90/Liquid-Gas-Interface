#ifndef vector2d_H
#define vector2d_H

#include "math.h"
#include <iostream>
using namespace std;

template <typename E, typename T, int N, int M>
class vector2dExpression
{
public:
    const T operator()(const int & i, const int & j) const{return static_cast<const E&>(*this)(i,j);}
    operator E&() { return static_cast<E&>(*this);}
    operator E const&() const { return static_cast<const E&>(*this);}
};

template<typename T ,int N, int M>
class vector2d: public vector2dExpression<vector2d<T,N,M>, T, N, M>
{
public:
    vector2d();
    vector2d(const vector2d &);
    ~vector2d();
    template<typename E> vector2d(const vector2dExpression<E,T,N,M> & m)
    {
        _data=new T * [N];
        for(int i=0;i<N;i++)
        {
            _data[i]=new T [M];
            for(int j=0;j<M;j++)
            {
                _data[i][j]=m(i,j);
            }
        }
    }

    T & operator()(const int & ,const int & );
    const T & operator()(const int & , const int & ) const;
    const vector2d<T,N,M> & operator=(const T &a)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _data[i][j]=a;
            }
        }
        return *this;
    }
    const vector2d<T,N,M> & operator=(const vector2d<T,N,M> & m)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _data[i][j]=m(i,j);
            }
        }
        return *this;
    }
    template<typename E> const vector2d<T,N,M> & operator=(const vector2dExpression<E,T,N,M> & m)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _data[i][j]=m(i,j);
            }
        }
        return *this;
    }
    template<typename E> const vector2d<T,N,M> & operator+=(const vector2dExpression<E,T,N,M> & m)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _data[i][j]+=m(i,j);
            }
        }
        return *this;
    }
    template<typename E> const vector2d<T,N,M> & operator-=(const vector2dExpression<E,T,N,M> & m)
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<M;j++)
            {
                _data[i][j]-=m(i,j);
            }
        }
        return *this;
    }
    const vector2d<T,N,M> & operator/=(const T & );
    const vector2d<T,N,M> & operator*=(const T & );
private:
    T ** _data;
};

template<typename T, int N, int M>
vector2d<T,N,M>::vector2d()
{
    _data=new T * [N];
    for(int i=0;i<N;i++)
    {
        _data[i]=new T [M];
        for(int j=0;j<M;j++)
        {
            _data[i][j]=0;
        }
    }
}

template<typename T, int N, int M>
vector2d<T,N,M>::vector2d(const vector2d<T,N,M> & m)
{
    _data=new T * [N];
    for(int i=0;i<N;i++)
    {
        _data[i]=new T [M];
        for(int j=0;j<M;j++)
        {
            _data[i][j]=m(i,j);
        }
    }
}

template<typename T, int N, int M>
vector2d<T,N,M>::~vector2d()
{
    for(int i=0;i<N;i++)
    {
        delete [] _data[i];
    }
    delete [] _data;
}

template<typename T, int N, int M>
T & vector2d<T,N,M>::operator()(const int & i, const int & j)
{
    return _data[i][j];
}

template<typename T, int N, int M>
const T & vector2d<T,N,M>::operator()(const int & i, const int & j) const
{
    return _data[i][j];
}

template<typename T, int N, int M>
const vector2d<T,N,M> & vector2d<T,N,M>::operator*=(const T & a)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            _data[i][j]*=a;
        }
    }
    return *this;
}

template<typename T, int N, int M>
const vector2d<T,N,M> & vector2d<T,N,M>::operator/=(const T & a)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            _data[i][j]/=a;
        }
    }
    return *this;
}

template <typename E1, typename E2, typename T, int N, int M>
class vector2dplus : public vector2dExpression<vector2dplus<E1, E2, T, N, M>, T, N, M>
{
    public:
        vector2dplus(vector2dExpression<E1, T, N, M > const& u, vector2dExpression<E2, T, N, M> const& v) : _u(u), _v(v){};
        const T operator()(const int & i, const int & j) const { return _u(i,j) + _v(i,j); }
    private:
        E1 const& _u;
        E2 const& _v;
};

template<typename E1, typename E2, typename T, int N, int M>
const vector2dplus<E1,E2,T,N,M> operator+(const vector2dExpression<E1,T,N,M> & u , const vector2dExpression<E2,T,N,M> & v)
{
    return vector2dplus<E1,E2,T,N,M>(u,v);
}

template <typename E1, typename E2, typename T, int N, int M>
class vector2ddifference : public vector2dExpression<vector2ddifference<E1, E2, T, N, M>, T, N, M>
{
    public:
        vector2ddifference(vector2dExpression<E1, T, N, M > const& u, vector2dExpression<E2, T, N, M> const& v) : _u(u), _v(v){};
        const T operator()(const int & i, const int & j) const { return _u(i,j) - _v(i,j); }
    private:
        E1 const& _u;
        E2 const& _v;
};

template<typename E1, typename E2, typename T, int N, int M>
const vector2ddifference<E1,E2,T,N,M> operator-(const vector2dExpression<E1,T,N,M> & u , const vector2dExpression<E2,T,N,M> & v)
{
    return vector2ddifference<E1,E2,T,N,M>(u,v);
}

template<typename E, typename T, int N, int M>
class vector2dScalarProduct: public vector2dExpression<vector2dScalarProduct<E,T,N,M>,T,N,M>
{
public:
    vector2dScalarProduct(vector2dExpression<E,T,N,M> const & u, const T & a):_u(u),_a(a){};
    const T operator()(const int & i, const int & j) const {return _u(i,j)*_a;}
private:
    E const & _u;
    T const & _a;
};

template<typename E, typename T, int N, int M>
const vector2dScalarProduct<E,T,N,M> operator*(const vector2dExpression<E,T,N,M> & u, const T & a)
{
    return vector2dScalarProduct<E,T,N,M>(u,a);
}

template<typename E, typename T, int N, int M>
const vector2dScalarProduct<E,T,N,M> operator*(const T & a, const vector2dExpression<E,T,N,M> & u)
{
    return vector2dScalarProduct<E,T,N,M>(u,a);
}

template<typename E, typename T, int N, int M>
class vector2dScalarDivision: public vector2dExpression<vector2dScalarDivision<E,T,N,M>,T,N,M>
{
public:
    vector2dScalarDivision(vector2dExpression<E,T,N,M> const & u, const T & a):_u(u),_a(a){};
    const T operator()(const int & i, const int & j) const {return _u(i,j)/_a;}
private:
    E const & _u;
    T const & _a;
};

template<typename E, typename T, int N, int M>
const vector2dScalarDivision<E,T,N,M> operator/(const vector2dExpression<E,T,N,M> & u, const T & a)
{
    return vector2dScalarDivision<E,T,N,M>(u,a);
}



template<typename E1, typename E2, typename T, int N, int M>
const T operator*(const vector2dExpression<E1,T,N,M> & u, const vector2dExpression<E2,T,N,M> & v)
{
    T temp=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            temp+=u(i,j)*v(i,j);
        }
    }
    return temp/(M*N);
}

template<typename E, typename T, int N, int M>
const T norm(const vector2dExpression<E,T,N,M> & u)
{
    T temp=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            temp+=u(i,j)*u(i,j);
        }
    }
    return sqrt(temp/(M*N));
}

template<typename E, typename T, int N, int M>
const T normmax(const vector2dExpression<E,T,N,M> & u)
{
    T temp=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            temp=((fabs(u(i,j))>temp)?fabs(u(i,j)):temp);
        }
    }
    return temp;
}

template<typename E, typename T, int N, int M>
ostream & operator<<(ostream & cout, const vector2dExpression<E,T,N,M> & u)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            cout<<u(i,j)<<"  ";
        }
        cout<<endl;
    }
    return cout;
}

#endif // vector2d_H
