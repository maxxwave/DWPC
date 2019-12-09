// Simple multi-dimensional arrays in C++
// Author: Matt Ellis

#ifndef __ARRAY_H__
#define __ARRAY_H__

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif


#include <vector>
#include <iostream>
#include <cassert>

template <int _D, typename _T>
class array_t
{
    public:

        array_t ()
        {
            _Nelements = 0;
            for( int i = 0; i < _D; i++) _size[i] = 0;
        }

        ~array_t()
        {};

        array_t ( const int* n, const _T value = _T(0) )
        {
            _size[0] = n[0];
            _Nelements = _size[0];
            for( int i = 1; i < _D; i++) {
                _size[i] = n[i];
                _Nelements *= _size[i];
            }
            _data.assign( _Nelements, value);
        }

        size_t size( const int dim)
        {
            assert(dim < _D);
            return _size[dim];
        }

        void clear ()
        {
            _Nelements = 0;
            for( int i = 0; i < _D; i++)
                _size[i] = 0;
            _data.clear();
        }

        inline _T& RESTRICT operator () ( const int* n)
        {
            assert( n[0] < _size[0]);
            size_t index = n[0];
            for( int i = 1; i < _D; i++) {
                assert( n[i] < _size[i] );
                index = n[i] + _size[i] * index;
            }
            return _data[index];
        }



    private:
        size_t              _Nelements;
        size_t              _size[_D];
        std::vector<_T>     _data;


};

template < typename _T >
class array_t<1, _T>
{

    public:

        array_t()
            : _Nelements(0), _size(0), _data()
        {};

        ~array_t()
        {};

        array_t ( const int n, const _T value )
        {
            _size = n;
            _Nelements = n;
            _data.assign( n, value);
        }

        size_t size()
        {
            return _size;
        }

        size_t size( const int dim)
        {
            assert(dim < 1);
            return _size;
        }

        void clear ()
        {
            _Nelements = 0;
            _size = 0;
            _data.clear();
        }

        void assign( const int n0, const _T value )
        {
            clear();
            _size = n0;
            _Nelements = n0;
            _data.assign( _Nelements, value);
        }

        inline _T& RESTRICT operator () ( const int n)
        {
            assert( n < _size);
            return _data[n];
        }

    private:

        size_t _Nelements;
        size_t _size;
        std::vector<_T> _data;
};


template < typename _T >
class array_t<2, _T>
{

    public:

        array_t()
            : _Nelements(0), _size(), _data()
        {};

        ~array_t()
        {};

        array_t ( const int n0, const int n1, const _T value )
        {
            _size[0] = n0;
            _size[1] = n1;
            _Nelements = n0*n1;
            _data.assign( _Nelements, value);
        }

        size_t size( const int dim)
        {
            assert(dim < 2);
            return _size[dim];
        }

        void clear ()
        {
            _Nelements = 0;
            _size[0] = 0;
            _size[1] = 0;
            _data.clear();
        }

        void assign( const int n0, const int n1, const _T value )
        {
            clear();
            _size[0] = n0;
            _size[1] = n1;
            _Nelements = n0*n1;
            _data.assign( _Nelements, value);
        }

        inline _T& RESTRICT operator () ( const int n0, const int n1)
        {
            assert( n0 < _size[0]);
            assert( n1 < _size[1]);
            return _data[ n1 + _size[1] * n0];
        }

    private:

        size_t _Nelements;
        size_t _size[2];
        std::vector<_T> _data;
};


template < typename _T >
class array_t<3, _T>
{

    public:

        array_t()
            : _Nelements(0), _size(), _data()
        {};

        ~array_t()
        {};

        array_t ( const int n0, const int n1, const int n2)
        {
            _size[0] = n0;
            _size[1] = n1;
            _size[2] = n2;
            _Nelements = n0*n1*n2;
            _data.assign( _Nelements, _T(0));
        }

        size_t size( const int dim)
        {
            assert(dim < 3);
            return _size[dim];
        }

        array_t ( const int n0, const int n1, const int n2, const _T value )
        {
            _size[0] = n0;
            _size[1] = n1;
            _size[2] = n2;
            _Nelements = n0*n1*n2;
            _data.assign( _Nelements, value);
        }

        void assign( const int n0, const int n1, const int n2, const _T value)
        {
            _size[0] = n0;
            _size[1] = n1;
            _size[2] = n2;
            _Nelements = n0*n1*n2;
            _data.assign( _Nelements, value);
        }


        void clear ()
        {
            _Nelements = 0;
            _size[0] = 0;
            _size[1] = 0;
            _size[2] = 0;
            _data.clear();
        }

        inline _T& RESTRICT operator () ( const int n0, const int n1, const int n2)
        {
            assert( n0 < _size[0]);
            assert( n1 < _size[1]);
            assert( n2 < _size[2]);
            return _data[ n2 + _size[2]*(n1 + _size[1] * n0) ];
        }

    private:

        size_t _Nelements;
        size_t _size[3];
        std::vector<_T> _data;
};

template < typename _T >
class array_t<4, _T>
{

    public:

        array_t()
            : _Nelements(0), _size(), _data()
        {};

        ~array_t()
        {};

        array_t ( const int n0, const int n1, const int n2, const int n3)
        {
            _size[0] = n0;
            _size[1] = n1;
            _size[2] = n2;
            _size[3] = n3;
            _Nelements = n0*n1*n2*n3;
            _data.assign( _Nelements, _T(0));
        }

        array_t ( const int n0, const int n1, const int n2, const int n3, const _T value )
        {
            _size[0] = n0;
            _size[1] = n1;
            _size[2] = n2;
            _size[3] = n3;
            _Nelements = n0*n1*n2*n3;
            _data.assign( _Nelements, value);
        }

        size_t size( const int dim)
        {
            assert(dim < 4);
            return _size[dim];
        }

        void clear ()
        {
            _Nelements = 0;
            _size[0] = 0;
            _size[1] = 0;
            _size[2] = 0;
            _size[3] = 0;
            _data.clear();
        }

        inline _T& RESTRICT operator () ( const int n0, const int n1, const int n2, const int n3)
        {
            assert( n0 < _size[0]);
            assert( n1 < _size[1]);
            assert( n2 < _size[2]);
            assert( n3 < _size[3]);
            return _data[ n3 + _size[3]*(n2 + _size[2]*(n1 + _size[1] * n0)) ];
        }

    private:

        size_t _Nelements;
        size_t _size[4];
        std::vector<_T> _data;
};


#endif //__ARRAY_H__
