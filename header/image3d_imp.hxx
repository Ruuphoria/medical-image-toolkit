/*
=========================================================================================
3D Image Toolkit

Copyright 2013:
Marco Fedele, Luca Barbarotta, Francesco Cremonesi, Elena Faggiano.
All rights reserved.
Read details of the license in the file license.txt provided with the library.
=========================================================================================
*/

/*!
  \file image3d_imp.hxx

   \brief file with the implementation of all the methods of class \ref im3d::image3d
*/

#ifndef IMAGE3D_IMP_HPP_INCLUDED
#define IMAGE3D_IMP_HPP_INCLUDED


#include "interface.hxx"



//CONSTRUCTORS

template <typename T>
im3d::image3d<T>::image3d () :
    rawimage (0), dimx (0), dimy (0), dimz (1),
    hx (0), hy (0), hz (1)
{   }


template <typename T>
im3d::image3d<T>::image3d (uint const& x, uint const& y, uint const& z,
                           double const& hx, double const& hy, double const& hz) :
    dimx (x), dimy (y), dimz (z), hx (hx), hy (hy), hz (hz)
{
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
}


template <typename T>
im3d::image3d<T>::image3d (image3d const& tocopy)
{
    this->setdim (tocopy.getdimx(), tocopy.getdimy(), tocopy.getdimz() );
    this->seth (tocopy.gethx(), tocopy.gethy(), tocopy.gethz() );

    #pragma omp parallel for
    for (uint i = 0; i < dimx; ++i)
        for (uint j = 0; j < dimy; ++j)
            for (uint k = 0; k < dimz; ++k)
            {
                (*this) (i, j, k) = tocopy (i, j, k);
            }
}



//MEMBERS TO GET AND SET PRIVATE PARAMETERS

template <typename T>
void im3d::image3d<T>::setdimx (uint const& x)
{
    this->dimx = x;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdimy (uint const& y)
{
    this->dimy = y;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdimz (uint const& z)
{
    this->dimz = z;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz);
    return;
}


template <typename T>
void im3d::image3d<T>::setdim (uint const& x, uint const& y, uint const& z,
                               T const& value)
{
    this->dimx = x;
    this->dimy = y;
    this->dimz = z;
    this->rawimage.resize (0);
    this->rawimage.resize ( this->dimx * this->dimy * this->dimz, value);
    return;
}


template <typename T>
void im3d::image3d<T>::sethx (double const& hx)
{
    this->hx = hx;
    return;
}


template <typename T>
void im3d::image3d<T>::sethy (double const& hy)
{
    this->hy = hy;
    return;
}


template <typename T>
void im3d::image3d<T>::sethz (double const& hz)
{
    this->hz = hz;
    return;
}


template <typename T>
void im3d::image3d<T>::seth (double const& hx, double const& hy, double const& hz)
{
    this->hx = hx;
    this->hy = hy;
    this->hz = hz;
    return;
}

//USEFULL OVERLOADING OF VAROIUS OPERATORS

template <typename T>
T im3d::image3d<T>::operator() (uint const& i, uint const& j, uint const& k) const
{
    return rawimage[ this->dimy * this->dimz * i + this->dimz * j + k ];
}


template <typename T>
T& im3d::image3d<T>::operator() (uint const& i, uint const& j, uint const& k)
{
    return rawimage[ this->dimy * this->dimz * i + this->dimz * j + k ];
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator= (S const& toassign)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] = static_cast<T> (toassign);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator= (image3d<S> const& toassign)
{
    this->setdim (toassign.getdimx(), toassign.getdimy(), toassign.getdimz() );
    this->seth (toassign.gethx(), toassign.gethy(), toassign.gethz() );

    #pragma omp parallel for
    for (uint i = 0; i < this->dimx; ++i)
        for (uint j = 0; j < this->dimy; ++j)
            for (uint k = 0; k < this->dimz; ++k)
            {
                (*this) (i, j, k) = static_cast<T> (toassign (i, j, k) );
            }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator += (S const& addend)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] += static_cast<T> (addend);
    }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator += (image3d<S> const& addend)
{
    if (this->dimx == addend.getdimx() &&
            this->dimy == addend.getdimy() &&
            this->dimz == addend.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) += static_cast<T> (addend (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator+= : dimensions must agree" << std::endl;
    }

    return *this;
}




template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator -= (S const& addend)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] -= static_cast<T> (addend);
    }

    return *this;
}



template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator -= (image3d<S> const& addend)
{
    if (this->dimx == addend.getdimx() &&
            this->dimy == addend.getdimy() &&
            this->dimz == addend.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) -= static_cast<T> (addend (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator-= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator *= (S const& factor)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] *= static_cast<T> (factor);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator *= (image3d<S> const& factor)
{
    if (this->dimx == factor.getdimx() &&
            this->dimy == factor.getdimy() &&
            this->dimz == factor.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) *= static_cast<T> (factor (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator*= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator /= (S const& factor)
{
    #pragma omp parallel for
    for (uint i = 0; i < this->dimx * this->dimy * this->dimz; ++i)
    {
        this->rawimage[i] /= static_cast<T> (factor);
    }

    return *this;
}


template <typename T>
template <typename S>
im3d::image3d<T>& im3d::image3d<T>::operator /= (image3d<S> const& factor)
{
    if (this->dimx == factor.getdimx() &&
            this->dimy == factor.getdimy() &&
            this->dimz == factor.getdimz() )
    {
        #pragma omp parallel for
        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    (*this) (i, j, k) /= static_cast<T> (factor (i, j, k) );
                }
    }
    else
    {
        std::cout << "WARNING::image3d::operator/= : dimensions must agree" << std::endl;
    }

    return *this;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator + (image3d<S> const& addend1, R const& addend2)
{
    image3d<S> sum (addend1);
    return sum += addend2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator - (image3d<S> const& addend1, R const& addend2)
{
    image3d<S> sum (addend1);
    return sum -= addend2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator * (image3d<S> const& factor1, R const& factor2)
{
    image3d<S> prod (factor1);
    return prod *= factor2;
}


template <typename S, typename R>
im3d::image3d<S> const im3d::operator / (image3d<S> const& factor1, R const& factor2)
{
    image3d<S> prod (factor1);
    return prod /= factor2;
}



//MEMBERS WORKING ON IMAGE'S VALUES


template <typename S, typename R>
S const im3d::scalarprod ( image3d<S> const& factor1, image3d<R> const& factor2 )
{
    if (factor2.dimx == factor1.dimx &&
            factor2.dimy == factor1.dimy &&
            factor2.dimz == factor1.dimz )
    {
        S result (0);

        #pragma omp parallel for reduction (+:result)
        for (uint i = 0; i < factor1.dimx * factor1.dimy * factor1.dimz; ++i)
        {
            result += factor1.rawimage[i] * static_cast<S> (factor2.rawimage[i]);
        }