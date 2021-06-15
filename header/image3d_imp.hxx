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

        return result;
    }
    else
    {
        std::cout << "WARNING::scalarprod: dimensions must agree, 0 is returned" << std::endl;
        return 0;
    }
}



template <typename S, typename R>
S const im3d::scalarprod_L2 ( image3d<S> const& factor1, image3d<R> const& factor2 )
{

    if (factor2.hx == factor1.hx &&
            factor2.hy == factor1.hy &&
            factor2.hz == factor1.hz )
    {
        return factor1.hx * factor1.hy * factor1.hz * scalarprod (factor1, factor2);
    }
    else
    {
        std::cout << "WARNING::scalarprod_L2: spacing must agree, 0 is returned" << std::endl;
        return 0;
    }
}



template <typename T>
void im3d::image3d<T>::crop (image3d<T>& res,
                             uint const& XSTART, uint const& YSTART, uint const& ZSTART,
                             uint const& XEND, uint const& YEND, uint const& ZEND) const
{

    uint xstart (XSTART), ystart (YSTART), zstart (ZSTART), xend (XEND), yend (YEND), zend (ZEND);

    // check if coordinates are inside the image
    if (xend > this->dimx)
    {
        xend = this->dimx - 1;
    }
    if (xstart > this->dimx)
    {
        xstart = this->dimx - 1;
    }

    if (ystart > this->dimy)
    {
        ystart = this->dimy - 1;
    }
    if (yend > this->dimy)
    {
        yend = this->dimy - 1;
    }


    // check if starts are smaller than ends
    if (xstart > xend)
    {
        uint aux = xstart;
        xstart = xend;
        xend = aux;
    }

    if (ystart > yend)
    {
        uint aux = ystart;
        ystart = yend;
        yend = aux;
    }

    // if 2d, set zstart=zend=0
    if (this->dimz == 1)
    {
        zstart = 0;
        zend = 0;
    }
    // check if z coordinates are inside the image and start is smaller than end
    else
    {
        if (zstart > this->dimz)
        {
            zstart = this->dimz - 1;
        }
        if (zend > this->dimz)
        {
            zend = this->dimz - 1;
        }
        if (zstart > zend)
        {
            uint aux = zstart;
            zstart = zend;
            zend = aux;
        }
    }

    int newdimx = xend - xstart + 1;
    int newdimy = yend - ystart + 1;
    int newdimz = zend - zstart + 1;

    res.setdim ( newdimx, newdimy, newdimz);
    res.seth ( this->hx, this->hy, this->hz);

    #pragma omp parallel for
    for (uint i = 0; i < res.getdimx(); ++i)
        for (uint j = 0; j < res.getdimy(); ++j)
            for (uint k = 0; k < res.getdimz(); ++k)
            {
                res (i, j, k) = (*this) (i + xstart, j + ystart, k + zstart) ;
            }

    return;
}



template <typename T>
void im3d::image3d<T>::crop (image3d<T>& res,
                             double const& XSTART, double const& YSTART, double const& ZSTART,
                             double const& XEND, double const& YEND, double const& ZEND) const
{
    double xstart (XSTART), ystart (YSTART), zstart (ZSTART), xend (XEND), yend (YEND), zend (ZEND);

    // if negative set it to zero
    if (xstart < 0)
    {
        xstart = 0;
    }
    if (ystart < 0)
    {
        ystart = 0;
    }
    if (zstart < 0)
    {
        zstart = 0;
    }
    if (xend < 0)
    {
        xend = 0;
    }
    if (yend < 0)
    {
        yend = 0;
    }
    if (zend < 0)
    {
        zend = 0;
    }

    //  uint xs = floor(dimx*xstart), ys = floor(dimy*ystart), zs = floor(dimz*zstart);
    //  uint xe = ceil(dimx*xend), ye = ceil(dimy*yend), ze = ceil(dimz*zend);

    uint xs = floor (xstart / hx), ys = floor (ystart / hy), zs = floor (zstart / hz);
    uint xe = ceil (xend / hx), ye = ceil (yend / hy), ze = ceil (zend / hz);

    // in crop will be checked the indexes
    this->crop (res, xs, ys, zs, xe, ye, ze);

    return;
}



template <typename T>
void im3d::image3d<T>::change_resolution (image3d<T>& res, uint ratio, bool increase) const
{

    // decreasing resolution
    if (!increase)
    {
        res.setdim (ceil ( static_cast<T> (this->dimx) / static_cast<T> (ratio) ),
                    ceil ( static_cast<T> (this->dimy) / static_cast<T> (ratio) ),
                    ceil ( static_cast<T> (this->dimz) / static_cast<T> (ratio) ) );

        res.seth ( this->hx * ratio, this->hy * ratio, this->hz * ratio );

        #pragma omp parallel for
        for (uint i = 0; i < res.getdimx(); ++i)
            for (uint j = 0; j < res.getdimy(); ++j)
                for (uint k = 0; k < res.getdimz(); ++k)
                {
                    res (i, j, k) = (*this) (i * ratio, j * ratio, k * ratio) ;
                }
    }

    // increasing resolution
    else
    {
        // make ratio a power of 2
        while ( (ratio % 2) != 0 )
        {
            ++ratio;
        }

        image3d<T> resold (*this);

        for (uint r = 2; r <= ratio; r *= 2)
        {
            res.setdim ( (this->dimx - 1) *r + 1, (this->dimy - 1) *r + 1, (this->dimz - 1) *r + 1 );
            res.seth (this->hx / static_cast<double> (r),
                      this->hy / static_cast<double> (r),
                      this->hz / static_cast<double> (r) );

            // copying resold elements every 2 pixels
            #pragma omp parallel for
            for (uint i = 0; i < res.getdimx(); i += 2)
                for (uint j = 0; j < res.getdimy(); j += 2)
                    for (uint k = 0; k < res.getdimz(); k += 2)
                    {
                        res (i, j, k) = resold (i / 2, j / 2, k / 2) ;
                    }


            // 3d case
            if (dimz != 1)
            {

                #pragma omp parallel for
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                        for (uint k = 0; k < res.getdimz() - 2; k += 2)
                        {

                            res (i, j, k + 1) = ( res (i, j, k) + res (i, j, k + 2) ) / 2.;

                            res (i + 1, j, k) = ( res (i, j, k) + res (i + 2, j, k) ) / 2.;

                            res (i, j + 1, k) = ( res (i, j, k) + res (i, j + 2, k) ) / 2.;

                            res (i + 1, j + 1, k + 1) = ( res (i, j, k) + res (i, j, k + 2) +
                                                          res (i, j + 2, k) + res (i, j + 2, k + 2) +
                                                          res (i + 2, j, k) + res (i + 2, j, k + 2) +
                                                          res (i + 2, j + 2, k) + res (i + 2, j + 2, k + 2) ) / 8.;

                            res (i + 1, j, k + 1) = ( res (i, j, k) + res (i, j, k + 2) +
                                                      res (i + 2, j, k) + res (i + 2, j, k + 2) ) / 4.;

                            res (i + 1, j + 1, k) = ( res (i, j, k) + res (i + 2, j, k) +
                                                      res (i, j + 2, k) + res (i + 2, j + 2, k) ) / 4.;

                            res (i, j + 1, k + 1) = ( res (i, j, k) + res (i, j + 2, k) +
                                                      res (i, j, k + 2) + res (i, j + 2, k + 2) ) / 4.;
                        }


                //  complete interpolation on the last three faces
                //  of the domain with an even index
                uint X (res.getdimx() - 1), Y (res.getdimy() - 1), Z (res.getdimz() - 1);

                // face i=X
                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    for (uint k = 0; k < res.getdimz() - 2; k += 2)
                    {
                        res (X, j, k + 1) = ( res (X, j, k) + res (X, j, k + 2) ) / 2.;
                        res (X, j + 1, k) = ( res (X, j, k) + res (X, j + 2, k) ) / 2.;
                        res (X, j + 1, k + 1) = ( res (X, j, k) + res (X, j + 2, k) +
                                                  res (X, j, k + 2) + res (X, j + 2, k + 2) ) / 4.;
                    }

                // face j=Y
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint k = 0; k < res.getdimz() - 2; k += 2)
                    {
                        res (i, Y, k + 1) = ( res (i, Y, k) + res (i, Y, k + 2) ) / 2.;
                        res (i + 1, Y, k) = ( res (i, Y, k) + res (i + 2, Y, k) ) / 2.;
                        res (i + 1, Y, k + 1) = ( res (i, Y, k) + res (i, Y, k + 2) +
                                                  res (i + 2, Y, k) + res (i + 2, Y, k + 2) ) / 4.;
                    }

                // face k=Z
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    {
                        res (i + 1, j, Z) = ( res (i, j, Z) + res (i + 2, j, Z) ) / 2.;
                        res (i, j + 1, Z) = ( res (i, j, Z) + res (i, j + 2, Z) ) / 2.;
                        res (i + 1, j + 1, Z) = ( res (i, j, Z) + res (i + 2, j, Z) +
                                                  res (i, j + 2, Z) + res (i + 2, j + 2, Z) ) / 4.;
                    }

                // three common edges
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                {
                    res (i + 1, Y, Z) = ( res (i, Y, Z) + res (i + 2, Y, Z) ) / 2.;
                }

                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                {
                    res (X, j + 1, Z) = ( res (X, j, Z) + res (X, j + 2, Z) ) / 2.;
                }

                for (uint k = 0; k < res.getdimz() - 2; k += 2)
                {
                    res (X, Y, k + 1) = ( res (X, Y, k) + res (X, Y, k + 2) ) / 2.;
                }

            }// end 3d case

            // 2d case
            else
            {
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                    for (uint j = 0; j < res.getdimy() - 2; j += 2)
                    {

                        res (i + 1, j, 0) = ( res (i, j, 0) + res (i + 2, j, 0) ) / 2.;

                        res (i, j + 1, 0) = ( res (i, j, 0) + res (i, j + 2, 0) ) / 2.;

                        res (i + 1, j + 1, 0) = ( res (i, j, 0) + res (i + 2, j, 0) +
                                                  res (i, j + 2, 0) + res (i + 2, j + 2, 0) ) / 4.;
                    }

                //  complete interpolation on the last two edges
                //  of the domain with an even index
                uint X (res.getdimx() - 1), Y (res.getdimy() - 1);

                //  edge i=X
                for (uint j = 0; j < res.getdimy() - 2; j += 2)
                {
                    res (X, j + 1, 0) = ( res (X, j, 0) + res (X, j + 2, 0) ) / 2.;
                }

                //  edge j=Y
                for (uint i = 0; i < res.getdimx() - 2; i += 2)
                {
                    res (i + 1, Y, 0) = ( res (i, Y, 0) + res (i + 2, Y, 0) ) / 2.;
                }

            }// end 2d case

            resold = res;

        }// end for on r = 2...ratio


    }// end increasing resolution
    return;
}



template <typename T>
template <typename S>
void im3d::image3d<T>::grad (std::vector<image3d<S> >& res) const
{

    // CASE OF 3D IMAGES
    if (dimz > 1)
    {

        if (res.size() != 3)
        {
            res.resize (3);
        }
        for (uint i = 0; i < 3; ++i)
        {
            res[i].setdim (this->dimx, this->dimy, this->dimz);
            res[i].seth (this->hx, this->hy, this->hz);
        }

        #pragma omp parallel
        {
            // starting parallel section

            #pragma omp for
            for (uint i = 1; i < this->dimx - 1; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    for (uint k = 0; k < this->dimz; ++k)
                        res[0] (i, j, k) =
                            (static_cast<S> ( (*this) (i + 1, j, k) ) -
                             static_cast<S> ( (*this) (i - 1, j, k) ) ) / (2.*hx) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 1; j < this->dimy - 1; ++j)
                    for (uint k = 0; k < this->dimz; ++k)
                        res[1] (i, j, k) =
                            (static_cast<S> ( (*this) (i, j + 1, k) ) -
                             static_cast<S> ( (*this) (i, j - 1, k) ) ) / (2.*hy) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    for (uint k = 1; k < this->dimz - 1; ++k)
                        res[2] (i, j, k) =
                            ( static_cast<S> ( (*this) (i, j, k + 1) ) -
                              static_cast<S> ( (*this) (i, j, k - 1) ) ) / (2.*hz) ;

        }   // ending parallel section

        // boundary values

        for (uint j = 0; j < this->dimy; ++j)
            for (uint k = 0; k < this->dimz; ++k)
            {
                res[0] (0, j, k) = ( 4 * static_cast<S> ( (*this) (1, j, k) ) -
                                     3 * static_cast<S> ( (*this) (0, j, k) ) -
                                     static_cast<S> ( (*this) (2, j, k) ) ) / (2.*hx) ;
                res[0] (this->dimx - 1, j, k) = (3 * static_cast<S> ( (*this) (dimx - 1, j, k) ) -
                                                 4 * static_cast<S> ( (*this) (dimx - 2, j, k) ) +
                                                 static_cast<S> ( (*this) (dimx - 3, j, k) ) ) / (2.*hx);
            }

        for (uint i = 0; i < this->dimx; ++i)
            for (uint k = 0; k < this->dimz; ++k)
            {
                res[1] (i, 0, k) = ( 4 * static_cast<S> ( (*this) (i, 1, k) ) -
                                     3 * static_cast<S> ( (*this) (i, 0, k) ) -
                                     static_cast<S> ( (*this) (i, 2, k) ) ) / (2.*hy) ;
                res[1] (i, this->dimy - 1, k) = (3 * static_cast<S> ( (*this) (i, dimy - 1, k) ) -
                                                 4 * static_cast<S> ( (*this) (i, dimy - 2, k) ) +
                                                 static_cast<S> ( (*this) (i, dimy - 3, k) ) ) / (2.*hy);
            }

        for (uint i = 0; i < this->dimx; ++i)
            for (uint j = 0; j < this->dimy; ++j)
            {
                res[2] (i, j, 0) = ( 4 * static_cast<S> ( (*this) (i, j, 1) ) -
                                     3 * static_cast<S> ( (*this) (i, j, 0) ) -
                                     static_cast<S> ( (*this) (i, j, 2) ) ) / (2.*hz) ;
                res[2] (i, j, this->dimz - 1) = (3 * static_cast<S> ( (*this) (i, j, dimz - 1) ) -
                                                 4 * static_cast<S> ( (*this) (i, j, dimz - 2) ) +
                                                 static_cast<S> ( (*this) (i, j, dimz - 3) ) ) / (2.*hz);
            }
    } // end if(dimz>1)

    // CASE OF 2D IMAGES
    else
    {

        if (res.size() != 2)
        {
            res.resize (2);
        }
        for (uint i = 0; i < 2; ++i)
        {
            res[i].setdim (this->dimx, this->dimy, this->dimz);
            res[i].seth (this->hx, this->hy, this->hz);
        }

        #pragma omp parallel
        {
            // starting parallel section

            #pragma omp for
            for (uint i = 1; i < this->dimx - 1; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                    res[0] (i, j, 0) =
                        ( static_cast<S> ( (*this) (i + 1, j, 0) ) -
                          static_cast<S> ( (*this) (i - 1, j, 0) ) ) / (2.*hx) ;

            #pragma omp for
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 1; j < this->dimy - 1; ++j)
                    res[1] (i, j, 0) =
                        ( static_cast<S> ( (*this) (i, j + 1, 0) ) -
                          static_cast<S> ( (*this) (i, j - 1, 0) ) ) / (2.*hy) ;

        }   // ending parallel section


        // boundary values

        for (uint j = 0; j < this->dimy; ++j)
        {
            res[0] (0, j, 0) = ( 4 * static_cast<S> ( (*this) (1, j, 0) ) -
                                 3 * static_cast<S> ( (*this) (0, j, 0) ) -
                                 static_cast<S> ( (*this) (2, j, 0) ) ) / (2.*hx) ;
            res[0] (this->dimx - 1, j, 0) = ( 3 * static_cast<S> ( (*this) (dimx - 1, j, 0) ) -
                                              4 * static_cast<S> ( (*this) (dimx - 2, j, 0) ) +
                                              static_cast<S> ( (*this) (dimx - 3, j, 0 ) ) ) / (2.*hx);
        }

        for (uint i = 0; i < this->dimx; ++i)
        {
            res[1] (i, 0, 0) = ( 4 * static_cast<S> ( (*this) (i, 1, 0) ) -
                                 3 * static_cast<S> ( (*this) (i, 0, 0) ) -
                                 static_cast<S> ( (*this) (i, 2, 0) ) ) / (2.*hy) ;
            res[1] (i, this->dimy - 1, 0) = ( 3 * static_cast<S> ( (*this) (i, dimy - 1, 0) ) -
                                              4 * static_cast<S> ( (*this) (i, dimy - 2, 0) ) +
                                              static_cast<S> ( (*this) (i, dimy - 3, 0) ) ) / (2.*hy);
        }

    }// end else ( CASE OF 2D IMAGES)

    return;
}



template <typename S, typename R>
void im3d::vector_abs (image3d<S>& res, std::vector<image3d<R> > const& fun)
{
    // check dimensions of input vector
    if (fun.size() != 2 && fun.size() != 3)
    {
        std::cout << "In image3d::vector_abs: input vector could be only 2d or 3d, ";
        std::cout << fun.size() << " is a not allowed dimension." << std::endl;
        return;
    }

    res.setdim (fun[0].dimx, fun[0].dimy, fun[0].dimz );
    res.seth (fun[0].hx, fun[0].hy, fun[0].hz );

    // 3d case
    if (fun.size() == 3 &&
            fun[0].dimx == fun[1].dimx && fun[1].dimx == fun[2].dimx  &&
            fun[0].dimy == fun[1].dimy && fun[1].dimy == fun[2].dimy  &&
            fun[0].dimz == fun[1].dimz && fun[1].dimz == fun[2].dimz )
    {
        #pragma omp parallel for
        for (uint i = 0; i < res.dimx * res.dimy * res.dimz; ++i)
            res.rawimage[i] =
                static_cast<S> (sqrt ( fun[0].rawimage[i] * fun[0].rawimage[i] +
                                       fun[1].rawimage[i] * fun[1].rawimage[i] +
                                       fun[2].rawimage[i] * fun[2].rawimage[i] ) );
        //note: sqrt works only with float, double or long double
    }

    // 2d case
    else if (fun.size() == 2 &&
             fun[0].dimx == fun[1].dimx &&
             fun[0].dimy == fun[1].dimy &&
             fun[0].dimz == fun[1].dimz )
    {
        #pragma omp parallel for
        for (uint i = 0; i < res.dimx * res.dimy * res.dimz; ++i)
            res.rawimage[i] =
                static_cast<S> (sqrt ( fun[0].rawimage[i] * fun[0].rawimage[i] +
                                       fun[1].rawimage[i] * fun[1].rawimage[i] ) );
        //note: sqrt works only with float, double or long double
    }

    else
    {
        std::cout << "In image3d::vector_abs: dimensions of image3d " <<
                  "inside input vector must agree" << std::endl;
        return;
    }

    return;
}



template <typename S, typename R>
void im3d::div (image3d<S>& res, std::vector<image3d<R> > const& fun)
{

    // check dimensions of input vector
    if (fun.size() != 2 && fun.size() != 3)
    {
        std::cout << "In image3d::div: input vector could be only 2d or 3d, ";
        std::cout << fun.size() << " is a not allowed dimension." << std::endl;
        return;
    }

    res.setdim (fun[0].getdimx(), fun[0].getdimy(), fun[0].getdimz() );
    res.seth (fun[0].gethx(), fun[0].gethy(), fun[0].gethz() );

    double hx = res.gethx(), hy = res.gethy(), hz = res.gethz();
    uint X = res.getdimx(), Y = res.getdimy(), Z = res.getdimz();

    // 2d case
    if (fun.size() == 2 &&
            fun[0].getdimx() == fun[1].getdimx() &&
            fun[0].getdimy() == fun[1].getdimy() &&
            fun[0].getdimz() == fun[1].getdimz() )
    {

        #pragma omp parallel for
        for (uint i = 1; i < X - 1; ++i)
            for (uint j = 1; j < Y - 1; ++j)
                res (i, j, 0) =
                    ( (static_cast<S> (fun[0] (i + 1, j, 0) ) -
                       static_cast<S> (fun[0] (i - 1, j, 0) ) ) / (2.*hx) +
                      (static_cast<S> (fun[1] (i, j + 1, 0) ) -
                       static_cast<S> (fun[1] (i, j - 1, 0) ) ) / (2.*hy) );

        //edges i-variabili
        for (uint i = 1; i < X - 1; ++i)
        {
            //j=0, k=0
            res (i, 0, 0) =
                ( (static_cast<S> (fun[0] (i + 1, 0, 0) ) -
                   static_cast<S> (fun[0] (i - 1, 0, 0) ) ) / (2.*hx) +
                  (4 * static_cast<S> (fun[1] (i, 1, 0) ) -
                   3 * static_cast<S> (fun[1] (i, 0, 0) ) -
                   static_cast<S> (fun[1] (i, 2, 0) ) ) / (2.*hy) );

            //j=Y-1, k=0
            res (i, Y - 1, 0) =
                ( (static_cast<S> (fun[0] (i + 1, Y - 1, 0) ) -
                   static_cast<S> (fun[0] (i - 1, Y - 1, 0) ) ) / (2.*hx) +
                  (3 * static_cast<S> (fun[1] (i, Y - 1, 0) ) -
                   4 * static_cast<S> (fun[1] (i, Y - 2, 0) ) +
                   static_cast<S> (fun[1] (i, Y - 3, 0) ) ) / (2.*hy) );
        }

        //edges j variabili
        for (uint j = 1; j < Y - 1; ++j)
        {
            //i=0, k=0
            res (0, j, 0) =
                ( (4 * static_cast<S> (fun[0] (1, j, 0) ) -
                   3 * static_cast<S> (fun[0] (0, j, 0) ) -
                   static_cast<S> (fun[0] (2, j, 0) ) ) / (2.*hx) +
                  (static_cast<S> (fun[1] (0, j + 1, 0) ) -
                   static_cast<S> (fun[1] (0, j - 1, 0) ) ) / (2.*hy) );

            //i=X-1, k=0
            res (X - 1, j, 0) =
                ( (3 * static_cast<S> (fun[0] (X - 1, j, 0) ) -
                   4 * static_cast<S> (fun[0] (X - 2, j, 0) ) +
                   static_cast<S> (fun[0] (X - 3, j, 0) ) ) / (2.*hx) +
      