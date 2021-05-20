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
im3d::image