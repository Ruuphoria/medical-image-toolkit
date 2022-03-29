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
  \file poisson_imp.hxx

 \brief file with the implementation of all the methods and functions of namespace \ref lapl
*/

#ifndef POISSON_IMP_HXX_INCLUDED
#define POISSON_IMP_HXX_INCLUDED



template <typename T>
void lapl::NeuGaussSeidel (im3d::image3d<T>& res, im3d::image3d<T> const& b,
                           T const& dt, T const& bc, im3d::image3d<T> const& input)
{

    // we impose external normal derivative equal to bc value on the border

    uint const X = res.getdimx(), Y = res.getdimy(), Z = res.getdimz();
    double const hx = res.gethx(), hy = res.gethy(), hz = res.gethz();

    if (Z == 1)
    {
        double htilde = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy);
        htilde = 1. / htilde;

        double htildei = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy);
        htildei = 1. / htildei;

        double htildej = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy);
        htildej = 1. / htildej;

        double htildeij = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy);
        htildeij = 1. / htildeij;

        // COMPUTE THE OBLIQUE NORMAL COMPONENTS
        T nxy (1.);
        nxy /= sqrt (2.);

        // the computation must be done in the correct order, in res there is the previous
        // step result Gauss-Seidel method is semi-implicit, it uses both current step
        // value and previous step values we read the previous step value in a pixel,
        // we compute the new value and write it in the same pixel the computation
        // involves pixels near to the one on which we write the result
        // (the central pixel) so the pixels with indexes smaller than the indexes of
        // the central pixel will be treated as current step values because at current
        // step we have already computed them the pixels with indexes greater than the
        // indexes of the central pixel will be treated as previous step values because
        // we have not computed t