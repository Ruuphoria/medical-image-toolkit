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
        // we have not computed them yet we start calculating the result from the pixel
        // (0,0,0)

        // corner i=0, j=0
        res (0, 0, 0) = ( ( res (1, 0, 0) - nxy * bc * hx ) / (hx * hx) +
                          ( res (0, 1, 0) - nxy * bc * hy ) / (hy * hy) +
                          b (0, 0, 0) + res (0, 0, 0) / dt ) * htildeij;

        // edge j=0 , i>0
        for (uint i = 1; i < X - 2; ++i)
            res (i, 0, 0) = ( (res (i + 1, 0, 0) + res (i - 1, 0, 0) ) / (hx * hx) +
                              (res (i, 1, 0) - bc * hy ) / (hy * hy) +
                              b (i, 0, 0) + res (i, 0, 0) / dt  ) * htildej;

        // corner i=X-1 j=0
        res (X - 1, 0, 0) = ( ( res (X - 2, 0, 0) + nxy * bc * hx ) / (hx * hx) +
                              ( res (X - 1, 1, 0) - nxy * bc * hy ) / (hy * hy) +
                              b (X - 1, 0, 0) + res (X - 1, 0, 0) / dt ) * htildeij;

        // INTERNAL NODES
        for (uint j = 1; j < Y - 2; ++j)
        {
            // edge i=0 j>0
            res (0, j, 0) =
                ( ( res (1, j, 0) - bc * hx ) / (hx * hx) + (res (0, j + 1, 0) + res (0, j - 1, 0) ) * (hy * hy)
                  + b (0, j, 0) + res (0, j, 0) / dt   ) * htildei;

            // internal nodes i>0, j>0
            for (uint i = 1; i < X - 1; ++i)
                res (i, j, 0) =
                    ( (res (i + 1, j, 0) + res (i - 1, j, 0) ) / (hx * hx) + (res (i, j + 1, 0) + res (i, j - 1, 0) ) /
                      (hy * hy) + b (i, j, 0) + res (i, j, 0) / dt  ) * htilde;

            // edge i=max, j>0
            res (X - 1, j, 0) =
                ( ( res (X - 2, j, 0) + bc * hx ) / (hx * hx) + ( res (X - 1, j + 1, 0) + res (X - 1, j - 1, 0) ) /
                  (hy * hy) + b (X - 1, j, 0) + res (X - 1, j, 0) / dt ) * htildei;
        }
        // corner i=0, j=Y-1
        res (0, Y - 1, 0) =
            ( ( res (1, Y - 1, 0) - nxy * bc * hx ) / (hx * hx) + ( res (0, Y - 2, 0) + nxy * bc * hy ) /
              (hy * hy) + b (0, Y - 1, 0) + res (0, Y - 1, 0) / dt ) * htildeij;

        //edge i>0, j=Y-1
        for (uint i = 1; i < X - 2; ++i)
            res (i, Y - 1, 0) =
                ( (res (i + 1, Y - 1, 0) + res (i - 1, Y - 1, 0) ) * (hx * hx) + ( res (i, Y - 2, 0) + bc * hy ) /
                  (hy * hy) + b (i, Y - 1, 0) + res (i, Y - 1, 0) / dt ) * htildej;

        // corner i=X-1, j=Y-1
        res (X - 1, Y - 1, 0) =
            ( ( res (X - 2, Y - 1, 0) + nxy * bc * hx ) / (hx * hx) + ( res (X - 1, Y - 2, 0) + nxy * bc * hy ) /
              (hy * hy) + b (X - 1, Y - 1, 0) + res (X - 1, Y - 1, 0) / dt ) * htildeij;

    }
    else
    {
        double htilde = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz);
        htilde = 1. / htilde;

        double htildei = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz);
        htildei = 1. / htildei;

        double htildej = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy) + 2. / (hz * hz);
        htildej = 1. / htildej;

        double htildek = 1. / dt + 2. / (hx * hx) + 2. / (hy * hy) + 1. / (hz * hz);
        htildek = 1. / htildek;

        double htildejk = 1. / dt + 2. / (hx * hx) + 1. / (hy * hy) + 1. / (hz * hz);
        htildejk = 1. / htildejk;

        double htildeik = 1. / dt + 1. / (hx * hx) + 2. / (hy * hy) + 1. / (hz * hz);
        htildeik = 1. / htildeik;

        double htildeij = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy) + 2. / (hz * hz);
        htildeij = 1. / htildeij;

        double htildeijk = 1. / dt + 1. / (hx * hx) + 1. / (hy * hy) + 1. / (hz * hz);
        htildeijk = 1. / htildeijk;

        // the computation must be done in the correct order, in res there is the previous
        // step result Gauss-Seidel method is semi-implicit, it uses both curent step
        // value and previous step values we read the previous step value in a pixel,
        // we compute the new value and write it in the same pixel the computation
        // involves pixels near to the one on which we write the result
        // (the central pixel) so the pixels with indexes smaller than the indexes of the
        // central pixel will be treated as current step values because at current step we
        // have already computed them the pixels with indexes greater than the indexes of
        // the central pixel will be treated as previous step values because we have not
        // computed them yet

        // OBLIQUE NORMAL COMPUTATION
        // edge  [(1,1,0) (1,0,1) (0,1,1)]*(1/sqrt(2))
        // angle (1,1,1)*(1/sqrt(3))
        T n_edge (1.), n_angle (1.);
        n_edge /= sqrt (2.);
        n_angle /= sqrt (3.);

        // COMPUTE FACE k=0
        // corner i=0, j=0, k=0
        res (0, 0, 0) = ( (res (1, 0, 0) - n_angle * bc * hx) / (hx * hx) +
                          (res (0, 1, 0) - n_angle * bc * hy) / (hy * hy) +
                          (res (0, 0, 1) - n_angle * bc * hz) / (hz * hz) +
                          b (0, 0, 0) + res (0, 0, 0) / dt ) * htildeijk;
        // edge i>0, j=0, k=0
        for (uint i = 1; i < X - 1; ++i)
            res (i, 0, 0) =  ( (res (i + 