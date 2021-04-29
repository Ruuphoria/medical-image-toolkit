
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
  \file convolution_imp.hxx

 \brief file with the implementation of all the methods and functions of namespace \ref conv
*/

#ifndef CONVOLUTION_IMP_HXX_INCLUDED
#define CONVOLUTION_IMP_HXX_INCLUDED


template <typename T>
void conv::filtering<T>::build_gaussian_filter (im3d::image3d<T> const& f, T const& sigma,
                                                double const& radius, bool pixel_approach)
{
    // get spacing from factor
    double hx (f.gethx() ), hy (f.gethy() ), hz (f.gethz() );

    if (pixel_approach)
    {
        hx = hy = hz = 1.;
    }

    int Nhx, Nhy, Nhz;  //time_filter dimensions

    // Set time_filter dimensions to make it about zero in the bounds.
    Nhx = ceil (2 * radius * sigma / hx);
    Nhy = ceil (2 * radius * sigma / hy);
    if (f.getdimz() == 1) // in case of bidimensional images
    {
        Nhz = 1;
    }
    else
    {
        Nhz = ceil (2 * radius * sigma / hz);
    }

    // Make them odds.
    if (Nhx % 2 == 0)
    {
        Nhx++;
    }
    if (Nhy % 2 == 0)
    {
        Nhy++;
    }
    if (Nhz % 2 == 0)
    {
        Nhz++;
    }

    // assigning dimensions to time_filter
    filtering<T>::time_filter.setdim (Nhx, Nhy, Nhz);
    filtering<T>::time_filter.seth (f.gethx(), f.gethy(), f.gethz() );


    //VALUES COMPUTATION

    // compute position of the maximal value of the gaussian filter in every direction
    int const shiftx = floor (static_cast<double> (Nhx) / 2.);
    int const shifty = floor (static_cast<double> (Nhy) / 2.);
    int const shiftz = floor (static_cast<double> (Nhz) / 2.);

    // compute time_filter
    #pragma omp parallel for
    for (int i = 0; i < Nhx; ++i)
        for (int j = 0; j < Nhy; ++j)
            for (int k = 0; k < Nhz; ++k)
                filtering<T>::time_filter (i, j, k) =
                    exp ( (- (i - shiftx) * (i - shiftx) * hx * hx
                           - (j - shifty) * (j - shifty) * hy * hy
                           - (k - shiftz) * (k - shiftz) * hz * hz) / (2 * sigma * sigma) );

    // make normL1 of time_filter equal to one
    filtering<T>::time_filter /= filtering<T>::time_filter.normL1();

    // compute corresponding cyclic and acyclic filters
    filtering<T>::build_freq_filters (f);

    return;
}



template <typename T>
void conv::filtering<T>::build_high_pass_filter (im3d::image3d<T> const& f, double const& cnc,
                                                 bool pixel_approach)
{
    if (cnc <= 0)
    {
        std::cout << "WARNING::build_high_pass_filter: central_node_coeff " <<
                  "must be positive" << std::endl;
        return;
    }

    // get spacing from factor
    double hx (f.gethx() ), hy (f.gethy() ), hz (f.gethz() );

    if (pixel_approach)
    {
        hx = hy = hz = 1.;
    }

    // 2d case
    if (f.getdimz() == 1)
    {
        // assigning dimensions an spacing to time_filter
        filtering<T>::time_filter.setdim (3, 3, 1);
        filtering<T>::time_filter.seth (hx, hy, hz);

        // compute time_filter
        filtering<T>::time_filter = 0.;
        filtering<T>::time_filter (0, 1, 0) = filtering<T>::time_filter (2, 1, 0) = -1. / (hy * hy);
        filtering<T>::time_filter (1, 0, 0) = filtering<T>::time_filter (1, 2, 0) = -1. / (hx * hx);
        filtering<T>::time_filter (1, 1, 0) = 2. / (hx * hx) + 2. / (hy * hy) + 2.*cnc / (hx * hx + hy * hy);

        // make sum of elements of time_filter equal to one
        double sum =
            2 * filtering<T>::time_filter (0, 1, 0) +
            2 * filtering<T>::time_filter (1, 0, 0) +
            filtering<T>::time_filter (1, 1, 0);

        filtering<T>::time_filter /= sum;
    }

    // 3d case
    else
    {
        // assigning dimensions to time_filter
        filtering<T>::time_filter.setdim (3, 3, 3);
        filtering<T>::time_filter.seth (hx, hy, hz);

        // compute time_filter
        filtering<T>::time_filter = 0.;
        filtering<T>::time_filter (0, 1, 1) = filtering<T>::time_filter (2, 1, 1) = -1. / (hy * hy);
        filtering<T>::time_filter (1, 0, 1) = filtering<T>::time_filter (1, 2, 1) = -1. / (hx * hx);
        filtering<T>::time_filter (1, 1, 0) = filtering<T>::time_filter (1, 1, 2) = -1. / (hz * hz);
        filtering<T>::time_filter (1, 1, 1) =
            2. / (hx * hx) + 2. / (hy * hy) + 2. / (hz * hz) + 3.*cnc / (hx * hx + hy * hy + hz * hz);

        // make sum of elements of time_filter equal to one
        double sum =
            2 * filtering<T>::time_filter (0, 1, 1) +
            2 * filtering<T>::time_filter (1, 0, 1) +
            2 * filtering<T>::time_filter (1, 1, 0) +
            filtering<T>::time_filter (1, 1, 1);

        filtering<T>::time_filter /= sum;
    }

    // compute corresponding cyclic and acyclic filters
    filtering<T>::build_freq_filters (f);

    return;
}



template <typename T>
void conv::filtering<T>::build_laplacian_filter (im3d::image3d<T> const& f,
                                                 bool pixel_approach)
{
    // get spacing from factor
    double hx (f.gethx() ), hy (f.gethy() ), hz (f.gethz() );

    if (pixel_approach)
    {
        hx = hy = hz = 1.;
    }

    // 2d case
    if (f.getdimz() == 1)
    {
        // assigning dimensions an spacing to time_filter
        filtering<T>::time_filter.setdim (3, 3, 1);
        filtering<T>::time_filter.seth (hx, hy, hz);

        // compute time_filter
        filtering<T>::time_filter = 0.;
        filtering<T>::time_filter (0, 1, 0) = filtering<T>::time_filter (2, 1, 0) = 1. / (hy * hy);
        filtering<T>::time_filter (1, 0, 0) = filtering<T>::time_filter (1, 2, 0) = 1. / (hx * hx);
        filtering<T>::time_filter (1, 1, 0) = -2. / (hx * hx) - 2. / (hy * hy);
    }

    // 3d case
    else
    {
        // assigning dimensions to time_filter
        filtering<T>::time_filter.setdim (3, 3, 3);
        filtering<T>::time_filter.seth (hx, hy, hz);

        // compute time_filter
        filtering<T>::time_filter = 0.;
        filtering<T>::time_filter (0, 1, 1) = filtering<T>::time_filter (2, 1, 1) = 1. / (hy * hy);
        filtering<T>::time_filter (1, 0, 1) = filtering<T>::time_filter (1, 2, 1) = 1. / (hx * hx);
        filtering<T>::time_filter (1, 1, 0) = filtering<T>::time_filter (1, 1, 2) = 1. / (hz * hz);
        filtering<T>::time_filter (1, 1, 1) = -2. / (hx * hx) - 2. / (hy * hy) - 2. / (hz * hz);
    }

    // compute corresponding cyclic and acyclic filters
    filtering<T>::build_freq_filters (f);

    return;
}



template <typename T>
void conv::filtering<T>::build_freq_filters (im3d::image3d<T> const& f)
{
    // allowing usage of fftw3_omp library to improve performances
    if (omp_get_max_threads() > 1)
    {
        int good = fftw_init_threads();
        if (good)
        {
            fftw_plan_with_nthreads (omp_get_max_threads() );
        }
        else
        {
            std::clog << "WARNING::build_freq_filters: thread creation error" << std::endl;
        }
    }

    // getting time_filter dimensions
    int const Nhx = filtering<T>::time_filter.getdimx();
    int const Nhy = filtering<T>::time_filter.getdimy();
    int const Nhz = filtering<T>::time_filter.getdimz();


    // SET FREQUENCY FILTERS DIMENSIONS

    // Set correct acyclic_freq_filter dimensions.
    int const Nx = f.getdimx() + Nhx - 1;
    int const Ny = f.getdimy() + Nhy - 1;
    int const Nz = f.getdimz() + Nhz - 1;

    // Set correct cyclic_freq_filter dimensions.
    int Ncx, Ncy, Ncz;
    Ncx = f.getdimx();
    Ncy = f.getdimy();
    Ncz = f.getdimz();
    while (Ncx < Nhx)
    {
        Ncx += f.getdimx();
    }
    while (Ncy < Nhy)
    {
        Ncy += f.getdimy();
    }
    while (Ncz < Nhz)
    {
        Ncz += f.getdimz();
    }

    // assigning dimensions to frequency filters
    filtering<T>::acycl_freq_filter.setdim (Nx, Ny, Nz);
    filtering<T>::cycl_freq_filter.setdim (Ncx, Ncy, Ncz);


    // VALUES COMPUTATION

    // compute position of the maximal value of the time_filter in every direction
    int const shiftx = floor (static_cast<double> (Nhx) / 2.);
    int const shifty = floor (static_cast<double> (Nhy) / 2.);
    int const shiftz = floor (static_cast<double> (Nhz) / 2.);

    // ***** acyclic_freq_filter *****

    fftw_complex* in, *out;
    in = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * (Nx * Ny * Nz) );
    out = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * Nx * Ny * Nz);

    // copying time_filter values in the correct way
    // to obtain a coherent Fourier transform using fftw

    #pragma omp parallel
    {
        // starting parallel section

        // initializing
        #pragma omp for
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                for (int k = 0; k < Nz; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] = 0.;
                    in[Ny * Nz * i + Nz * j + k][1] = 0.;
                }

        //blocco 1: x=0...X/2, y=0...Y/2, z=0...Z/2
        #pragma omp for
        for (int i = 0; i < shiftx + 1; ++i)
            for (int j = 0; j < shifty + 1; ++j)
                for (int k = 0; k < shiftz + 1; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (shiftx + i, shifty + j, shiftz + k);
                }

        //blocco 2: x=0...X/2, y=0...Y/2, z=Z/2...Z
        #pragma omp for
        for (int i = 0; i < shiftx + 1; ++i)
            for (int j = 0; j < shifty + 1; ++j)
                for (int k = Nz - shiftz; k < Nz ; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (shiftx + i, shifty + j, k - Nz + shiftz);
                }

        //blocco 3: x=0...X/2, y=Y/2...Y, z=0...Z/2
        #pragma omp for
        for (int i = 0; i < shiftx + 1; ++i)
            for (int j = Ny - shifty; j < Ny ; ++j)
                for (int k = 0; k < shiftz + 1; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (shiftx + i, j - Ny + shifty, shiftz + k);
                }

        //blocco 4: x=X/2...X, y=0...Y/2, z=0...Z/2
        #pragma omp for
        for (int i = Nx - shiftx; i < Nx; ++i)
            for (int j = 0; j < shifty + 1; ++j)
                for (int k = 0; k < shiftz + 1; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (i - Nx + shiftx, shifty + j, shiftz + k);
                }

        //blocco 5: x=0...X/2, y=Y/2...Y, z=Z/2...Z
        #pragma omp for
        for (int i = 0; i < shiftx + 1; ++i)
            for (int j = Ny - shifty; j < Ny; ++j)
                for (int k = Nz - shiftz; k < Nz; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (shiftx + i, j - Ny + shifty, k - Nz + shiftz);
                }

        //blocco 6: x=X/2...X, y=0...Y/2, z=Z/2...Z
        #pragma omp for
        for (int i = Nx - shiftx; i < Nx; ++i)
            for (int j = 0; j < shifty + 1; ++j)
                for (int k = Nz - shiftz; k < Nz; ++k)
                {
                    in[Ny * Nz * i + Nz * j + k][0] =
                        filtering<T>::time_filter (i - Nx + shiftx, shifty + j, k - Nz + shiftz);
                }

        //blocco 7: x=X/2...X, y=Y/2...Y, z=0...Z/2