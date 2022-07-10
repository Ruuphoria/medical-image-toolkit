
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
  \file  rsfe_splitbregman_imp.hxx

  \brief file with the implementation of all the methods of class \ref segm::rsfe_splitbregman
*/


#ifndef SEGMENTATION_IMP_HXX_INCLUDED
#define SEGMENTATION_IMP_HXX_INCLUDED




template <typename T>
segm::rsfe_splitbregman<T>::rsfe_splitbregman() :
    segmentation<T>(), getpotfile ("data"), verbosity (true), maxiter (100),
    showfrequency (0), dumpfrequency (0),  outputname ("result"), logfilename (""),
    onthego (true), save_current (false), end_now (false),
    sigma (3), lamda (5.e-2), lamda1 (1.e-5), lamda2 (1.e-5),
    beta (100.), dt (1), ls_steps (5), epsilon (1.e-3), a0 (0.), b0 (1.),
    auto_extract_conn_comp (false), cc_init_pixel_x (-1), cc_init_pixel_y (-1),
    cc_init_pixel_z (-1), cc_binary_output (true), cc_init_variables (true), auto_extract_conn_comp_freq (0), bc (NEU),
    ls_tol (0), edge_detector (GRAD), edge_detector_sigma (0), nu (1.),
    onestep_poisson (lapl::DirGaussSeidel), cv (conv::acyclic_fftw_convolution),
    gaussian_pixel_approach (false), tol (1.e-5),
    ix (-1), iy (-1), iz (-1), fx (-1), fy (-1), fz (-1), dimx (0), dimy (0), dimz (1),
    space_dim (3), hx (1), hy (1), hz (1), init_variables (false)
{
    this->alpha = ( this->a0 + this->b0 ) / 2.;
}




// PRIVATE MEMBER IMPLEMENTATION

template <typename T>
void segm::rsfe_splitbregman<T>::initialize_phi_with_cube ()
{
    this->phi.setdim (this->dimx, this->dimy, this->dimz);

    if (this->ix == -1)
    {
        double frac = 0.33;

        this->ix = this->dimx * frac;
        this->iy = this->dimy * frac;
        this->iz = this->dimz * frac;
        this->fx = this->dimx * (1 - frac);
        this->fy = this->dimy * (1 - frac);
        this->fz = this->dimz * (1 - frac);
    }

    // initializing phi with the external value a0
    this->phi = this->a0;


    // modifying values of internal nodes

    // 2d case
    if (this->dimz == 1)
    {
        #pragma omp parallel for
        for (int i = this->ix; i < this->fx; ++i)
            for (int j = this->iy; j < this->fy; ++j)
            {
                this->phi (i, j, 0) = this->b0;
            }
    }
    // 3d case
    else
    {

        #pragma omp parallel for
        for (int i = this->ix; i < this->fx; ++i)
            for (int j = this->iy; j < this->fy; ++j)
                for (int k = this->iz; k < this->fz; ++k)
                {
                    this->phi (i, j, k) = this->b0;
                }
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::initialize_phi_with_init (im3d::image3d<T> const& init)
{
    uint X (init.getdimx() ), Y (init.getdimy() ), Z (init.getdimz() );

    // 1st CASE: init has the same dimensions of the image to segment
    // just assign it to phi
    if ( X == this->dimx && Y == this->dimy && Z == this->dimz )
    {
        this->phi = init;
        this->phi.change_range_of_intensity (this->b0, this->a0);
    }

    // 2nd CASE: in at least one direction init is a pixel smaller
    // assign init to phi and copy results of last but one faces on the last ones
    // (this allows to pass algorithm an interpolation of a result with a lower resolution)
    else if ( (X == this->dimx || X == this->dimx - 1) &&
              (Y == this->dimy || Y == this->dimy - 1) &&
              (Z == this->dimz || Z == this->dimz - 1) )
    {
        this->phi.setdim (this->dimx, this->dimy, this->dimz);

        #pragma omp parallel for
        for (uint i = 0; i < X; ++i)
            for (uint j = 0; j < Y; ++j)
                for (uint k = 0; k < Z; ++k)
                {
                    this->phi (i, j, k) = init (i, j, k);
                }

        // copying result of last but one faces on last threes in case of different
        // dimensions between phi and init
        if (this->dimx != X)
        {
            for (uint j = 0; j < this->dimy; ++j)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    this->phi (X, j, k) = this->phi (X - 1, j, k);
                }
        }

        if (this->dimy != Y)
        {
            for (uint i = 0; i < this->dimx; ++i)
                for (uint k = 0; k < this->dimz; ++k)
                {
                    this->phi (i, Y, k) = this->phi (i, Y - 1, k);
                }
        }

        if (this->dimz != Z)
        {
            for (uint i = 0; i < this->dimx; ++i)
                for (uint j = 0; j < this->dimy; ++j)
                {
                    this->phi (i, j, Z) = this->phi (i, j, Z - 1);
                }
        }

        this->phi.change_range_of_intensity (this->b0, this->a0);
    }

    // 3rd CASE: init has incompatible dimensions to initialize phi, cube is used
    else
    {
        this->initialize_phi_with_cube();
        std::clog << "WARNING::apply: wrong dimensions of init, cube initial " <<
                  "contour is used instead." << std::endl;
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::set_heaviside (im3d::image3d<T>& Heps) const
{
    Heps.setdim (dimx, dimy, dimz);
    Heps.seth (hx, hy, hz);

    T center = (this->a0 + this->b0) / 2.;

    #pragma omp parallel for
    for (uint i = 0; i < dimx; ++i)
        for (uint j = 0; j < dimy; ++j)
            for (uint k = 0; k < dimz; ++k)
                Heps (i, j, k) =
                    0.5 * (1. + 2. / M_PI *
                           atan ( (this->phi (i, j, k) - center)
                                  / (this->epsilon/**std::abs(this->b0-this->a0)*/) ) );

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::edge_detect (im3d::image3d<T>& res, im3d::image3d<T> const& f) const
{
    im3d::image3d<T> aux (f);
    char choice = 'n';
    im3d::interface<T> shower;

    res.setdim (this->dimx, this->dimy, this->dimz);
    res.seth (this->hx, this->hy, this->hz);

    if (this->verbosity == 1)
    {
        std::clog << "- Setting Edge Detector ..." << std::endl;
    }

    if (this->edge_detector == GRAD || this->edge_detector == GRAD_THRESHOLD || this->dimz != 1)
    {
        std::vector<im3d::image3d<T> > aux_vec;

        if (this->edge_detector_sigma != 0)
        {
            // Build filters of the correct dimensions
            conv::filtering<T>::build_gaussian_filter (aux,
                                                       this->edge_detector_sigma,
                                                       3,
                                                       this->gaussian_pixel_approach);
            cv (aux, f);

            conv::filtering<T>::build_gaussian_filter (aux,
                                                       this->sigma,
                                                       3,
                                                       this->gaussian_pixel_approach);
        }

        aux.grad (aux_vec);

        // case GRAD
        if (this->edge_detector == GRAD)
        {
            // aux = (Modulus of aux_vec)^2
            aux = aux_vec[0] * aux_vec[0];
            aux += aux_vec[1] * aux_vec[1];

            if (this->space_dim == 3)
            {
                aux += aux_vec[2] * aux_vec[2];
            }
        }

        // case GRAD_THRESHOLD
        else
        {
            // aux = Modulus of aux_vec
            im3d::vector_abs (res, aux_vec);

            shower.convertfromimage3d (res);
            std::cout << "showing modulus of gradient to help you to set " <<
                      "the correct threshold ..." << std::endl;
            shower.show_image();

            double t;

            while (choice == 'n')
            {
                std::cout << "Choose a threshold for modulus of gradient:" << std::endl;
                std::cout << "\tthreshold: ";
                cin >> t;
                res.im_to_black_and_white (aux, t);

                shower.convertfromimage3d (aux);
                std::cout << "showing output of chosen edge detector ..." << std::endl;
                shower.show_image();

                std::cout << "Do you accept this edge detector ";
                std::cout << "(if not you can recompute it)? [y/n]\t" << std::endl;
                cin >> choice;
            }
        }

    }// end GRAD/GRAD_THRESHOLD edge detector


    else if (this->edge_detector == LBP)
    {
        T constant;
        double t1, t2;

        while (choice == 'n')
        {
            std::cout << "Choose parameters of the local binary pattern" << std::endl;
            std::cout << "\tconstant T: ";
            cin >> constant;
            std::cout << "\tthreshold t1: ";
            cin >> t1;
            std::cout << "\tthreshold t2 (>t1): ";
            cin >> t2;
            // aux = Local Bynary Pattern Edge Detector Output
            f.local_binary_pattern_edge_detector (aux, constant, t1, t2);

            shower.convertfromimage3d (aux);
            std::cout << "showing output of chosen edge detector ..." << std::endl;

            shower.show_image();

            std::cout << "Do you accept this edge detector ";
            std::cout << "(if not you can recompute it)? [y/n]\t" << std::endl;
            cin >> choice;
        }

    }// end LBP edge detector

    // Build gforshrink
    aux *= this->beta;
    aux += static_cast<T> (1.);
    res = static_cast<T> (1.);
    res /= aux;
    res /= this->lamda;
    res *= this->nu;

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::cut_phi ()
{
    for (uint k = 0; k < this->dimz; ++k)
        for (uint j = 0; j < this->dimy; ++j)
            for (uint i = 0; i < this->dimx; ++i)
            {
                if (this->phi (i, j, k) < this->a0)
                {
                    this->phi (i, j, k) = this->a0;
                }

                else if (this->phi (i, j, k) > this->b0)
                {
                    this->phi (i, j, k) = this->b0;
                }
            }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::shrink
( std::vector<im3d::image3d<T> >& res, std::vector<im3d::image3d<T> > const& f1, im3d::image3d<T> const& f2) const
{

    for (uint y = 0; y < this->space_dim; ++y)
        #pragma omp parallel for
        for (uint i = 0; i < dimx; ++i)
            for (uint j = 0; j < dimy; ++j)
                for (uint k = 0; k < dimz; ++k)
                {

                    //try{
                    res[y] (i, j, k) =
                        static_cast<T> ( static_cast<int> (f1[y] (i, j, k) > 0.) - static_cast<int> (f1[y] (i, j, k) < 0.) ) * // sign(f1[y])
                        max ( std::abs (f1[y] (i, j, k) ) - f2 (i, j, k), static_cast<T> (0.) ) ;
                    //test_fpe_exception();
                    //}
                    //catch (BadFPOper & x){std::clog << x.what() << std::endl;}
                    //catch (BadDivision & x){std::clog << x.what() << std::endl;}
                }
    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::internal_show ()
{
    if (this->space_dim == 2)
    {
        this->levelset.show_contour_with_background_image (this->image, this->alpha);
        //this->levelset.show_contour(this->alpha);
        this->levelset.show_image_and_contour (this->alpha);
    }
    else
    {
        this->levelset.show_contour_with_background_image (this->image, this->alpha);
        //this->levelset.show_contour(this->alpha);
        //this->levelset.show_image();
        this->levelset.show_image_and_contour (this->alpha);
    }

    return;
}



template <typename T>
void segm::rsfe_splitbregman<T>::extract_connected_component ()
{

    im3d::image3d<T> connected_component;

    // first case: manual extraction with user interaction
    if (!auto_extract_conn_comp)
    {
        char choice;

        std::cout << "Do you want to see a connected component of the level set? (y/n): ";
        cin >> choice;

        if (choice == 'y')
        {
            std::cout << "Choose pixel of desired connected component: " << std::endl;

            uint i, j, k;

            this->phi.im_to_black_and_white (connected_component, (this->alpha - this->a0) /
                                             (this->b0 - this->a0) );

            im3d::interface<T> helper (connected_component);

            helper.get_coordinates (i, j, k);

            this->phi.connected_component (connected_component, i, j, k);

            helper.convertfromimage3d (connected_component);
            helper.show_image_and_contour (0.5);

            std::cout << "Do you want to set this connected component " <<
                      "as current level set? (y/n): ";
            cin >> choice;

            if (choice == 'y')
            {
                this->phi = connected_component;
                this->phi.change_range_of_intensity (this->b0, this->a0);
                if (verbosity)
                {
                    std::clog << "-- Updating current levelset to a chosen connected component" << std::endl;
                }

                std::cout << "Do you want to initialize all other variables of the algorithm? (y/n): ";
                cin >> choice;
                if (choice == 'y')
                {
                    this->init_variables = true;
                }

            }
        }

    }

    // second case: automatic extraction using private pixel
    else
    {
        int i = floor ( this->cc_init_pixel_x / this->hx );

        int j = floor ( this->cc_init_pixel_y / this->hy );

        int k = floor ( this->cc_init_pixel_z / this->hz );

        if ( i < 0 || j < 0 || (this->dimz != 1 && k < 0) )
        {
            std::clog << "WARNING::apply::extract_connected_component: " <<
                      "coordinates outside allowed range";
            std::clog << std::endl;
            return;
        }

        uint I, J, K;
        I = static_cast<uint> (i);
        J = static_cast<uint> (j);
        K = static_cast<uint> (k);

        if ( I > this->dimx - 1 || J > this->dimy - 1 || (this->dimz != 1 && K > this->dimz - 1) )
        {
            std::clog << "WARNING::apply::extract_connected_component: " <<
                      "coordinates outside allowed range";
            std::clog  << std::endl;
            return;
        }

        this->phi.im_to_black_and_white (connected_component, (this->alpha - this->a0) / (this->b0 - this->a0) );

        // continue only if private pixel is a white pixel after b&w conversion
        if ( connected_component (I, J, K) == 0 )
        {
            return;
        }

        this->phi.connected_component (connected_component, I, J, K, 0.5, false, this->cc_binary_output);

        // first subcase: automatic setting of connected_component as current levelset
        if ( this->auto_extract_conn_comp_freq != 0 )