
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
 \file image3d.hxx

 \brief header file of class \ref im3d::image3d
 */

#ifndef IMAGE3D_HPP_INCLUDED
#define IMAGE3D_HPP_INCLUDED

#include "names.hxx"

#include <algorithm>

namespace im3d
{

/*!
 \brief A template class conceived as a container of data and typical operations
 of a 3D image.

 This class is a container for data types that can be reassembled to a 3D image
 (mainly in size), making it possible to manipulate them, for instance in an
 algorithm.
 In particular it overloads arithmetic operators and implements various methods
 to access the data in both reading and writing mode.
 The main usage of this class is to represent a three dimensional image or an
 associate function defined on its uniform structured grid of pixels.
 It could be also used to represent a more general 3D signals, for instance a 3d audio
 signal or a video using the third dimension as time.
 It works good also with 2D images or signals if \ref image3d::dimz is set to 1.

 \tparam T is the type of data contained in the class, it can be any kind of number
 variable (int, float, double...) to represent black and white images.

 \warning This class is not conceived to represent RGB images. To make it work with
 RGB images user has to implement an appropriate son class reimplementing most of
 members.
 */
template <typename T = double>
class image3d
{

protected:
    //! a template vector containing the data of type T of the image.
    std::vector<T> rawimage;
    //! number of pixels in X direction.
    uint dimx;
    //! number of pixels in Y direction.
    uint dimy;
    //! number of pixels in Z direction.
    uint dimz;
    //! spacing between every pixel in X direction.
    double hx;
    //! spacing between every pixel in Y direction.
    double hy;
    //! spacing between every pixel in Z direction.
    double hz;

    /*!
     \brief private member used by both version of public member
     \ref connected_component.
     */
    template <typename S>
    void connected_component (image3d<S>& res, image3d<S>& bw,
                              bool const& full_connected) const;

    /*!
     \brief private member used by both version of public member
     \ref connected_component.
     */
    template <typename S>
    void cc_not_binary_output (image3d<S>& res, double threshold) const;

public:

    //CONSTRUCTORS

    //! Default Constructor to build an empty \ref image3d.
    image3d ();

    /*!
     \brief Variable Constructor to build a black (all zeros) \ref image3d of
     specified size and spacing.

     \param x set X dimension \ref dimx.
     \param y set Y dimension \ref dimy.
     \param z set Z direction \ref dimz.
     \param hx set X spacing \ref hx.
     \param hy set Y spacing \ref hy.
     \param hz set Z spacing \ref hz.
     */
    image3d (uint const& x, uint const& y, uint const& z,
             double const& hx = 1., double const& hy = 1., double const& hz = 1.);

    /*!
     \brief Copy Constructor.

     \param tocopy is the \ref image3d type object to be copied.
     */
    image3d (image3d const& tocopy);

    //!Default Denstructor
    virtual ~image3d() {};


    // MEMBERS TO GET AND SET PRIVATE PARAMETERS

    // Spacing

    //! Member to get X dimension \ref dimx.
    inline uint getdimx () const
    {
        return this->dimx;
    };

    //! Member to get Y dimension \ref dimy.
    inline uint getdimy () const
    {
        return this->dimy;
    };

    //! Member to get Z dimension \ref dimz.
    inline uint getdimz () const
    {
        return this->dimz;
    };

    /*!
     \brief Member to set X dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param x is the desired \ref dimx.
     */
    inline void setdimx (uint const& x);

    /*!
     \brief Member to set Y dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param y is the desired \ref dimy.
     */
    inline void setdimy (uint const& y);

    /*!
     \brief Member to set Z dimension.

     \warning Changing dimensions reinitializes this \ref image3d as a black
     (all zeros) image. To change dimensions without modify values of the image use
     \ref crop member.

     \param z is the desired \ref dimz.
     */
    inline void setdimz (uint const& z);