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
 \file   segmentation.hxx

   \brief header file of class \ref segm::segmentation
 */

#ifndef SEGMENTATION_HXX_INCLUDED
#define SEGMENTATION_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"
#include "interface.hxx"

using namespace std;

/*!
 \brief A namespace in which are stored all classes or function related to a segmentation
 algorithm.
 */
namespace segm
{

/*!
 \brief A template class conceived as a container for the fundamental tools and
 properties related to the segmentation of an image using a level set method.

 This class implements some tools needed to manipulate an image during a segmentation
 algorithm. In particular it contains members that allow to visualize the image as
 well as members that are able to perform the segmentation algorithm and all the
 related computations and finally also the level set function itself.
 This class does not implement a particular algorithm on purpose, leaving the user
 the liberty to exploit it in complete freedom.
 Any particular algorithm should be implemented in a child
 class, e.g. \ref rsfe_splitbregman.

 \tparam T should solely be used to set the precision of the pixel values: float,
 double...
 */
template <typename T>
class segmentation
{

protected:

    /*!
     \brief a number representing the level of the levelset in which contour is
     expected to be.
     */
    T alpha;