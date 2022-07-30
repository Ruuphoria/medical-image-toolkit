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

 This class impl