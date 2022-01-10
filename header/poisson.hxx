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
 \file   poisson.hxx

  \brief file with functor \ref lapl::unsteady_poisson_functor and the linked functions
 */

#ifndef POISSON_HXX_INCLUDED
#define POISSON_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"


/*!
 \brief A namespace grouping all object or function linked to solve the Poisson equation.
 */
namespace lapl
{
/*!
 \brief A template functor conceived to solve one step of the unsteady Poisson problem
 using finite difference.

 The functor contains a function pointer and a method to change the function pointed
 by the pointer.\n
 The function pointed by the pointer takes four arguments and return a void
 (the result is written in the first argument), it