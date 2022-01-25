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
 (the result is written in the first argument), it aims to solve an elliptic PDE as:
 \f[
    \left\{
        \begin{array}{l}
        - \Delta u=f \\
        + \mbox{ }bc
    \end{array}
    \right.
 \f]
 The boundary condition depends on the function pointed.\n
 Within the namespace there are two already implemented functions
 \ref lapl::NeuGaussSeidel and \ref lapl::DirGaussSeidel
 which solve the problem with the relative boundary condition type.\n
 For our purposes the value on the boundary is always constant so the class as it
 is can solve the problem above only in case of constant boundary value.

 */
template <typename T>
class unsteady_poisson_functor
{

private:
    //! Private member storing a function pointer
    void (*f) (im3d::image3d<T>& res, im3d::image3d<T> const& b, T const& dt, T const& bc,
               im3d::image3d<T> const& input);

public:

    /*!
     \brief Constructor taking a function pointer to assign to the private function
     pointer
     */
    unsteady_poisson_functor (void (*f) (im3d::image3d<T>& res, im3d::image3d<T> const& b,
                                         T const& dt, T const& bc,
                                         im3d::image3d<T> const& input) ) : f (f) {};

    /*!
     \brief This operator makes class \ref unsteady_poisson_functor a functor.

     Operator() calls the private function f with the same arguments

     \param res is the first argument of the private function where result is written.
     \param b is the second argument of the private function and is the known term of
     the equation.
     \param dt is the third argument of the private function and is the time spacing
     between two steps.
     \param bc is the fourth argument of the private function and is the value of the
     boundary condition
     \param input is the input image to be update. Default is an empty image because
     if you use an 'in-place' method to solve linear system (like
     \ref NeuGaussSeidel) you would prefer to save input image directly in parameter
     res to limit memory usage.
     */
    void operator () (im3d::image3d<T>& res, im3d::image3d<T> const& b,
                      T const& dt = 1., T const& bc = 0.,
                      im3d::image3d<T> const& input = im3d::image3d<T>() )
    {
        this->f (res, b, dt, bc, input);
        return;
    };
    /*!
     \brief This function reassign the private function pointer with the new one
     passed as argument.

     \param f is a void returning function pointer taking three parameters:
     \ref im3d::image3d res, \ref im3d::image3d b and