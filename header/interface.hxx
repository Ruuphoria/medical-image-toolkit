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
 \file  interface.hxx

 \brief header file of class \ref im3d::interface
 */

#ifndef INTERFACE_HXX_INCLUDED
#define INTERFACE_HXX_INCLUDED

#include "names.hxx"
#include "image3d.hxx"
#include "abstractinterface.hxx"

#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkImageViewer.h>
#include <vtkMetaImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkJPEGReader.h>
#include <vtkBMPReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>
#include <vtkMarchingCubes.h>
#include <vtkContourFilter.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkCamera.h>
#include <vtkImageCast.h>



/*!
 * \brief A namespace in which all objects linked to
 * an image 3d are defined.
 */
namespace im3d
{

// STATIC VARIABLES
static vtkSmartPointer<vtkRenderer> RENDERER = NULL;
static vtkSmartPointer<vtkRenderWindowInteractor> RENWININT = NULL;

/*!
 * \brief A class inheriting from \ref abstractinterface conceived as a linker
 * between an \ref image3d and an \ref imageptr of type vtkSmartPointer<vtkImageData>.
 *
 * This class follows the guideline of \ref abstractinterface providing the user with
 * an interface class specialized in communication between an \ref image3d and
 * elements from the vtk library.
 * The vtk library is an open source project for 3D computer graphics, see www.vtk.org
 * for more details.
 * The use of this library makes the visualization process easy and aesthetically
 * satisfying, with planes the allow to explore the image as well as the possibility
 * to freely rotate it in the 3-dimensional space.
 * To make the program work without vtk libraries the user has to implement another
 * child class using an alternative way to implement abstract members.
 */
template <typename S>
class interface : public abstractinterface < vtkSmartPointer<vtkImageData>, S >
{

protected:

    //! member to define the background color of the visualization.
    double colour[3];
    //! member to define the opacity of the image during visualization.
    double opacity;

public:

    // CONSTRUCTORS

    /*!
     \brief default constructor
     */
    interface ();

    /*!
     \brief constructor from an image name.

     Will construct the class automatically calling the load function with the desired
     filename.

     \param imagename a string containing the path to the image to be processed,
     including the file's extension.
     Supported extensions are .mhd, .dcm, .jpg and .bmp.
     */
    interface (std::string const& imagename);

    /*!
     \brief copy constructor
     */
    interface (image3d<S> const& myim);

    /*!
     \brief default destructor
     */
    virtual ~interface () {};

    // CONVERSION

    /*!
     \brief function to convert an image3d into a T image
     */
    void convertfromimage3d (image3d<S> const& myim);

    /*!
     \brief function to convert T image into an image3d
     */
    void convert2image3d (image3d<S>& myim) const;

    // SHOWING

    /*!
     \brief member to set background color.

     the color is specified in RGB values.
     \param red the Red value.
     \param green the Green value.
     \param blue the Blue value.
     */
    void setcolour (double const& red = 1., double const& green = 1., double const& blue = 0.)
    {
        colour[0] = red;
        colour[1] = green;
        colour[2] = blue;
        return;
    };

    /*!
     \brief member to set opacity

     \param opacity should be a value between 0 and 1.
     */
    inline void setopacity (double const& opacity = 0.5)
    {
        this->opacity = opacity;
        return;
    };

    /*!
     \brief member to show only a level of a T \ref imageptr.

     This member will produce an image where only the pixels with a value equal to
     level are highlighted. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations.

     \param level the value of the level to be shown.
     */
    void show_contour (double const& level = 0);

    /*!
     \brief Member to show an image using three orthogonal planes.

     This member will show the image in a box made of three movable planes.
     These planes may be directed using the mouse, by pressing both buttons at the
     same time and moving it. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations and the addition of planes.
     */
    void show_image ();

    /*!
     \brief Member to show an image using three orthogonal planes and a selected level.

     This function is a combination of \ref show_contour
     and \ref show_image. An extensive use of the vtk library is made, concerning
     in particular rendering and/or shading operations and the addition of planes.
     \param level the value of the level to be shown.
     */
    void show_image_and_contour (double const& level = 0);

    /*!
     \brief