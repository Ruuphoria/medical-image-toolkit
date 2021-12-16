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
 * to freely rotate it in the 3-dimensional spa