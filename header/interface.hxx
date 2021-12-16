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
 * \brief A namespace in which all objects link