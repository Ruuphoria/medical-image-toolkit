
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
  \file  interface_imp.hxx

   \brief file with the implementation of all the methods of class \ref im3d::interface
*/

#ifndef INTERFACE_IMP_HXX_INCLUDED
#define INTERFACE_IMP_HXX_INCLUDED



template <typename S>
im3d::interface<S>::interface ()
{
    setcolour();
    setopacity();
    this->imageptr = vtkImageData::New();
}



template <typename S>
im3d::interface<S>::interface (std::string const& imagename)
{
    setcolour();
    setopacity();

    // Read a specified .mhd file.
    if ( imagename.compare (imagename.size() - 4, 4, ".mhd") == 0 )
    {
        vtkSmartPointer<vtkMetaImageReader> metaReader = vtkMetaImageReader::New();
        metaReader->SetFileName (imagename.c_str() );
        metaReader->Update();
        this->imageptr = metaReader->GetOutput();
    }

    // Read a specified .dcm file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".dcm") == 0 )
    {
        vtkSmartPointer<vtkDICOMImageReader> dicom_reader = vtkDICOMImageReader::New();
        dicom_reader->SetFileName (imagename.c_str() );
        dicom_reader->Update();
        this->imageptr = dicom_reader->GetOutput();
    }

    // Read a specified .jpg file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".jpg") == 0 ||
              imagename.compare (imagename.size() - 4, 4, ".JPG") == 0 ||
              imagename.compare (imagename.size() - 5, 5, ".jpeg") == 0 ||
              imagename.compare (imagename.size() - 5, 5, ".JPEG") == 0 )
    {
        vtkSmartPointer<vtkJPEGReader> jpeg_reader = vtkJPEGReader::New();
        jpeg_reader->SetFileName (imagename.c_str() );
        jpeg_reader->Update();
        this->imageptr = jpeg_reader->GetOutput();
    }

    // Read a specified .bmp file.
    else if ( imagename.compare (imagename.size() - 4, 4, ".bmp") == 0 )
    {
        vtkSmartPointer<vtkBMPReader> bmp_reader = vtkBMPReader::New();
        bmp_reader->SetFileName (imagename.c_str() );
        bmp_reader->Update();
        this->imageptr = bmp_reader->GetOutput();
    }

    // Read all the DICOM files in the specified directory.
    else
    {
        vtkSmartPointer<vtkDICOMImageReader> dicom_reader = vtkDICOMImageReader::New();
        dicom_reader->SetDirectoryName (imagename.c_str() );
        dicom_reader->Update();
        this->imageptr = dicom_reader->GetOutput();
    }

    return;
}



template <typename S>
im3d::interface<S>::interface (image3d<S> const& myim)
{
    setcolour();
    setopacity();
    this->imageptr = vtkImageData::New();
    this->convertfromimage3d (myim);
}



template <typename S>
void im3d::interface<S>::convertfromimage3d (image3d<S> const& myim)
{
    this->imageptr->SetSpacing (static_cast<double> (myim.gethx() ),
                                static_cast<double> (myim.gethy() ),
                                static_cast<double> (myim.gethz() ) );
    this->imageptr->SetExtent (0, myim.getdimx() - 1, 0, myim.getdimy() - 1, 0, myim.getdimz() - 1);
#if VTK_MAJOR_VERSION <= 5
    this->imageptr->SetNumberOfScalarComponents (1);
    this->imageptr->SetScalarTypeToDouble();
#else
    this->imageptr->AllocateScalars (VTK_DOUBLE, 1);
#endif
    for (uint i = 0; i < myim.getdimx(); ++i)
        for (uint j = 0; j < myim.getdimy(); ++j)
            for (uint k = 0; k < myim.getdimz(); ++k)
                this->imageptr->SetScalarComponentFromDouble
                (i, j, k, 0, static_cast<double> (myim (i, j, k) ) );
    return;
}



template <typename S>
void im3d::interface<S>::convert2image3d (image3d<S>& myim) const
{
    int dims[3];
    double h[3];

    this->imageptr->GetDimensions (dims);

    myim.setdim (dims[0], dims[1], dims[2]);

    this->imageptr->GetSpacing (h);

    myim.sethx ( static_cast<double> (h[0]) );
    myim.sethy ( static_cast<double> (h[1]) );
    myim.sethz ( static_cast<double> (h[2]) );

    #pragma omp parallel for
    for (uint i = 0; i < myim.getdimx(); ++i)
        for (uint j = 0; j < myim.getdimy(); ++j)
            for (uint k = 0; k < myim.getdimz(); ++k)
                myim (i, j, k) =
                    static_cast<S> (this->imageptr->GetScalarComponentAsDouble (i, j, k, 0) );

    return;
}




template <typename S>
void im3d::interface<S>::show_contour (double const& level)
{
    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    setopacity();
    setcolour();
    addcontour2renderer (renderer, level);

    startshowing (renderer, renWinInt);

    deleterendering (renderer, renWinInt);

    return;
}



template <typename S>
void im3d::interface<S>::show_image ()
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];


    vtkSmartPointer<vtkRenderer> renderer = NULL;
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInt = NULL;

    newrendering (renderer, renWinInt);

    vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();

    renWin->SetInteractor (renWinInt);
    renWin->AddRenderer (renderer);

    vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetZ = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
    imagePlaneWidgetZ->SetInputData (this->imageptr);
#else
    imagePlaneWidgetZ->SetInput (this->imageptr);
#endif
    imagePlaneWidgetZ->SetInteractor (renWinInt);
    imagePlaneWidgetZ->SetPlaneOrientationToZAxes();
    imagePlaneWidgetZ->DisplayTextOn();
    imagePlaneWidgetZ->UseContinuousCursorOn();
    imagePlaneWidgetZ->On();

    if (dimz != 1)
    {
        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetX = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetX->SetInputData (this->imageptr);
#else
        imagePlaneWidgetX->SetInput (this->imageptr);
#endif
        imagePlaneWidgetX->SetInteractor (renWinInt);
        imagePlaneWidgetX->SetPlaneOrientationToXAxes();
        imagePlaneWidgetX->DisplayTextOn();
        imagePlaneWidgetX->UseContinuousCursorOn();
        imagePlaneWidgetX->On();

        vtkSmartPointer<vtkImagePlaneWidget> imagePlaneWidgetY = vtkImagePlaneWidget::New();
#if VTK_VERSION_MAJOR >= 6
        imagePlaneWidgetY->SetInputData (this->imageptr);
#else
        imagePlaneWidgetY->SetInput (this->imageptr);
#endif
        imagePlaneWidgetY->SetInteractor (renWinInt);
        imagePlaneWidgetY->SetPlaneOrientationToYAxes();
        imagePlaneWidgetY->DisplayTextOn();
        imagePlaneWidgetY->UseContinuousCursorOn();
        imagePlaneWidgetY->On();
    }

    interface::startshowing (renderer, renWinInt);

    setopacity();
    setcolour();

    deleterendering (renderer, renWinInt);
}



template <typename S>
void im3d::interface<S>::get_coordinates (double& x, double& y, double& z)
{
    int dims[3], dimz;

    this->imageptr->GetDimensions (dims);
    dimz = dims[2];

    std::cout << "Showing your image ..." << std::endl;
    std::cout << "Click with mouse on the image to know coordinates." << std::endl;
    std::cout << "When you have chosen a pixel write ";
    std::cout << "coordinates (positive) in a piece of paper and close the image.";
    std::cout << std::endl << std::endl;

    this->show_image();

    // 2d case
    if (dimz == 1)
    {
        std::cout << "Write coordinates of the pixel (x,y):" << std::endl;

        do
        {
            std::cout << "\tx: ";
            cin >> x;
        }
        while (x < 0);

        do
        {
            std::cout << "\ty: ";
            cin >> y;
        }
        while (y < 0);

        z = 0;
    }
    // 3d case
    else
    {
        std::cout << "Write coordinates of the pixel (x,y,z):" << std::endl;

        do
        {
            std::cout << "\tx: ";
            cin >> x;
        }
        while (x < 0);

        do
        {
            std::cout << "\ty: ";
            cin >> y;
        }
        while (y < 0);

        do
        {
            std::cout << "\tz: ";
            cin >> z;
        }
        while (z < 0);
    }


    return;
}



template <typename S>
void im3d::interface<S>::get_coordinates (uint& i, uint& j, uint& k)
{
    double x, y, z;
    int dims[3];
    double h[3];

    this->imageptr->GetDimensions (dims);
    this->imageptr->GetSpacing (h);


    this->get_coordinates (x, y, z);

    i = floor (x / h[0]);
    j = floor (y / h[1]);

    if (dims[2] == 1)
    {
        k = 0;
    }
    else
    {
        k = floor (z / h[2]);
    }

    if ( static_cast<int> (i) >= dims[0] )
    {
        i = dims[0] - 1;
    }
    if ( static_cast<int> (j) >= dims[1] )