# MedicalImageToolkit: A comprehensive toolkit for Medical Image Segmentation and Filtering
A C++ library to perform core image filtering and a semi-automatic segmentation through the Connected Component Localization of the Region-Scalable Fitting Energy.

The theory supporting this segmentation procedure is detailed in the following [publication](https://doi.org/10.1109/ISPA.2015.7306035):

- Fedele, M., Faggiano, E., Barbarotta, L., Cremonesi, F., Formaggia, L., & Perotto, S. (2015, September). Semi-automatic three-dimensional vessel segmentation using a connected component localization of the Region-Scalable Fitting Energy. In *Image and Signal Processing and Analysis (ISPA), 2015 9th International Symposium on (pp. 72-77). IEEE.*

If you exploit this work, please kindly attribute this paper.

## License
Copyright © 2013-2016 Luca Barbarotta, Francesco Cremonesi, Elena Faggiano, Marco Fedele. All Rights Reserved.

This library is dispensed with the BSD license, see file [license.md](./license.md) for further details.

## Authors and Contacts
- **Ruuphoria**

The original implementation of the segmentation algorithm was supervised by [Prof. Luca Formaggia](https://mox.polimi.it/people-detail/?id=142) and [Prof. Simona Perotto](https://mox.polimi.it/people-detail/?id=117).

![Alt text](./images/logo.jpg)

## Installation
This library is validated on Linux and harnesses the following libraries:

- [cmake](http://www.cmake.org) to construct the MakeFile
- [VTK](http://www.vtk.org) to visualize, to import and to drop 3d images
- [fftw](http://www.fftw.org) to efficiently compute DFT
- [GetPot](http://getpot.sourceforge.net/) to import parameters from a data file

Prior to initiating the installation, all these libraries must be installed.

Note:

- The library requires the fftw version in tandem with openmp (if you compile it, employ the option --enable-shared)
- You can download GetPot from the website and copy it in your path, for example `/usr/include/`
- In debian grounded Linux distribution all the requisite packages (except GetPot) should be acquirable with apt-get

Upon installing all these packages, you can exploit cmake to compile all the executables:

```
mkdir build_directory
cd path_to_build_directory
cmake path_to_source_directory
make
```

## Documentaion
To shape documentation you need to install [doxygen](http://www.doxygen.org) and [graphviz](http://www.graphviz.org) and then type:

```
cmake path_to_source_directory -DBUILD_DOCUMENTATION=ON
make doc
```

The html guide is at `build_directory/doc/html/`

## Usage
Launch an executable sans any option to witness its helper.