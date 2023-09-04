# DIRT - Dartmouth Introductory Ray Tracer

Dirt is a simple ray tracer written in C++ and originally was used for the Rendering Algorithms course taught at Dartmouth College. It provides basic functionalities that are required to complete the assignments for Physics-Based Rendering (CS 15-468) taught at Carnegie Mellon.

### Scene Files
The scene files used in assignments must be downloaded separately. Go to the [v0 release](https://github.com/cmu-15-468/dirt-s23/releases/tag/v0) and download the `scenes.tar.gz` bundle. You can unzip this bundle using `tar -xf scenes.tar.gz` on the command line. Place the unzipped folder inside of your repo (it should be on the same level as your src and include directories). Note that you can ignore all of the other bundles on the release page.

### High Level Overview
***********************

The Dirt repository contains source code that are a part of the actual ray tracer as well as other external libraries that dirt depends upon.

The external dependencies are located within the `ext` subdirectory, and include:

* `filesystem` - for manipulating paths on Linux/Windows/Mac OS
* `json` - for handling json input file formats of scenes
* `stb_image` - for loading and storing images
* `tinyformat` - pretty formatting and printing of strings
* `pcg32` - a small, self-contained cross-platform random number generator

The major components of the dirt can be split into the following.

* Utility
    + argparse.h - class to help parse command-line arguments
    + common.h - utility code used throughout dirt. 
    + image.h - manipulates image data. 
    + obj.h - loading geometry data from OBJ files.
    + parser.h - parse dirt scenes and data structuresa from a JSON file. 
    + progress.h - manage/print a progress bar on the command-line
    + timer.h - 
* Math
    + transform.h - contains code that implements 3-D transformations.
    + vec.h - functionality related to N-D vectors (including colors), 4x4 matrices, rays, and axis-aligned bounding boxes
* Raytracing
    + camera.h - contains the camera class that represents the virtual camera placed inside a scene.
    + material.h - contains the material class that determines the shaded color of an object.
    + surface.h - contains the base surface class that represents any geometry within the scene.
    + surfacegroup.h - contains the base class for aggregate/grouped geometry (an acceleration structure)
    + scene.h - contains the scene class that stores all geometry, materials and camera.
    + testscenes.cpp - a few hard-coded scenes

### Building Instructions
*************************

Building from the command-line generally looks like this:

```
mkdir build
cd build
cmake ..
cmake --build . -j 4
```

More detailed instructions for building dirt are provided in the course website.

### License Information
***********************

Check the LICENSE file
