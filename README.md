# Make It Stand

In this program, we implemented a methodology to modify 3D objects in order to balance them
on defined balance points. The methodology is described in detail in the paper:
_Make it stand: balancing shapes for 3D fabrication_ [[1]](#1).

## Program

A 3D object printed out in reality may not balance as we want it to. Our program 
helps a user create a balanced 3D object given any point of balance
and orientation of the object.

## Implementation



## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies.

## Run

From within the `build` directory just issue:

    ./make-it-stand

A glfw app should launch displaying the default bunny mesh with menu options to select the balancing options.

## References
<a id="1">[1]</a> 
Romain Prévost, Emily Whiting, Sylvain Lefebvre, and Olga Sorkine-Hornung. 2013. 
Make it stand: balancing shapes for 3D fabrication. 
ACM Trans. Graph. 32, 4, Article 81 (July 2013), 10 pages. 
DOI: https://doi.org/10.1145/2461912.2461957