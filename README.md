# Make It Balance

In this program, we implemented a methodology to modify 3D objects in order to balance them
on defined balance points. The methodology is described in detail in the paper:
_Make it stand: balancing shapes for 3D fabrication_ [[1]](#1).

## Program

A 3D object printed out in reality may not balance as we want it to. Our program 
helps a user create a balanced 3D object given a specified point-of-balance
and an orientation of the object. To balance the object, our program selects regions 
within the object that it decides will be empty, thus adjusting the center of mass
for the object. We assume that the objects we are working with solid objects (not hollow).

Once the program is running, we can begin choosing our balancing options using the
**Make It Stand** menu on the left. First, we configure
the point of balance by clicking on object. The selected balance point
will appear as a yellow dot.

<img src="images/select-point.png" alt="drawing" width="600">
<br/>
<br/>


Next, we adjust the rotation of the object. Rotation is defined in terms of yaw, pitch, 
and roll. Note that the plane shown on the screen represents the ground.

<img src="images/rotate.png" alt="drawing" width="600">
<br/>
<br/>


Now that we defined how our object will look statically. We can begin the carving process
to make a balanced object. At first, nothing is carved, meaning that the inside of the
object is full/solid.

<img src="images/begin-carving.png" alt="drawing" width="600">
<br/>
<br/>


Click the **Carve** button run one iteration of carving. Click the  **Finish Carving** button
to complete the carving of the image.

<img src="images/finish-carving.png" alt="drawing" width="600">


## Implementation of Balancing
The first step to balancing an object is to discretize the object into voxels — commonly known as
voxelization. Converting the object into voxels simplies the problem so that we can carve 
an object by removing a voxel. An important constraint to follow is that voxels within a certain
distance from the border of the mesh should not be removed. This maintains the object's 
integrity when printed.

Next, we want to carve the object to adjust its center of mass. Having a center of mass
directly above the point-of-balance would mean the object is balanced. To carve the object, we 
draw a plane on the point-of-balance. This plane is parallel to the direction of gravity,
and perpendicular to the point-of-balance to center-of-mass vector. An illustration is shown in the paper.

One side of this plane contains the object's current center of mass. This is the side we want to
remove voxels from because removing a voxel will potentially bring the center of mass closer to the plane 
and balance point, and make the object more balanced. Removing a voxel from the other 
side of the plane will always make the shape more unbalanced.

## Code Summary
`main.cpp`
* Displays mesh objects and balancing options
* Controls the workflow of balance properties selection

`inner_void_mesh.h` and `inner_void_mesh.cpp`
* `InnerVoidMesh` class holds the current state of carving; e.g. holds a list of voxels, whether they
are removed, etc.
* `Voxel` holds the definition for a voxel
* Voxalizes the mesh
* Decides which voxel to carve

`list3d.h` and `list3d.cpp`
* 3D vector wrapper. Used to hold the 3D list of voxels

`center_of_mass.h` and `center_of_mass.cpp`
* Computes the center of mass of an object

`grid_util.h` and `grid_util.cpp`
* Helper functions to help transform (translate, scale, rotate) and voxelize objects

## Next Step
This implementation of 3D balancing could be improved by introducing deformation 
to balance the object, allowing the use of points of balance, and using polygons
as bases of balance. These were features described in the paper and were left out
in consideration of time.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

Make sure to use the `--recursive` flag if you're cloning this project. This will clone the submodules for this project — specifically, `LIBIGL`.

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