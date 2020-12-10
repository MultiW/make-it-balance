#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <string>

#include "make_it_stand.h"


const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

int main(int argc, char *argv[])
{
	// Outer and inner meshes
	Eigen::MatrixXd V, MiV, planeV;
	Eigen::MatrixXi F, MiF, planeF;

  	igl::read_triangle_mesh(argc > 1 ? argv[1] : DEFAULT_MESH_FILE, V, F);
	igl::readOBJ("../data/plane.obj", planeV, planeF);

	make_it_stand(V, F, MiV, MiF);

	// TODO: REMOVE. Temporary test to display meshes
	//igl::read_triangle_mesh("../data/humpty_outer.stl", V, F);
	//igl::read_triangle_mesh("../data/humpty_inner.stl", MiV, MiF);

	// TODO DELETE. code to rotate shape
	//Eigen::Matrix3d rotate;
	//Eigen::Vector3d ang(-5, 5, 5);
	//Eigen::Vector3d v0 = planeV.row(1).transpose();
	//Eigen::Vector3d v1 = planeV.row(1).transpose() + ang;
	//rotate = igl::rotation_matrix_from_directions(v0, v1);
	//std::cout << rotate << std::endl;
	//planeV *= rotate.transpose();
	//planeV.col(0).array() += 0.5;
	//planeV.col(1).array() += 1.0;
	//planeV.col(2).array() += 0.5;


	// == Visualization ==
	igl::opengl::glfw::Viewer viewer;

	// display main object
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);

	// display inner voxelized mesh (representing the empty space)
	viewer.append_mesh();
	viewer.data().set_mesh(MiV, MiF);
	viewer.data().set_face_based(true);

	// display plane (representing the ground)
	viewer.append_mesh();
	viewer.data().set_mesh(planeV, planeF);
	viewer.data().set_face_based(true);

	viewer.selected_data_index = 0;

	viewer.launch();
}
