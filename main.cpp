#include <igl/opengl/glfw/Viewer.h>
#include <igl/fast_winding_number.h>
#include <igl/winding_number.h>
#include "make_it_stand.h"

const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

int main(int argc, char *argv[])
{
	// Outer and inner meshes
	Eigen::MatrixXd V, MiV;
	Eigen::MatrixXi F, MiF;

  	igl::read_triangle_mesh(argc > 1 ? argv[1] : DEFAULT_MESH_FILE, V, F);

	make_it_stand(V, F, MiV, MiF);

  	// == Visualization ==
  	igl::opengl::glfw::Viewer viewer;

	// display main object
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.data().show_faces = false;

	// display inner voxelized mesh (representing the empty space)
	viewer.append_mesh();
	viewer.data().set_mesh(MiV, MiF);
	viewer.data().set_face_based(true);

	viewer.launch();
}
