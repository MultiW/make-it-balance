#include <iostream>
#include <vector>
#include <string>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>

#include "make_it_stand.h"


const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);

struct State
{
	bool selectingBalancePoint = true;
	Eigen::MatrixXd balancePoints;
} state;

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

	// Select balance point
	// SOURCE: from deformation assignment starter code
	Eigen::RowVector3f last_mouse; 
	viewer.callback_mouse_down =
		[&](igl::opengl::glfw::Viewer&, int, int) -> bool
	{
		last_mouse = Eigen::RowVector3f(
			viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
		if (state.selectingBalancePoint)
		{
			// Find closest point on mesh to mouse position
			int fid;
			Eigen::Vector3f bary;
			if (igl::unproject_onto_mesh(
				last_mouse.head(2), 
				viewer.core().view,
				viewer.core().proj, 
				viewer.core().viewport, 
				V, F, 
				fid, bary))
			{
				long c;
				bary.maxCoeff(&c);
				Eigen::RowVector3d selectedPoint = V.row(F(fid, c));

				// Snap to closest vertex on hit face
				viewer.data().clear_points();
				viewer.data().set_points(selectedPoint, blue);
				return true;
			}
		} 
	};

	viewer.launch();
}
