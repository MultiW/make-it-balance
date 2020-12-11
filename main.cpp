#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include <Eigen/Geometry>

#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/bounding_box.h>
#include <igl/grid.h>

#include "include/make_it_stand.h"
#include "include/grid_util.h"


#define PI 3.14159265


const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

igl::opengl::glfw::Viewer viewer;

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
const Eigen::RowVector3d black(0.0,0.0,0.0);

struct State
{
	// object current position/orientation
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd planeV;
	Eigen::MatrixXi planeF;

	bool selectingBalancePoint = true;
	Eigen::RowVector3d balancePoint;

	// gravity unit vector
	Eigen::Vector3d gravity;

	// object orientation 
	float yaw;
	float pitch;
	float roll;
} state;

double sin_deg(double angle)
{
	return std::sin(angle * PI / 180);
}

double cos_deg(double angle)
{
	return std::cos(angle * PI / 180);
}

void updateGravity()
{
	// Convert yaw, pitch, roll to vector
	// - source: https://stackoverflow.com/questions/1568568/how-to-convert-euler-angles-to-directional-vector
	state.gravity(0) = -cos_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) * -sin_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(1) = -sin_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) + cos_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(2) = cos_deg(state.pitch) * sin_deg(state.roll);

	// TODO: REMOVE gravity vector
	viewer.selected_data_index = 3;
	viewer.data().clear_points();
	viewer.data().add_points(Eigen::RowVector3d(0.0,0.0,0.0), black);
	viewer.data().add_points(state.gravity.transpose() / 4.0, black);
}

void updateObject()
{
	// TODO change object position
	// TODO shift balance point
	
}

bool findLowestPoint(const Eigen::MatrixXd &V, Eigen::RowVector3d &lowest)
{
	double minZ = V.col(2).minCoeff();
	int minIdx = -1;
	for (int i = 0; i < V.rows(); i++) {
		if (V(i, 2) == minZ) {
			minIdx = i;
			break;
		}
	}

	if (minIdx == -1) {
		return false;
	}

	lowest = V.row(minIdx);

	return true;
}

void initState(int argc, char *argv[])
{
	// Load object
  	igl::read_triangle_mesh(argc > 1 ? argv[1] : DEFAULT_MESH_FILE, state.V, state.F);

	// Align object to axis planes
	Eigen::MatrixXd unalignedV = state.V;
	alignToAxis(unalignedV, state.V);

	// Create plane from bottom of the object's bounding box
	Eigen::AlignedBox3d boundingBox;
	createAlignedBox(state.V, boundingBox);
	state.planeV.resize(4,3);
	state.planeV <<
		boundingBox.corner(boundingBox.BottomLeftCeil).transpose(),
		boundingBox.corner(boundingBox.BottomRightCeil).transpose(),
		boundingBox.corner(boundingBox.BottomLeftFloor).transpose(),
		boundingBox.corner(boundingBox.BottomRightFloor).transpose();
	state.planeF.resize(2,3);
	state.planeF <<
		0, 1, 2,
		2, 1, 3;

	// Gravity points down by default
	state.gravity << 0.0, -1.0, 0.0;
	
	// set lowest point as balance point
	if (!findLowestPoint(state.V, state.balancePoint)) {
		std::cerr << "Failed to find index of lowest mesh point";
		exit(-1);
	}
}

int main(int argc, char *argv[])
{
	initState(argc, argv);

	// Outer and inner meshes
	Eigen::MatrixXd MiV;
	Eigen::MatrixXi MiF;

	make_it_stand(state.V, state.F, MiV, MiF);

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


	// == GUI ==
	// Select balance point
	// SOURCE: from deformation assignment starter code
	Eigen::RowVector3f last_mouse; 
	viewer.callback_mouse_down =
		[&](igl::opengl::glfw::Viewer&, int, int) -> bool
	{
		last_mouse = Eigen::RowVector3f(
			viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
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
				state.V, state.F, 
				fid, bary))
			{
				long c;
				bary.maxCoeff(&c);
				state.balancePoint = state.V.row(state.F(fid, c));

				// Snap to closest vertex on hit face
				viewer.data().clear_points();
				viewer.data().set_points(state.balancePoint, blue);
				return true;
			}
		} 
	};

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		if (ImGui::CollapsingHeader("Balance", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::SliderFloat("Yaw", &(state.yaw), 0.0, 360.0);
			if (ImGui::IsItemActive())
			{
				updateGravity();
				updateObject();
			}
			ImGui::SliderFloat("Pitch", &(state.pitch), 0.0, 360.0);
			if (ImGui::IsItemActive())
			{
				updateGravity();
				updateObject();
			}
			ImGui::SliderFloat("Roll", &(state.roll), 0.0, 360.0);
			if (ImGui::IsItemActive())
			{
				updateGravity();
				updateObject();
			}
		};
	};

	// == Display Meshes ==
	// display main object
	viewer.data().set_mesh(state.V, state.F);
	viewer.data().set_face_based(true);

	// display inner voxelized mesh (representing the empty space)
	viewer.append_mesh();
	viewer.data().set_mesh(MiV, MiF);
	viewer.data().set_face_based(true);

	// display plane (representing the ground)
	viewer.append_mesh();
	viewer.data().set_mesh(state.planeV, state.planeF);
	viewer.data().set_face_based(true);
	
	// TODO: REMOVE gravity vector
	// show gravity vector
	viewer.append_mesh();
	viewer.data().clear_points();
	viewer.data().add_points(Eigen::RowVector3d(0.0,0.0,0.0), black);
	viewer.data().add_points(state.gravity.transpose(), black);

	// currently selected mesh is the input object
	viewer.selected_data_index = 0;

	viewer.launch();
}
