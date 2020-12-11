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
#include <igl/centroid.h>

#include "include/make_it_stand.h"
#include "include/grid_util.h"


#define PI 3.14159265

const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

// Viewer data indices
const int MESH_DATA_IDX = 0;
const int INNER_DATA_MESH_IDX = 1;
const int PLANE_DATA_IDX = 2;
const int GRAVITY_DATA_IDX = 3;

const Eigen::Vector3d DEFAULT_GRAVITY(0.0,-1.0,0.0);

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
const Eigen::RowVector3d black(0.0,0.0,0.0);

igl::opengl::glfw::Viewer viewer;


struct State
{
	// object current position/orientation
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd planeV;
	Eigen::MatrixXi planeF;
	Eigen::RowVector3d planeCenter;

	bool selectingBalancePoint = true;
	int balancePointIdx;

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

Eigen::RowVector3d getBalancePoint()
{
	return state.V.row(state.balancePointIdx);
}

void updateBalancePoint()
{
	viewer.selected_data_index = MESH_DATA_IDX;
	viewer.data().clear_points();
	viewer.data().set_points(getBalancePoint(), yellow);
}

void updateGravity()
{
	// Convert yaw, pitch, roll to vector
	// - source: https://stackoverflow.com/questions/1568568/how-to-convert-euler-angles-to-directional-vector
	state.gravity(0) = -cos_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) * -sin_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(1) = -sin_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) + cos_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(2) = cos_deg(state.pitch) * sin_deg(state.roll);

	if (state.gravity.isZero()) {
		state.gravity = DEFAULT_GRAVITY;
	}

	// TODO: REMOVE gravity vector
	viewer.selected_data_index = GRAVITY_DATA_IDX;
	viewer.data().clear_points();
	viewer.data().add_points(Eigen::RowVector3d(0.0,0.0,0.0), black);
	viewer.data().add_points(state.gravity.transpose() / 10.0, black);
}

void updateMesh()
{
	// Rotate mesh (in reverse direction from gravity)
	Eigen::MatrixXd newOrientationPoints;
	newOrientationPoints.resize(2,3);
	newOrientationPoints <<
		Eigen::RowVector3d::Zero(),
		getBalancePoint();
	Eigen::AlignedBox3d newOrientation;
	createAlignedBox(newOrientationPoints, newOrientation);

	Eigen::MatrixXd oldV = state.V;
	transformVertices(oldV, newOrientation, true, state.V);

	//Eigen::Matrix3d R;
	//R <<
	//	cos_deg(state.yaw)*cos_deg(state.pitch), -cos_deg(state.yaw)*sin_deg(state.pitch)*sin_deg(state.roll) - sin_deg(state.yaw)*cos_deg(state.roll), -cos_deg(state.yaw)*sin_deg(state.pitch)*cos_deg(state.roll) + sin_deg(state.yaw)*sin_deg(state.roll),
	//	sin_deg(state.yaw)*cos_deg(state.pitch), -sin_deg(state.yaw)*sin_deg(state.pitch)*sin_deg(state.roll) + cos_deg(state.yaw)*cos_deg(state.roll), -sin_deg(state.yaw)*sin_deg(state.pitch)*cos_deg(state.roll) - cos_deg(state.yaw)*sin_deg(state.roll),
	//	sin_deg(state.pitch), cos_deg(state.pitch)*sin_deg(state.roll), cos_deg(state.pitch)*sin_deg(state.roll);
	//state.V *= R.transpose();

	// Translate mesh to middle of plane
	Eigen::RowVector3d translate = state.planeCenter - getBalancePoint();
	for (int i = 0; i < state.V.rows(); i++) {
		state.V.row(i) += translate;
	}

	// Update viewer
	viewer.selected_data_index = MESH_DATA_IDX;
	viewer.data().set_mesh(state.V, state.F);
}

bool findLowestPointIdx(const Eigen::MatrixXd &V, int &lowestIdx)
{
	int heightCol = 1;
	double minZ = V.col(heightCol).minCoeff();
	lowestIdx = -1;
	for (int i = 0; i < V.rows(); i++) {
		if (V(i, heightCol) == minZ) {
			lowestIdx = i;
			break;
		}
	}

	if (lowestIdx == -1) {
		return false;
	}

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

	state.planeCenter = state.planeV.colwise().sum() / 4.0;

	// Gravity points down by default
	state.yaw = 0.0;
	state.pitch = 0.0;
	state.roll = 0.0;
	state.gravity = DEFAULT_GRAVITY;
	
	// set lowest point as balance point
	if (!findLowestPointIdx(state.V, state.balancePointIdx)) {
		std::cerr << "Failed to find index of lowest mesh point";
		exit(-1);
	}
	updateBalancePoint();
}

int main(int argc, char *argv[])
{
	initState(argc, argv);

	// Outer and inner meshes
	Eigen::MatrixXd MiV;
	Eigen::MatrixXi MiF;

	make_it_stand(state.V, state.F, MiV, MiF);

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
				state.balancePointIdx = state.F(fid, c);

				// Snap to closest vertex on hit face
				updateBalancePoint();
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
				updateMesh();
			}
			ImGui::SliderFloat("Pitch", &(state.pitch), 0.0, 360.0);
			if (ImGui::IsItemActive())
			{
				updateGravity();
				updateMesh();
			}
			ImGui::SliderFloat("Roll", &(state.roll), 0.0, 360.0);
			if (ImGui::IsItemActive())
			{
				updateGravity();
				updateMesh();
			}
		};
	};

	// == Display Meshes ==
	// 0: display main mesh 
	viewer.selected_data_index = MESH_DATA_IDX;
	viewer.data().set_face_based(true);
	viewer.data().set_mesh(state.V, state.F);

	// 1: display inner voxelized mesh (representing the empty space)
	viewer.append_mesh();
	viewer.data().set_face_based(true);
	viewer.data().set_mesh(MiV, MiF);

	// 2: display plane (representing the ground)
	viewer.append_mesh();
	viewer.data().set_face_based(true);
	viewer.data().set_mesh(state.planeV, state.planeF);
	
	// TODO: REMOVE gravity vector
	// 3: show gravity vector
	viewer.append_mesh();
	viewer.data().add_points(Eigen::RowVector3d(0.0,0.0,0.0), black);
	viewer.data().add_points(state.gravity.transpose() / 10.0, black);

	// currently selected mesh is the input object
	viewer.selected_data_index = MESH_DATA_IDX;

	viewer.launch();
}
