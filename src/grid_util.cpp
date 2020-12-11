#include "../include/grid_util.h"

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>

const static int GRID_LEN = 15; // number of voxels on each side of the grid

void createAlignedBox(const Eigen::MatrixXd &V, Eigen::AlignedBox3d &box) {
	for (int i = 0; i < V.rows(); i++) {
		box.extend(V.row(i).transpose());
	}
}

void transformGrid(const Eigen::MatrixXd &V, const Eigen::AlignedBox3d &destination, Eigen::MatrixXd &Vout) {
	Vout.resize(V.rows(), 3);

	Eigen::Vector3d corner1 = destination.corner(destination.BottomLeftFloor);
	Eigen::Vector3d corner2 = destination.corner(destination.TopRightCeil);

	// Compute diagonals of input and output grids
	Eigen::AlignedBox<double, 3> inBox;
	createAlignedBox(V, inBox);
	Eigen::Vector3d inBLF = inBox.corner(inBox.BottomLeftFloor);
	Eigen::Vector3d inTRC = inBox.corner(inBox.TopRightCeil);
	Eigen::Vector3d inDiag = inBLF - inTRC;

	Eigen::Vector3d outDiag = corner1 - corner2;

	// Scale V to desired size
	Eigen::Vector3d scale = outDiag.cwiseQuotient(inDiag);

	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) = V.row(i).cwiseProduct(scale.transpose());
	}

	// Move V to desired location
	Eigen::Vector3d translate = corner1 - inBLF.cwiseProduct(scale);
	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) += translate.transpose();
	}
}

void createVoxelGrid(const Eigen::MatrixXd &V, Eigen::MatrixXd &centers, Eigen::MatrixXd &corners) {
	// Compute bounding box
	Eigen::MatrixXd BV, BF;
	igl::bounding_box(V, BV, BF);
	Eigen::AlignedBox3d boundBox;
	for (int i = 0; i < V.rows(); i++) {
		boundBox.extend(V.row(i).transpose());
	}

	// Create voxel grid (defined by center vertices of voxels)
	Eigen::RowVector3i dimOut; // column vector of dimentions
	igl::voxel_grid(boundBox, GRID_LEN, 0, centers, dimOut);
	Eigen::Vector3i dimVoxels = dimOut.transpose();
	Eigen::AlignedBox3d centersBox;
	createAlignedBox(centers, centersBox);

	Eigen::Vector3i one(1.0, 1.0, 1.0);
	Eigen::Vector3d voxelLen = centersBox.sizes().cwiseQuotient((dimVoxels - one).cast<double>());

	// Compute corner vertices of voxel grid
	Eigen::MatrixXd defaultGrid;
	Eigen::Vector3d dimVertices = (dimVoxels + one).cast<double>();
	igl::grid(dimVertices, defaultGrid);

	// scale default grid to voxel grid
	Eigen::Vector3d outBLF = centersBox.corner(centersBox.BottomLeftFloor) - (voxelLen / 2.0);
	Eigen::Vector3d outTRC = centersBox.corner(centersBox.TopRightCeil) + (voxelLen / 2.0);
	Eigen::AlignedBox3d outBox(outBLF, outTRC);
	transformGrid(defaultGrid, outBox, corners);
}