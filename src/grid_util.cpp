#include "../include/grid_util.h"

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
#include <igl/rotation_matrix_from_directions.h>

const static int GRID_LEN = 15; // number of voxels on each side of the grid

void createAlignedBox(const Eigen::MatrixXd &V, Eigen::AlignedBox3d &box) {
	for (int i = 0; i < V.rows(); i++) {
		box.extend(V.row(i).transpose());
	}
}

void transformVertices(const Eigen::MatrixXd &V, const Eigen::AlignedBox3d &newOrientation, bool keepSize, Eigen::MatrixXd &Vout) {
	Vout.resize(V.rows(), 3);

	Eigen::Vector3d corner1 = newOrientation.corner(newOrientation.BottomLeftFloor);
	Eigen::Vector3d corner2 = newOrientation.corner(newOrientation.TopRightCeil);

	// Compute diagonals of input and output grids
	Eigen::AlignedBox<double, 3> inBox;
	createAlignedBox(V, inBox);
	Eigen::Vector3d inBLF = inBox.corner(inBox.BottomLeftFloor);
	Eigen::Vector3d inTRC = inBox.corner(inBox.TopRightCeil);
	Eigen::Vector3d inDiag = inBLF - inTRC;

	Eigen::Vector3d outDiag = corner1 - corner2;

	// Scale V to desired size
	Eigen::Matrix3d scale = Eigen::Matrix3d::Identity();
	if (!keepSize) {
		scale.diagonal() << outDiag.cwiseQuotient(inDiag).array();
	}
	Vout = V * scale;

	// Rotate
	Eigen::Matrix3d R = igl::rotation_matrix_from_directions(inDiag.normalized(), outDiag.normalized());
	Vout *= R.transpose();

	// Move V to desired location
	Eigen::Vector3d translate = corner1 - scale * inBLF;
	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) += translate.transpose();
	}
}

void transformVertices(const Eigen::MatrixXd& V, const Eigen::AlignedBox3d& newOrientation, Eigen::MatrixXd& Vout) {
	transformVertices(V, newOrientation, false, Vout);
}

void createVoxelGrid(
	const Eigen::MatrixXd &V, 
	Eigen::MatrixXd &centers, 
	Eigen::MatrixXd &corners, 
	Eigen::Vector3i &dimensions, 
	Eigen::Vector3d &voxelSize
) {
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
	dimensions = dimOut.transpose();
	Eigen::AlignedBox3d centersBox;
	createAlignedBox(centers, centersBox);

	Eigen::Vector3i one(1.0, 1.0, 1.0);
	voxelSize = centersBox.sizes().cwiseQuotient((dimensions - one).cast<double>());

	// Compute corner vertices of voxel grid
	Eigen::MatrixXd defaultGrid;
	Eigen::Vector3d dimCorners = (dimensions + one).cast<double>();
	igl::grid(dimCorners, defaultGrid);

	// scale default grid to voxel grid
	Eigen::Vector3d outBLF = centersBox.corner(centersBox.BottomLeftFloor) - (voxelSize / 2.0);
	Eigen::Vector3d outTRC = centersBox.corner(centersBox.TopRightCeil) + (voxelSize / 2.0);
	Eigen::AlignedBox3d outBox(outBLF, outTRC);
	transformVertices(defaultGrid, outBox, corners);
}

void alignToAxis(const Eigen::MatrixXd& V, Eigen::MatrixXd& alignedV) {
	Eigen::AlignedBox3d inputBox;
	createAlignedBox(V, inputBox);

	Eigen::AlignedBox3d alignedBox;
	Eigen::MatrixXd alignedGrid;
	Eigen::Vector3d dim = Eigen::Vector3d::Constant(10.0);
	igl::grid(dim, alignedGrid);
	createAlignedBox(alignedGrid, alignedBox);

	transformVertices(V, alignedBox, true, alignedV);
}