#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <limits>

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>

#include <iostream>

const static int GRID_LEN = 15; // number of voxels on each side of the grid

class Voxel {
	Eigen::Vector3d center;
	bool isFilled;
	// is voxel inside the given mesh
	bool isInner;
public:
	Voxel(Eigen::Vector3d center, bool isInner) {
		center = center;
		isInner = isInner;
		isFilled = false;
	}
};

class InnerVoidMesh {
	std::vector<Voxel *> *innerMesh;
public:
	InnerVoidMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	void convertToMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
};

void createAlignedBox(const Eigen::MatrixXd &V, Eigen::AlignedBox<double, 3> &box) {
	for (int i = 0; i < V.rows(); i++) {
		box.extend(V.row(i).transpose());
	}
}

void transformGrid(const Eigen::MatrixXd &V, const Eigen::Vector3d &outBLF, const Eigen::Vector3d &outTRC, Eigen::MatrixXd &Vout) {
	Vout.resize(V.rows(), 3);

	// Compute diagonals of input and output grids
	Eigen::AlignedBox<double, 3> inBox;
	createAlignedBox(V, inBox);
	Eigen::Vector3d inBLF = inBox.corner(inBox.BottomLeftFloor);
	Eigen::Vector3d inTRC = inBox.corner(inBox.TopRightCeil);
	Eigen::Vector3d inDiag = inBLF - inTRC;

	Eigen::Vector3d outDiag = outBLF - outTRC;

	// Scale V to desired size
	Eigen::Vector3d scale = outDiag.cwiseQuotient(inDiag);

	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) = V.row(i).cwiseProduct(scale.transpose());
	}

	// Move V to desired location
	Eigen::Vector3d translate = outBLF - inBLF.cwiseProduct(scale);
	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) += translate.transpose();
	}
}

void createVoxelGrid(const Eigen::MatrixXd &V, Eigen::MatrixXd &centers, Eigen::MatrixXd &corners) {
	// Compute bounding box
	Eigen::MatrixXd BV, BF;
	igl::bounding_box(V, BV, BF);
	Eigen::AlignedBox<double, 3> boundBox;
	for (int i = 0; i < V.rows(); i++) {
		boundBox.extend(V.row(i).transpose());
	}

	// Create voxel grid (defined by center vertices of voxels)
	Eigen::RowVector3i dimOut; // column vector of dimentions
	igl::voxel_grid(boundBox, GRID_LEN, 0, centers, dimOut);
	Eigen::Vector3i dimVoxels = dimOut.transpose();
	Eigen::AlignedBox<double, 3> centersBox;
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
	transformGrid(defaultGrid, outBLF, outTRC, corners);
}

InnerVoidMesh::InnerVoidMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
	Eigen::MatrixXd centers, corners;
	createVoxelGrid(V, centers, corners);

	// Compute voxel distances to mesh
	Eigen::VectorXd distances, closestFaces;
	Eigen::MatrixXd closestPoints, closestNormals;
	igl::signed_distance(centers, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, 
		std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), distances, 
		closestFaces, closestPoints, closestNormals);

	// TODO: Use 3D array
	// TODO: Initialize inner void mesh structure
	this->innerMesh = new std::vector<Voxel *>();
	for (int i = 0; i < centers.rows(); i++) {
		bool isInner = distances(i) < 0;

		Voxel* voxel = new Voxel(centers.row(i).transpose(), isInner);
		innerMesh->push_back(voxel);
	}
}

void InnerVoidMesh::convertToMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {

}
