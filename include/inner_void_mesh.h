#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <limits>

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>

#include <iostream>

#include "voxel_grid.h"

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
