#include <Eigen/Core>
#include <Eigen/Geometry>

#include <set>

#include "list3d.h"

class Voxel {
public:
	Eigen::Vector3d center;
	double distanceToMesh; // distance to closest face on mesh
	Eigen::Vector3d closestNormal;

	int xIdx, yIdx, zIdx; // location of voxel within the grid

	std::set<int> cornerIndices; // index to corner vertices of voxel

	bool isFilled; // for carving algorithm
	bool isInner; // is voxel inside the given mesh

	Voxel();
};

class InnerVoidMesh {
	List3d<Voxel> innerMesh;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::Vector3d gravity;
	Eigen::Vector3d balancePoint;
	Eigen::Vector3i dimensions; // in voxels
	Eigen::MatrixXd corners; // corners of all voxels
public:
	InnerVoidMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::Vector3d &gravity, const Eigen::Vector3d &balancePoint);
	void convertToMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
};