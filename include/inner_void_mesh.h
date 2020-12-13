#include <Eigen/Core>
#include <Eigen/Geometry>

#include <set>

#include "list3d.h"

class Voxel {
public:
	Eigen::Vector3d center;
	double distanceToMesh; // distance to closest face on mesh

	int xIdx, yIdx, zIdx; // location of voxel within the grid

	std::set<int> cornerIndices; // index to corner vertices of voxel

	bool isFilled; // for carving algorithm
	bool isInner; // is voxel inside the given mesh

	Voxel() : center(Eigen::Vector3d::Zero()), distanceToMesh(0), 
		xIdx(-1), yIdx(-1), zIdx(-1), isFilled(true), isInner(false) {}
};

class InnerVoidMesh {
	List3d<Voxel> innerMesh;
	Eigen::MatrixXd V; // object
	Eigen::MatrixXi F;
	Eigen::MatrixXd planeV; // plane
	Eigen::MatrixXi planeF;
	Eigen::Vector3d gravity;
	Eigen::Vector3d balancePoint;

	Eigen::Vector3i dimensions; // in voxels
	Eigen::Vector3d voxelSize;
	Eigen::MatrixXd corners; // corners of all voxels

	Eigen::Vector3d targetCOM; // optimal center of mass for object to balance
	Eigen::Vector3d currCOM; // center of mass after carving the object
	double mass;
public:
	InnerVoidMesh(
		const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
		const Eigen::MatrixXd &planeV, const Eigen::MatrixXi &planeF, 
		const Eigen::Vector3d &gravity, 
		const Eigen::Vector3d &balancePoint);
	void convertToMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void carveInnerMesh(Eigen::MatrixXd &carvePlaneV, Eigen::MatrixXi &carvePlaneF, Eigen::Vector3d& com);
	double getCoMDistance();
	bool isOptimized();
};