#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>

class Voxel {
	Eigen::Vector3d center;
	bool isFilled;
	// is voxel inside the given mesh
	bool isInner;
public:
	Voxel(Eigen::Vector3d center, bool isInner);
};

class InnerVoidMesh {
	std::vector<Voxel *> *innerMesh;
public:
	InnerVoidMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	void convertToMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F);
};