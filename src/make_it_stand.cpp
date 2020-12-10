#include "../include/make_it_stand.h"
#include "../include/inner_void_mesh.h"

void make_it_stand(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &MiV, Eigen::MatrixXi &MiF)
{
	MiV.resize(V.rows(), V.cols());
	MiF.resize(F.rows(), F.cols());

	// TODO
	MiV = V;
	MiF = F;

	// Voxelize mesh
	InnerVoidMesh voxelMesh(V, F);
	voxelMesh.convertToMesh(MiV, MiF);
}