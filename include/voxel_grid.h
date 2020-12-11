#include <Eigen/Core>

//
// Computes a voxel grid enclosing a given mesh
//
// Inputs:
//   V  #vertices by 3 matrix of vertex positions
// Outputs:
//   centers  #x * #y * #z by 3 matrix of voxel centers of the voxels
//   corners  #x+1 * #y+1 * #z+1 by 3 matrix of voxel centers of the voxels
void createVoxelGrid(const Eigen::MatrixXd& V, Eigen::MatrixXd& centers, Eigen::MatrixXd& corners);
