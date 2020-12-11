#include <Eigen/Core>
#include <Eigen/Geometry>


//
// Inputs:
//   V    #vertices by 3 matrix of vertex representing a box
// Outputs:
//   box  AlignedBox constructed from V
void createAlignedBox(const Eigen::MatrixXd& V, Eigen::AlignedBox3d& box);


//
// Scales and translates the given set of vertices
//
// Transform the vertices based on the transformation between the vertices' bounding
// box and the given destination bounding box
//
// Inputs:
//   V            #vertices by 3 matrix of vertices
//   destination  desired bounding box that transformed V should be in
//   corner1      first corner of expected bounding box
//   corner2      opposite corner of corner1
//   keepSize    boolean indicating whether Vout should be the same size as V
// Outputs:
//   Vout     scaled and transformed v
void transformGrid(const Eigen::MatrixXd& V, const Eigen::AlignedBox3d &destination, Eigen::MatrixXd& Vout, bool keepSize);

//
// Computes a voxel grid enclosing a given mesh
//
// Inputs:
//   V  #vertices by 3 matrix of vertex positions
// Outputs:
//   centers  #x * #y * #z by 3 matrix of voxel centers of the voxels
//   corners  #x+1 * #y+1 * #z+1 by 3 matrix of voxel centers of the voxels
void createVoxelGrid(const Eigen::MatrixXd& V, Eigen::MatrixXd& centers, Eigen::MatrixXd& corners);

//
// Aligns given mesh to the axis planes. Keeps the size the same
//
// Inputs:
//   V         #vertices by 3 matrix of vertices
// Outputs:
//   alignedV  #vertices by 3 matrix of vertices transformed to align to axis planes
void alignToAxis(const Eigen::MatrixXd& V, Eigen::MatrixXd& alignedV);
