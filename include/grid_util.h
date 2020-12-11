#include <Eigen/Core>
#include <Eigen/Geometry>

/*
  libigl style helper functions.
*/


//
// Inputs:
//   V    #V by 3 matrix of vertex representing a box
// Outputs:
//   box  AlignedBox constructed from V
void createAlignedBox(const Eigen::MatrixXd& V, Eigen::AlignedBox3d& box);


//
// Scale, rotate, and translate the given mesh cloud
// - Based on the transformation between the vertices' bounding
//   box and the given destination bounding box
// - Translation is based on the BottomLeftFloor (Eigen::AlignedBox corner)
//
// Inputs:
//   V               #V by 3 matrix of vertices
//   newOrientation  desired bounding box that transformed V should be in
//   keepSize        boolean indicating whether Vout should be the same size as V
// Outputs:
//   Vout     scaled and transformed v
void transformVertices(const Eigen::MatrixXd& V, const Eigen::AlignedBox3d &newOrientation,  bool keepSize, Eigen::MatrixXd& Vout);

// keepSize = false by default
void transformVertices(const Eigen::MatrixXd& V, const Eigen::AlignedBox3d &newOrientation, Eigen::MatrixXd& Vout);

//
// Computes a voxel grid enclosing a given mesh
//
// Inputs:
//   V  #V by 3 matrix of vertex positions
// Outputs:
//   centers  x*y*z by 3 matrix of voxel centers
//   corners  (x+1)*(y+1)*(z+1) by 3 matrix of voxel centers
void createVoxelGrid(const Eigen::MatrixXd& V, Eigen::MatrixXd& centers, Eigen::MatrixXd& corners);

//
// Aligns given mesh to the axis planes. Keeps the size the same
//
// Inputs:
//   V         #vertices by 3 matrix of vertices
// Outputs:
//   alignedV  #vertices by 3 matrix of vertices transformed to align to axis planes
void alignToAxis(const Eigen::MatrixXd& V, Eigen::MatrixXd& alignedV);
