#include "inner_void_mesh.h"
#include "grid_util.h"
#include "center_of_mass.h"

#include <vector>
#include <set>
#include <limits>
#include <iostream>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
#include <igl/per_face_normals.h>
#include <igl/centroid.h>

const int MAX_ITER = 20;
const int MAX_CARVE = 1000; // max number of voxels to carve per iteration
const double MAX_EMPTY = 0.80; // stop carving when mesh is this empty

int getBinIdx(double min, double max, double binSize, double currLocation)
{
	assert(currLocation >= min && currLocation <= max);

	return std::floor((currLocation - min) / binSize);
}

void getNeighborsIdx(double min, double max, double binSize, double currBorder, double epsilon, std::vector<int> &out)
{
	assert(currBorder >= min - epsilon && currBorder <= max + epsilon);

	int borderIdx = std::round((currBorder - min) / binSize);

	assert(std::abs(currBorder - (min + (borderIdx) * binSize)) < epsilon);

	int maxIdx = std::round((max - min) / binSize);
	if (borderIdx == maxIdx) 
	{
		out[0] = borderIdx - 1;
		out[1] = -1;
	}
	else if (borderIdx == 0)
	{
		out[0] = 0;
		out[1] = -1;
	}
	else
	{
		out[0] = borderIdx - 1;
		out[1] = borderIdx;
	}
}

InnerVoidMesh::InnerVoidMesh(
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
	const Eigen::MatrixXd &planeV, const Eigen::MatrixXi &planeF, 
	const Eigen::Vector3d &gravity, 
	const Eigen::Vector3d &balancePoint): 
	planeV(planeV), planeF(planeF), gravity(gravity), targetCOM(balancePoint.transpose()), balancePoint(balancePoint), innerMesh(List3d<Voxel>()), 
	voxelsDisplaying(0), iterations(0)
{
	// Compute center of mass of initial un-carved object
	this->mass = center_of_mass(V, F, this->currCOM);

	// Define voxel grid
	// - corners is the vertices of the new voxelized mesh
	Eigen::MatrixXd centers;
	createVoxelGrid(V, centers, this->corners, this->dimensions, this->voxelSize);

	Eigen::AlignedBox3d gridBounds;
	createAlignedBox(corners, gridBounds);
	double minX = corners.col(0).minCoeff();
	double maxX = corners.col(0).maxCoeff();
	double minY = corners.col(1).minCoeff();
	double maxY = corners.col(1).maxCoeff();
	double minZ = corners.col(2).minCoeff();
	double maxZ = corners.col(2).maxCoeff();

	// Compute voxel distances to mesh
	Eigen::VectorXd distances, closestFaces;
	Eigen::MatrixXd closestPoints, closestNormals;
	igl::signed_distance(centers, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, 
		std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), distances, 
		closestFaces, closestPoints, closestNormals);

	// Initialize inner mesh list
	innerMesh.resize(this->dimensions(0), this->dimensions(1), this->dimensions(2));

	// iterate centers
	Eigen::RowVector3d currCenter;
	int xIdx, yIdx, zIdx;
	for (int i = 0; i < centers.rows(); i++)
	{
		currCenter = centers.row(i);
		xIdx = getBinIdx(minX, maxX, this->voxelSize(0), currCenter(0));
		yIdx = getBinIdx(minY, maxY, this->voxelSize(1), currCenter(1));
		zIdx = getBinIdx(minZ, maxZ, this->voxelSize(2), currCenter(2));

		innerMesh(xIdx, yIdx, zIdx).center = currCenter;
		innerMesh(xIdx, yIdx, zIdx).xIdx = xIdx;
		innerMesh(xIdx, yIdx, zIdx).yIdx = yIdx;
		innerMesh(xIdx, yIdx, zIdx).zIdx = zIdx;
		innerMesh(xIdx, yIdx, zIdx).distanceToMesh = distances(i);

		// label inner voxels (inside of mesh + padding)
		if (distances(i) < 0.0 && std::abs(distances(i)) >= voxelSize(0))
		{
			// make border voxels don't count as inner voxels
			if (xIdx > 0 && xIdx < dimensions(0) - 1 &&
				yIdx > 0 && yIdx < dimensions(1) - 1 &&
				zIdx > 0 && zIdx < dimensions(2) - 1)
			{
				innerMesh(xIdx, yIdx, zIdx).isInner = true;
			}
		}
	}

	// iterate corners
	Eigen::RowVector3d currCorner;
	double epsilon = voxelSize(0) / 10.0;
	std::vector<int> xNeighs(2);
	std::vector<int> yNeighs(2);
	std::vector<int> zNeighs(2);
	for (int rowIdx = 0; rowIdx < corners.rows(); rowIdx++)
	{
		currCorner.setZero();
		currCorner = corners.row(rowIdx);
		getNeighborsIdx(minX, maxX, voxelSize(0), currCorner(0), epsilon, xNeighs);
		getNeighborsIdx(minY, maxY, voxelSize(1), currCorner(1), epsilon, yNeighs);
		getNeighborsIdx(minZ, maxZ, voxelSize(2), currCorner(2), epsilon, zNeighs);

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					if (xNeighs[i] != -1 && yNeighs[j] != -1 && zNeighs[k] != -1)
					{
						innerMesh(xNeighs[i], yNeighs[j], zNeighs[k]).cornerIndices.insert(rowIdx);
					}
				}
			}
		}
	}
}

// Should voxels show on the screen (represents the void space)
bool shouldDisplayVoxel(const Voxel voxel)
{
	// isFilled is false means this voxel is a void space.
	// Remember voxels represents void space
	return voxel.isInner && !voxel.isFilled;
}

bool areFacingSameDirection(const Eigen::RowVector3d& v1, const Eigen::RowVector3d& v2)
{
	return v1.dot(v2) > 0;
}

// Given two neighboring voxels, triangulate their common side and append to the list of faces "outF"
void appendCommonSide(
	const Voxel &voxel1, const Voxel &voxel2, 
	const Eigen::MatrixXd &V, Eigen::MatrixXi &outF, 
	bool onlyNecessary=true)
{
	Voxel insideVoxel;
	if (onlyNecessary)
	{
		// No need to display side of two displaying voxels
		// No need to display side of two empty voxels
		if ((shouldDisplayVoxel(voxel1) && shouldDisplayVoxel(voxel2)) ||
			(!shouldDisplayVoxel(voxel1) && !shouldDisplayVoxel(voxel2)))
		{
			return;
		}
		insideVoxel = shouldDisplayVoxel(voxel1) ? voxel1 : voxel2;
	}
	else
	{
		insideVoxel = voxel1;
	}

	// Find common points
	std::vector<int> intersect(0);
	std::set_intersection(
		voxel1.cornerIndices.begin(), voxel1.cornerIndices.end(), 
		voxel2.cornerIndices.begin(), voxel2.cornerIndices.end(),
		std::inserter(intersect, intersect.begin()));

	assert(intersect.size() == 4);

	// Triangulate
	Eigen::MatrixXi sideF;
	sideF.resize(2, 3);
	sideF <<
		Eigen::RowVector3i(intersect[0], intersect[1], intersect[2]),
		Eigen::RowVector3i(intersect[2], intersect[1], intersect[3]);

	// assertion that triangles face the same direction
	Eigen::MatrixXd N;
	igl::per_face_normals(V, sideF, N);
	assert(areFacingSameDirection(N.row(0), N.row(1)));

	// Fix direction of faces. Voxels should be inside out (to represent void space)
	Eigen::RowVector3d vectToSide;
	vectToSide = V.row(sideF(0)) - insideVoxel.center.transpose();
	if (areFacingSameDirection(N.row(0), vectToSide))
	{
		sideF << 
			Eigen::RowVector3i(intersect[2], intersect[1], intersect[0]), 
			Eigen::RowVector3i(intersect[3], intersect[1], intersect[2]);
	}

	outF.conservativeResize(outF.rows() + 2, outF.cols());
	outF.block(outF.rows() - 2, 0, 2, 3) = sideF;
}

void InnerVoidMesh::convertToMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) 
{
	V.resize(corners.rows(), corners.cols());
	V = corners;

	F.resize(0, 3);
	Eigen::MatrixXi currF;
	Voxel currVoxel;
	for (int i = 1; i < dimensions(0)-1; i++)
	{
		for (int j = 1; j < dimensions(1)-1; j++)
		{
			for (int k = 1; k < dimensions(2)-1; k++)
			{
				currF.resize(0, 3);

				currVoxel = innerMesh(i, j, k);

				if (!shouldDisplayVoxel(currVoxel))
				{
					continue;
				}

				appendCommonSide(innerMesh(i - 1, j, k), currVoxel, corners, currF);
				appendCommonSide(innerMesh(i + 1, j, k), currVoxel, corners, currF);
				appendCommonSide(innerMesh(i, j - 1, k), currVoxel, corners, currF);
				appendCommonSide(innerMesh(i, j + 1, k), currVoxel, corners, currF);
				appendCommonSide(innerMesh(i, j, k - 1), currVoxel, corners, currF);
				appendCommonSide(innerMesh(i, j, k + 1), currVoxel, corners, currF);
				
				// Add faces to output
				int newFacesCount = currF.rows();
				F.conservativeResize(F.rows() + newFacesCount, F.cols());
				F.block(F.rows() - newFacesCount, 0, newFacesCount, currF.cols()) = currF;
			}
		}
	}
}

// Energy function to reduce
// Find the distance between target and current center of mass projected onto the ground
double InnerVoidMesh::getCoMDistance()
{
	Eigen::MatrixXd comPoints;
	comPoints.resize(2, 3);
	comPoints <<
		this->targetCOM.transpose(),
		this->currCOM.transpose();

	// get distance from each COM to ground
	Eigen::VectorXd distances, closestFaces;
	Eigen::MatrixXd closestPoints, closestNormals;
	igl::signed_distance(comPoints, this->planeV, this->planeF, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, 
		std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), distances, 
		closestFaces, closestPoints, closestNormals);

	// distance between two CoM points
	double comDistance = (comPoints.row(0) - comPoints.row(1)).squaredNorm();
	// difference in height between CoM points
	double comHeightDiff = std::abs(distances(1) - distances(0));

	// projected distance using pythagorean theorem
	double term1 = std::pow(comDistance, 2);
	double term2 = std::pow(comHeightDiff, 2);
	return std::sqrt(std::max(term1, term2) - std::min(term1, term2));
}

bool InnerVoidMesh::isOptimized()
{
	// Stop carving when:
	// - CoM is within a voxel away from the point of balance
	// - too many voxels were removed
	// - too many iterations have passed
	double percentDisplaying = this->voxelsDisplaying / (double) this->dimensions.prod();
	return getCoMDistance() < this->voxelSize.maxCoeff() || percentDisplaying > MAX_EMPTY || this->iterations > MAX_ITER;
}

void flipFaces(Eigen::MatrixXi& F)
{
	for (int i = 0; i < F.rows(); i++)
	{
		int x = F(i, 0);
		F(i, 0) = F(i, 2);
		F(i, 2) = x;
	}
}

// Add or remove the side between the two voxels to the global com
void updateCenterOfMass(
	const Voxel &insertVoxel, const Voxel &neighbor, const Eigen::MatrixXd &V, double mass, 
	Eigen::Vector3d &globalCom)
{
	// Get common side mesh
	Eigen::MatrixXi F;
	F.resize(0, 3);
	appendCommonSide(insertVoxel, neighbor, V, F, true);

	// Get center of mass for current side
	Eigen::Vector3d CoM;
	center_of_mass(V, F, mass, CoM);

	// Case 1: remove this side
	if (shouldDisplayVoxel(neighbor))
	{
		globalCom -= CoM;
	}
	// Case 2: add new side
	else
	{
		globalCom += CoM;
	}
}

void InnerVoidMesh::carveInnerMeshStep(Eigen::MatrixXd &carvePlaneV, Eigen::MatrixXi &carvePlaneF, Eigen::Vector3d& com)
{
	this->iterations++;

	// plane perpendicular to the cutting plane
	Eigen::MatrixXd cuttingPlaneTV(3,3);
	Eigen::MatrixXi cuttingPlaneTF(1,3);
	cuttingPlaneTV <<
		this->targetCOM.transpose(),
		this->currCOM.transpose(),
		(this->targetCOM + this->gravity).transpose();
	cuttingPlaneTF << 0, 1, 2;

	Eigen::MatrixXd N;
	igl::per_face_normals(cuttingPlaneTV, cuttingPlaneTF, N);

	// Define cutting plane (make large enough to cover object)
	double gravityMultiplier = (this->dimensions.maxCoeff() * this->voxelSize.maxCoeff() * 10) / this->gravity.squaredNorm();
	double otherVectorMultiplier = (this->dimensions.maxCoeff() * this->voxelSize.maxCoeff() * 10) / N.row(0).squaredNorm();

	carvePlaneV.resize(5, 3);
	carvePlaneV <<
		this->targetCOM.transpose(),
		(this->targetCOM + this->gravity * gravityMultiplier).transpose(),
		this->targetCOM.transpose() + N.row(0) * otherVectorMultiplier,
		(this->targetCOM - this->gravity * gravityMultiplier).transpose(),
		this->targetCOM.transpose() - N.row(0) * otherVectorMultiplier;
	carvePlaneF.resize(4, 3);
	carvePlaneF <<
		0, 1, 2,
		2, 3, 0,
		0, 3, 4,
		4, 1, 0;

	N.resize(0, 0);
	igl::per_face_normals(carvePlaneV, carvePlaneF, N);
	assert(areFacingSameDirection(N.row(0), N.row(1)));
	assert(areFacingSameDirection(N.row(1), N.row(2)));
	assert(areFacingSameDirection(N.row(2), N.row(3)));

	// make sure plane faces the current center of mass
	Eigen::VectorXd distances, closestFaces;
	Eigen::MatrixXd closestPoints, closestNormals;
	igl::signed_distance(this->currCOM.transpose(), carvePlaneV, carvePlaneF, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, 
		std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), distances, 
		closestFaces, closestPoints, closestNormals);
	if (distances(0) < 0.0)
	{
		flipFaces(carvePlaneF);
	}

	// Carve inner mesh
	// sort voxels by distance to cutting plane
	std::vector<std::pair<Voxel*, double>> voxelList; // Voxel, distance-to-plane pair
	Voxel *currVoxel;
	for (int i = 1; i < dimensions(0) - 1; i++)
	{
		for (int j = 1; j < dimensions(1) - 1; j++)
		{
			for (int k = 1; k < dimensions(2) - 1; k++)
			{
				currVoxel = &innerMesh(i, j, k);
				if (!currVoxel->isInner || !currVoxel->isFilled)
				{
					continue;
				}

				Eigen::VectorXd distances, closestFaces;
				Eigen::MatrixXd closestPoints, closestNormals;
				igl::signed_distance(currVoxel->center.transpose(), carvePlaneV, carvePlaneF, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, 
					std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), distances, 
					closestFaces, closestPoints, closestNormals);

				if (distances(0) > 0)
				{
					voxelList.push_back(std::make_pair(currVoxel, distances(0)));
				}
			}
		}
	}

	// sort by distance
	struct compareDistance
	{
		inline bool operator() (const std::pair<Voxel*, double> pair1, const std::pair<Voxel*, double> pair2)
		{
			return (pair1.first < pair2.first);
		}
	};
	std::sort(voxelList.begin(), voxelList.end(), compareDistance());

	// carve inner mesh (starting with the furthest from the plane)
	for (int i = voxelList.size() - 1; i >= std::max(0, (int) voxelList.size() - MAX_CARVE); i--)
	{
		Voxel* currVoxel = voxelList[i].first;
		int x = currVoxel->xIdx;
		int y = currVoxel->yIdx;
		int z = currVoxel->zIdx;

		// Add voxel
		currVoxel->isFilled = false;
		this->voxelsDisplaying++;

		// Update current center of mass
		updateCenterOfMass(*currVoxel, this->innerMesh(x + 1, y, z), this->corners, this->mass, this->currCOM);
		updateCenterOfMass(*currVoxel, this->innerMesh(x - 1, y, z), this->corners, this->mass, this->currCOM);
		updateCenterOfMass(*currVoxel, this->innerMesh(x, y + 1, z), this->corners, this->mass, this->currCOM);
		updateCenterOfMass(*currVoxel, this->innerMesh(x, y - 1, z), this->corners, this->mass, this->currCOM);
		updateCenterOfMass(*currVoxel, this->innerMesh(x, y, z + 1), this->corners, this->mass, this->currCOM);
		updateCenterOfMass(*currVoxel, this->innerMesh(x, y, z - 1), this->corners, this->mass, this->currCOM);
	}

	// return current center of mass
	com = this->currCOM;
}


void InnerVoidMesh::carveInnerMesh()
{
	while(!isOptimized())
	{
		Eigen::MatrixXd tempV;
		Eigen::MatrixXi tempF;
		Eigen::Vector3d tempCom;
		carveInnerMeshStep(tempV, tempF, tempCom);
	}
}
