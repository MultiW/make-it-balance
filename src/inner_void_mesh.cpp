#include "inner_void_mesh.h"
#include "grid_util.h"

#include <vector>
#include <set>
#include <limits>
#include <iostream>
#include <tuple>
#include <math.h>
#include <algorithm>

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>



Voxel::Voxel(): isFilled(true), isInner(false) {}

int getBinIdx(double min, double max, double binSize, double currLocation)
{
	assert(currLocation >= min && currLocation <= max);

	return std::floor((currLocation - min) / binSize);
}

void getNeighborsIdx(double min, double max, double binSize, double currBorder, double epsilon, std::vector<int> &out)
{
	assert(currBorder >= min - epsilon && currBorder <= max + epsilon);

	int borderIdx = std::round((currBorder - min) / binSize);
	assert(std::abs(currBorder - (min + borderIdx * binSize)) < epsilon);

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
	const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &F, 
	const Eigen::Vector3d &gravity, 
	const Eigen::Vector3d &balancePoint): V(V), F(F), gravity(gravity), balancePoint(balancePoint), innerMesh(List3d<Voxel>())
{
	// Define voxel grid
	// - corners is the vertices of the new voxelized mesh
	Eigen::MatrixXd centers;
	Eigen::Vector3d voxelSize;
	createVoxelGrid(V, centers, corners, dimensions, voxelSize);

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
	innerMesh.resize(dimensions(0), dimensions(1), dimensions(2));

	// iterate centers
	Eigen::RowVector3d currCenter;
	int xIdx, yIdx, zIdx;
	for (int i = 0; i < centers.rows(); i++)
	{
		currCenter = centers.row(i);
		xIdx = getBinIdx(minX, maxX, voxelSize(0), currCenter(0));
		yIdx = getBinIdx(minY, maxY, voxelSize(1), currCenter(1));
		zIdx = getBinIdx(minZ, maxZ, voxelSize(2), currCenter(2));

		innerMesh(xIdx, yIdx, zIdx).center = currCenter;
		innerMesh(xIdx, yIdx, zIdx).xIdx = xIdx;
		innerMesh(xIdx, yIdx, zIdx).yIdx = yIdx;
		innerMesh(xIdx, yIdx, zIdx).zIdx = zIdx;
		innerMesh(xIdx, yIdx, zIdx).distanceToMesh = distances(i);
		innerMesh(xIdx, yIdx, zIdx).closestNormal = closestNormals.row(i).transpose();

		// label inner voxels (inside of mesh + padding)
		if (distances(i) < 0 && std::abs(distances(i)) >= voxelSize(0))
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
	double epsilon = voxelSize(0) / 100.0;
	std::vector<int> xNeighs(2);
	std::vector<int> yNeighs(2);
	std::vector<int> zNeighs(2);
	for (int i = 0; i < corners.rows(); i++)
	{
		currCorner = corners.row(i);
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
						innerMesh(xNeighs[i], yNeighs[j], zNeighs[k]).cornerIndices.insert(i);
					}
				}
			}
		}
	}
}

// Should voxels show on the screen (represents the void space)
bool shouldDisplayVoxel(const Voxel &voxel)
{
	return voxel.isInner && voxel.isFilled;
}

void addCommonSide(const Voxel &voxel1, const Voxel &voxel2, const Eigen::MatrixXd &V, Eigen::MatrixXi &outF)
{
	// No need to display side of two filled voxels
	if (shouldDisplayVoxel(voxel1) && shouldDisplayVoxel(voxel2))
	{
		return;
	}

	// Find common points
	std::vector<int> intersect;
	std::set_intersection(
		voxel1.cornerIndices.begin(), voxel1.cornerIndices.end(), 
		voxel2.cornerIndices.begin(), voxel2.cornerIndices.begin(),
		std::inserter(intersect, intersect.begin()));

	assert(intersect.size() == 4);

	outF.conservativeResize(outF.rows() + 2, outF.cols());
	outF.row(outF.rows() - 2) << intersect[0], intersect[1], intersect[2];
	outF.row(outF.rows() - 1) << intersect[2], intersect[1], intersect[3];
}

void InnerVoidMesh::convertToMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) 
{
	V = corners;

	F.resize(0, 3);
	Eigen::MatrixXi currF;
	for (int i = 1; i < dimensions(0)-1; i++)
	{
		for (int j = 1; j < dimensions(1)-1; j++)
		{
			for (int k = 1; k < dimensions(k)-1; k++)
			{
				currF.resize(0, 3);

				Voxel currVoxel = innerMesh(i, j, k);
				if (!shouldDisplayVoxel(currVoxel))
				{
					continue;
				}

				addCommonSide(innerMesh(i - 1, j, k), currVoxel, corners, currF);
				addCommonSide(innerMesh(i + 1, j, k), currVoxel, corners, currF);
				addCommonSide(innerMesh(i, j - 1, k), currVoxel, corners, currF);
				addCommonSide(innerMesh(i, j + 1, k), currVoxel, corners, currF);
				addCommonSide(innerMesh(i, j, k - 1), currVoxel, corners, currF);
				addCommonSide(innerMesh(i, j, k + 1), currVoxel, corners, currF);
				
				// Add faces to output
				int newFacesCount = currF.rows();
				F.conservativeResize(F.rows() + newFacesCount, F.cols());
				F.block(F.rows() - newFacesCount, 0, newFacesCount, 3) = currF;
			}
		}
	}
}
