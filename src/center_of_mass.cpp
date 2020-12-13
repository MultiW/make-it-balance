#include "center_of_mass.h"

double center_of_mass(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::Vector3d& center)
{
	int n = V.rows();
	int m = F.rows();

	double mass = 0;
	center.setZero();

	Eigen::RowVector3i triangle;
	Eigen::RowVector3d X0;
	Eigen::RowVector3d X1;
	Eigen::RowVector3d X2;
	Eigen::RowVector3d A;
	Eigen::RowVector3d B;
	Eigen::RowVector3d N;
	double a;
	double b;
	double c;
	double s;
	double area;

	// Compute mass
	for (int i = 0; i < m; i++) // iterate triangles
	{
		triangle = F.row(i);
		X0 = V.row(triangle(0));
		X1 = V.row(triangle(1));
		X2 = V.row(triangle(2));

		// compute area (Heron's formula)
		A = (X0 - X1); // triangle sides
		B = (X0 - X2);
		a = (X0 - X1).norm(); // triangle side length
		b = (X0 - X2).norm();
		c = (X1 - X2).norm();
		s = (a + b + c) / 2.0; // semi-perimeter
		area = std::sqrt(s * (s - a) * (s - b) * (s - c));

		// compute unit normal
		N = A.cross(B).normalized();

		// compute mass
		mass += 2 * area * N(0) * (X0(0) + X1(0) + X2(0)) / 6.0;
	}

	// Compute center of mass
	center_of_mass(V, F, mass, center);

	return mass;
}

void center_of_mass(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double mass, Eigen::Vector3d& center)
{
	int n = V.rows();
	int m = F.rows();

	center.setZero();

	Eigen::RowVector3i triangle;
	Eigen::RowVector3d X0;
	Eigen::RowVector3d X1;
	Eigen::RowVector3d X2;
	Eigen::RowVector3d A;
	Eigen::RowVector3d B;
	Eigen::RowVector3d N;
	double a;
	double b;
	double c;
	double s;
	double area;
	for (int i = 0; i < m; i++) // iterate triangles
	{
		triangle = F.row(i);
		X0 = V.row(triangle(0));
		X1 = V.row(triangle(1));
		X2 = V.row(triangle(2));

		// compute area (Heron's formula)
		A = (X0 - X1); // triangle sides
		B = (X0 - X2);
		a = (X0 - X1).norm(); // triangle side length
		b = (X0 - X2).norm();
		c = (X1 - X2).norm();
		s = (a + b + c) / 2.0; // semi-perimeter
		area = std::sqrt(s * (s - a) * (s - b) * (s - c));

		// compute unit normal
		N = A.cross(B).normalized();

		// compute center of mass
		double N1 = N(0);
		double N2 = N(1);
		double N3 = N(2);
		double X01 = X0(0);
		double X02 = X0(1);
		double X03 = X0(2);
		double X11 = X1(0);
		double X12 = X1(1);
		double X13 = X1(2);
		double X21 = X2(0);
		double X22 = X2(1);
		double X23 = X2(2);
		center(0) += 2 * area * (N1 * (X01 * X11 + X01 * X21 + X11 * X21 + X01 * X01 + X11 * X11 + X21 * X21)) / 2.4E+1;
		center(1) += 2 * area * (N2 * (X02 * X12 + X02 * X22 + X12 * X22 + X02 * X02 + X12 * X12 + X22 * X22)) / 2.4E+1;
		center(2) += 2 * area * (N3 * (X03 * X13 + X03 * X23 + X13 * X23 + X03 * X03 + X13 * X13 + X23 * X23)) / 2.4E+1;
	}

	center /= mass;
}