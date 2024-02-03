#include "PlanarGeometry.h"

using namespace Eigen;


std::string PlanarGeometry::quads_to_string()
{
	// Big array of quads
	std::stringstream out;
	out << std::fixed;
	out << "{";

	// Iterate over each polygon
	for(Index i = 0; i < quads.cols(); i++)
	{
		out << "{";


		// Put each 3D vertex
		for(Index j = 0; j < 4; j++)
		{
			out << "{";

			Eigen::Vector3d vert = transform * vertices.col(quads(j,i));


			out << vert(0);
			out << ",";
			out << vert(1);
			out << ",";
			out << vert(2);

			out << "}";
			if(j != 3)
				out << ",";
		}

		out << "}";
		if(i != quads.cols() - 1)
			out << ",";
	}

	out << "}";
	return out.str();
}


std::string PlanarGeometry::normals_to_string(double scl)
{
	std::stringstream out;
	out << std::fixed;
	out << "{";

	// Iterate over each polygon
	for(Index i = 0; i < quads.cols(); i++)
	{
		// Put polygon center point
		Vector3d center = transform * vertices.col(quads(0, i)).matrix();
		center += transform * vertices.col(quads(1, i)).matrix();
		center += transform * vertices.col(quads(2, i)).matrix();
		center += transform * vertices.col(quads(3, i)).matrix();
		center /= 4.0;

		// Put normal end point
		Vector3d endpoint = normals.col(i) * scl;
		endpoint += center;

		out << "Arrow[{{" << center(0) << "," << center(1) << "," << center(2) << "}, ";

		out << "{" << endpoint(0) << "," << endpoint(1) << "," << endpoint(2) << "}}]";

		if(i != quads.cols() - 1)
			out << ",";
	}

	out << "}";
	return out.str();
}

void PlanarGeometry::generate_normals()
{
	normals = Array3Xd(3, quads.cols());
	for(Index i = 0; i < quads.cols(); i++)
	{
		Vector3d v1 = transform * vertices.col(quads(0, i));
		Vector3d v2 = transform * vertices.col(quads(1, i));
		Vector3d v3 = transform * vertices.col(quads(2, i));
		Vector3d v4 = transform * vertices.col(quads(3, i));


		normals.col(i) = (v3 - v1).cross(v2 - v4).normalized();
		Vector3d nrm = v4-v1;
		Vector3d nrm2 = v3-v2;
		nrm = nrm;
	}

}

PlanarGeometry::PlanarGeometry()
{
	transform.setIdentity();

}
