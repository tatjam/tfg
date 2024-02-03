#pragma once
#include <Eigen>

class Geometry
{
public:
	// 3D position of each vertex
	Eigen::Matrix3Xd vertices;
	// Index list into the previous matrix for each polygon
	// ordered clockwise as seen from above (positive y)
	Eigen::Matrix4Xi edges;


	// Export Mathematica representable string
	std::string to_string();


};
