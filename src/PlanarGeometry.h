#pragma once
#include <Eigen>
#include <memory>

// PlanarGeometry for a wing or wake, which is infinitely thin
class PlanarGeometry
{
public:
	// NOTE: Eigen uses column-major by default so 3X is optimal instead of X3 (for large X)
	// 3D position of each vertex
	Eigen::Array3Xd vertices;
	// Index list into the previous matrix for each quadrilateral
	// ordered clockwise as seen from above (positive y)
	Eigen::Array4Xi quads;

	Eigen::Array3Xd normals;

	// Transform matrix
	// WARNING: Make sure to recalculate normals if this changes!
	Eigen::Affine3d transform;

	// Neighbor for each quad, in arbitrary directions
	// (We define neighbor as sharing a vertex)
	// Maximum of 8 neighbors per element, empty neighbors are marked (-1)
	Eigen::ArrayXXi neighbors;


	// Export Mathematica list of vertices that can be readily
	// parsed by Polygon to display in 3D
	// (Use ToExpression[Import["...", "String"]] to import it)
	std::string quads_to_string(bool internal = false);

	// Export Mathematica list of points thta can be readily
	// parsed by Arrow to display in 3D
	// (Use ToExpression[Import["...", "String"]] to import it)
	std::string normals_to_string(double nrm_scale = 0.5);

	void generate_normals();

	PlanarGeometry();

};
