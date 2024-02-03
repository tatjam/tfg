#pragma once
#include <Eigen>
#include <memory>

// Geometry for a wing or wake
class Geometry
{
public:
	// NOTE: Eigen uses column-major by default so 3X is optimal instead of X3 (for large X)
	// 3D position of each vertex
	Eigen::Array3Xd vertices;
	// Index list into the previous matrix for each polygon
	// ordered clockwise as seen from above (positive y)
	Eigen::Array4Xi edges;


	// Export Mathematica list of vertices that can be readily
	// parsed by Polygon to display in 3D
	// (Use ToExpression[Import["...", "String"]] to import it)
	std::string to_string();

	// Generates a surface using cosine sampling
	// chord_fx(span_pos: [-1, 1]) -> pair(chord_norm, chord_line_norm)
	// camber_fx(span_pos: [-1, 1], chord_pos: [0, 1]) -> y offset
	// Coordinates are normalized, so that span positions range from -1.0 to 1.0
	// but chord_norm lengths must be greater than 0, while chord_line_norm may take on any value
	// (chord_line_norm is the location of the middle chord line relative to x = 0)
	// returned chord_norm is scaled by chord for the final geometry
	//
	// span is the wingspan (so both sides), and chord_scale is the scaling factor
	// for returned chord lengths
	static std::shared_ptr<Geometry> from_chord_and_camber(std::function<std::pair<double, double>(double)> chord_fx,
														   std::function<double(double, double)> camber_fx,
							   size_t num_chordwise, size_t num_spanwise, double span, double chord_scale);


};
