#pragma once
#include "SheetGeometry.h"
#include <Geometry>

// A thin wing situated at a relatively small angle of attack with respect to the incoming flow
// (The angle must be small enough as for the "Kutta condition" to hold)
class ThinWing : public SheetGeometry
{
public:

	size_t num_chorwise;
	size_t num_spanwise;

	// Indices of the vertices which build the trailing edge where the wake is shed
	Eigen::VectorXi trailing_edge;
	// Indices of the quads at the trailing edge
	Eigen::VectorXi trailing_panels;

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
	using ChordFx = std::function<std::pair<double, double>(double)>;
	using CamberFx = std::function<double(double, double)>;
	static std::shared_ptr<ThinWing> generate(ChordFx chord_fx, CamberFx camber_fx,
											  size_t num_chordwise, size_t num_spanwise,
											  double span, double chord_scale);

};
