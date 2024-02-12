#include "ThinWing.h"

using namespace Eigen;

std::shared_ptr<ThinWing>
ThinWing::from_chord_and_camber(std::function<std::pair<double, double>(double)> chord_fx,
									  std::function<double(double, double)> camber_fx,
									  size_t num_chordwise, size_t num_spanwise, double span, double chord_scale)
{
	auto out = std::make_shared<ThinWing>();

	assert(num_spanwise > 0);
	assert(num_chordwise > 0);

	out->vertices = Array3Xd(3, num_chordwise * num_spanwise);
	out->trailing_edge = VectorXi(num_spanwise);
	out->trailing_panels = VectorXi(num_spanwise - 1);

	Index trail = 0;
	// Place vertices
	for(size_t xi = 0; xi < num_spanwise; xi++)
	{
		double xprog = M_PI * ((double)xi / (double)(num_spanwise - 1));
		double xcos = (1.0 - std::cos(xprog)) / 2.0;
		double x = (xcos * 2.0 - 1.0);
		double uns_x = x;
		auto[chord_norm, chord_pos] = chord_fx(x);
		x *= span * 0.5;
		double chord = chord_norm * chord_scale;
		double chord_center = chord_pos * chord_scale;

		for(size_t zi = 0; zi < num_chordwise; zi++)
		{
			double zprog = M_PI * ((double)zi / (double)(num_chordwise - 1));
			double zcos = (1.0 - std::cos(zprog)) / 2.0;
			double y = camber_fx(uns_x, zcos);

			// This ranges from -0.5->0.5
			double z = zcos - 0.5;

			// And scale by chord
			z *= chord;

			// Make its center be chord_center
			z += chord_center;

			// Put the vertices
			out->vertices(0, xi * num_chordwise + zi) = x;
			out->vertices(1, xi * num_chordwise + zi) = y;
			out->vertices(2, xi * num_chordwise + zi) = z;

			if(zi == num_chordwise - 1)
			{
				out->trailing_edge(trail) = xi * num_chordwise + zi;
				trail++;
			}

		}
	}

	// Generate the rectangles
	out->quads = Array4Xi(4, (num_spanwise - 1) * (num_chordwise - 1));
	out->neighbors = ArrayXXi(8, (num_spanwise - 1) * (num_chordwise - 1));
	out->neighbors.setConstant(-1);

	trail = 0;

	for(size_t xi = 0; xi < num_spanwise - 1; xi++)
	{
		for(size_t zi = 0; zi < num_chordwise - 1; zi++)
		{
			Index idx = xi * (num_chordwise - 1) + zi;
			// We consider x going right and z going down for nomenclature
			// Top left
			out->quads(0, idx) = xi * num_chordwise + zi;
			// Top right
			out->quads(1, idx) = (xi + 1) * num_chordwise + zi;
			// Bottom right
			out->quads(2, idx) = (xi + 1) * num_chordwise + zi + 1;
			// Bottom left
			out->quads(3, idx) = xi * num_chordwise + zi + 1;

			if(zi == num_chordwise - 2)
			{
				out->trailing_panels(trail) = idx;
				trail++;
			}

			// Neighbor assignment
			// Direct neighbors (share an edge) (row < 4)
			if(xi > 0)
				out->neighbors(0, idx) = (xi - 1) * (num_chordwise - 1) + zi;
			if(zi > 0)
				out->neighbors(1, idx) = xi * (num_chordwise - 1) + zi - 1;
			if(xi < num_spanwise - 2)
				out->neighbors(2, idx) = (xi + 1) * (num_chordwise - 1) + zi;
			if(zi < num_chordwise - 2)
				out->neighbors(3, idx) = xi * (num_chordwise - 1) + zi + 1;
			// Diagonal neighbors (share only a vertex) (row >= 4)
			if(xi > 0 && zi > 0)
				out->neighbors(4, idx) = (xi - 1) * (num_chordwise - 1) + zi - 1;
			if(xi > 0 && zi < num_chordwise - 2)
				out->neighbors(5, idx) = (xi - 1) * (num_chordwise - 1) + zi + 1;
			if(zi > 0 && xi < num_spanwise - 2)
				out->neighbors(6, idx) = (xi + 1) * (num_chordwise - 1) + zi - 1;
			if(xi < num_spanwise - 2 && zi < num_chordwise - 2)
				out->neighbors(7, idx) = (xi + 1) * (num_chordwise - 1) + zi + 1;


		}
	}

	return out;
}
