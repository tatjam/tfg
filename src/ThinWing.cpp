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

	// Place vertices
	for(size_t xi = 0; xi < num_spanwise; xi++)
	{
		double xprog = M_PI * ((double)xi / (double)(num_spanwise - 1));
		double xcos = (1.0 - std::cos(xprog)) / 2.0;
		double x = (xcos * 2.0 - 1.0);
		auto[chord_norm, chord_pos] = chord_fx(x);
		x *= span * 0.5;
		double chord = chord_norm * chord_scale;
		double chord_center = chord_pos * chord_scale;

		for(size_t zi = 0; zi < num_chordwise; zi++)
		{
			double zprog = M_PI * ((double)zi / (double)(num_chordwise - 1));
			double zcos = (1.0 - std::cos(zprog)) / 2.0;
			double y = camber_fx(x, zcos);

			// This ranges from -0.5->0.5
			double z = zcos - 0.5;

			// Make its center be chord_center
			z += chord_center;
			// And scale by chord
			z *= chord;

			// Put the vertices
			out->vertices(0, xi * num_chordwise + zi) = x;
			out->vertices(1, xi * num_chordwise + zi) = y;
			out->vertices(2, xi * num_chordwise + zi) = z;

		}
	}

	// Generate the rectangles
	out->quads = Array4Xi(4, (num_spanwise - 1) * (num_chordwise - 1));

	for(size_t xi = 0; xi < num_spanwise - 1; xi++)
	{
		for(size_t zi = 0; zi < num_chordwise - 1; zi++)
		{
			// We consider x going right and z going down for nomenclature
			// Top left
			out->quads(0, xi * (num_spanwise - 1) + zi) = xi * num_chordwise + zi;
			// Top right
			out->quads(1, xi * (num_chordwise - 1) + zi) = (xi + 1) * num_chordwise + zi;
			// Bottom right
			out->quads(2, xi * (num_chordwise - 1) + zi) = (xi + 1) * num_chordwise + zi + 1;
			// Bottom left
			out->quads(3, xi * (num_chordwise - 1) + zi) = xi * num_chordwise + zi + 1;
		}
	}

	return out;
}
