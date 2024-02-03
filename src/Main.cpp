#include <iostream>
#include <chrono>
#include "Util.h"
#include "Geometry.h"


int main()
{
	auto chord_fx = [](double span_pos)
	{
		double chord_norm = 1.0;
		double chord_line = 0.0;
		return std::make_pair(chord_norm, chord_line);
	};

	auto camber_fx = [](double span_pos, double chord_pos)
	{
		return 0.0;
	};

	auto geom = Geometry::from_chord_and_camber(chord_fx, camber_fx, 20, 20, 10.0, 4.0);
	write_string_to_file("workdir/test.dat", geom->to_string());
}
