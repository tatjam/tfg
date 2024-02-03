#include <iostream>
#include <chrono>
#include "Util.h"
#include "PanelMethod.h"


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

	auto geom = ThinWing::from_chord_and_camber(chord_fx, camber_fx, 20, 20, 10.0, 4.0);
	geom->generate_normals();
	write_string_to_file("workdir/wing.dat", geom->quads_to_string());
	write_string_to_file("workdir/wing_nrm.dat", geom->normals_to_string());

	PanelMethod panels;

	panels.thin_wings.push_back(geom);

	panels.build_geometry_matrix();


}
