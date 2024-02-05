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

	auto geom = ThinWing::from_chord_and_camber(chord_fx, camber_fx, 10, 10, 10.0, 4.0);
	//geom->transform.translate(Eigen::Vector3d(-5, 0, 0));
	geom->generate_normals();
	write_string_to_file("workdir/wing.dat", geom->quads_to_string());
	write_string_to_file("workdir/wing_nrm.dat", geom->normals_to_string());

	PanelMethod panels;
	//panels.body_vel = Eigen::Vector3d(0, 1, -10);
	panels.body_vel = Eigen::Vector3d(2, 3, -10);
	panels.omega = Eigen::Vector3d(1, 1, 1);

	panels.thin_wings.push_back(geom);
	panels.shed_initial_wake(100, 0.1);

	panels.build_geometry_matrix();

	panels.build_dynamic();

	panels.solve();

	write_string_to_file("workdir/mat.dat", panels.geometry_matrix_to_string());
	write_string_to_file("workdir/sln.dat", panels.solution_to_string(0));
	write_string_to_file("workdir/wake.dat", panels.wake_geom_to_string(0));

}
