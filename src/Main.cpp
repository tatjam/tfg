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
		chord_line = 0.5*std::abs(span_pos);
		chord_norm = 1.0 - 0.5 * std::abs(span_pos);
		return std::make_pair(1.0, 0.0);
		return std::make_pair(chord_norm, chord_line);
	};

	auto camber_fx = [](double span_pos, double chord_pos)
	{
		double x = chord_pos;
		double m = 0.4;
		double p = 0.3;
		double sc = 0.5;
		return 0.0;
		if(std::abs(span_pos) > 0.9)
		{
			return -5.0 * (std::abs(span_pos) - 0.9);
		}
		return 0.0;
		// Camber line
		if(x <= p)
		{
			return -sc * m / (p * p) * (2.0 * p * x - x * x);
		}
		else
		{
			return -sc * m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0*p) + 2.0*p*x - x * x);
		}
	};

	auto geom = ThinWing::from_chord_and_camber(chord_fx, camber_fx, 40, 40, 10.0, 4.0);
	//geom->transform.rotate(Eigen::AngleAxisd(-0.06, Eigen::Vector3d(1, 0, 0))); // Angle of attack
	geom->generate_normals();
	write_string_to_file("workdir/wing.dat", geom->quads_to_string());
	write_string_to_file("workdir/wing_nrm.dat", geom->normals_to_string());

	PanelMethod panels;
	double AoA = 4.0 * M_PI / 180.0;
	double fvel = 10.0;
	panels.body_vel = Eigen::Vector3d(0, fvel * std::sin(AoA), -fvel * std::cos(AoA));
	//panels.body_vel = Eigen::Vector3d(0, 0, -10);
	panels.omega = Eigen::Vector3d(0, 0, 0);

	panels.thin_wings.push_back(geom);
	panels.shed_initial_wake(10, 0.25);

	panels.build_geometry_matrix();

	panels.build_dynamic_steady();

	panels.solve();

	//panels.compute_cps(0.7);
	panels.compute_cps_smart();
	std::cout << panels.compute_aero_force() << std::endl;
	//write_string_to_file("workdir/mat.dat", panels.geometry_matrix_to_string());
	//write_string_to_file("workdir/dyn_mat.dat", panels.dynamic_matrix_to_string());
	write_string_to_file("workdir/sln.dat", panels.solution_to_string(0));
	write_string_to_file("workdir/wake.dat", panels.wake_geom_to_string(0));
	write_string_to_file("workdir/cps.dat", panels.cps_to_string(0));
	write_string_to_file("workdir/flowmap.dat", panels.sample_flow_field_to_string(
			Eigen::Vector3d(0.0, -0.5, -2.5),
			Eigen::Vector3d(0, 0, 1),
			Eigen::Vector3d(0, 1, 0),
			150,150
			));


}
