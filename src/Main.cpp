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
		/*if(std::abs(span_pos) > 0.5)
		{
			return -5.0 * (std::abs(span_pos) - 0.9);
		}*/
		//return 0.0;
		// Camber line
		/*if(x <= p)
		{
			return -sc * m / (p * p) * (2.0 * p * x - x * x);
		}
		else
		{
			return -sc * m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0*p) + 2.0*p*x - x * x);
		}
		 */
	};

	auto geom = ThinWing::from_chord_and_camber(chord_fx, camber_fx, 16, 16, 1.0, 0.25);
	//geom->transform.rotate(Eigen::AngleAxisd(-0.06, Eigen::Vector3d(1, 0, 0))); // Angle of attack
	geom->generate_normals();
	write_string_to_file("workdir/wing.dat", geom->quads_to_string());
	write_string_to_file("workdir/wing_nrm.dat", geom->normals_to_string());

	PanelMethod panels;
	PanelMethod panels2;
	double AoA = 1.0 * M_PI / 180.0;
	double fvel = 0.001;

	size_t NPANELS = 200;
	double TSTEP = 0.01;

	Eigen::Vector3d bvel = Eigen::Vector3d(0, fvel * std::sin(AoA), -fvel * std::cos(AoA));
	Eigen::Vector3d omega = Eigen::Vector3d(0, 0, 0);

	panels.thin_wings.push_back(geom);
	panels.backward_difference_order = 1;
	panels.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
	panels.build_geometry_matrix();

	panels2.thin_wings.push_back(geom);
	panels2.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
	panels2.build_geometry_matrix();
	panels2.build_dynamic(true);
	panels2.solve(true);

	// Initial solution that will be convected in the wake
	panels.build_dynamic(true);
	panels.solve(true);

	panels.transfer_solution_to_wake();

	double t = 0.0;
	int MAX_IT = NPANELS;
	Eigen::ArrayXd forces = Eigen::ArrayXd(MAX_IT);
	Eigen::ArrayXd steady_forces = Eigen::ArrayXd(MAX_IT);
	for(size_t i = 0; i < MAX_IT; i++)
	{
		AoA = 1.0 * M_PI / 180.0;
		fvel = 1.0;
		bvel = Eigen::Vector3d(0, fvel * std::sin(AoA), -fvel * std::cos(AoA));
		omega = Eigen::Vector3d(0, 0, 0);
		panels.timestep(bvel, omega);
		std::cout << "It: " << i + 1 << " / " << MAX_IT << std::endl;
		panels.compute_cps(false);
		forces(i) = panels.compute_aero_force(true)(1);

		panels2.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
		panels2.build_dynamic(true);
		panels2.solve(true);
		panels2.compute_cps(true);
		steady_forces(i) = panels2.compute_aero_force(true)(1);
		t += TSTEP;

	}
	std::stringstream s;
	s << forces.transpose();
	write_string_to_file("workdir/lift_over_time.dat", s.str());

	std::stringstream s2;
	s2 << steady_forces.transpose();
	write_string_to_file("workdir/lift_over_time_steady.dat", s2.str());

	//write_string_to_file("workdir/mat.dat", panels.geometry_matrix_to_string());
	//write_string_to_file("workdir/dyn_mat.dat", panels.dynamic_matrix_to_string());
	write_string_to_file("workdir/sln.dat", panels.solution_to_string(0));
	write_string_to_file("workdir/wake_sln.dat", panels.wake_solution_to_string(0));
	write_string_to_file("workdir/wake.dat", panels.wake_geom_to_string(0));
	write_string_to_file("workdir/cps.dat", panels.cps_to_string(0));
	/*write_string_to_file("workdir/flowmap.dat", panels.sample_flow_field_to_string(
			Eigen::Vector3d(0.0, -0.5, -2.5),
			Eigen::Vector3d(0, 0, 1),
			Eigen::Vector3d(0, 1, 0),
			150,150
			));*/


}
