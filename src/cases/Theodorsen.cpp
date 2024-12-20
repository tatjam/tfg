#include <iostream>
#include <chrono>
#include "../Util.h"
#include "../PanelMethod.h"

using namespace Eigen;

const double ADVANCE_VEL = 1.0;
const size_t LONG_TERM_NPANELS = 200;
const double LONG_TERM_TSTEP = 0.1;
const double LONG_TERM_AOA = 0.1;

const size_t NPANELS = 50;
const double TSTEP = 0.03;
const size_t MAX_IT = NPANELS * 2;

const double HEAVE_AMPL = 0.01;
const double HEAVE_OMEGA = 4.3;

const double CHORD = 1.0;

const bool USE_CENTERLINE = false;

void write_params()
{
	std::stringstream s;
	s << ADVANCE_VEL << std::endl;
	s << TSTEP << std::endl;
	s << MAX_IT << std::endl;
	s << HEAVE_AMPL << std::endl;
	s << HEAVE_OMEGA << std::endl;
	s << CHORD << std::endl;
	write_string_to_file("workdir/theodorsen_params.dat", s.str());
}

void long_term_steady_case(std::shared_ptr<ThinWing> wing)
{
	Vector3d bvel = Vector3d(0, ADVANCE_VEL * std::sin(LONG_TERM_AOA), -ADVANCE_VEL * std::cos(LONG_TERM_AOA));
	Vector3d omega = Vector3d(0, 0, 0);

	PanelMethod pan;
	pan.thin_wings.push_back(wing);
	pan.shed_initial_wake(LONG_TERM_NPANELS, LONG_TERM_TSTEP, bvel, omega);
	pan.build_geometry_matrix();
	pan.build_dynamic(true);
	pan.solve(true);
	pan.compute_cps(true);

	Vector3d aero_force = pan.compute_aero_force(USE_CENTERLINE);

	{
		std::stringstream s;
		s << aero_force.transpose();
		write_string_to_file("workdir/theodorsen_long_term_forces.dat", s.str());
	}

	std::cout << "Wrote long term steady solution for correction factor" << std::endl;
}

void prepare_oscillating_case(std::shared_ptr<ThinWing> wing, PanelMethod& steady, PanelMethod& dynamic)
{
	dynamic.thin_wings.push_back(wing);
	dynamic.backward_difference_order = 1;

	Vector3d bvel = Vector3d(0, 0.0, -ADVANCE_VEL);
	Vector3d omega = Vector3d(0, 0, 0);

	// Start of method timstep with steady solve
	dynamic.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
	dynamic.build_geometry_matrix();
	dynamic.build_dynamic(true);
	dynamic.solve(true);
	dynamic.transfer_solution_to_wake();
	// Initial phi computation so first time-step can compute cps
	dynamic.compute_phis();

	steady.thin_wings.push_back(wing);
	steady.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
	steady.build_geometry_matrix();

}

void iterate_oscillating_case(PanelMethod& steady, PanelMethod& dynamic, bool write_results, double& t,
							  std::stringstream& bvels)
{
	Eigen::ArrayXd forces = Eigen::ArrayXd(MAX_IT);
	Eigen::ArrayXd steady_forces = Eigen::ArrayXd(MAX_IT);
	Eigen::ArrayXd alpha = Eigen::ArrayXd(MAX_IT);


	for (size_t i = 0; i < MAX_IT; i++)
	{
		double heave = -HEAVE_AMPL * ADVANCE_VEL * std::cos(HEAVE_OMEGA * t);
		Vector3d bvel = Eigen::Vector3d(0, heave, -ADVANCE_VEL);
		Vector3d omega = Eigen::Vector3d(0.0, 0.0, 0.0);
		dynamic.timestep(bvel, omega);
		bvels << "{{" << -bvel(2) << ", 0.0, " << bvel(1) <<"}, {0, 0, 0}}" << std::endl;
		std::cout << "It: " << i + 1 << " / " << MAX_IT << " (t = " << t << ")" << std::endl;
		if(write_results)
		{
			dynamic.compute_cps(false);
			forces(i) = dynamic.compute_aero_force(USE_CENTERLINE)(1);

			// Long steady wake!
			steady.shed_initial_wake(NPANELS, TSTEP * 8.0, bvel, omega);
			steady.build_dynamic(true);
			steady.solve(true);
			steady.compute_cps(true);
			steady_forces(i) = steady.compute_aero_force(USE_CENTERLINE)(1);

			alpha(i) = std::atan2(heave, ADVANCE_VEL);
		}

		t += TSTEP;
	}

	{
	}

	if(write_results)
	{
		{
			std::stringstream s;
			s << forces.transpose();
			write_string_to_file("workdir/theodorsen_lift_over_time.dat", s.str());
		}

		{
			std::stringstream s;
			s << alpha.transpose();
			write_string_to_file("workdir/theodorsen_alpha_over_time.dat", s.str());
		}

		{
			std::stringstream s;
			s << steady_forces.transpose();
			write_string_to_file("workdir/theodorsen_lift_over_time_steady.dat", s.str());
		}

		std::cout << "Results written to disk" << std::endl;
	}
}

int main()
{
	auto chord_fx = [](double span_pos)
	{
		return std::make_pair(1.0, 0.0);
	};

	auto camber_fx = [](double span_pos, double chord_pos)
	{
		return 0.0;
	};

	auto geom = ThinWing::generate(chord_fx, camber_fx, 24, 24, 50.0, CHORD);
	geom->generate_normals();
	write_string_to_file("workdir/geom.dat", geom->quads_to_string());

	write_params();
	long_term_steady_case(geom);

	PanelMethod steady;
	PanelMethod dynamic;
	double t = 0.0;

	std::stringstream bvels;
	prepare_oscillating_case(geom, steady, dynamic);
	auto start = std::chrono::system_clock::now();
	iterate_oscillating_case(steady, dynamic, true, t, bvels);
	auto end = std::chrono::system_clock::now();
	std::cout << "Pre-computation took: " <<
			  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	// Iterativo + sparse: 23929ms
	// Iterativo: 19817ms
	// Normal: 20896
	// Conclusion: bottle-neck is matrix building code (dynamic matrix)
	iterate_oscillating_case(steady, dynamic, true, t, bvels);

	write_string_to_file("workdir/theodorsen_bvel.dat", bvels.str());

	write_string_to_file("workdir/mat.dat", steady.geometry_matrix_to_string());
	write_string_to_file("workdir/phis_above.dat", dynamic.phis_to_string(0, false, true));
	write_string_to_file("workdir/phis_below.dat", dynamic.phis_to_string(0, false, false));

}
