#include <iostream>
#include <chrono>
#include "../Util.h"
#include "../PanelMethod.h"

using namespace Eigen;

const double ADVANCE_VEL = 1.0;
const size_t LONG_TERM_NPANELS = 200;
const double LONG_TERM_AOA = 0.1;

const size_t NPANELS = 50;
const double TSTEP = 0.025;
const size_t MAX_IT = NPANELS * 4;

const bool USE_CENTERLINE = true;

void write_params()
{
	std::stringstream s;
	s << ADVANCE_VEL << std::endl;
	s << LONG_TERM_AOA << std::endl;
	s << TSTEP << std::endl;
	s << MAX_IT << std::endl;
	write_string_to_file("workdir/wagner_params.dat", s.str());
}

void long_term_steady_case(std::shared_ptr<ThinWing> wing)
{
	Vector3d bvel = Vector3d(0, ADVANCE_VEL * std::sin(LONG_TERM_AOA), -ADVANCE_VEL * std::cos(LONG_TERM_AOA));
	Vector3d omega = Vector3d(0, 0, 0);

	PanelMethod pan;
	pan.thin_wings.push_back(wing);
	pan.shed_initial_wake(LONG_TERM_NPANELS, TSTEP, bvel, omega);
	pan.build_geometry_matrix();
	pan.build_dynamic(true);
	pan.solve(true);
	pan.compute_cps(true);

	Vector3d aero_force = pan.compute_aero_force(USE_CENTERLINE);

	{
		std::stringstream s;
		s << aero_force.transpose();
		write_string_to_file("workdir/wagner_long_term_forces.dat", s.str());
	}

	std::cout << "Wrote long term steady solution for correction factor" << std::endl;
}

void prepare_wagner_case(std::shared_ptr<ThinWing> wing, PanelMethod& dynamic)
{
	dynamic.thin_wings.push_back(wing);
	dynamic.backward_difference_order = 1;

	Vector3d bvel = Vector3d(0, ADVANCE_VEL * std::sin(LONG_TERM_AOA), -ADVANCE_VEL * std::cos(LONG_TERM_AOA));
	Vector3d omega = Vector3d(0, 0, 0);

	// Start of method timestep with steady solve
	dynamic.shed_initial_wake(NPANELS, TSTEP, bvel, omega);
	dynamic.build_geometry_matrix();
	dynamic.build_dynamic(true);
	dynamic.solve(true);
	dynamic.transfer_solution_to_wake();
	// Initial phi computation so first time-step can compute cps
	dynamic.compute_phis();
}

void iterate_wagner_case(PanelMethod& dynamic, double& t)
{
	Eigen::ArrayXd forces = Eigen::ArrayXd(MAX_IT);

	for (size_t i = 0; i < MAX_IT; i++)
	{
		Vector3d bvel = Vector3d(0, ADVANCE_VEL * std::sin(LONG_TERM_AOA), -ADVANCE_VEL * std::cos(LONG_TERM_AOA));
		Vector3d omega = Eigen::Vector3d(0.0, 0.0, 0.0);
		dynamic.timestep(bvel, omega);
		std::cout << "It: " << i + 1 << " / " << MAX_IT << " (t = " << t << ")" << std::endl;
		dynamic.compute_cps(false);
		forces(i) = dynamic.compute_aero_force(USE_CENTERLINE)(1);
		t += TSTEP;
	}

	{
		std::stringstream s;
		s << forces.transpose();
		write_string_to_file("workdir/wagner_lift_over_time.dat", s.str());
	}

	std::cout << "Results written to disk" << std::endl;
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

	auto geom = ThinWing::generate(chord_fx, camber_fx, 16, 32, 10.0, 0.25);
	geom->generate_normals();
	write_string_to_file("workdir/geom.dat", geom->quads_to_string());

	write_params();
	long_term_steady_case(geom);

	PanelMethod dynamic;
	double t = 0.0;

	prepare_wagner_case(geom, dynamic);
	auto start = std::chrono::system_clock::now();
	iterate_wagner_case(dynamic, t);
	auto end = std::chrono::system_clock::now();
	std::cout << "Computation took: " <<
			  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}
