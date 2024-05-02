#include <iostream>
#include <chrono>
#include "../Util.h"
#include "../PanelMethod.h"
#include <filesystem>

using namespace Eigen;

const double ADVANCE_VEL = 1.0;
const size_t LONG_TERM_NPANELS = 200;
const double LONG_TERM_TIMESTEP = 0.15;
const double SPAN = 10.0;
const double CHORD = 1.0;


const size_t NPANELS = 200;

std::string make_fname(double AoA, const std::string& sub)
{
	if(std::abs(AoA) < 0.00001)
	{
		AoA = 0.0;
	}
	std::stringstream fname_base;
	fname_base << std::fixed << std::setprecision(2);
	fname_base << "workdir/clalpha/" << sub << "_" << AoA << "_.dat";
	return fname_base.str();
}

void write_params(double AoA)
{
	std::stringstream s;
	s << ADVANCE_VEL << std::endl;
	s << AoA << std::endl;
	write_string_to_file(make_fname(AoA, "params"), s.str());
}

void long_term_steady_case(std::shared_ptr<ThinWing> wing, double AoA)
{
	write_params(AoA);
	Vector3d bvel = Vector3d(0, ADVANCE_VEL * std::sin(AoA), -ADVANCE_VEL * std::cos(AoA));
	Vector3d omega = Vector3d(0, 0, 0);

	PanelMethod pan;
	pan.thin_wings.push_back(wing);
	pan.shed_initial_wake(LONG_TERM_NPANELS, LONG_TERM_TIMESTEP, bvel, omega);
	pan.build_geometry_matrix();
	pan.build_dynamic(true);
	pan.solve(true);
	pan.compute_cps(true);

	Vector3d aero_force = pan.compute_aero_force(false);

	{
		std::stringstream s;
		s << aero_force.transpose();
		write_string_to_file(make_fname(AoA, "aero_force_3d"), s.str());
	}

	aero_force = pan.compute_aero_force(true);

	{
		std::stringstream s;
		s << aero_force.transpose();
		write_string_to_file(make_fname(AoA, "aero_force_2d"), s.str());
	}
}


int main()
{
	std::filesystem::create_directory("./workdir/clalpha/");
	auto chord_fx = [](double span_pos)
	{
		return std::make_pair(1.0, 0.0);
	};

	auto camber_fx = [](double span_pos, double chord_pos)
	{
		return 0.0;
	};

	auto geom = ThinWing::generate(chord_fx, camber_fx, 32, 32, SPAN, CHORD);
	geom->generate_normals();
	write_string_to_file("workdir/geom.dat", geom->quads_to_string());

	for(double AoA = -0.4; AoA < 0.4; AoA += 0.01)
	{
		long_term_steady_case(geom, AoA);
	}


}
