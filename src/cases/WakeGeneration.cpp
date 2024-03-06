#include <iostream>
#include <chrono>
#include "../Util.h"
#include "../PanelMethod.h"

using namespace Eigen;

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

	auto geom = ThinWing::generate(chord_fx, camber_fx, 8, 12, 1.0, 0.25);
	geom->generate_normals();

	PanelMethod dynamic;

	Vector3d vel = Vector3d(0, 1, -10.0);
	Vector3d omg = Vector3d(0, 0, 0);

	dynamic.thin_wings.push_back(geom);
	dynamic.build_geometry_matrix();
	dynamic.shed_initial_wake(100, 0.01, vel, omg);

	// Initial steady solution
	dynamic.build_dynamic(true);
	dynamic.solve(true);
	dynamic.transfer_solution_to_wake();

	double t = 0.0;
	for(size_t i = 0; i < 100; i++)
	{
		vel = Vector3d(0, std::sin(t * 10.0) * 2.0, -10.0 - std::cos(t * 10.0) * 2.0);
		omg = Vector3d(t * 10.0 + 1.0, std::cos(t * 10.0) * 4.0, std::sin(t * 15.0) * 5.0 + t * 10.0);
		dynamic.timestep(vel, omg);
		t += 0.01;
	}


	write_string_to_file("workdir/wake_generation_geom.dat", dynamic.thin_wings[0]->quads_to_string());
	write_string_to_file("workdir/wake_generation_wake_geom.dat", dynamic.wake_geom_to_string(0));
	write_string_to_file("workdir/wake_generation_sln.dat", dynamic.solution_to_string(0));
	write_string_to_file("workdir/wake_generation_wake_sln.dat", dynamic.wake_solution_to_string(0));


}
