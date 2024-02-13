#pragma once
#include "PlanarGeometry.h"
#include "ThinWing.h"
#include <Geometry>
#include <vector>
#include <deque>

class Wake : public PlanarGeometry
{
public:

	std::vector<Eigen::Vector3d> pos_history;
	std::vector<Eigen::Isometry3d> orient_history;

	std::deque<Eigen::Vector3d> vel_history;
	std::deque<Eigen::Vector3d> angvel_history;

	// Maps each wake panel to a panel from the original wing
	Eigen::ArrayXi from_panel;

	Eigen::Index num_edges;

	void integrate_velocities(double wake_scale);
	void build_from_history(const ThinWing& wing);

	void shed_from(const ThinWing& wing, double wake_scale, Eigen::Index num_edges, const Eigen::Vector3d& body_vel,
				   const Eigen::Vector3d& omega);

	void timestep(const ThinWing& wing, double wake_scale, const Eigen::Vector3d& body_vel,
				  const Eigen::Vector3d& omega);

};
