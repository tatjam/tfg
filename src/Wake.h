#pragma once
#include "PlanarGeometry.h"
#include "ThinWing.h"
#include <Geometry>

class Wake : public PlanarGeometry
{
public:
	// Maps each wake panel to a panel from the original wing
	Eigen::ArrayXi from_panel;

	Eigen::Index num_edges;

	void shed_from(const ThinWing& wing, double wake_scale, Eigen::Index num_edges, const Eigen::Vector3d& body_vel,
				   const Eigen::Vector3d& omega);

};
