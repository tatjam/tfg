#pragma once
#include "PlanarGeometry.h"
#include "ThinWing.h"
#include <Geometry>
#include <vector>
#include <deque>

class PanelMethod;

class Wake : public PlanarGeometry
{
public:
	// Maps each wake panel to a panel from the original wing
	Eigen::ArrayXi from_panel;

	void build_from_history(const ThinWing& wing, const PanelMethod& method);

	// Doesn't actually lay out the vertices, but creates them
	void build_initial_geometry(const ThinWing& wing, const PanelMethod& method);
	void timestep(const ThinWing& wing, const PanelMethod& method);

};
