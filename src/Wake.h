#pragma once
#include "SheetGeometry.h"
#include "ThinWing.h"
#include <Geometry>
#include <vector>
#include <deque>

class PanelMethod;

class Wake : public SheetGeometry
{
public:
	// Maps each wake panel to a panel from the original wing
	Eigen::ArrayXi from_panel;

	// mu value at trailing edge panel (history)
	// (Used for wake convection)
	std::deque<Eigen::ArrayXd> mu_history;

	// Fixed part of mu value at each wake panel
	Eigen::ArrayXd mus;
	// "Variable" part of mu value at each wake panel (proportional to trailing edge)
	// Influence at each wake panel from its trailing edge, due to interpolation
	// Mostly 0s except for panels near trailing edge or very fast movement
	Eigen::ArrayXd influences;

	void build_from_history(const ThinWing& wing, const PanelMethod& method);

	// Doesn't actually lay out the vertices, but creates them
	void build_initial_geometry(const ThinWing& wing, const PanelMethod& method);

	void mu_convect(size_t wing_idx, const PanelMethod& method);

	void inherit_solution(size_t wake_idx, const PanelMethod& method);
	void transfer_unsteady_solution(size_t wake_idx, const PanelMethod& method);

};
