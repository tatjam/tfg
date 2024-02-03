#pragma once
#include "ThinWing.h"

class PanelMethod
{
private:

	// This matrix doesn't include influence of dynamic elements (wakes for now)
	//  An element of the matrix (effect, cause) represents
	// 	the influence of panel 'cause' on panel 'effect' in the sense of
	//  induced velocity (dot) surface normal
	// Panels are indexed in the order of the geometry, but because multiple geometries may be present
	// we will offset them:
	//  index = index_within_geometry + geom_sizes[geometry]
	Eigen::MatrixXd geometry_matrix;

	// normal . freestream, regenerated as needed
	Eigen::MatrixXd freestream_rhs;

	// This matrix includes the influence of dynamic elements
	Eigen::MatrixXd dynamic_matrix;

	std::vector<size_t> geom_sizes;

	// Return the normal projected induced velocity vector at effect panel caused by cause panel, which would be scaled
	// by "mu" (the doublet strength) to obtain the real induced normal velocity
	double induced_norm_vel(const ThinWing& cause, Eigen::Index cause_panel, const ThinWing& effect, size_t effect_panel);


public:

	Eigen::Vector3d induced_vel(const ThinWing& cause, Eigen::Index cause_panel, const Eigen::Vector3d& pos);
	// Warning: Do not modify after calling build_geometry_matrix()
	std::vector<std::shared_ptr<ThinWing>> thin_wings;

	void build_geometry_matrix();
	void build_dynamic_matrix();

	void solve();

};
