#pragma once
#include "ThinWing.h"
#include "Wake.h"

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

	Eigen::ArrayXd solution;

	Eigen::ArrayXd cps;

	// - normal . freestream, regenerated as needed
	Eigen::MatrixXd rhs;

	// This matrix includes the influence of dynamic elements
	Eigen::MatrixXd dynamic_matrix;

	std::vector<size_t> geom_sizes;

	// Return the normal projected induced velocity vector at effect panel caused by cause panel, which would be scaled
	// by "mu" (the doublet strength) to obtain the real induced normal velocity
	Eigen::Vector3d induced_vel(const ThinWing& cause, Eigen::Index cause_panel, const ThinWing& effect, Eigen::Index effect_panel);
	Eigen::Vector3d induced_vel_wake(const Wake& wake, Eigen::Index cause_trailing, const ThinWing& effect, Eigen::Index effect_panel);
	double induced_norm_vel(const ThinWing& cause, Eigen::Index cause_panel, const ThinWing& effect, Eigen::Index effect_panel);
	double induced_norm_vel_wake(const Wake& wake, Eigen::Index cause_trailing, const ThinWing& effect, Eigen::Index effect_panel);

	void build_dynamic_matrix_steady();
	void build_rhs();

	std::vector<Wake> wakes;

	Eigen::Vector3d get_center(const ThinWing& wing, Eigen::Index panel);

public:

	// In body coordinates, which means that incoming airflow
	// is constant regardless of omega
	Eigen::Vector3d body_vel;
	// Omega is assumed to be rotation around the origin, in body axes!
	Eigen::Vector3d omega;

	Eigen::Vector3d induced_vel(const ThinWing& cause, Eigen::Index cause_panel, const Eigen::Vector3d& pos);
	Eigen::Vector3d induced_vel_wake(const Wake& wake, Eigen::Index cause_trailing, const Eigen::Vector3d& pos);
	Eigen::Vector3d induced_vel_verts(const Eigen::Matrix<double, 3, 4>& vertices,
									  const Eigen::Vector3d& nrm,
									  const Eigen::Vector3d& pos);
	// Warning: Do not modify after calling build_geometry_matrix()
	std::vector<std::shared_ptr<ThinWing>> thin_wings;

	void build_geometry_matrix();

	void shed_initial_wake(Eigen::Index num_wake_panels, double wake_distance);
	// Builds the dynamic matrix and RHS matrix
	// Uses wake convection
	void build_dynamic();
	// Same as before but assumes steady condition and thus fixed wake
	void build_dynamic_steady();

	std::string geometry_matrix_to_string();
	std::string dynamic_matrix_to_string();
	// (Use Import["...", "Table"][[1]] to import it in Mathematica)
	std::string solution_to_string(size_t for_geom);
	// Same as before but for pressure coefficients
	std::string cps_to_string(size_t for_geom);

	void compute_cps(double epsilon);
	// Requires cps to be computed
	Eigen::Vector3d compute_aero_force();

	std::string wake_geom_to_string(size_t for_geom);

	// Obtains flow field in a rectangle defined by a corner and two sides (x_axis and y_axis) using given sample points
	std::string sample_flow_field_to_string(Eigen::Vector3d corner, Eigen::Vector3d x_axis, Eigen::Vector3d y_axis, size_t num_x, size_t num_y);

	void solve();

	PanelMethod();

};
