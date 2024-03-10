#pragma once
#include "ThinWing.h"
#include "Wake.h"

class PanelMethod
{
	friend class Wake;
protected:

	// This matrix doesn't include influence of dynamic elements (wakes for now)
	//  An element of the matrix (effect, cause) represents
	// 	the influence of panel 'cause' on panel 'effect' in the sense of
	//  induced velocity (dot) surface normal
	// Panels are indexed in the order of the geometry, but because multiple geometries may be present
	// we will offset them:
	//  index = index_within_geometry + geom_sizes[geometry]
	Eigen::MatrixXd geometry_matrix;

	std::deque<Eigen::ArrayXd> sln_hist;

	std::deque<Eigen::ArrayXd> phi_hist_above;
	std::deque<Eigen::ArrayXd> phi_hist_below;

	Eigen::ArrayXd cps;

	// - normal . freestream, regenerated as needed
	Eigen::MatrixXd rhs;

	// This matrix includes the influence of dynamic elements
	Eigen::MatrixXd dynamic_matrix;

	std::vector<size_t> geom_sizes;

	// Return the normal projected induced velocity vector at effect panel caused by cause panel, which would be scaled
	// by "mu" (the doublet strength) to obtain the real induced normal velocity
	Eigen::Vector3d induced_vel(const ThinWing& cause, Eigen::Index cause_panel, const ThinWing& effect, Eigen::Index effect_panel);
	// MODE = 0 -> Returns induced velocity for steady mode
	// MODE = 1 -> Returns induced velocity multiplied by wake.influence
	// MODE = 2 -> Returns induced velocity multiplied by wake.mus
	Eigen::Vector3d induced_vel_wake(
			const Wake& wake, Eigen::Index cause_trailing, const ThinWing& effect, Eigen::Index effect_panel,
			const int MODE);
	double induced_norm_vel(const ThinWing& cause, Eigen::Index cause_panel, const ThinWing& effect, Eigen::Index effect_panel);

	double induced_norm_vel_wake(const Wake& wake, Eigen::Index cause_trailing, const ThinWing& effect,
								 Eigen::Index effect_panel, const bool STEADY);

	// Uses latest solution!
	double induced_phi(size_t cause_geom, Eigen::Index cause_panel, const ThinWing& effect, Eigen::Index effect_panel,
					   bool above);
	// Uses latest solution!
	double induced_phi_wake(size_t cause_geom, Eigen::Index cause_trailing, const ThinWing& effect, Eigen::Index effect_panel,
							bool above);


	void build_dynamic_matrix(const bool STEADY);
	void build_rhs(const bool STEADY);

	Eigen::Vector3d get_center(const ThinWing& wing, Eigen::Index panel);
	double get_area(const ThinWing& wing, Eigen::Index panel);

	void integrate_velocities(double wake_scale);
	void build_wakes_from_history();

public:
	std::vector<Wake> wakes;

	// We will store backward_difference_order + 1 solutions
	int backward_difference_order;

	double wake_scale;
	// NOTE: Includes trailing edge panel (which is always fixed!)
	Eigen::Index num_wake_edges;

	std::vector<Eigen::Vector3d> pos_history;
	std::vector<Eigen::Isometry3d> orient_history;

	std::deque<Eigen::Vector3d> vel_history;
	std::deque<Eigen::Vector3d> angvel_history;

	Eigen::Vector3d induced_vel(const ThinWing& cause, Eigen::Index cause_panel, const Eigen::Vector3d& pos);
	Eigen::Vector3d induced_vel_wake(const Wake& wake, Eigen::Index cause_trailing, const Eigen::Vector3d& pos,
									 const int MODE);
	Eigen::Vector3d induced_vel_verts(const Eigen::Matrix<double, 3, 4>& vertices,
									  const Eigen::Vector3d& nrm,
									  const Eigen::Vector3d& pos);

	double induced_phi(const ThinWing& cause, Eigen::Index cause_panel, const Eigen::Vector3d& pos,
					   bool above);
	double induced_phi_wake(size_t cause_geom, Eigen::Index cause_trailing, const Eigen::Vector3d& pos,
						bool above);
	double induced_phi_verts(const Eigen::Matrix<double, 3, 4>& vertices,
							 	const Eigen::Vector3d& nrm,
								const Eigen::Vector3d& pos, bool above);

	// Warning: Do not modify after calling build_geometry_matrix()
	std::vector<std::shared_ptr<ThinWing>> thin_wings;

	void build_geometry_matrix();

	// Pushes front into phis_history
	void compute_phis();


	// This "resets" all wakes to stationary condition,
	// potentially resizing vertex counts (all wakes must have same panel count)
	// num_wake_edges includes trailing edge!
	void shed_initial_wake(Eigen::Index num_wake_edges, double wake_scale,
						   const Eigen::Vector3d& body_vel, const Eigen::Vector3d& body_angvel);

	// Meant to be used after steady solution to have an starting point for
	// wake convection
	void transfer_solution_to_wake();

	// Builds the dynamic matrix and RHS matrix
	void build_dynamic(const bool STEADY);

	std::string geometry_matrix_to_string();
	std::string dynamic_matrix_to_string();
	// (Use Import["...", "Table"][[1]] to import it in Mathematica)
	std::string solution_to_string(size_t for_geom);
	// Same as before but for pressure coefficients
	std::string cps_to_string(size_t for_geom);
	// Same as before but for the mus over the wake of a geometry
	// (It's a 1D array that maps to each polygon the solution, order of polygons
	//  set as in wake_geom_to_string)
	std::string wake_solution_to_string(size_t for_geom);

	// (Use Import["...", "Table"][[1]] to import it in Mathematica)
	std::string phis_to_string(size_t for_geom, bool include_freestream, bool above);

	void compute_cps(bool STEADY);

	void timestep(const Eigen::Vector3d& cur_vel, const Eigen::Vector3d& cur_angvel);

	// Requires cps to be computed
	Eigen::Vector3d compute_aero_force(bool centerline);

	std::string wake_geom_to_string(size_t for_geom);

	// Obtains flow field in a rectangle defined by a corner and two sides (x_axis and y_axis) using given sample points
	std::string sample_flow_field_to_string(Eigen::Vector3d corner, Eigen::Vector3d x_axis, Eigen::Vector3d y_axis, size_t num_x, size_t num_y);

	// If STEADY is true, solutions are overwritten to simulate steady state
	void solve(bool STEADY);

	PanelMethod();

};
