#include "PanelMethod.h"

#include <IterativeLinearSolvers>
#include <iostream>

using namespace Eigen;

void PanelMethod::build_geometry_matrix()
{
	geom_sizes.clear();
	geom_sizes.reserve(thin_wings.size());
	size_t total_size = 0;
	for(const auto& ptr : thin_wings)
	{
		geom_sizes.push_back(total_size);
		total_size += ptr->quads.cols();
	}
	geom_sizes.push_back(total_size);

	geometry_matrix = MatrixXd(total_size, total_size);
	dynamic_matrix = MatrixXd(total_size, total_size);
	geometry_matrix.setZero();
	dynamic_matrix.setZero();
	rhs = MatrixXd(total_size, 1);
	cps = MatrixXd(total_size, 1);

	for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
	{
		for(size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
		{
			for(Index effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
			{
				for(Index cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
				{
					double induced = induced_norm_vel(*thin_wings[cause_geom], cause, *thin_wings[effect_geom], effect);
					auto effect_idx = (Index)(effect + geom_sizes[effect_geom]);
					auto cause_idx = (Index)(cause + geom_sizes[cause_geom]);

					geometry_matrix(effect_idx, cause_idx) = induced;
				}
			}
		}
	}

}

Eigen::Vector3d PanelMethod::induced_vel(const ThinWing &cause, Eigen::Index cause_panel, const ThinWing &effect,
										 Eigen::Index effect_panel)
{
	Vector3d center = get_center(effect, effect_panel);
	return induced_vel(cause, cause_panel, center);
}

double PanelMethod::induced_phi(size_t cause_geom, Eigen::Index cause_panel, const ThinWing &effect,
								Eigen::Index effect_panel, bool above)
{
	Vector3d center = get_center(effect, effect_panel);
	auto cause_idx = (Index)(cause_panel + geom_sizes[cause_geom]);

	return induced_phi(*thin_wings[cause_geom], cause_panel, center, above) * sln_hist.front()(cause_idx);
}

Eigen::Vector3d PanelMethod::induced_vel_wake(const Wake &wake, Eigen::Index cause_trailing, const ThinWing &effect,
											  Eigen::Index effect_panel, const int MODE)
{
	Vector3d center = get_center(effect, effect_panel);
	return induced_vel_wake(wake, cause_trailing, center, MODE);
}

double PanelMethod::induced_phi_wake(size_t cause_geom, Eigen::Index cause_trailing, const ThinWing &effect,
									 Eigen::Index effect_panel, bool above)
{
	Vector3d center = get_center(effect, effect_panel);
	return induced_phi_wake(cause_geom, cause_trailing, center, above);
}

Eigen::Vector3d PanelMethod::induced_vel_wake(const Wake &wake, Eigen::Index cause_trailing, const Vector3d &pos,
											  const int MODE)
{
	Vector3d acc = Vector3d(0, 0, 0);
	for(Index i = 0; i < wake.num_shed; i++)
	{
		Matrix<double, 3, 4> verts;
		Index cause_panel = cause_trailing * (num_wake_edges - 1) + i;
		verts.col(0) = wake.vertices.col(wake.quads(0, cause_panel));
		verts.col(1) = wake.vertices.col(wake.quads(1, cause_panel));
		verts.col(2) = wake.vertices.col(wake.quads(2, cause_panel));
		verts.col(3) = wake.vertices.col(wake.quads(3, cause_panel));

		Vector3d nrm = wake.normals.col(cause_panel);
		Vector3d val = induced_vel_verts(verts, nrm, pos);

		Index panel_idx = cause_trailing * (num_wake_edges - 1) + i;
		if(MODE == 1)
		{
			val *= wake.influences(panel_idx);
		}
		else if(MODE == 2)
		{
			double infl = wake.influences(panel_idx);
			val *= wake.mus(panel_idx) * (1.0 - infl);
		}

		acc += val;
	}
	return acc;
}


double PanelMethod::induced_phi_wake(size_t cause_geom, Eigen::Index cause_trailing, const Vector3d &pos, bool above)
{
	double acc = 0.0;
	const Wake& wake = wakes[cause_geom];
	for(Index i = 0; i < wake.num_shed; i++)
	{
		Index cause = thin_wings[cause_geom]->trailing_panels(cause_trailing);
		auto cause_idx = (Index)(cause + geom_sizes[cause_geom]);
		Index panel_idx = cause_trailing * (num_wake_edges - 1) + i;


		// Compute mu
		double infl = wake.influences(panel_idx);
		double mu = wake.mus(panel_idx) * (1.0 - infl);
		mu += sln_hist.front()(cause_idx) * infl;

		Matrix<double, 3, 4> verts;
		Index cause_panel = cause_trailing * (num_wake_edges - 1) + i;
		verts.col(0) = wake.vertices.col(wake.quads(0, cause_panel));
		verts.col(1) = wake.vertices.col(wake.quads(1, cause_panel));
		verts.col(2) = wake.vertices.col(wake.quads(2, cause_panel));
		verts.col(3) = wake.vertices.col(wake.quads(3, cause_panel));

		Vector3d nrm = wake.normals.col(cause_panel);
		double val = induced_phi_verts(verts, nrm, pos, above);

		acc += val * mu;
	}
	return acc;
}

double
PanelMethod::induced_norm_vel(const ThinWing& cause, Index cause_panel, const ThinWing& effect, Index effect_panel)
{
	Vector3d ivel = induced_vel(cause, cause_panel, effect, effect_panel);
	Vector3d normal = effect.normals.col(effect_panel).matrix();
	return ivel.dot(normal);
}

double PanelMethod::induced_norm_vel_wake(const Wake &wake, Eigen::Index cause_trailing, const ThinWing &effect,
										  Eigen::Index effect_panel, const bool STEADY)
{
	// 0 = STEADY mode, 1 = Multiply by influence
	Vector3d acc = induced_vel_wake(wake, cause_trailing, effect, effect_panel, STEADY ? 0 : 1);
	Vector3d normal = effect.normals.col(effect_panel).matrix();
	return acc.dot(normal);
}

// Velocity induced by an (assumed planar) doublet sheet, doublets oriented such that they point in the
// direction of the panel normal.
// (Developed from Low Speed Aerodynamics (...) by Joseph Katz and Allen Plotkin)
// This function is singular on the quadrilateral edges (technically the flat projection)
Eigen::Vector3d PanelMethod::induced_vel(const ThinWing &cause, Index cause_panel, const Vector3d &pos)
{
	Matrix<double, 3, 4> verts;
	verts.col(0) = cause.transform * cause.vertices.col(cause.quads(0, cause_panel));
	verts.col(1) = cause.transform * cause.vertices.col(cause.quads(1, cause_panel));
	verts.col(2) = cause.transform * cause.vertices.col(cause.quads(2, cause_panel));
	verts.col(3) = cause.transform * cause.vertices.col(cause.quads(3, cause_panel));
	Vector3d nrm = cause.normals.col(cause_panel);
	return induced_vel_verts(verts, nrm, pos);
}

double PanelMethod::induced_phi(const ThinWing &cause, Eigen::Index cause_panel, const Vector3d &pos, bool above)
{
	Matrix<double, 3, 4> verts;
	verts.col(0) = cause.transform * cause.vertices.col(cause.quads(0, cause_panel));
	verts.col(1) = cause.transform * cause.vertices.col(cause.quads(1, cause_panel));
	verts.col(2) = cause.transform * cause.vertices.col(cause.quads(2, cause_panel));
	verts.col(3) = cause.transform * cause.vertices.col(cause.quads(3, cause_panel));
	Vector3d nrm = cause.normals.col(cause_panel);
	return induced_phi_verts(verts, nrm, pos, above);
}

void PanelMethod::build_rhs(const bool STEADY)
{
	rhs.setZero();
	for(size_t geom = 0; geom < thin_wings.size(); geom++)
	{
		for(size_t panel = 0; panel < thin_wings[geom]->quads.cols(); panel++)
		{
			Vector3d norm = thin_wings[geom]->normals.col(panel);
			Vector3d center = get_center(*thin_wings[geom], panel);
			Vector3d freestream = -vel_history.front();
			freestream -= angvel_history.front().cross(center);
			rhs(panel + geom_sizes[geom], 0) = -norm.dot(freestream);

			if(!STEADY)
			{
				Vector3d induced_wakes; induced_wakes.setZero();
				// Include fixed wake mus
				for(size_t wake_geom = 0; wake_geom < thin_wings.size(); wake_geom++)
				{
					for(Index trail_pan = 0; trail_pan < thin_wings[geom]->trailing_panels.rows(); trail_pan++)
					{
						induced_wakes += induced_vel_wake(wakes[wake_geom], trail_pan, *thin_wings[geom], panel, 2);
					}
				}
				rhs(panel + geom_sizes[geom], 0) -= norm.dot(induced_wakes);
			}

		}
	}
}

void PanelMethod::build_dynamic_matrix(const bool STEADY)
{
	// We simply prescribe the same circulation to the wake than in the trailing
	// edge panels, thus the wake become "integrated" into the contribution
	// from the last panel.

	// VERY IMPORTANT, otherwise some untouched panels will get uninitialized values
	// or values from previous steps
	dynamic_matrix.setZero();

	for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
	{
		for(size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
		{
			for(Index effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
			{
				for(Index cause_trail = 0; cause_trail < thin_wings[cause_geom]->trailing_panels.rows(); cause_trail++)
				{
					Index cause = thin_wings[cause_geom]->trailing_panels(cause_trail);

					double induced = induced_norm_vel_wake(wakes[cause_geom], cause_trail,
																   *thin_wings[effect_geom], effect, STEADY);
					auto effect_idx = (Index)(effect + geom_sizes[effect_geom]);
					auto cause_idx = (Index)(cause + geom_sizes[cause_geom]);

					dynamic_matrix(effect_idx, cause_idx) += induced;
				}
			}
		}
	}
}

void PanelMethod::build_dynamic(const bool STEADY)
{
	build_dynamic_matrix(STEADY);
	build_rhs(STEADY);
}

std::string PanelMethod::geometry_matrix_to_string()
{
	std::stringstream o;
	o << geometry_matrix;
	return o.str();
}

std::string PanelMethod::dynamic_matrix_to_string()
{
	std::stringstream o;
	o << dynamic_matrix;
	return o.str();
}

std::string PanelMethod::solution_to_string(size_t for_geom)
{
	std::stringstream o;
	o << sln_hist.front().block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
	return o.str();
}

std::string PanelMethod::wake_solution_to_string(size_t for_geom)
{
	std::stringstream o;
	o << wakes[for_geom].mus.transpose();
	return o.str();
}

std::string PanelMethod::cps_to_string(size_t for_geom)
{
	std::stringstream o;
	o << cps.block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
	return o.str();
}

PanelMethod::PanelMethod()
{
	// TODO: It only works well for = 1
	backward_difference_order = 1;
}

void PanelMethod::solve(bool STEADY)
{
	Eigen::ArrayXd guess;
	if (sln_hist.size() == backward_difference_order + 1)
	{
		sln_hist.pop_back();
	}

	if(phi_hist_above.size() == backward_difference_order + 1)
	{
		phi_hist_above.pop_back();
		phi_hist_below.pop_back();
	}

	if(sln_hist.size() > 0)
	{
		guess = sln_hist.back();
	}
	else
	{
		guess = Eigen::ArrayXd(rhs.size());
		guess.setZero();
	}

	MatrixXd mat = geometry_matrix + dynamic_matrix;
	/*
	BiCGSTAB<MatrixXd> solver;
	solver.compute(mat);
	VectorXd sln = solver.solveWithGuess(rhs, guess);
	*/
	VectorXd sln = mat.householderQr().solve(rhs);
	sln_hist.push_front(sln);
	double relative_error = (mat*sln - rhs).norm() / rhs.norm();
	std::cout << "Solve relative error is " << relative_error << std::endl;

	if(!STEADY)
	{
		//compute_phis();
	}

	if(STEADY)
	{
		if(sln_hist.size() != backward_difference_order + 1)
		{
			sln_hist.resize(backward_difference_order + 1);
		}
		for(size_t i = 1; i < backward_difference_order + 1; i++)
		{
			sln_hist[i] = sln_hist[0];
		}
	}
}

void PanelMethod::shed_initial_wake(Eigen::Index num_wake_edges, double wake_scale,
									const Vector3d& body_vel, const Vector3d& omega)
{
	this->num_wake_edges = num_wake_edges;
	this->wake_scale = wake_scale;
	// State vector
	vel_history.clear();
	angvel_history.clear();

	pos_history.clear();
	orient_history.clear();

	pos_history.reserve(num_wake_edges - 1);
	orient_history.reserve(num_wake_edges - 1);

	for(Index zi = 1; zi < num_wake_edges; zi++)
	{
		vel_history.push_back(body_vel);
		angvel_history.push_back(omega);
	}
	integrate_velocities(wake_scale);

	// (Re)create wakes
	wakes.clear();
	for(const auto& wing : thin_wings)
	{
		Wake& w = wakes.emplace_back();
		w.build_initial_geometry(*wing, *this);
		w.build_from_history(*wing, *this);
	}


}

std::string PanelMethod::wake_geom_to_string(size_t for_geom)
{
	return wakes[for_geom].quads_to_string();
}


Eigen::Vector3d
PanelMethod::induced_vel_verts(const Matrix<double, 3, 4> &vertices, const Vector3d& normal, const Vector3d &pos)
{
	// Ring approach
	// (Velocity induced by ring vortex)
	Vector3d out; out.setZero();

	for(Index i = 0; i < 4; i++)
	{
		Index ni = i + 1;
		if(ni == 4)
			ni = 0;
		Vector3d pi = vertices.col(i);
		Vector3d pj = vertices.col(ni);
		Vector3d ri = pos - pi;
		Vector3d rj = pos - pj;
		Vector3d rij = pj - pi;

		Vector3d cr = ri.cross(rj);

		out += cr * rij.dot(ri.normalized() - rj.normalized()) / cr.squaredNorm();
	}

	out /= 4.0 * M_PI;

	return out;

}

double PanelMethod::induced_phi_verts(const Matrix<double, 3, 4> &vertices, const Vector3d &nrm, const Vector3d &pos,
									  bool above)
{
	// The doublet potential "happens" to match the induced normal velocity for a constant source sheet!
	// (Normal to the source sheet)
	// So we implement a "brute-force" way, by transforming into panel coordinates and back
	// Taken from https://github.com/byuflowlab/FLOWPanel.jl
	Vector3d center = vertices.col(0).matrix();
	center += vertices.col(1).matrix();
	center += vertices.col(2).matrix();
	center += vertices.col(3).matrix();
	center /= 4.0;

	Quaterniond quat; quat.setFromTwoVectors(nrm, Vector3d(0, 0, 1));
	Isometry3d transform; transform.setIdentity();

	transform = transform * quat * Translation3d(-center);

	// Transform the points into local coordinates
	Vector3d p = transform * pos;

	if(std::abs(p(2)) < 0.001)
	{
		bool all_pos = true;

		for(Index i = 0; i < 4; i++)
		{
			Vector3d pi = transform * vertices.col(i).matrix();
			Index next_i = i + 1;
			if(next_i == 4) next_i = 0;
			Vector3d pj = transform * vertices.col(i + 1).matrix();

			double dij = (pj - pi).norm();
			double ri = (p - pi).norm();
			double rj = (p - pj).norm();

			double Qij = std::log((ri + rj + dij)/(ri + rj - dij + 0.001));

			double Sij = (pj(1) - pi(1)) / dij;
			double Cij = (pj(0) - pi(0)) / dij;

			double siji = (pi(0) - p(0)) * Cij + (pi(1) - p(1)) * Sij;
			double sijj = (pj(0) - p(0)) * Cij + (pj(1) - p(1)) * Sij;
			double Rij = (p(0) - pi(0)) * Sij - (p(1) - pi(1)) * Cij;

			if(Rij < 0)
			{
				all_pos = false;
				break;
			}
		}

		if(all_pos)
		{
			// Inside
			return above ? 0.5 : -0.5;
		}
		else
		{
			// Outside
			return 0.0;
		}

	}

	double w = 0.0;
	double dtheta = 2.0 * M_PI;
	for(Index i = 0; i < 4; i++)
	{
		Vector3d pi = transform * vertices.col(i).matrix();
		Index next_i = i + 1;
		if(next_i == 4) next_i = 0;
		Vector3d pj = transform * vertices.col(i + 1).matrix();

		double dij = (pj - pi).norm();
		double ri = (p - pi).norm();
		double rj = (p - pj).norm();

		double Qij = std::log((ri + rj + dij)/(ri + rj - dij + 0.001));

		double Sij = (pj(1) - pi(1)) / dij;
		double Cij = (pj(0) - pi(0)) / dij;

		double siji = (pi(0) - p(0)) * Cij + (pi(1) - p(1)) * Sij;
		double sijj = (pj(0) - p(0)) * Cij + (pj(1) - p(1)) * Sij;
		double Rij = (p(0) - pi(0)) * Sij - (p(1) - pi(1)) * Cij;

		double Jij = std::atan2( Rij * std::abs(p(2)) * ( ri * sijj - rj * siji),
								 ri * rj * Rij * Rij + p(2) * p(2) * sijj * siji);

		w -= 1.0 / (4.0 * M_PI) * (Rij * Qij + std::abs(p(2)) * Jij);
		dtheta *= Rij >= 0.0 ? 1.0 : 0.0;
	}

	w += (1.0) / (4.0 * M_PI) * std::abs(p(2)) * dtheta;

	return w;
}

Eigen::Vector3d PanelMethod::get_center(const ThinWing &wing, Eigen::Index wing_panel)
{
	Vector3d center = wing.transform * wing.vertices.col(wing.quads(0, wing_panel)).matrix();
	center += wing.transform * wing.vertices.col(wing.quads(1, wing_panel)).matrix();
	center += wing.transform * wing.vertices.col(wing.quads(2, wing_panel)).matrix();
	center += wing.transform * wing.vertices.col(wing.quads(3, wing_panel)).matrix();
	center /= 4.0;
	return center;
}

double PanelMethod::get_area(const ThinWing &wing, Eigen::Index panel)
{
	// We approximate the area of the cuadrilateral as two triangles
	Vector3d v1 = wing.vertices.col(wing.quads(0, panel));
	Vector3d v2 = wing.vertices.col(wing.quads(1, panel));
	Vector3d v3 = wing.vertices.col(wing.quads(2, panel));
	Vector3d v4 = wing.vertices.col(wing.quads(3, panel));
	return 0.5 * (v2 - v1).cross(v3 - v1).norm() + 0.5 * (v4 - v2).cross(v4 - v3).norm();
}

	Eigen::Vector3d PanelMethod::compute_aero_force(bool centerline, int center_line_pos)
	{
		Vector3d acc;
		acc.setZero();

		if (centerline)
		{
			if(center_line_pos < 0 || center_line_pos >= thin_wings[0]->num_spanwise - 1)
			{
				center_line_pos = thin_wings[0]->num_spanwise / 2;
			}

			for (Index i = 0; i < thin_wings[0]->num_chorwise - 1; i++)
			{
				Index idx = center_line_pos * (thin_wings[0]->num_chorwise - 1) + i;
				Vector3d nrm = thin_wings[0]->normals.col(idx);
				double area = get_area(*thin_wings[0], idx);
				Vector3d a = thin_wings[0]->vertices.col(thin_wings[0]->quads(1, idx));
				Vector3d b = thin_wings[0]->vertices.col(thin_wings[0]->quads(2, idx));
				double length = (b - a).norm();
				Vector3d force = nrm * cps(idx);
				acc += force * length;
			}
		}
		else
		{
			for(size_t geom = 0; geom < thin_wings.size(); geom++)
			{
				for(Index panel = 0; panel < thin_wings[geom]->quads.cols(); panel++)
				{
					Vector3d nrm = thin_wings[geom]->normals.col(panel);
					double area = get_area(*thin_wings[geom], panel);

					Vector3d force = nrm * cps(geom_sizes[geom] + panel);
					acc += force * area;
				}
			}
		}
		return acc;
	}

	std::string
	PanelMethod::sample_flow_field_to_string(Vector3d corner, Vector3d x_axis, Vector3d y_axis,
											 size_t num_x, size_t num_y)
	{
		std::stringstream o;
		o << std::fixed;
		o << "{";

		for(size_t y = 0; y < num_y; y++)
		{
			o << "{";
			for(size_t x = 0; x < num_x; x++)
			{
				double xpos = (double)x / ((double)num_x - 1.0);
				double ypos = (double)y / ((double)num_y - 1.0);
				Vector3d pos = corner + x_axis * xpos + y_axis * ypos;

				Vector3d total_ind; total_ind.setZero();
				Vector3d freestream = -vel_history.front();
				freestream -= angvel_history.front().cross(pos);

				for (size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
				{
					// Panel effect
					for (Index cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
					{
						Vector3d ind = induced_vel(*thin_wings[cause_geom], cause, pos);
						ind *= sln_hist.front()(geom_sizes[cause_geom] + cause);
						total_ind += ind;
					}
					// Wake effect
					for(Index trail_cause = 0; trail_cause < thin_wings[cause_geom]->trailing_panels.rows(); trail_cause++)
					{
						Vector3d ind = induced_vel_wake(wakes[cause_geom], trail_cause, pos, 2);
						// NOTE: Mode 2 already premultiplies by the correct mu value
						total_ind += ind;
					}
				}

				Vector3d total = freestream + total_ind;

				o << "{";
				o << total(0);
				o << ",";
				o << total(1);
				o << ",";
				o << total(2);
				o << "}";
				if(x != num_x - 1)
					o << ",";
			}
			o << "}";
			if(y != num_y - 1)
				o << ",";
		}

		o << "}";
		return o.str();
	}

	// Source: https://stackoverflow.com/a/44719219
	constexpr inline size_t binom(size_t n, size_t k) noexcept
	{
		return
				(        k> n  )? 0 :          // out of range
				(k==0 || k==n  )? 1 :          // edge
				(k==1 || k==n-1)? n :          // first
				binom(n - 1, k - 1) * n / k;   // recursive
	}


	void PanelMethod::compute_cps(bool STEADY)
	{
		// Square of reference velocity
		const double Vref2 = 1.0;
		assert(sln_hist.size() > 1 || STEADY);

		Eigen::ArrayXd& solution = sln_hist.front();
		bool message_put = false;

		// TODO: This could be improved much by caching weighted averages
		for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
		{
			const ThinWing& tw = *thin_wings[effect_geom];
			for (Index effect = 0; effect < tw.quads.cols(); effect++)
			{
				Vector3d center = get_center(tw, effect);
				Vector3d normal = tw.normals.col(effect);

				// TODO: Signs
				Vector3d freestream = -vel_history.front();
				// Non-dimensional freestream doesn't include angular velocity component!
				Vector3d ref_freestream = freestream;
				freestream -= angvel_history.front().cross(center);

				Array4d areas; areas.setZero();
				Array4d mus; mus.setZero();

				double self_area = get_area(tw, effect);
				double self_sol = solution(geom_sizes[effect_geom] + effect);

				// Compute area-weighted average of mu at each node of the quadrilateral,
				// exploting the mesh structure
				for(Index vert = 0; vert < 4; vert++)
				{
					// Self contribution
					areas(vert) += self_area;
					mus(vert) += self_sol * self_area;

					Index vert_idx = tw.quads(vert, effect);

					for(Index nbor = 0; nbor < 8; nbor++)
					{
						Index nbor_idx = tw.neighbors(nbor, effect);
						if(nbor_idx < 0)
							continue;

						// Check that the neighbor also contains this vertex
						// TODO: This is technically not needed in such a structured mesh!
						bool has_vert = false;
						for(Index svert = 0; svert < 4; svert++)
						{
							Index vidx = tw.quads(svert, nbor_idx);
							if(vidx == vert_idx)
							{
								has_vert = true;
								break;
							}
						}

						if(!has_vert)
							continue;

						double area = get_area(tw, nbor_idx);
						areas(vert) += area;
						mus(vert) += area * solution(geom_sizes[effect_geom] + nbor_idx);

					}

					mus(vert) /= areas(vert);

					// this neighbor is always in -z dir, so if it's not present
					// we are a leading edge
					bool is_leading_edge =
							tw.neighbors(1, effect) < 0;
					if(is_leading_edge && vert <= 1)
					{
						// If that's the case, the area approximation is horrible, as it
						// assumes that there's a "wake" coming from the freestream
						// We instead set the mu to a negative value:
						// 0 would make the most sense, but because we are using
						// a small number of panels a "exageration" factor can be introduced.
						const double exag_factor = -0.0;
						// wowaw
						mus(vert) *= exag_factor;
					}
				}


				// We know mu at each of the quadrilateral corners, and at its center
				// For convenience, we will project the quadrilateral into a flat plane,
				// (assuming it's flat, small error) such that the normal vector points
				// in the z direction, and the quadrilateral is contained in the x/y plane

				// Because we don't really care about the orientation of x and y, we can simply
				// use the axis-to-axis rotation quaternion
				Quaterniond quat; quat.setFromTwoVectors(normal, Vector3d(0, 0, 1));
				Isometry3d transform; transform.setIdentity();

				Isometry3d tinv = (transform * quat).inverse();
				// Note that this implies that quat is applied before the translation!
				transform = transform * quat * Translation3d(-center);

				// Now the center of the quadrilateral is at its origin, and its vertices
				// can be assumed to lay flat on the plane (up to a good approximation)
				Matrix<double, 3, 5> ps;
				ps.col(0) = transform * tw.transform * tw.vertices.col(tw.quads(0, effect));
				ps.col(1) = transform * tw.transform * tw.vertices.col(tw.quads(1, effect));
				ps.col(2) = transform * tw.transform * tw.vertices.col(tw.quads(2, effect));
				ps.col(3) = transform * tw.transform * tw.vertices.col(tw.quads(3, effect));
				// Center point
				ps.col(4) = Vector3d(0, 0, 0);

				// We find the best fit for the following conic
				// phi = Ax + By + C
				// Which has gradient
				// (A, B)
				// Thus at the origin (x = 0, y = 0) the gradient is simply
				// (A, B)
				// Which can be converted back into 3D by inverse transformation
				// We use the least squares method to find such fit

				MatrixXd lsq = MatrixXd(5 /*data points*/, 3 /* parameters */);
				for(Index p = 0; p < 5; p++)
				{
					lsq(p, 0) = ps(0, p); // Ax
					lsq(p, 1) = ps(1, p); // By
					lsq(p, 2) = 1.0; // C
				}

				// Solve the least squares problem (see "Introducción a los métodos numéricos", José Antonio Ezquerro)
				Vector<double, 5> rhs;
				for(Index p = 0; p < 4; p++)
				{
					rhs(p) = mus(p);
				}
				rhs(4) = self_sol;

				// 3x5 * 5x1 = 3x1
				Vector<double, 3> rhs_t = lsq.transpose() * rhs;
				// 3x5 * 5x3 = 3x3
				Matrix<double, 3, 3> A = lsq.transpose() * lsq;
				Vector<double, 3> parameters = A.fullPivLu().solve(rhs_t);

				// grad phi
				Vector3d grad = Vector3d(parameters(0), parameters(1), 0.0);
				grad = tinv * grad;

				Index panel_idx = geom_sizes[effect_geom] + effect;

				cps(panel_idx) = -2.0 * (freestream.dot(grad)) / Vref2 * Vref2;
				// Corrección debido a la perdida del gran gradiente en los extremos
				//cps(panel_idx) *= 1.25;

				if(!STEADY)
				{
					// backward differences for phi evolution
					// which is roughly equal to mu evolution (see tfg writeup)
					double partial_mu = 0.0;
					double sgn = 1.0;
					if(sln_hist.size() == 0 && !message_put)
					{
						std::cout << "WARNING: No sln_history!" << std::endl;
						message_put = true;
					}
					for(size_t i = 0; i < sln_hist.size(); i++)
					{
						partial_mu += sgn * (double)binom(sln_hist.size() - 1, i) * sln_hist[i](panel_idx);
						sgn *= -1.0;
					}

					// wake_scale is our timestep in essence
					partial_mu /= std::pow(wake_scale, sln_hist.size() - 1);

					// TODO: + or - (- seems to work better)
					// Also we may have to introduce the varition in phi due to change of boundary
					// conditions!!!!
					cps(panel_idx) -= 2.0 * partial_mu / Vref2;
				}

			}

		}

	}

	void PanelMethod::integrate_velocities(double wake_scale)
	{
		pos_history.clear();
		orient_history.clear();

		Vector3d body_pos; body_pos.setZero();
		Isometry3d body_orient; body_orient.setIdentity();

		// Consider an observer situated in the ground, seeing the wing move.
		// The (unperturbed) wake is the set of points that the trailing edge of the
		// wing touches as it moves.
		// the set of those points are generated by integrating:
		// dx/dt (inertial) = dx/dt (rotating) + Omega * x
		// Where dx / dt (rotating) = body_vel  by definition
		// While this problem has an analytical solution, it's way simpler
		// to numerically integrate the equations
		// (We in fact integrate backwards in time)
		// Note that first edge is trailing edge!
		for(Index zi = 1; zi < num_wake_edges; zi++)
		{
			double progress = (double)zi / (double)(num_wake_edges - 1);
			double arr_progress = (double)(zi - 1) / (double)(num_wake_edges - 1);

			double pos = 1.0 - std::cos(M_PI * arr_progress * 0.5);

			double speed = std::sin(M_PI * progress * 0.5);

			// Interpolated sampling for history
			// TODO: Maybe higher than linear is good?
			double array_prog = pos * (double)(num_wake_edges - 1);
			double step = 1.0 / (double)(num_wake_edges - 1);
			Index prog_sup = (size_t)(std::ceil(array_prog));
			Index prog_inf = (size_t)(std::floor(array_prog));
			double fac_sup =  array_prog - (double)prog_inf;
			double fac_inf = 1.0 - fac_sup;

			Vector3d vel = vel_history[prog_sup] * fac_sup + vel_history[prog_inf] * fac_inf;
			Vector3d omega = angvel_history[prog_sup] * fac_sup + angvel_history[prog_inf] * fac_inf;

			// Integrate velocity
			Vector3d inertial_vel = body_orient * vel;
			body_pos -= inertial_vel * speed * wake_scale;

			double omega_length = omega.norm();
			Vector3d omega_norm = omega / omega_length;

			// Integrate rotation
			if(omega_length > 0.0)
			{
				body_orient = body_orient.rotate(AngleAxisd(-omega_length * speed * wake_scale, omega_norm));
			}

			pos_history.push_back(body_pos);
			orient_history.push_back(body_orient);
		}
	}

	void PanelMethod::timestep(const Vector3d &cur_vel, const Vector3d &cur_angvel)
	{
		// Rebuild wake

		// Pop-out oldest velocity profile
		vel_history.pop_back();
		angvel_history.pop_back();

		// Push new velocity profile
		vel_history.push_front(cur_vel);
		angvel_history.push_front(cur_angvel);

		integrate_velocities(wake_scale);

		for(size_t i = 0; i < thin_wings.size(); i++)
		{
			// Transfer new wake geometry
			wakes[i].build_from_history(*thin_wings[i], *this);
			// Convect mus and build new influence coefficients
			wakes[i].mu_convect(i, *this);
		}

		build_dynamic(false);
		solve(false);

		for(size_t i = 0; i < thin_wings.size(); i++)
		{
			wakes[i].transfer_unsteady_solution(i, *this);
		}
	}

	void PanelMethod::transfer_solution_to_wake()
	{
		for(size_t i = 0; i < wakes.size(); i++)
		{
			wakes[i].inherit_solution(i, *this);
		}

	}

	void PanelMethod::compute_phis()
	{
		Eigen::ArrayXd phis_above(rhs.rows());
		Eigen::ArrayXd phis_below(rhs.rows());

		for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
		{
			for (size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
			{
				for (Index effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
				{
					double phi_above = 0.0;
					double phi_below = 0.0;
					auto effect_idx = (Index)(effect + geom_sizes[effect_geom]);
					// Effect of main panels
					for (Index cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
					{
						phi_above += induced_phi(cause_geom, cause, *thin_wings[effect_geom], effect, true);
						phi_below += induced_phi(cause_geom, cause, *thin_wings[effect_geom], effect, false);
					}
					// Effect of wake panels
					for(Index wake_cause = 0; wake_cause < thin_wings[cause_geom]->trailing_panels.cols(); wake_cause++)
					{
						phi_above += induced_phi_wake(cause_geom, wake_cause, *thin_wings[effect_geom], effect, true);
						phi_below += induced_phi_wake(cause_geom, wake_cause, *thin_wings[effect_geom], effect, false);
					}

					phis_above(effect_idx) = phi_above;
					phis_below(effect_idx) = phi_below;
				}
			}
	}

	phi_hist_above.push_front(phis_above);
	phi_hist_below.push_front(phis_below);

}

std::string PanelMethod::phis_to_string(size_t for_geom, bool include_freestream, bool above)
{
	std::stringstream o;
	if(include_freestream)
	{
		// TODO: Implement
		abort();
	}

	if(above)
	{
		o << phi_hist_above.front().block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
	}
	else
	{
		o << phi_hist_below.front().block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
	}
	return o.str();
}





