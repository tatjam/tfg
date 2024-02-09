#include "PanelMethod.h"

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

Eigen::Vector3d PanelMethod::induced_vel_wake(const Wake &wake, Eigen::Index cause_trailing, const ThinWing &effect,
											  Eigen::Index effect_panel)
{
	Vector3d center = get_center(effect, effect_panel);
	return induced_vel_wake(wake, cause_trailing, center);
}

Eigen::Vector3d PanelMethod::induced_vel_wake(const Wake &wake, Eigen::Index cause_trailing, const Vector3d &pos)
{
	Vector3d acc = Vector3d(0, 0, 0);
	for(Index i = 0; i < wake.num_edges - 1; i++)
	{
		Matrix<double, 3, 4> verts;
		Index cause_panel = cause_trailing * (wake.num_edges - 1) + i;
		verts.col(0) = wake.vertices.col(wake.quads(0, cause_panel));
		verts.col(1) = wake.vertices.col(wake.quads(1, cause_panel));
		verts.col(2) = wake.vertices.col(wake.quads(2, cause_panel));
		verts.col(3) = wake.vertices.col(wake.quads(3, cause_panel));

		Vector3d nrm = wake.normals.col(cause_panel);
		acc += induced_vel_verts(verts, nrm, pos);
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
										  Eigen::Index effect_panel)
{
	Vector3d acc = induced_vel_wake(wake, cause_trailing, effect, effect_panel);
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

void PanelMethod::build_rhs()
{
	for(size_t geom = 0; geom < thin_wings.size(); geom++)
	{
		for(size_t panel = 0; panel < thin_wings[geom]->quads.cols(); panel++)
		{
			Vector3d norm = thin_wings[geom]->normals.col(panel);
			rhs(panel + geom_sizes[geom], 0) = norm.dot(body_vel);
		}
	}
}

void PanelMethod::build_dynamic_matrix_steady()
{
	// We simply prescribe the same circulation to the wake than in the trailing
	// edge panels, thus the wake become "integrated" into the contribution
	// from the last panel.

	for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
	{
		for(size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
		{
			for(Index effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
			{
				for(Index cause_trail = 0; cause_trail < thin_wings[cause_geom]->trailing_panels.rows(); cause_trail++)
				{
					Index cause = thin_wings[cause_geom]->trailing_panels(cause_trail);

					double induced = induced_norm_vel_wake(wakes[cause_geom], cause_trail, *thin_wings[effect_geom], effect);
					auto effect_idx = (Index)(effect + geom_sizes[effect_geom]);
					auto cause_idx = (Index)(cause + geom_sizes[cause_geom]);

					dynamic_matrix(effect_idx, cause_idx) = induced;
				}
			}
		}
	}




}

void PanelMethod::build_dynamic_steady()
{
	build_rhs();
	build_dynamic_matrix_steady();

}


void PanelMethod::build_dynamic()
{

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
	o << solution.block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
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
	body_vel = Vector3d(0, 0, 1);
	omega = Vector3d(0, 0, 0);

}

void PanelMethod::solve()
{
	solution = (geometry_matrix + dynamic_matrix).colPivHouseholderQr().solve(rhs);

}

void PanelMethod::shed_initial_wake(Eigen::Index num_wake_panels, double wake_distance)
{
	wakes.clear();

	for(const auto& wing : thin_wings)
	{
		Wake& w = wakes.emplace_back();
		w.shed_from(*wing, wake_distance, num_wake_panels, body_vel, omega);
	}

}

std::string PanelMethod::wake_geom_to_string(size_t for_geom)
{
	return wakes[for_geom].quads_to_string();
}

Eigen::Vector3d
PanelMethod::induced_vel_verts(const Matrix<double, 3, 4> &vertices, const Vector3d& normal, const Vector3d &pos)
{
	// We first of all create the transform matrix that goes from global coordinate
	// to coordinates in which the z axis is the panel normal, and the x y axes are arbitrary
	Vector3d center = vertices.col(0).matrix();
	center += vertices.col(1).matrix();
	center += vertices.col(2).matrix();
	center += vertices.col(3).matrix();
	center /= 4.0;

	// Because we don't really care about the orientation of x and y, we can simply
	// use the axis-to-axis rotation quaternion
	Quaterniond quat; quat.setFromTwoVectors(normal, Vector3d(0, 0, 1));
	Isometry3d transform; transform.setIdentity();

	// Note that this implies that quat is applied before the translation!
	transform = transform * quat * Translation3d(-center);

	Vector3d p = transform * pos;

	// Vertices in new coordinate system
	Vector3d p1 = transform * vertices.col(0);
	Vector3d p2 = transform * vertices.col(1);
	Vector3d p3 = transform * vertices.col(2);
	Vector3d p4 = transform * vertices.col(3);

	double r1 = (p1 - p).norm();
	double r2 = (p2 - p).norm();
	double r3 = (p3 - p).norm();
	double r4 = (p4 - p).norm();

	Vector3d out;

	// WARNING: There's an errata in the bibliography, which states that
	// the denominators contain a substraction instead of addition (ie, r1 * r2 - (...)) instead of what's given!
	double den1 =
			(r1 * r2 * (r1 * r2 + ((p(0) - p1(0)) * (p(0) - p2(0)) + (p(1) - p1(1)) * (p(1) - p2(1)) + p(2) * p(2))));
	double den2 =
			(r2 * r3 * (r2 * r3 + ((p(0) - p2(0)) * (p(0) - p3(0)) + (p(1) - p2(1)) * (p(1) - p3(1)) + p(2) * p(2))));
	double den3 =
			(r3 * r4 * (r3 * r4 + ((p(0) - p3(0)) * (p(0) - p4(0)) + (p(1) - p3(1)) * (p(1) - p4(1)) + p(2) * p(2))));
	double den4 =
			(r4 * r1 * (r4 * r1 + ((p(0) - p4(0)) * (p(0) - p1(0)) + (p(1) - p4(1)) * (p(1) - p1(1)) + p(2) * p(2))));

	out(0) = p(2) * (p1(1) - p2(1)) * (r1 + r2) / den1;
	out(0) += p(2) * (p2(1) - p3(1)) * (r2 + r3) / den2;
	out(0) += p(2) * (p3(1) - p4(1)) * (r3 + r4) / den3;
	out(0) += p(2) * (p4(1) - p1(1)) * (r4 + r1) / den4;
	out(0) /= 4.0 * M_PI;

	out(1) = -p(2) * (p1(0) - p2(0)) * (r1 + r2) / den1;
	out(1) -= p(2) * (p2(0) - p3(0)) * (r2 + r3) / den2;
	out(1) -= p(2) * (p3(0) - p4(0)) * (r3 + r4) / den3;
	out(1) -= p(2) * (p4(0) - p1(0)) * (r4 + r1) / den4;
	out(1) /= 4.0 * M_PI;

	out(2) = ((p(0) - p2(0)) * (p(1) - p1(1)) - (p(0) - p1(0)) * (p(1) - p2(1))) * (r1 + r2) / den1;
	out(2) += ((p(0) - p3(0)) * (p(1) - p2(1)) - (p(0) - p2(0)) * (p(1) - p3(1))) * (r2 + r3) / den2;
	out(2) += ((p(0) - p4(0)) * (p(1) - p3(1)) - (p(0) - p3(0)) * (p(1) - p4(1))) * (r3 + r4) / den3;
	out(2) += ((p(0) - p1(0)) * (p(1) - p4(1)) - (p(0) - p4(0)) * (p(1) - p1(1))) * (r4 + r1) / den4;
	out(2) /= 4.0 * M_PI;

	Isometry3d itform; itform.setIdentity();
	itform = itform * quat;
	// Now untransform it, only the rotational component of course
	out = itform.inverse() * out;

	return out;
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

void PanelMethod::compute_cps(double epsilon)
{
	// This value requires some fine-tuning for certain geometries
	// TODO: Better method instead of this brute-force one

	// Panel-on-panel effect
	for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
	{
			for (Index effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
			{
				Vector3d center = get_center(*thin_wings[effect_geom], effect);
				// TODO: Signs
				Vector3d freestream = -body_vel;
				freestream -= omega.cross(center);
				Vector3d total_ind_top; total_ind_top.setZero();
				Vector3d total_ind_bottom; total_ind_bottom.setZero();

				Vector3d normal = thin_wings[effect_geom]->normals.col(effect);
				Vector3d above = center + normal * epsilon;
				Vector3d below = center - normal * epsilon;

				for (size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
				{
					// Panel-on-panel
					for (Index cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
					{
						Vector3d ind_top = induced_vel(*thin_wings[cause_geom], cause, above);
						ind_top *= solution(geom_sizes[cause_geom] + cause);
						Vector3d ind_bottom = induced_vel(*thin_wings[cause_geom], cause, below);
						ind_bottom *= solution(geom_sizes[cause_geom] + cause);
						total_ind_top += ind_top;
						total_ind_bottom += ind_bottom;
					}
					// Wake-on-panel effect
					for(Index trail_cause = 0; trail_cause < thin_wings[cause_geom]->trailing_panels.rows(); trail_cause++)
					{
						Vector3d tangential_top = induced_vel_wake(wakes[cause_geom], trail_cause, above);
						tangential_top *= solution(geom_sizes[cause_geom] + thin_wings[cause_geom]->trailing_panels(trail_cause));
						Vector3d tangential_bottom = induced_vel_wake(wakes[cause_geom], trail_cause, below);
						tangential_bottom *= solution(geom_sizes[cause_geom] + thin_wings[cause_geom]->trailing_panels(trail_cause));
						total_ind_top += tangential_top;
						total_ind_bottom += tangential_bottom;
					}
				}


				// Thus by applying Bernouilli we arrive at the following expression
				// Thus 2 * (P_top - P_bottom) / rho = v_bottom^2 - v_top^2
				// By the definition of c_P
				// c_P = 2 * (P_top - P_bottom) / rho / v_infty^2
				//     = (v_bottom^2 - v_top^2) / v_infty^2
				// Where the bottom and superior terms also include the freestream term,
				// which we approximate as the same for both
				Vector3d v_top = total_ind_top + freestream;
				Vector3d v_bottom = total_ind_bottom + freestream;
				double cp = (v_bottom.squaredNorm() - v_top.squaredNorm()) / freestream.squaredNorm();
				cps(geom_sizes[effect_geom] + effect) = cp;
		}
	}

}

Eigen::Vector3d PanelMethod::compute_aero_force()
{
	Vector3d acc; acc.setZero();

	for(size_t geom = 0; geom < thin_wings.size(); geom++)
	{
		for(Index panel = 0; panel < thin_wings[geom]->quads.cols(); panel++)
		{
			Vector3d nrm = thin_wings[geom]->normals.col(panel);
			// We approximate the area of the cuadrilateral as two triangles
			Vector3d v1 = thin_wings[geom]->vertices.col(thin_wings[geom]->quads(0, panel));
			Vector3d v2 = thin_wings[geom]->vertices.col(thin_wings[geom]->quads(1, panel));
			Vector3d v3 = thin_wings[geom]->vertices.col(thin_wings[geom]->quads(2, panel));
			Vector3d v4 = thin_wings[geom]->vertices.col(thin_wings[geom]->quads(3, panel));
			double area = 0.5 * (v2 - v1).cross(v3 - v1).norm() + 0.5 * (v4 - v2).cross(v4 - v3).norm();

			Vector3d force = nrm * cps(geom_sizes[geom] + panel);
			acc += force * area;
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
			Vector3d freestream = -body_vel;
			freestream -= omega.cross(pos);

			for (size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
			{
				// Panel effect
				for (Index cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
				{
					Vector3d ind = induced_vel(*thin_wings[cause_geom], cause, pos);
					ind *= solution(geom_sizes[cause_geom] + cause);
					total_ind += ind;
				}
				// Wake effect
				for(Index trail_cause = 0; trail_cause < thin_wings[cause_geom]->trailing_panels.rows(); trail_cause++)
				{
					Vector3d ind = induced_vel_wake(wakes[cause_geom], trail_cause, pos);
					ind *= solution(geom_sizes[cause_geom] + thin_wings[cause_geom]->trailing_panels(trail_cause));
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








