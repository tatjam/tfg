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
	for(Index i = 0; i < num_wake_edges - 1; i++)
	{
		Matrix<double, 3, 4> verts;
		Index cause_panel = cause_trailing * (num_wake_edges - 1) + i;
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
			rhs(panel + geom_sizes[geom], 0) = norm.dot(vel_history.front());
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
}

void PanelMethod::solve()
{
	solution = (geometry_matrix + dynamic_matrix).colPivHouseholderQr().solve(rhs);

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
	/*
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

	return out;*/
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

Eigen::Vector3d PanelMethod::compute_aero_force()
{
	Vector3d acc; acc.setZero();


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


	// Center line lift

	/*for(Index i = 0; i < thin_wings[0]->num_chorwise - 1; i++)
	{
		Index idx = thin_wings[0]->num_spanwise / 2 * (thin_wings[0]->num_chorwise - 1) + i;
		Vector3d nrm =thin_wings[0]->normals.col(idx);
		double area = get_area(*thin_wings[0], idx);
		Vector3d a = thin_wings[0]->vertices.col(thin_wings[0]->quads(0, idx));
		Vector3d b = thin_wings[0]->vertices.col(thin_wings[0]->quads(3, idx));
		double length = (b - a).norm();
		Vector3d force = nrm * cps(idx);
		acc += force * length;
	}*/

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
void PanelMethod::compute_cps()
{
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
			// phi = Ax^2 + Bxy + Cy^2 + Dx + Ey + F
			// Which has gradient
			// (2Ax + By + D, Bx + 2Cy + E)
			// Thus at the origin (x = 0, y = 0) the gradient is simply
			// (D, E)
			// Which can be converted back into 3D by inverse transformation
			// We use the least squares method to find such fit

			MatrixXd lsq = MatrixXd(5 /*data points*/, 6 /* parameters */);
			for(Index p = 0; p < 5; p++)
			{
				lsq(p, 0) = ps(0, p)*ps(0, p); // Ax^2
				lsq(p, 1) = ps(0, p)*ps(1, p); // Bxy
				lsq(p, 2) = ps(1, p)*ps(1, p); // Cy^2
				lsq(p, 3) = ps(0, p); //Dx
				lsq(p, 4) = ps(1, p); //Ey
				lsq(p, 5) = 1.0; // F
			}

			// Solve the least squares problem (see "Introducción a los métodos numéricos", José Antonio Ezquerro)
			Vector<double, 5> rhs;
			for(Index p = 0; p < 4; p++)
			{
				rhs(p) = mus(p);
			}
			rhs(4) = self_sol;

			// 6x5 * 5x1 = 6x1
			Vector<double, 6> rhs_t = lsq.transpose() * rhs;
			// 6x5 * 5x6 = 6x6
			Matrix<double, 6, 6> A = lsq.transpose() * lsq;
			Vector<double, 6> parameters = A.fullPivLu().solve(rhs_t);

			//
			Vector3d grad = Vector3d(parameters(3), parameters(4), 0.0);
			grad = tinv * grad;

			// Now,
			cps(geom_sizes[effect_geom] + effect) = -2.0 * (freestream.dot(grad)) / freestream.squaredNorm();
		}
	}

}

void PanelMethod::build_wakes_from_history()
{

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
	// Pop-out oldest velocity profile
	vel_history.pop_back();
	angvel_history.pop_back();

	// Push new velocity profile
	vel_history.push_front(cur_vel);
	angvel_history.push_front(cur_angvel);

	integrate_velocities(wake_scale);

	for(size_t i = 0; i < thin_wings.size(); i++)
	{
		wakes[i].timestep(*thin_wings[i], *this);
	}

}








