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
	rhs = MatrixXd(total_size, 1);

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

double
PanelMethod::induced_norm_vel(const ThinWing& cause, Index cause_panel, const ThinWing& effect, Index effect_panel)
{
	Vector3d center = effect.transform * effect.vertices.col(effect.quads(0, effect_panel)).matrix();
	center += effect.transform * effect.vertices.col(effect.quads(1, effect_panel)).matrix();
	center += effect.transform * effect.vertices.col(effect.quads(2, effect_panel)).matrix();
	center += effect.transform * effect.vertices.col(effect.quads(3, effect_panel)).matrix();
	center /= 4.0;

	Vector3d ivel = induced_vel(cause, cause_panel, center);

	Vector3d normal = effect.normals.col(effect_panel).matrix();

	return ivel.dot(normal);
}

// Velocity induced by an (assumed planar) doublet sheet, doublets oriented such that they point in the
// direction of the panel normal.
// (Developed from Low Speed Aerodynamics (...) by Joseph Katz and Allen Plotkin)
// This function is singular on the quadrilateral edges (technically the flat projection)
Eigen::Vector3d PanelMethod::induced_vel(const ThinWing &cause, Index cause_panel, const Vector3d &pos)
{
	// We first of all create the transform matrix that goes from global coordinate
	// to coordinates in which the z axis is the panel normal, and the x y axes are arbitrary
	Vector3d normal = cause.normals.col(cause_panel);
	Vector3d center = cause.vertices.col(cause.quads(0, cause_panel)).matrix();
	center += cause.vertices.col(cause.quads(1, cause_panel)).matrix();
	center += cause.vertices.col(cause.quads(2, cause_panel)).matrix();
	center += cause.vertices.col(cause.quads(3, cause_panel)).matrix();
	center /= 4.0;

	// Because we don't really care about the orientation of x and y, we can simply
	// use the axis-to-axis rotation quaternion
	Quaterniond quat; quat.setFromTwoVectors(normal, Vector3d(0, 0, 1));
	Isometry3d transform; transform.setIdentity();

	// Note that this implies that quat is applied before the translation!
	transform = transform * quat * Translation3d(-center);

	Vector3d p = transform * pos;

	// Vertices in new coordinate system
	Vector3d p1 = transform * cause.transform * cause.vertices.col(cause.quads(0, cause_panel)).matrix();
	Vector3d p2 = transform * cause.transform * cause.vertices.col(cause.quads(1, cause_panel)).matrix();
	Vector3d p3 = transform * cause.transform * cause.vertices.col(cause.quads(2, cause_panel)).matrix();
	Vector3d p4 = transform * cause.transform * cause.vertices.col(cause.quads(3, cause_panel)).matrix();

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

void PanelMethod::build_dynamic()
{
	build_rhs();

}

std::string PanelMethod::geometry_matrix_to_string()
{
	std::stringstream o;
	o << geometry_matrix;
	return o.str();
}

std::string PanelMethod::solution_to_string(size_t for_geom)
{
	std::stringstream o;
	o << solution.block(geom_sizes[for_geom], 0, geom_sizes[for_geom + 1], 1).transpose();
	return o.str();
}

PanelMethod::PanelMethod()
{
	body_vel = Vector3d(0, 0, 1);
	omega = Vector3d(0, 0, 0);

}

void PanelMethod::solve()
{
	solution = geometry_matrix.colPivHouseholderQr().solve(rhs);

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




