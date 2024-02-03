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

	geometry_matrix = MatrixXd(total_size, total_size);

	for(size_t effect_geom = 0; effect_geom < thin_wings.size(); effect_geom++)
	{
		for(size_t cause_geom = 0; cause_geom < thin_wings.size(); cause_geom++)
		{
			for(size_t effect = 0; effect < thin_wings[effect_geom]->quads.cols(); effect++)
			{
				for(size_t cause = 0; cause < thin_wings[cause_geom]->quads.cols(); cause++)
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
PanelMethod::induced_norm_vel(const ThinWing& cause, Index cause_panel, const ThinWing& effect, size_t effect_panel)
{

	return 0.0;

}

// Velocity induced by an (assumed planar) doublet sheet, doublets oriented such that they point in the
// direction of the panel normal.
// (Developed from Low Speed Aerodynamics (...) by Joseph Katz and Allen Plotkin)
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
	Quaterniond quat; quat.setFromTwoVectors(Vector3d(0, 0, 1), normal);
	Isometry3d transform; transform.setIdentity();

	// Note that this implies that quat is applied before the translation!
	transform = transform * quat * Translation3d(-center);

	Vector3d pos_trans = transform * pos;

	// Vertices in new coordinate system
	Vector3d p1 = transform * cause.transform * cause.vertices.col(cause.quads(0, cause_panel)).matrix();
	Vector3d p2 = transform * cause.transform * cause.vertices.col(cause.quads(1, cause_panel)).matrix();
	Vector3d p3 = transform * cause.transform * cause.vertices.col(cause.quads(2, cause_panel)).matrix();
	Vector3d p4 = transform * cause.transform * cause.vertices.col(cause.quads(3, cause_panel)).matrix();

	double r1 = (p1 - pos_trans).norm();
	double r2 = (p2 - pos_trans).norm();
	double r3 = (p3 - pos_trans).norm();
	double r4 = (p4 - pos_trans).norm();






}
