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
PanelMethod::induced_norm_vel(const ThinWing& cause, size_t cause_panel, const ThinWing& effect, size_t effect_panel)
{

	return 0.0;

}

Eigen::Vector3d PanelMethod::induced_vel(const ThinWing &cause, size_t cause_panel, const Vector3d &pos)
{
	return Eigen::Vector3d();
}
