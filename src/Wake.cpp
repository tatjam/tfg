#include "Wake.h"
#include "PanelMethod.h"

using namespace Eigen;

void Wake::build_initial_geometry(const ThinWing &wing, const PanelMethod &method)
{
	// Create vertices, starting from the trailing edge
	// We index similarly to wings, (xi * num_edges + zi)
	// where xi represents moving along the trialing edge, and zi moving back into the wake
	size_t num_trailing_edge = wing.trailing_edge.rows();
	vertices = Array3Xd(3, num_trailing_edge * method.num_wake_edges);
	for(Index i = 0; i < wing.trailing_edge.rows(); i++)
	{
		vertices.col(i * method.num_wake_edges + 0) = wing.transform * wing.vertices.col(wing.trailing_edge(i));
	}

	// Generate the rectangles, same as for the wing
	quads = Array4Xi(4, (num_trailing_edge - 1) * (method.num_wake_edges - 1));

	for(int xi = 0; xi < num_trailing_edge - 1; xi++)
	{
		for(int zi = 0; zi < method.num_wake_edges - 1; zi++)
		{
			quads(0, xi * (method.num_wake_edges - 1) + zi) = xi * method.num_wake_edges + zi;
			quads(1, xi * (method.num_wake_edges - 1) + zi) = (xi + 1) * method.num_wake_edges + zi;
			quads(2, xi * (method.num_wake_edges - 1) + zi) = (xi + 1) * method.num_wake_edges + zi + 1;
			quads(3, xi * (method.num_wake_edges - 1) + zi) = xi * method.num_wake_edges + zi + 1;
		}
	}

}

void Wake::timestep(const ThinWing &wing, const PanelMethod& panels)
{
	build_from_history(wing, panels);
	generate_normals();
}

void Wake::build_from_history(const ThinWing &wing, const PanelMethod& panels)
{
	size_t num_trailing_edge = wing.trailing_edge.rows();
	// We start at the last stored time-step and progress forwards
	for(Index zi = 1; zi < panels.num_wake_edges; zi++)
	{
		// This is kind of "velocity" cosine sampling
		double progress = (double)zi / (double)(panels.num_wake_edges - 1);
		//double pos = 1.0 - std::cos(M_PI * progress * 0.5);
		// This is kind of like a "time-step", that increases further away from the body
		double speed = std::sin(M_PI * progress * 0.5);

		Vector3d body_pos = panels.pos_history[zi - 1];
		Isometry3d body_orient = panels.orient_history[zi - 1];

		for(Index xi = 0; xi < num_trailing_edge; xi++)
		{
			Vector3d vertex_pos = wing.transform * wing.vertices.col(wing.trailing_edge(xi)).matrix();
			Vector3d rotated = body_orient * vertex_pos;
			Vector3d rel_pos = body_pos + rotated;
			vertices.col(xi * panels.num_wake_edges + zi) = rel_pos;
		}
	}

}


