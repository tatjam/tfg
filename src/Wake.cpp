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
	from_panel = ArrayXi((num_trailing_edge - 1) * (method.num_wake_edges - 1));

	for(int xi = 0; xi < num_trailing_edge - 1; xi++)
	{
		for(int zi = 0; zi < method.num_wake_edges - 1; zi++)
		{
			quads(0, xi * (method.num_wake_edges - 1) + zi) = xi * method.num_wake_edges + zi;
			quads(1, xi * (method.num_wake_edges - 1) + zi) = (xi + 1) * method.num_wake_edges + zi;
			quads(2, xi * (method.num_wake_edges - 1) + zi) = (xi + 1) * method.num_wake_edges + zi + 1;
			quads(3, xi * (method.num_wake_edges - 1) + zi) = xi * method.num_wake_edges + zi + 1;

			Index panel_idx = wing.trailing_panels(xi);
			from_panel(xi * (method.num_wake_edges - 1) + zi) = panel_idx;
		}
	}

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

	generate_normals();

}

void Wake::inherit_solution(size_t wake_idx, const PanelMethod &method)
{
	size_t num_trailing_edge = method.thin_wings[wake_idx]->trailing_panels.rows();
	mus = ArrayXd(quads.cols());
	influences = ArrayXd(quads.cols());

	ArrayXd stat_sample(num_trailing_edge);

	// Initial values don't need any advanced computation
	for(size_t xi = 0; xi < num_trailing_edge; xi++)
	{
		Index trail_panel = method.thin_wings[wake_idx]->trailing_panels(xi);
		double val = method.sln_hist.front()(method.geom_sizes[wake_idx] + trail_panel);
		stat_sample(xi) = val;
		for(size_t zi = 0; zi < method.num_wake_edges - 1; zi++)
		{
			mus(xi * (method.num_wake_edges - 1) + zi) = val;
			influences(xi * (method.num_wake_edges - 1) + zi) = 0.0;
		}
	}

	// Initial flat history
	for(size_t zi = 1; zi < method.num_wake_edges; zi++)
	{
		mu_history.push_back(stat_sample);
	}

}

void Wake::mu_convect(size_t wake_idx, const PanelMethod &method)
{
	size_t num_trailing_edge = method.thin_wings[wake_idx]->trailing_panels.rows();
	// New solution, pushed as 0s for now
	ArrayXd zeros = ArrayXd(num_trailing_edge);
	zeros.setZero();
	mu_history.pop_back();
	mu_history.push_front(zeros);
	for(Index zi = 1; zi < method.num_wake_edges; zi++)
	{
		double arr_progress = (double)(zi - 1) / (double)(method.num_wake_edges - 1);

		double pos = 1.0 - std::cos(M_PI * arr_progress * 0.5);

		// TODO: This could REALLY benefit from smoother interp
		double array_prog = pos * (double)(method.num_wake_edges - 1);
		Index prog_sup = (size_t)(std::ceil(array_prog));
		Index prog_inf = (size_t)(std::floor(array_prog));
		double fac_sup =  array_prog - (double)prog_inf;
		double fac_inf = 1.0 - fac_sup;

		for(Index xi = 0; xi < num_trailing_edge; xi++)
		{
			double mu = 0.0;

			double influence = 0.0;
			if(prog_inf == 0)
			{
				influence += fac_inf;
				// !!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!
				// Note that this value is used for later interpolation
				// so it's interesting to set it to the superior mu
				// (As otherwise we would have to interpolate again!)
				// TODO: THIS IS VERY IMPORTANT, DOCUMENT IT WELL!
				mu = mu_history[prog_sup](xi);
			}
			else
			{
				mu = mu_history[prog_sup](xi) * fac_sup + mu_history[prog_inf](xi) * fac_inf;
			}

			mus(xi * (method.num_wake_edges - 1) + zi - 1) = mu;
			influences(xi * (method.num_wake_edges - 1) + zi - 1) = influence;
		}
	}

}

void Wake::transfer_unsteady_solution(size_t wake_idx, const PanelMethod &method)
{
	size_t num_trailing_edge = method.thin_wings[wake_idx]->trailing_panels.rows();
	for(Index i = 0; i < mus.rows(); i++)
	{
		Index panel = from_panel(i);
		double sln = method.sln_hist.front()(method.geom_sizes[wake_idx] + panel);
		double infl = influences(i);
		mus(i) = mus(i) * (1.0 - infl) + sln * infl;
	}

	// Overwrite history too for proper convection
	ArrayXd overwrite_hist = ArrayXd(num_trailing_edge);
	for(Index i = 0; i < num_trailing_edge; i++)
	{
		Index panel = method.thin_wings[wake_idx]->trailing_panels(i);
		overwrite_hist(i) = method.sln_hist.front()(method.geom_sizes[wake_idx] + panel);
	}
	mu_history[0] = overwrite_hist;
}


