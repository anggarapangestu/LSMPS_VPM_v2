// #include "remeshing.hpp"

// double remeshing::basic(const double xgi, const double ygi, const double hG, const vector<int> &neighbor,
// 						const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &gpz)
// {
// 	double _gpz_grid = 0.0e0;
// 	double dex, dey, Wx, Wy;
// 	static int kernel_type = Pars::opt_remesh_W;

// 	for (int i = 0; i < neighbor.size(); i++)
// 	{
// 		dex = std::abs((xp[neighbor[i]] - xgi) / hG);
// 		dey = std::abs((yp[neighbor[i]] - ygi) / hG);
// 		d_base_remeshing.kernel(dex, dey, kernel_type, Wx, Wy);
// 		_gpz_grid += gpz[neighbor[i]] * Wx * Wy;
// 	}

// 	return _gpz_grid;
// }