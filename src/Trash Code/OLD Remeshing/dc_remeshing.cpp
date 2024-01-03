// #include "remeshing.hpp"
// #include "../Eigen/Dense" // linear algebra library
// using namespace Eigen;	// for linear algebra operation

// // =============================================================
// // =============================================================
// std::vector<double> remeshing::dc_remeshing(const Particle &p_scatter, const Particle &p_grid, const int l,
// 											const std::vector<std::vector<int>> neighborList)
// {
// 	// // -- define local variables pointing to particle class
// 	// const std::vector<double> &_x_grd = p_grid.x;
// 	// const std::vector<double> &_y_grd = p_grid.y;
// 	// const std::vector<double> &_s_grd = p_grid.s;
// 	// const std::vector<double> &_f_grd = p_grid.gz;
// 	// const int _n_grd = p_grid.num;

// 	// const std::vector<double> &_x_par = p_scatter.x;
// 	// const std::vector<double> &_y_par = p_scatter.y;
// 	// const std::vector<double> &_s_par = p_scatter.s;
// 	// const std::vector<double> &_f_par = p_scatter.gz;
// 	// const int _n_par = p_scatter.num;

// 	// -- internal variable
// 	std::vector<double> _Q00(p_grid.x.size());

// 	/*looping started ...*/
// 	for (size_t i = 0; i < p_grid.x.size(); i++)
// 	{
// 		std::vector<int> _neighbor = neighborList[i];
// 		// -- perform basic remeshing for particle outside body
// 		if (/*sign_par[i] == true ||*/ sign_par[i] == false)
// 		{
// 			_Q00[i] = basic(p_grid.x[i], p_grid.y[i], p_grid.s[i], _neighbor, p_scatter.x, p_scatter.y, p_scatter.gz);
// 			// continue;
// 		}
// 		// -- perform DC remeshing for particle inside body
// 		else
// 		{
// 			int n_neighbor = _neighbor.size();
// 			double xi = p_grid.x[i];
// 			double yi = p_grid.y[i];

// 			/*arrange matrix E(xp) and V(xp)*/
// 			Eigen::MatrixXd phi_e = Eigen::MatrixXd::Zero(n_neighbor, n_neighbor); // or use sparse matrix SpMat<double> from Armadillo
// 			Eigen::MatrixXd E = Eigen::MatrixXd::Zero(n_neighbor, n_neighbor);	 // or use sparse matrix SpMat<double> from Armadillo
// 			Eigen::MatrixXd V(n_neighbor, l);									   // vandermonde matrix
// 			Eigen::RowVectorXd df(n_neighbor);									   // to store difference of intensity

// 			/*arrange matrix b(xp)*/
// 			Eigen::VectorXd b00 = Eigen::VectorXd::Zero(l);
// 			b00(0) = 1;
// 			Eigen::VectorXd eta_a = Eigen::VectorXd::Zero(n_neighbor); // correction function

// 			// #pragma omp parallel for schedule(dynamic)
// 			for (size_t j = 0; j < n_neighbor; j++)
// 			{
// 				int idx = _neighbor[j];
// 				double xj = p_scatter.x[idx];
// 				double yj = p_scatter.y[idx];
// 				double dx = (xi - xj) / p_grid.s[i];
// 				double dy = (yi - yj) / p_grid.s[i];
// 				double R = std::sqrt(std::pow(dx * p_grid.s[i], 2) + std::pow(dy * p_grid.s[i], 2));

// 				df(j) = p_scatter.gz[idx] - 0.0e0;

// 				E(j, j) = std::exp(-std::pow(R / p_grid.s[i], 2) * 0.5);	 // sqrt of Gaussian kernel
// 				phi_e(j, j) = std::exp(-std::pow(R / p_grid.s[i], 2) * 1.0); // Gaussian kernel
// 				// #pragma omp parallel for
// 				for (size_t k = 0; k < l; k++)
// 				{
// 					V(j, k) = d_base_dc.Vandermonde_zero(dx, dy, p_grid.s[i], k);
// 				}

// 				// if (static_cast<float>(R) <= 1.0e-6) // determine the same particle
// 				// {
// 				// 	// b00(0) = 0.0e0;
// 				// 	// eta_a(j) = 1;
// 				// 	// for (size_t k = 1; k < l; k++)
// 				// 	// {
// 				// 	// 	b00(k) = -d_base_dc.Vandermonde_zero(dx, dy, p_grid.s[i], k);
// 				// 	// }
// 				// 	// continue;
// 				// }
// 			}

// 			/*arrange matrix B(xp)*/
// 			Eigen::MatrixXd B = E.transpose() * V;
// 			/*arrange matrix A(xp)*/
// 			Eigen::MatrixXd A = B.transpose() * B;
// 			/*calculate a(xp) */
// 			Eigen::VectorXd a00 = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b00); // least square solver
// 			// Eigen::VectorXd a00 = A.fullPivHouseholderQr().solve(b00);

// 			Eigen::VectorXd K_e_00 = V * a00;		   // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
// 			Eigen::VectorXd eta_e_00 = phi_e * K_e_00; // (nxn)(nx1) = (nx1), DC kernel
// 			VectorXd sumQ = df * (eta_e_00 + eta_a);   // (1xn)(nx1) = (1x1)
// 			_Q00[i] = sumQ.sum();
// 		}
// 	}

// 	return _Q00;
// }

// // // =============================================================
// // // =============================================================
// // double remeshing::quartic_spline(const double r, const double eps)
// // {
// // 	double phi_tilde;
// // 	double s = r / (eps * 3.5 / 2); // eps*3.5 = cutoff radius !!!
// // 	if (s <= 1)
// // 	{
// // 		phi_tilde = -3 * pow(s, 4) + 8 * pow(s, 3) - 6 * pow(s, 2) + 1;
// // 	}
// // 	else
// // 	{
// // 		phi_tilde = 0.0e0;
// // 	}
// // 	return phi_tilde;
// // }
