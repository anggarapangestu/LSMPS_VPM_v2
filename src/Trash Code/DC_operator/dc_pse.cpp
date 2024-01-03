#include "../Eigen/Dense"
#include "base_dc.hpp"
// #include <armadillo>
using namespace Eigen;

std::vector<double> base_dc::Q22(const int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp,
								 const std::vector<double> &fp, const std::vector<std::vector<int>> &neighborList, const int l)
{
	// internal variables
	std::vector<double> Q22_dc(np); // returning values
	double Rc;						// core size multiplied by a constant[3,4,5]
	int n_neighbor;					// number of neighbor for each particle

	// // std::vector<std::vector<int>> neighbor_list(np, std::vector<int>());
	// // clock_t t1 = clock();
	// // for (size_t n = 0; n < pair_i.size(); n++)
	// // {
	// // 	int i = pair_i[n];
	// // 	int j = pair_j[n];
	// // 	neighbor_list[i].push_back(j);
	// // 	neighbor_list[j].push_back(i);
	// // }
	// // t1 = clock()-t1;
	// // printf("%f\n", (double)t1/CLOCKS_PER_SEC);

	// --> arrange matrix b(xp)
	// arma::mat b(l,1,arma::fill::zeros);
	// b.row(2) = 2; b.row(3) = 2;
	// VectorXd b = VectorXd::Zero(l);
	VectorXd b = VectorXd::Zero(l);
	b(2) = 2;
	b(3) = 2;

	// --> looping started <--
	for (size_t i = 0; i < np; i++)
	{
		std::vector<int> neigh = neighborList[i];
		n_neighbor = neigh.size();
		double xi = xp[i];
		double yi = yp[i];

		// --> arrange matrix E(xp) and V(xp)
		// arma::mat E(n_neighbor,n_neighbor,arma::fill::zeros); // or use sParse matrix SpMat<double>
		// arma::mat V(n_neighbor,l);
		// arma::mat dfp(1,n_neighbor);
		MatrixXd E = MatrixXd::Zero(n_neighbor, n_neighbor); // or use sParse matrix SpMat<double>
		MatrixXd V(n_neighbor, l);
		RowVectorXd dfp(n_neighbor);
		MatrixXd phi_e = MatrixXd::Zero(n_neighbor, n_neighbor); // (nxn), Gaussian kernel

		// #pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < n_neighbor; j++)
		{
			int idx = neigh[j];
			double xj = xp[idx];
			double yj = yp[idx];
			double dx = (xi - xj) / sp[i];
			double dy = (yi - yj) / sp[i];
			double R = std::sqrt(std::pow(dx * sp[i], 2) + std::pow(dy * sp[i], 2));
			dfp(j) = fp[idx] - fp[i];

			E(j, j) = (std::exp(-std::pow(R / sp[i], 2) * 0.5));	 // sqrt of Gaussian kernel
			phi_e(j, j) = (std::exp(-std::pow(R / sp[i], 2) * 1.0)); // Gaussian kernel
			// #pragma opm parallel for
			for (int k = 0; k < l; k++)
			{
				V(j, k) = Vandermonde(dx, dy, sp[i], k);
			}
		}

		// --> arrange matrix B(xp)
		MatrixXd B = E.transpose() * V;
		// --> arrange matrix A(xp)
		MatrixXd A = B.transpose() * B;
		// --> calculate a(xp)
		// arma::mat L  = chol(A,"lower"); // use Choleski decomposition
		// arma::mat Li = solve(L,b);
		// arma::mat a  = solve(L.t(),Li);
		VectorXd a = A.colPivHouseholderQr().solve(b);
		// A.reset(); B.reset();

		// VectorXd phi_e  = arma::pow(E,2); // (nxn), Gaussian kernel
		VectorXd K_e = V * a;		   // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		VectorXd eta_e = phi_e * K_e;  // (nxn)(nx1) = (nx1), DC kernel
		VectorXd sumQ22 = dfp * eta_e; // (1xn)(nx1) = (1x1)

		// Q22_dc[i] = Pars::vis * arma::accu(sumQ22) / std::pow(sp[i], 2)
		Q22_dc[i] = Pars::vis * sumQ22.sum() / std::pow(sp[i], 2);
	}

	// -- delete internal variables
	// // for (size_t i = 0; i < np; i++)
	// // {
	// // 	neighbor_list[i].clear();
	// // }
	// // neighbor_list.clear();

	return Q22_dc;
}
