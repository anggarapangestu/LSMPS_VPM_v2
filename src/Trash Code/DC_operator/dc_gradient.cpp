#include "../Eigen/Dense"
#include "base_dc.hpp"
// #include <armadillo>

std::vector<std::vector<double>> base_dc::Q11(const int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp,
											  const std::vector<double> &fp, const std::vector<std::vector<int>> &neighborList, const int l)
{
	std::vector<std::vector<double>> Q11(np, std::vector<double>(2));
	// internal variables
	double Rc;		// core size multiplied by a constant[3,4,5]
	int n_neighbor; // number of neighbor for each particle

	// // std::vector<std::vector<int>> neighborList;
	// // for (int i = 0; i < np; i++)
	// // {
	// // 	neighborList.push_back(vector<int>());
	// // }
	// // for (int n = 0; n < pair_i.size(); n++)
	// // {
	// // 	int i = pair_i[n];
	// // 	int j = pair_j[n];
	// // 	neighborList[i].push_back(j);
	// // 	neighborList[j].push_back(i);
	// // }

	/*arrange matrix b(xp)*/
	Eigen::VectorXd b10 = Eigen::VectorXd::Zero(l);
	b10(1) = -1;
	Eigen::VectorXd b01 = Eigen::VectorXd::Zero(l);
	b01(2) = -1;

	// /*looping started ...*/
	for (int i = 0; i < np; i++)
	{
		std::vector<int> neigh = neighborList[i];
		n_neighbor = neigh.size();
		double xi = xp[i];
		double yi = yp[i];

		/*arrange matrix E(xp) and V(xp)*/
		Eigen::MatrixXd E = Eigen::MatrixXd::Zero(n_neighbor, n_neighbor); // or use sParse matrix SpMat<double>
		Eigen::MatrixXd V(n_neighbor, l);
		Eigen::RowVectorXd df(n_neighbor);
		Eigen::MatrixXd phi_e = Eigen::MatrixXd::Zero(n_neighbor, n_neighbor); // (nxn), Gaussian kernel

		// #pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < n_neighbor; j++)
		{
			int idx = neigh[j];
			double xj = xp[idx];
			double yj = yp[idx];
			double dx = (xi - xj) / sp[i];
			double dy = (yi - yj) / sp[i];
			double R = std::sqrt(std::pow(dx * sp[i], 2) + std::pow(dy * sp[i], 2));

			// df_x(0,j)		= du(idx)+du(i);
			// df_y(0,j)		= dv(idx)+dv(i);
			df(0, j) = fp[idx] + fp[i];

			E(j, j) = std::exp(-std::pow(R / sp[i], 2) * 0.5);	 // sqrt of Gaussian kernel
			phi_e(j, j) = std::exp(-std::pow(R / sp[i], 2) * 1.0); // Gaussian kernel

			for (int k = 0; k < l; k++)
			{
				V(j, k) = Vandermonde_zero(dx, dy, sp[i], k);
			}
		}

		/*arrange matrix B(xp)*/
		Eigen::MatrixXd B = E.transpose() * V;
		/*arrange matrix A(xp)*/
		Eigen::MatrixXd A = B.transpose() * B;
		/*calculate a(xp) */
		Eigen::VectorXd a10 = A.colPivHouseholderQr().solve(b10);
		Eigen::VectorXd a01 = A.colPivHouseholderQr().solve(b01);

		Eigen::MatrixXd K_e_10 = V * a10;				// (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		Eigen::MatrixXd eta_e_10 = phi_e * K_e_10;		// (nxn)(nx1) = (nx1), DC kernel
		Eigen::MatrixXd sumQ_10 = df /*_y*/ * eta_e_10; // (1xn)(nx1) = (1x1)
		Eigen::MatrixXd K_e_01 = V * a01;				// (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		Eigen::MatrixXd eta_e_01 = phi_e * K_e_01;		// (nxn)(nx1) = (nx1), DC kernel
		Eigen::MatrixXd sumQ_01 = df /*_x*/ * eta_e_01; // (1xn)(nx1) = (1x1)

		Q11[i][0] = sumQ_10.sum() / pow(sp[i], 1);
		Q11[i][1] = sumQ_01.sum() / pow(sp[i], 1);
	}

	return Q11;
}
