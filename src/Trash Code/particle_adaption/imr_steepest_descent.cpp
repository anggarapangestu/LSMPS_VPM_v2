#include "particle_adaption.hpp"

// ===================
// ===================
double particle_adaption::V1r(const double rpq, const double Dpq)
{
	double r = rpq / Dpq;
	if (r <= 0.5)
	{
		return 0.0e0;
	}
	else
	{
		return 0.8 * std::pow(2.5, 1 - (5 * r)) - std::pow(2.5, -4 * r);
	}
}
// ===================
// ===================
double particle_adaption::V2r(const double rpq, const double Dpq)
{
	double r = rpq / Dpq;
	if (r <= 0.5)
	{
		return 0.0e0;
	}
	else
	{
		return (pow(r, -2) / 2) + (pow(r, -6) / 6);
	}
}
// =================================================================================================
// =================================================================================================
void particle_adaption::steepestDescent()
{
	// == internal variables
	double minp = Pars::D_0 * Pars::r_scale; // for convergence criteria
	double minq = Pars::D_0 * Pars::r_scale; // for convergence criteria
	std::vector<double> Dpqmin(2);			 // for convergence criteria
	int minIndexp = 0;						 // for convergence criteria
	int minIndexq = 0;						 // for convergence criteria
	double xpmin, ypmin, xqmin, yqmin;		 // for convergence criteria

	double minDpq = 0.0e0;
	double xj, yj, xi, yi, dx, dr, Vaksen, Vr;

	// double sptemp;
	// double delta_rstar = 1.0e0;

	std::vector<int> pair;
	std::vector<double> Dpq(2, 0.0e0);
	std::vector<int> countiac(this->_npnew, 0);
	// std::vector<double> Wp(_npnew, 0.0e0);
	// std::vector<std::vector<double>> wp(this->_npnew, std::vector<double>(2, 0.0e0));
	double Wp[_npnew];	// interparticle energy
	double wp[_npnew][2]; // displacement's vector

	int l = Pars::l_0_0_3; // number moment condition

	// == internal dynamic variables
	int n_neighbor, idx; // number of neighbor for each particle
	double deltax, deltay, R, deltaR;
	// arma::cube Aijk(l,l,this->_npnew);
	// arma::cube phie_ijk;
	std::vector<std::vector<int>> neighborlist(this->_npnew, std::vector<int>());

	// --> arrange matrix b(xp)
	arma::mat b(l, 1, arma::fill::zeros);
	b(0, 0) = 1;
	arma::mat b10(l, 1, arma::fill::zeros);
	b10(1, 0) = -1;
	arma::mat b01(l, 1, arma::fill::zeros);
	b01(2, 0) = -1;

	printf("performing steepest descent algorithm...\n");
	// == processsing
	for (int i = 0; i < _npnew; i++)
	{
		_directFindMRinterp(_xp[i], _yp[i], _sp[i], _npnew, _sp, _xp, _yp, pair);
		countiac[i] = pair.size();
		deltaR = 1; // reset cutoff multiplier
		if (pair.size() <= l)
		{
			_directFind(_xp[i], _yp[i], _sp[i] * deltaR, _npnew, _sp, _xp, _yp, pair);
			deltaR += 0.5;
		}

		n_neighbor = pair.size();
		neighborlist[i] = pair;

		// -- finding convergence
		if (minp >= _sp[i])
		{
			minp = _sp[i];
			minIndexp = i;
		}

		xi = this->_xp[i];
		yi = this->_yp[i];
		// --> arrange matrix E(xp) and V(xp)
		arma::sp_mat E(n_neighbor, n_neighbor);		// or use sParse matrix SpMat<double>
		arma::sp_mat phi_e(n_neighbor, n_neighbor); // or use sParse matrix SpMat<double>
		arma::mat V(n_neighbor, l);
		arma::sp_mat dField(1, n_neighbor);

		for (int j = 0; j < n_neighbor; j++)
		{
			dx = std::sqrt(std::pow(_xp[i] - _xp[pair[j]], 2) + std::pow(_yp[i] - _yp[pair[j]], 2));
			Dpq[0] = this->_sp[i] / Pars::r_scale;
			Dpq[1] = this->_sp[pair[j]] / Pars::r_scale;
			minDpq = *std::min_element(Dpq.begin(), Dpq.end());
			dr = dx / minDpq;

			idx = pair[j];
			xj = this->_xp[idx];
			yj = this->_yp[idx];
			// deltax = (xi - xj)/this->_sp[i];
			// deltay = (yi - yj)/this->_sp[i];
			deltax = (xi - xj) / (minDpq * Pars::r_scale);
			deltay = (yi - yj) / (minDpq * Pars::r_scale);
			dField(j) = std::pow(minDpq, 2);
			E(j, j) = std::sqrt(std::abs(V2r(dx, minDpq))); // sqrt of Potential kernel
			phi_e(j, j) = V2r(dx, minDpq);					// Potential kernel

			// if(dr > 0.5)
			// {
			for (int k = 0; k < l; k++)
			{
				V(j, k) = d_base_dc.Vandermonde_zero(deltax, deltay, minDpq, k);
			}
			// }
		}
		pair.clear(); // !!!!!

		// --> arrange matrix B(xp)
		arma::mat B = E.t() * V;
		// --> arrange matrix A(xp)
		arma::mat A = B.t() * B;
		// --> calculate a(xp)
		// arma::mat L  = chol(A,"lower"); // use Choleski decomposition
		// arma::mat Li = solve(L,b);
		// arma::mat a  = solve(L.t(),Li);
		arma::mat a = solve(A, b);

		arma::mat K_e = V * a;			   // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		arma::mat eta_e = phi_e * K_e;	 // (nxn)(nx1) = (nx1), DC kernel
		arma::mat sumQ00 = dField * eta_e; // (1xn)(nx1) = (1x1)

		Wp[i] = sumQ00(0);
	}

	// -----------------------------------------------------------------------------------------
	for (int dci = 0; dci < neighborlist[minIndexp].size(); dci++)
	{
		if (minq > _sp[neighborlist[minIndexp][dci]] && neighborlist[minIndexp][dci] != minIndexp)
		{
			minq = _sp[neighborlist[minIndexp][dci]];
			minIndexq = neighborlist[minIndexp][dci];
		}
	}

	// Dpqmin[0] = minp/Pars::r_scale;
	// Dpqmin[1] = minq/Pars::r_scale;
	xpmin = this->_xp[minIndexp];
	ypmin = this->_yp[minIndexp];
	xqmin = this->_xp[minIndexq];
	yqmin = this->_yp[minIndexq];
	// double Dminpminq = *std::min_element(Dpqmin.begin(), Dpqmin.end());
	double Dminpminq = minp < minq ? minp / Pars::r_scale : minq / Pars::r_scale;
	this->_dcmin = std::sqrt(std::pow(xpmin - xqmin, 2) + std::pow(ypmin - yqmin, 2)) / Dminpminq;
	this->_minneighbor = *std::min_element(countiac.begin(), countiac.end());
	printf("fewest neighbour = %d \t", this->_minneighbor);
	printf("dcmin = %f\n", this->_dcmin);
	// -----------------------------------------------------------------------------------------

	// == processsing
	for (int i = 0; i < _npnew; i++)
	{
		// _directFindMRinterp(_xp[i], _yp[i], _sp[i], _npnew, _sp, _xp, _yp, pair);
		pair = neighborlist[i];

		n_neighbor = pair.size() /* > 20 ? 20 : pair.size()*/;

		xi = this->_xp[i];
		yi = this->_yp[i];
		// --> arrange matrix E(xp) and V(xp)
		arma::sp_mat E(n_neighbor, n_neighbor);		// or use sParse matrix SpMat<double>
		arma::sp_mat phi_e(n_neighbor, n_neighbor); // or use sParse matrix SpMat<double>
		arma::mat V(n_neighbor, l);
		arma::mat dField(1, n_neighbor);

		for (int j = 0; j < n_neighbor; j++)
		{
			dx = std::sqrt(std::pow(_xp[i] - _xp[pair[j]], 2) + std::pow(_yp[i] - _yp[pair[j]], 2));

			idx = pair[j];
			xj = this->_xp[idx];
			yj = this->_yp[idx];
			deltax = (xi - xj) / this->_sp[i];
			deltay = (yi - yj) / this->_sp[i];
			dField(j) = Wp[idx] + Wp[i];
			E(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 0.5);	 // sqrt of Gaussian kernel
			phi_e(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 1.0); // Gaussian kernel

			// if(dr >= 0.5)
			// {
			for (int k = 0; k < l; k++)
			{
				V(j, k) = d_base_dc.Vandermonde_zero(deltax, deltay, this->_sp[i], k);
			}
			// }
		}
		pair.clear(); // !!!!!

		// --> arrange matrix B(xp)
		mat B = E.t() * V;
		// --> arrange matrix A(xp)
		mat A = B.t() * B;
		// --> calculate a(xp)
		arma::mat a10 = solve(A, b10);
		arma::mat a01 = solve(A, b01);

		arma::mat K_e_10 = V * a10;			  // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		arma::mat eta_e_10 = phi_e * K_e_10;  // (nxn)(nx1) = (nx1), DC kernel
		arma::mat sumQ10 = dField * eta_e_10; // (1xn)(nx1) = (1x1)

		arma::mat K_e_01 = V * a01;			  // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		arma::mat eta_e_01 = phi_e * K_e_01;  // (nxn)(nx1) = (nx1), DC kernel
		arma::mat sumQ01 = dField * eta_e_01; // (1xn)(nx1) = (1x1)

		if (_Dtildetemp[i] > 2 * _Dptemp[i] || _Dptemp[i] > 0.4 * Pars::D_0)
		{
			// printf("...is not moved\n");
			wp[i][0] = 0.0e0;
			wp[i][1] = 0.0e0;
		}
		else
		{
			wp[i][0] = sumQ10(0) / this->_sp[i];
			wp[i][1] = sumQ01(0) / this->_sp[i];
		}
	}

	for (int i = 0; i < _npnew; i++)
	{
		_xp[i] -= wp[i][0] * 4.0e-1; // 4.0e-1 is constant step size
		_yp[i] -= wp[i][1] * 4.0e-1; // 4.0e-1 is constant step size
	}

	// // == saving temporary data
	// std::ofstream outputs;
	// outputs.open("output/displacement.dat");
	// for (int i = 0; i < _npnew; i++)
	// {
	// 	outputs << _xp[i] << " " << _yp[i] << " " << -1.0e0 * wp[i][0] << " " << -1.0e0 * wp[i][1] << " " << _sp[i] << " " << Wp[i] << "\n";
	// }
	// outputs.close();
}
// =================================================================================================
// =================================================================================================
