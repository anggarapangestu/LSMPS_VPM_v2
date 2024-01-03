#include "particle_adaption.hpp"

// ===================
// ===================
void particle_adaption::resolutionInterpolation(const int np,
											   const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp)
{
	// =============================================================
	// input:
	// np = number of particles of old partocles distribution, particles that are going to be adapted
	// xp = x coordinate of particles of old partocles distribution
	// yp = y coordinate of particles of old partocles distribution
	// sp = core size  of particles of old partocles distribution
	// ==============================================================
	// == internal variables
	int niter;
	double deltaR;
	double dx, dy, Wx, Wy, finterp(0), dr(0), nomDp(0), nomDtilde(0), denom(0) /*, sptemp*/;
	std::vector<int> pair; // dynamic neihgbor index pairs
	int n_neighbor, idx;
	double xi, xj, yi, yj, deltax, deltay;
	printf("adapting particle's core size ... ");
	printf("%d -> %d\n", _np, _npnew);

	// == reset previous variables
	_Dptemp.clear();
	_Dtildetemp.clear();
	_Dptemp.resize(_npnew);
	_Dtildetemp.resize(_npnew);
	// clock_t t;
	// t = clock();

	// --> arrange matrix b(xp)
	int l = Pars::l_0_0_3;
	arma::mat b(l, 1, arma::fill::zeros);
	b(0, 0) = 1;

	for (int i = 0; i < _npnew; i++)
	{
		_directFindMRinterp(_xp[i], _yp[i], _sp[i], np, sp, xp, yp, pair);

		deltaR = 1;
		while (pair.size() <= l)
		{
			_directFind(_xp[i], _yp[i], _sp[i] * deltaR, np, sp, xp, yp, pair);
			deltaR += 0.5;
			// printf("-- sp = %f, n = %d\n", _sp[i], pair.size());
		}

		n_neighbor = pair.size();
		xi = this->_xp[i];
		yi = this->_yp[i];
		// --> arrange matrix E(xp) and V(xp)
		arma::sp_mat E(n_neighbor, n_neighbor);		// or use sParse matrix SpMat<double>
		arma::sp_mat phi_e(n_neighbor, n_neighbor); // or use sParse matrix SpMat<double>
		arma::mat V(n_neighbor, l);
		arma::mat dDp(1, n_neighbor);
		arma::mat dDtilde(1, n_neighbor);

		for (int j = 0; j < pair.size(); j++)
		{
			idx = pair[j];
			dx = std::sqrt(std::pow(this->_xp[i] - xp[idx], 2) + std::pow(this->_yp[i] - yp[idx], 2));
			xj = xp[idx];
			yj = yp[idx];
			deltax = (xi - xj) / this->_sp[i];
			deltay = (yi - yj) / this->_sp[i];
			dDp(j) = _Dp[idx];
			dDtilde(j) = _Dtilde[idx];
			E(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 0.5);	 // sqrt of Gaussian kernel
			phi_e(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 1.0); // Gaussian kernel

			for (int k = 0; k < l; k++)
			{
				V(j, k) = d_base_dc.Vandermonde_zero(deltax, deltay, this->_sp[i], k);
			}
		}

		// --> arrange matrix B(xp)
		mat B = E.t() * V;
		// --> arrange matrix A(xp)
		mat A = B.t() * B;
		// --> calculate a(xp)
		arma::mat a = solve(A, b);

		arma::mat K_e = V * a;					// (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
		arma::mat eta_e = phi_e * K_e;			// (nxn)(nx1) = (nx1), DC kernel
		arma::mat sumQDp = dDp * eta_e;			// (1xn)(nx1) = (1x1)
		arma::mat sumQDtilde = dDtilde * eta_e; // (1xn)(nx1) = (1x1)

		_Dptemp[i] = sumQDp(0);
		_Dtildetemp[i] = sumQDtilde(0);
		this->_sp[i] = _Dptemp[i] * Pars::r_scale;

		// if (_sp[i] < 0)
		// {
		// 	printf("sp[%d] = %f\n", i, _sp[i]);
		// }

		// // use dc-remeshing!!!!!!!!!!!!!!!
		// if (pair.size() > Pars::l_0_0_2)
		// {
		// 	for (int j = 0; j < pair.size(); j++)
		// 	{
		// 		dx = std::abs((this->_xp[i] - xp[pair[j]])/_sp[i]);
		// 		dy = std::abs((this->_yp[i] - yp[pair[j]])/_sp[i]);
		// 		// M4(dx,Wx); M4(dy,Wy);
		// 		// finterp += sp[pair[j]]*Wx*Wy;
		// 		// -- performing linear interpolation using weighting approach
		// 		dr = std::sqrt(std::pow(dx,2)+std::pow(dy,2));
		// 		// nom   += sp[pair[j]]*dr;
		// 		nomDp     += _Dp[pair[j]]*dr;
		// 		nomDtilde += _Dtilde[pair[j]]*dr;
		// 		denom += dr;
		// 	}
		// 	// _sp[i] = finterp;
		// 	// finterp = 0.0e0;
		// 	this->_Dptemp[i]     = nomDp/denom;
		// 	this->_Dtildetemp[i] = nomDtilde/denom;
		// 	// printf("%f %f\n",_Dptemp[i], _Dtildetemp[i] );
		// 	// if (_Dtildetemp[i] > 2*_Dptemp[i])
		// 	//  {
		// 	//  	printf("Dt > 2Dp, %d\n", i);
		// 	//  }
		// 	// this->_sp[i] = nom/denom;
		// 	this->_sp[i] = _Dptemp[i]*Pars::rstar;
		// 	// nom = 0.0e0;
		// 	nomDp = 0.0e0; nomDtilde = 0.0e0;
		// 	denom = 0.0e0;
		// }
		// pair.clear(); // clearing pair into null vector!
	}
	// t = clock() - t;
	// printf("%f\n", (float)t/CLOCKS_PER_SEC);

	// // -- saving temporary data
	// std::ofstream outputs;
	// outputs.open("output/interp.dat");
	// for (int i = 0; i < _npnew; i++)
	// 	outputs << _xp[i] << " " << _yp[i] << " " << _sp[i] << "\n";
	// outputs.close();
}
// ===================
// ===================
// void particle_adaption::M4(const double dr, double &W)
// {
// 	if (dr<=1)
// 		W = 1 - 2.5*std::pow(dr,2) + 1.5*std::pow(dr,3);
// 	else if(dr>1 && dr<=2)
// 		W = 0.5*std::pow((2-dr),2)*(1-dr);
// 	else
// 		W = 0.0e0;
// }
// ===================
// ===================
