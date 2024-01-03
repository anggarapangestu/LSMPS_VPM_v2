#include "particle_adaption.hpp"

// =================================================================================================
// =================================================================================================
void particle_adaption::iterativeAdaption(const Body &body, Particle &particle)
{
	// -- internal variables
	int &np = particle.num;				   // number of particles
	std::vector<double> &xp = particle.x;  // x position
	std::vector<double> &yp = particle.y;  // y position
	std::vector<double> &sp = particle.s;  // core size
	std::vector<double> &fp = particle.gz; // vorticity

	// -- initialization
	std::vector<double> _resfield(np, 0.0e0); // chi near body
	// d_base_penalization.kai_def(body, particle, _lambda, _resfield); // generate chi's value
	this->resolution_field(body, xp, yp, _resfield, sp); // !!! : change fp with CHI. so particle near bodies will be smaller, and outside body will be larger
	this->set_properties(np, xp, yp, sp, fp);			 // copy variables to local variables

	// -- iterative adaption
	// do
	// {
	this->removal();   // remove near particle
	this->insertion(); // add particle
	// this->resolutionInterpolation(np, xp, yp, sp); // interpolate field from old to new particle
	// 	this->steepestDescent();					   // move particle against gradient
	// } while (this->_dcmin < Pars::dcmin ||			   // conditional: ensure no particle is too close
	// 		 this->_minneighbor < Pars::l_2_0_2);	  // conditional: ensure enough neighbor

	// this->removal(); // remove once again

	_fptemp.clear();			   // clear internal variables
	_fptemp.resize(_npnew, 0.0e0); // interpolated intensities

#pragma region interpolate_intensities
	// double dx, dy, Wx, Wy;
	// double xi, yi, xj, yj;
	// std::vector<int> pair;
	// int l = Pars::l_0_0_3; // number moment condition
	// // == internal dynamic variables
	// int n_neighbor, idx; // number of neighbor for each particle
	// double deltax, deltay, R, deltaR;

	// // --> arrange matrix b(xp)
	// printf("Interpolating particles' intensities ...");
	// arma::mat b(l, 1, arma::fill::zeros);
	// b(0, 0) = 1;

	// for (int i = 0; i < _npnew; i++)
	// {
	// 	// _NeighborSearchAlgorithm.DirectFindMRinterp(_xp[i], _yp[i], _sp[i], np, sp, xp, yp, pair);
	// 	// _NeighborSearchAlgorithm.DirectFind(_xp[i], _yp[i], _sp[i], np, sp, xp, yp, pair);
	// 	_directFindMRinterp(_xp[i], _yp[i], _sp[i], np, sp, xp, yp, pair);

	// 	deltaR = 1; // reset cutoff radius multiplier to 1. this variable is used to ensure total neighbor particles
	// 	while (pair.size() < l)
	// 	{
	// 		// _NeighborSearchAlgorithm.DirectFind(_xp[i], _yp[i], _sp[i] * deltaR, np, sp, xp, yp, pair);
	// 		_directFind(_xp[i], _yp[i], _sp[i] * deltaR, np, sp, xp, yp, pair);
	// 		deltaR += 0.5;
	// 	}

	// 	n_neighbor = pair.size();

	// 	xi = this->_xp[i]; // x position of interpolated particle
	// 	yi = this->_yp[i]; // y position of interpolated particle
	// 	// --> arrange matrix E(xp) and V(xp)
	// 	arma::sp_mat E(n_neighbor, n_neighbor);		// matrix of sqrt gaussian kernel
	// 	arma::sp_mat phi_e(n_neighbor, n_neighbor); // matrix of gaussian kernel
	// 	arma::mat V(n_neighbor, l);					// vandermonde matrix
	// 	arma::mat dField(1, n_neighbor);			// difference of field

	// 	for (int j = 0; j < n_neighbor; j++)
	// 	{
	// 		idx = pair[j];
	// 		dx = std::sqrt(std::pow(this->_xp[i] - xp[idx], 2) + std::pow(this->_yp[i] - yp[idx], 2));
	// 		xj = xp[idx];
	// 		yj = yp[idx];
	// 		deltax = (xi - xj) / this->_sp[i];
	// 		deltay = (yi - yj) / this->_sp[i];
	// 		dField(j) = fp[idx];												   // assign nearby particle intensities
	// 		E(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 0.5 * 1.e-0);	 // sqrt of Gaussian kernel. c^2 = 0.04
	// 		phi_e(j, j) = std::exp(-std::pow(dx / this->_sp[i], 2) * 1.0 * 1.e-0); // Gaussian kernel. c^2 = 0.04

	// 		for (int k = 0; k < l; k++)
	// 		{
	// 			V(j, k) = d_base_dc.Vandermonde_zero(deltax, deltay, this->_sp[i], k);
	// 		}
	// 	}

	// 	// --> arrange matrix B(xp)
	// 	mat B = E.t() * V;
	// 	// --> arrange matrix A(xp)
	// 	mat A = B.t() * B;
	// 	// --> calculate a(xp)
	// 	arma::mat a = solve(A, b);

	// 	arma::mat K_e = V * a;			 // (nxl)(lx1) = (nx1), monomial vector multiplied by monomial basis
	// 	arma::mat eta_e = phi_e * K_e;   // (nxn)(nx1) = (nx1), DC kernel
	// 	arma::mat sumQ = dField * eta_e; // (1xn)(nx1) = (1x1)

	// 	_fptemp[i] = sumQ(0);
	// }
	// printf(" done\n");
#pragma endregion

	// -- resize particle core size
	for (size_t i = 0; i < _sp.size(); i++)
	{
		_sp[i] /= Pars::r_scale;
	}
	// for (auto &i : _sp)
	// {
	// 	_sp[i] /= Pars::r_scale;
	// }

	// == update particle's properties
	// xp = _xp;	 // assign new particle x position
	// yp = _yp;	 // assign new particle y position
	// sp = _sp;	 // assign new particle core size
	// fp = _fptemp; // assign new particle field intensities
	// np = _npnew;  // assign new particle number
	particle.x = this->_xp;
    particle.y = this->_yp;
    particle.s = this->_sp;
    particle.num = this->_npnew;
    particle.u.resize(particle.num, 0.0e0);
    particle.v.resize(particle.num, 0.0e0);
    particle.gz.resize(particle.num, 0.0e0);
    particle.neighbor.resize(particle.num, std::vector<int>());
    particle.isActive.resize(particle.num, 0.0e0);

	// // == saving temporary data
	// std::ofstream outputs;
	// outputs.open("output/xpnew.dat");
	// for (int i = 0; i < _npnew; i++)
	// {
	// 	outputs << _xp[i] << " " << _yp[i] << " " << _sp[i] << " " << _fptemp[i] << "\n";
	// }
	// outputs.close();
}
// =================================================================================================
// =================================================================================================