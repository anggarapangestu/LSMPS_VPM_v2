#include "particle_adaption.hpp"

//==================================
//==================================
void particle_adaption::resolution_field(const Body &b, const std::vector<double> &xp, const std::vector<double> &yp,
										 std::vector<double> &field, std::vector<double> &sp)
{
	int ntotal = xp.size(); // !!! : change fp with CHI. so particle near bodies will be smaller, and outside body will be larger
	// == reset variables
	_Dtilde.clear();
	_Dp.clear();
	// == assign new size of variables
	_Dtilde.resize(ntotal);
	_Dp.resize(ntotal);

#pragma region old_code
// == internal variables
// std::vector<int> pair_i;	   // first list of pair
// std::vector<int> pair_j;	   // second list of pair
// std::vector<double> dfield_dx; // df/dx
// std::vector<double> dfield_dy; // df/dy
// std::vector<double> _Dp_rstar(ntotal);

// // == neighbor searching
// _directFindMR(ntotal, sp, xp, yp, pair_i, pair_j); // get pair_i and pair_j values

// // == gradient calculation
// std::vector<double> dfield_dx = d_dc_operator.dc_gradient(ntotal, xp, yp, sp, field, pair_i, pair_j, Pars::l_1_0_2, 1);
// std::vector<double> dfield_dy = d_dc_operator.dc_gradient(ntotal, xp, yp, sp, field, pair_i, pair_j, Pars::l_1_0_2, 2);
// // dfield_dx = d_dc_operator.dc_gradient(p, field, 1);
// // dfield_dy = d_dc_operator.dc_gradient(p, field, 2);

// // == calculating field resolution
// double _norm_gradf, _denominator;
// for (int i = 0; i < ntotal; i++)
// {
// 	_norm_gradf = std::sqrt(std::pow(dfield_dx[i], 2) + std::pow(dfield_dy[i], 2));				  // |grad(f)|
// 	_denominator = std::sqrt(1 + std::pow(_norm_gradf, 2));										  // sqrt(1 + |grad(f)^2|)
// 	_Dtilde[i] = Pars::D_0 / _denominator > Pars::D_min ? Pars::D_0 / _denominator : Pars::D_min; // assign value
// 																								  // if (_denominator > 1)
// 																								  // 	std::cout << i << " " << ntotal << std::endl;
// }
// _Dp = _Dtilde; // _Dp = _Dtilde

// double minDp = *std::min_element(_Dp.begin(), _Dp.end());
// double maxDp = *std::max_element(_Dp.begin(), _Dp.end());
// printf("\nmin = %f \nmax = %f\n", minDp, maxDp);

// // -- cutoff radius for transient neighboring
// for (int i = 0; i < ntotal; i++)
// {
// 	_Dp_rstar[i] = _Dtilde[i] * Pars::r_scale; // or divide by 2
// }

// // -- assign _Dp from Eq. (2.3)
// int n_neighbor;
// std::vector<int> pair;
// for (int i = 0; i < ntotal; i++)
// {
// 	_directFindMRinterp(xp[i], yp[i], _Dp_rstar[i], ntotal, _Dp_rstar, xp, yp, pair);
// 	n_neighbor = pair.size();
// 	for (int j = 0; j < n_neighbor; j++)
// 	{
// 		int idj = pair[j];
// 		_Dp[idj] = _Dtilde[i] < _Dp[idj] ? _Dtilde[i] : _Dp[idj];
// 	}
// }
#pragma endregion

#pragma region new_code
	double *bodycenter = new double[2];
	bodycenter[0] = (*std::max_element(b.x.begin(), b.x.end()) + *std::min_element(b.x.begin(), b.x.end())) * 0.5;
	bodycenter[1] = (*std::max_element(b.y.begin(), b.y.end()) + *std::min_element(b.y.begin(), b.y.end())) * 0.5;

	double _rSigmaIn = Pars::ly * 0.5 + 10 * Pars::sigma; // !
	double _rSigmaOut = Pars::ly * 1;					  // !
	double _sigmaMin = Pars::sigma;						  // !
	double _sigmaMax = Pars::sigma * 2;					  // !
	for (int i = 0; i < ntotal; i++)
	{
		double r2bodycenter = std::sqrt(std::pow(xp[i] - bodycenter[0], 2) + std::pow(yp[i] - bodycenter[1], 2));
		if (r2bodycenter < _rSigmaIn)
		{
			_Dp[i] = _sigmaMin;
		}
		else if (r2bodycenter >= _rSigmaIn && r2bodycenter < _rSigmaOut)
		{
			double _r = r2bodycenter - _rSigmaIn;
			double _L = _rSigmaOut - _rSigmaIn;
			_Dp[i] = (_r / _L) * (_sigmaMax - _sigmaMin) + _sigmaMin;
		}
		else
		{
			_Dp[i] = _sigmaMax;
		}
	}
	delete[] bodycenter;

	// // int nmerge = 0;
	// for (size_t i = 0; i < ntotal; i++)
	// {
	// 	// double ratio = _Dp[i] / sp[i];
	// 	// if (ratio > 1.4)
	// 	// {
	// 	// 	nmerge++;
	// 	// 	printf("\n%d/%d = %f %f <-> %f", nmerge, ntotal, _Dp[i], sp[i], ratio);
	// 	// }
	// 	// _Dp[i] = sp[i]; // ! uncomment
	// 	field[i] = _Dp[i];
	// }

#pragma endregion
	// == assigning core size (core size = cutoff radius)
	for (int i = 0; i < ntotal; i++)
	{
		sp[i] = _Dp[i] * Pars::r_scale;
	}

	// std::ofstream outputs;
	// outputs.open("output/dtilde.dat");
	// for (int i = 0; i < ntotal; i++)
	// {
	// 	outputs << xp[i] << " " << yp[i] << " " << sp[i] << " " << _Dp[i] << "\n";
	// }
	// outputs.close();
}
