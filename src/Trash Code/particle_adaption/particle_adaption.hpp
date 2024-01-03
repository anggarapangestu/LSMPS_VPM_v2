#ifndef INCLUDED_PARTICLE_ADAPTION
#define INCLUDED_PARTICLE_ADAPTION

#ifndef INCLUDED_GLOBAL
#include "../../../global.hpp"
#endif

#ifndef INCLUDED_UTILS
#include "../../../Utils.hpp"
#endif

// #ifndef INCLUDED_DC_OPERATOR
// #include "../dc_operator/dc_operator.hpp"
// #endif

#ifndef INCLUDED_DC_BASE
#include "../../DC_operator/base_dc.hpp"
#endif

// #ifndef INCLUDED_PENALIZATION_BASE
// #include "../penalization/base_penalization.hpp"
// #endif

#include <armadillo>
using namespace arma;

// class dc_operator;
class base_dc;
// class base_penalization;
class particle_adaption
{
	// == creating class instances
	// dc_operator d_dc_operator;
	base_dc d_base_dc;
	// base_penalization d_base_penalization;

#pragma region internal_variables
	std::vector<double> _xcorner;
	std::vector<double> _ycorner;

	std::vector<bool> _isFused;
	std::vector<double> _Dtilde;
	std::vector<double> _Dp;
	std::vector<double> _Dtildetemp;
	std::vector<double> _Dptemp;

	int _np, _npnew;
	int _minneighbor;
	double _dcmin;
	std::vector<double> _xp;
	std::vector<double> _yp;
	std::vector<double> _sp;
	std::vector<double> _fp;
	std::vector<double> _xptemp;
	std::vector<double> _yptemp;
	std::vector<double> _sptemp;
	std::vector<double> _fptemp;
#pragma endregion

#pragma region searching_algorithm
	void _directFind(double xpi, double ypi, double hsmli,
					 const int ntotal, const std::vector<double> &hsml,
					 const std::vector<double> &xp, const std::vector<double> &yp, std::vector<int> &pair);
	void _directFindMR(int i, const int ntotal, const std::vector<double> &hsml,
					   const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<bool> &isFused,
					   std::vector<int> &pair);
	void _directFindMRinterp(double xpi, double ypi, double hsmli,
							 const int ntotal, const std::vector<double> &hsml,
							 const std::vector<double> &xp, const std::vector<double> &yp, std::vector<int> &pair);
	void _directFindMR(int ntotal, const std::vector<double> &hsml,
					   const std::vector<double> &xp, const std::vector<double> &yp,
					   std::vector<int> &pair_i, std::vector<int> &pair_j);
#pragma endregion

	// == internal functions
	// -- for inter-particle force calculation
	double V1r(const double rpq, const double Dpq);
	double V2r(const double rpq, const double Dpq);

	// for iterative adaption methods
	void removal();
	void insertion();
	void resolutionInterpolation(const int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp);

	void steepestDescent();

	// -- for initialization
	void set_properties(const int np, const std::vector<double> &xp, const std::vector<double> &yp,
						const std::vector<double> &sp, const std::vector<double> &fp);

	// void resolution_field(const std::vector<double> &xp, const std::vector<double> &yp,
	// 					  const std::vector<double> &field, std::vector<double> &sp);
	void resolution_field(const Body &b, const std::vector<double> &xp, const std::vector<double> &yp,
						  std::vector<double> &field, std::vector<double> &sp);

public:
	// constructor
	particle_adaption();
	// method
	// void iterativeAdaption(int &np, std::vector<double> &xp, std::vector<double> &yp,
	// 					   std::vector<double> &sp, std::vector<double> &fp);
	void iterativeAdaption(const Body &body, Particle &particle);
	// destructor
	~particle_adaption();
};

#endif