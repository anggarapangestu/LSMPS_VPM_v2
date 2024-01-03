#ifndef INCLUDED_DC_BASE
#define INCLUDED_DC_BASE

#ifndef INCLUDED_GLOBAL
#include "../../global.hpp"
#endif

class base_dc
{
public:
	// TODO: generate Vandermonde matrix component
	double Vandermonde(const double &dx, const double &dy, const double &eps, const int &j);
	double Vandermonde_zero(const double &dx, const double &dy, const double &eps, const int &j);

	// TODO: to calculate Laplacian operator
	std::vector<double> Q22(const int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp,
							const std::vector<double> &fp, const std::vector<std::vector<int>> &neighborList, const int l);
	// TODO: to calculate gradient(GRAD) operator
	std::vector<std::vector<double>> Q11(const int np, const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<double> &sp,
										 const std::vector<double> &fp, const std::vector<std::vector<int>> &neighborList, const int l);
};

#endif
