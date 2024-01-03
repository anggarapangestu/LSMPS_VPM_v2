#include "base_dc.hpp"

double base_dc::Vandermonde(const double &dx, const double &dy, const double &eps, const int &j)
{
	double pp;

	if (j == 0)
		pp = dx;
	else if (j == 1)
		pp = dy;
	else if (j == 2)
		pp = dx * dx;
	else if (j == 3)
		pp = dy * dy;
	else if (j == 4)
		pp = dx * dy;
	else if (j == 5)
		pp = (dx * dx) * dx;
	else if (j == 6)
		pp = (dy * dy) * dy;
	else if (j == 7)
		pp = (dx * dx) * dy;
	else if (j == 8)
		pp = dx * (dy * dy);
	else if (j == 9)
		pp = dy * dy * dy * dy;
	else if (j == 10)
		pp = dx * dy * dy * dy;
	else if (j == 11)
		pp = dx * dx * dy * dy;
	else if (j == 12)
		pp = dx * dx * dx * dy;
	else if (j == 13)
		pp = dx * dx * dx * dx;
	else if (j == 14)
		pp = dy * dy * dy * dy * dy;
	else if (j == 15)
		pp = dx * dy * dy * dy * dy;
	else if (j == 16)
		pp = dx * dx * dy * dy * dy;
	else if (j == 17)
		pp = dx * dx * dx * dy * dy;
	else if (j == 18)
		pp = dx * dx * dx * dx * dy;
	else if (j == 19)
		pp = dx * dx * dx * dx * dx;

	return pp;
}

double base_dc::Vandermonde_zero(const double &dx, const double &dy, const double &eps, const int &j)
{
	double pp;

	if (j == 0)
		pp = 1;
	else if (j == 1)
		pp = dx;
	else if (j == 2)
		pp = dy;
	else if (j == 3)
		pp = dx * dy;
	else if (j == 4)
		pp = dx * dx;
	else if (j == 5)
		pp = dy * dy;
	else if (j == 6)
		pp = dx * dx * dx;
	else if (j == 7)
		pp = dx * dx * dy;
	else if (j == 8)
		pp = dx * dy * dy;
	else if (j == 9)
		pp = dy * dy * dy;
	else if (j == 10)
		pp = dy * dy * dy * dy;
	else if (j == 11)
		pp = dx * dy * dy * dy;
	else if (j == 12)
		pp = dx * dx * dy * dy;
	else if (j == 13)
		pp = dx * dx * dx * dy;
	else if (j == 14)
		pp = dx * dx * dx * dx;

	return pp;
}
