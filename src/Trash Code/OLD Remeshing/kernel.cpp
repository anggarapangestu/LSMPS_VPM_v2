#include "base_remeshing.hpp"

// ** KERNELs FUNCTION

// * A popular third order M4' interpolation kernel
void base_remeshing::M4(double dex, double dey, double &Wx, double &Wy)
{
	Wx = 0.0e0;
	Wy = 0.0e0;
	if (dex > 2.0e0)
	{
		Wx = 0.0e0;
	}
	else if ((dex <= 2.0e0) && (dex > 1.0e0))
	{
		Wx = (0.5e0 * pow((2.0e0 - dex), 2)) * (1.0e0 - dex);
	}
	else if (dex <= 1.0e0)
	{
		Wx = 1.0e0 - (2.5e0 * pow(dex, 2)) + (1.5e0 * pow(dex, 3));
	}

	if (dey > 2.0e0)
	{
		Wy = 0.0e0;
	}
	else if ((dey <= 2.0e0) && (dey > 1.0e0))
	{
		Wy = (0.5e0 * pow((2.0e0 - dey), 2)) * (1.0e0 - dey);
	}
	else if (dey <= 1.0e0)
	{
		Wy = 1.0e0 - (2.5e0 * pow(dey, 2)) + (1.5e0 * pow(dey, 3));
	}
}

void base_remeshing::M04(double dex, double dey, double &Wx, double &Wy)
{
	// * M4 interpolation kernel, Monaghan 1985
	Wx = 0.0e0;
	Wy = 0.0e0;
	if (dex > 2.0e0)
	{
		Wx = 0.0e0;
	}
	else if ((dex <= 2.0e0) && (dex > 1.0e0))
	{
		Wx = 1.0e0 / 6.0e0 * pow((2.0e0 - dex), 3);
	}
	else if (dex <= 1.0e0)
	{
		Wx = 2.0e0 / 3.0e0 - pow(dex, 2) + 0.5e0 * pow(dex, 3);
	}

	if (dey > 2.0e0)
	{
		Wy = 0.0e0;
	}
	else if ((dey <= 2.0e0) && (dey > 1.0e0))
	{
		Wy = 1.0e0 / 6.0e0 * pow((2.0e0 - dey), 3);
	}
	else if (dey <= 1.0e0)
	{
		Wy = 2.0e0 / 3.0e0 - pow(dey, 2) + 0.5e0 * pow(dey, 3);
	}
}

void base_remeshing::M6(double dex, double dey, double &Wx, double &Wy)
{
	// * M*6 This kernel has a support of six points and has an error of O(h^4)
	Wx = 0.0e0;
	Wy = 0.0e0;

	if (dex > 3.0e0)
	{
		Wx = 0.0e0;
	}
	else if ((dex <= 3.0e0) && (dex > 2.0e0))
	{
		Wx = (-1.0e0 / 24.0e0) * (dex - 2.0e0) * pow((dex - 3.0e0), 3) * (5.0e0 * dex - 8.0e0);
	}
	else if ((dex <= 2.0e0) && (dex > 1.0e0))
	{
		Wx = (1.0e0 / 24.0e0) * (dex - 1.0e0) * (dex - 2.0e0) *
			 (25.0e0 * pow(dex, 3) - 114.0e0 * pow(dex, 2) + 153.0e0 * dex - 48.0e0);
	}
	else if (dex <= 1.0e0)
	{
		Wx = (-1.0e0 / 12.0e0) * (dex - 1.0e0) *
			 (25.0e0 * pow(dex, 4) - 38e0 * pow(dex, 3) - 3.0e0 * pow(dex, 2) +
			  12.0e0 * dex + 12.0e0);
	}

	if (dey > 3.0e0)
	{
		Wy = 0.0e0;
	}
	else if ((dey <= 3.0e0) && (dey > 2.0e0))
	{
		Wy = (-1.0e0 / 24.0e0) * (dey - 2.0e0) * pow((dey - 3.0e0), 3) * (5.0e0 * dey - 8.0e0);
	}
	else if ((dey <= 2.0e0) && (dey > 1.0e0))
	{
		Wy = (1.0e0 / 24.0e0) * (dey - 1.0e0) * (dey - 2.0e0) *
			 (25.0e0 * pow(dey, 3) - 114.0e0 * pow(dey, 2) + 153.0e0 * dey - 48.0e0);
	}
	else if (dey <= 1.0e0)
	{
		Wy = (-1.0e0 / 12.0e0) * (dey - 1.0e0) *
			 (25.0e0 * pow(dey, 4) - 38e0 * pow(dey, 3) - 3.0e0 * pow(dey, 2) +
			  12.0e0 * dey + 12.0e0);
	}
}

void base_remeshing::kernel(double dex, double dey, int kernel_type, double &Wx, double &Wy)
{
	// ----The Algorithms using kernel function for remeshing
	//      kernel_type = 1 : M'4 A popular third order M4' interpolation kernel
	//                    2 : M*6 This kernel has a support of six points and has an error of O(h^4)

	//      dex--   [input]
	//      dey--   [input]
	//      dez--   [input]
	//      Wx--    [output]
	//      Wy--    [output]
	//      Wz--    [output]
	double Wx4, Wy4, Wx04, Wy04, M4iso;
	Wx = 0;
	Wy = 0;
	// Choosing The kernel function type for Remeshing
	if (kernel_type == 1)
	{
		M4(dex, dey, Wx, Wy);
	}
	if (kernel_type == 2)
	{
		M6(dex, dey, Wx, Wy);
	}
	if (kernel_type == 3)
	{
		M04(dex, dey, Wx04, Wy04);
		M4(dex, dey, Wx4, Wy4);

		M4iso = Wx4 * Wy4 - (Wx4 - Wx04) * (Wy4 - Wy04);
		Wx = M4iso;
		Wy = 1.0e0;
	}
}
