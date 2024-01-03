#ifndef INCLUDED_BIOT_SAVART
#define INCLUDED_BIOT_SAVART

#ifndef INCLUDED_UTILS
#include "../../Utils.hpp"
#endif

// #ifndef PROGRESSBAR_PROGRESSBAR_HPP
// #include "../display/ProgressBar.hpp" // displaying percentage bar
// #endif
#include <sstream> 

class velocity_biot_savart
{
public:
	// regularization function 2d
	void gaussian_function(double rij, double sij, double &gaussij);
	void regul_func_2d		(double rij, double sij, double &q);

	// -- Biot-Savart direct
	void biotsavart_direct_2d(Particle &pi, Particle &pj);
	// -- Biot-Savart FMM
	void biotsavart_fmm_2d(Particle &pi, Particle &pj, const int icutoff, const int n_s, const int n_inter, const int ndp);

private:
	// -- subroutine all fmm 2d
	void par_loc(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi, 
		double xb_mini, double xb_maxi, double yb_mini, double yb_maxi, int &nip, std::vector<int> &ipp);

	void amount_inbox(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi, 
										int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, 
										double xb_mini, double xb_maxi, double yb_mini,double yb_maxi, int &npr_in, int n_inter);

	void hierarchy_mesh(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi, 
											int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, 
		int n_s, int n_inter, double xmin, double ymin, double xmax, double ymax, std::vector<int> &nb, 
		std::vector<std::vector<double> > &xb_min, std::vector<std::vector<double> > &yb_min,
		std::vector<std::vector<double> > &xb_max, std::vector<std::vector<double> > &yb_max,
		std::vector<std::vector<double> > &xb_cen, std::vector<std::vector<double> > &yb_cen,
		int &lev, std::vector<std::vector<int> > &nchild, std::vector<std::vector<std::vector<int> > > &ichild, 
		std::vector<std::vector<int> > &iparent, std::vector<std::vector<std::vector<int>>> &particleinbox );

	void list_one(int ib, int k, const std::vector<int> &nb, int lev, 
		const std::vector<std::vector<double> > &xb_min, const std::vector<std::vector<double> > &yb_min,
		const std::vector<std::vector<double> > &xb_max, const std::vector<std::vector<double> > &yb_max,
		const std::vector<std::vector<int> > &nchild, 
		int &nls1, std::vector<int> &ils1, std::vector<int> &kls1 );

	void list_two(int ib, int k, const std::vector<int> &nb,
		const std::vector<std::vector<double> > &xb_min, const std::vector<std::vector<double> > &yb_min,
		const std::vector<std::vector<double> > &xb_max, const std::vector<std::vector<double> > &yb_max,
		const std::vector<std::vector<int> > &iparent, const std::vector<std::vector<int> > &nchild, 
		const std::vector<std::vector<std::vector<int> > > &ichild, 
		int &nls2, std::vector<int> &ils2, std::vector<int> &kls2 );

	void list_three(int ib, int k, const std::vector<int> &nb, int lev, 
		const std::vector<std::vector<double> > &xb_min, const std::vector<std::vector<double> > &yb_min,
		const std::vector<std::vector<double> > &xb_max, const std::vector<std::vector<double> > &yb_max,
		const std::vector<std::vector<int> > &nchild, const std::vector<std::vector<std::vector<int> > > &ichild,
		int &nls3, std::vector<int> &ils3, std::vector<int> &kls3 );

	void list_four(int ib, int k, const std::vector<int> &nb, 
		const std::vector<std::vector<double> > &xb_min, const std::vector<std::vector<double> > &yb_min,
		const std::vector<std::vector<double> > &xb_max, const std::vector<std::vector<double> > &yb_max,
		const std::vector<std::vector<int> > &iparent, const std::vector<std::vector<int> > &nchild, 
		int &nls4, std::vector<int> &ils4, std::vector<int> &kls4 );

	void bico(int n, int &k, double &c);

	void direct_sum(int nip, int npi, const std::vector<int> &ipp, const std::vector<double> &xi, 
		const std::vector<double> &yi, const std::vector<double> &si, std::vector<double> &ui, std::vector<double> &vi, 
		int njp, int npj, const std::vector<int> &jpp, const std::vector<double> &xj, const std::vector<double> &yj, 
		const std::vector<double> &sj, const std::vector<double> &gj, int icutoff);

	void amount_inbox_new(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi,
						int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, double xb_mini, double xb_maxi,
						double yb_mini, double yb_maxi, int &npr_in, int n_inter, std::vector<std::vector<std::vector<int>>> &numberofparticle, int parentx, int parenty, int x, int y);

};

#endif
