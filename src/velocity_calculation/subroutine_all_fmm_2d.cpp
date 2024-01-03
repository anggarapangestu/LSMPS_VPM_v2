// !!!!============================================================
// !!!!==== subroutine for biotsavart-fmm 2D ======================
// !!!!============================================================
#include "velocity_biot_savart.hpp"
// #include "poisson.hpp"
// !!!!============================================================
// !!!!============================================================

void velocity_biot_savart::par_loc(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi,
						   double xb_mini, double xb_maxi, double yb_mini, double yb_maxi, int &nip, std::vector<int> &ipp)
{

	
	nip = 0;
	for (int i = n0; i <= n1; i++)
	{
		if ((xi[i] >= xb_mini) && (xi[i] < xb_maxi) &&
			(yi[i] >= yb_mini) && (yi[i] < yb_maxi))
		{
			nip = nip + 1;
			ipp[nip] = i;
		}
	}

	if (nip == 0)
	{
		ipp[1] = 0;
	}
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::amount_inbox(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi,
								int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, double xb_mini, double xb_maxi,
								double yb_mini, double yb_maxi, int &npr_in, int n_inter)
{
	// internal variables
	int isame;

	if (n_inter == 1)
	{
		npr_in = 0;
		//#pragma omp parallel for reduction(+:npr_in)
		for (int i = n0; i <= n1; i++)
		{
			if ((xi[i] >= xb_mini) && (xi[i] < xb_maxi) &&
				(yi[i] >= yb_mini) && (yi[i] < yb_maxi))
			{
				npr_in = npr_in + 1;
			}
		}
	}
	else if (n_inter == 2)
	{
		npr_in = 0;
		for (int j = n2; j <= n3; j++)
		{
			if ((xj[j] >= xb_mini) && (xj[j] < xb_maxi) &&
				(yj[j] >= yb_mini) && (yj[j] < yb_maxi))
			{
				npr_in = npr_in + 1;
			}
		}
	}
	else if (n_inter == 3)
	{
		npr_in = 0;
		for (int i = n0; i <= n1; i++)
		{
			if ((xi[i] >= xb_mini) && (xi[i] < xb_maxi) &&
				(yi[i] >= yb_mini) && (yi[i] < yb_maxi))
			{
				npr_in = npr_in + 1;
			}
		}
		for (int j = n2; j <= n3; j++)
		{
			for (int i = n1; i <= n1; i++)
			{
				isame = 0;
				if ((xj[j] == xi[i]) && (yj[j] == yi[i]))
				{
					isame = 1;
				}
			}
			if ((xj[j] >= xb_mini) && (xj[j] < xb_maxi) &&
				(yj[j] >= yb_mini) && (yj[j] < yb_maxi) && (isame == 0))
			{
				npr_in = npr_in + 1;
			}
		}
	}
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::hierarchy_mesh(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi,
								  int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, int n_s, int n_inter,
								  double xmin, double ymin, double xmax, double ymax, std::vector<int> &nb,
								  std::vector<std::vector<double>> &xb_min, std::vector<std::vector<double>> &yb_min,
								  std::vector<std::vector<double>> &xb_max, std::vector<std::vector<double>> &yb_max,
								  std::vector<std::vector<double>> &xb_cen, std::vector<std::vector<double>> &yb_cen,
								  int &lev, std::vector<std::vector<int>> &nchild, std::vector<std::vector<std::vector<int>>> &ichild,
								  std::vector<std::vector<int>> &iparent, std::vector<std::vector<std::vector<int>>> &particleinbox)
{

	//Create a quadtree for the basis of FMM
	
	particleinbox.resize(Pars::lmax + 1, std::vector<std::vector<int>>(Pars::nbmrl+1));
	int **nprin = new int *[Pars::nbmrl + 1];
	
	for (size_t i = 0; i < Pars::nbmrl + 1; i++)
	{
		nprin[i] = new int[Pars::lmax + 1];
	} 

	double cent, db, ratio, tlout, dr, rdx, rdy, xb_mint, xb_maxt, yb_mint, yb_maxt, dbx, dby;
	int nb1, nb2, nbt, ib, i_stop, k1, n_stop, npr_in;

	xmin = 1.0e0 / Pars::tol2;  //initiate very big value fo finding minimum value
	xmax = -1.0e0 / Pars::tol2; //initiate very small value for finding maximum value
	ymin = 1.0e0 / Pars::tol2;
	ymax = -1.0e0 / Pars::tol2;

	for (int i = n0; i <= n1; i++)
	{
		xmin = std::min(xmin, xi[i]);
		xmax = std::max(xmax, xi[i]);
		ymin = std::min(ymin, yi[i]);
		ymax = std::max(ymax, yi[i]);
	}

	//bisa dicut gak sih yang ini karena xi == xj 
	/*for (int j = n2; j <= n3; j++)
	{
		xmin = std::min(xmin, xj[j]);
		xmax = std::max(xmax, xj[j]);
		ymin = std::min(ymin, yj[j]);
		ymax = std::max(ymax, yj[j]);
	}*/
	//

	tlout = 1.0e-5;
	dr = std::min(std::abs(xmax - xmin), std::abs(ymax - ymin)) / ((double)(pow(2, Pars::lmax)));
	xmin = xmin - dr * tlout;
	xmax = xmax + dr * tlout;
	ymin = ymin - dr * tlout;
	ymax = ymax + dr * tlout;

	//determine which side is more longer than the other, after that make that into square domain
	rdx = std::abs(xmax - xmin); // x-side
	rdy = std::abs(ymax - ymin); // y-side
	if (rdy >= rdx)
	{
		cent = (xmin + xmax) / 2.0e0;
		xmin = cent - rdy / 2.0e0;
		xmax = cent + rdy / 2.0e0;
		db = rdy / 2.0e0;
		nb1 = 2;
		nb2 = 2;
	}
	else if (rdx > rdy)
	{
		cent = (ymin + ymax) / 2.0e0;
		ymin = cent - rdx / 2.0e0;
		ymax = cent + rdx / 2.0e0;
		db = rdx / 2.0e0;
		nb1 = 2;
		nb2 = 2;
	}

	//build first 4 quadrant of the tree
	ib = 0; //at start, there is no box

	for (int ib1 = 1; ib1 <= nb1; ib1++)
	{
		for (int ib2 = 1; ib2 <= nb2; ib2++)
		{
			xb_mint = xmin + (ib1 - 1) * db;
			xb_maxt = xmin + ib1 * db;
			yb_mint = ymin + (ib2 - 1) * db;
			yb_maxt = ymin + ib2 * db;
			
			//calculate total number of particles inside the box
			//amount_inbox(n0, n1, npi, xi, yi, n2, n3, npj, xj, yj, xb_mint, xb_maxt, yb_mint, yb_maxt, npr_in, n_inter);
			amount_inbox_new(n0, n1, npi, xi, yi, n2, n3, npj, xj, yj, xb_mint, xb_maxt, yb_mint, yb_maxt, npr_in, n_inter, particleinbox, -1, -1, 1, ib+1); 

			//If there are particles present in the box
			if (npr_in > 0)
			{
				//printf("%d\n", npr_in);
				ib = ib + 1; 								 //increase the number of box
				nb[1] = ib; 								 //update the total number of box at layer 1.. note: layer starts from 0 
				//nb[1] = ib_new[ib1][ib2];
				nprin[ib][1] = npr_in; 						 //update the total particles present in the box at layer 1
				iparent[ib][1] = 0; 						 //link the box to the parents: layer 0
				xb_min[ib][1] = xb_mint;					 //save the box minimum x coordinate
				xb_max[ib][1] = xb_maxt;            		 //save the box maximum x coordinate
				yb_min[ib][1] = yb_mint;            		 //save the box minimum y coordinate
				yb_max[ib][1] = yb_maxt;					 //save the box maximum y coordinate
				xb_cen[ib][1] = (xb_maxt + xb_mint) / 2.0e0; //center of box x-direction
				yb_cen[ib][1] = (yb_maxt + yb_mint) / 2.0e0; //center of box y-direction
				nchild[ib][1] = 0;  						 //number of child, at this moment still 0
			}
		}
	}

	//above loop provides us with 4 quadrant containing particles. Every quadrant were linked to their parent the overall domain. 
	
	//initiatize values for nb
	for (int i = 0; i < nb.size(); i++)
	{
		nb[i] = ib;
	}

	i_stop = 0;
	//start creating tree from layer 1
	for (int k = 1; k <= Pars::lmax - 1; k++)
	{
		k1 = k + 1;		//tree layer; starts at 2 until n-1 
		nb[k1] = 0;		//initialize the layer number of box
		n_stop = 0;     
		for (int ib = 1; ib <= nb[k]; ib++)  //check all the boxes of previous layer
		{
			nchild[ib][k] = 0; //initialize children for boxes from previous layer
			if (nprin[ib][k] > n_s) //check the total particles present in the box
			{
				n_stop = n_stop + 1; 							 
				dbx = (xb_max[ib][k] - xb_min[ib][k]) / 2.0e0;	//divide the selected box equally in x direction
				dby = (yb_max[ib][k] - yb_min[ib][k]) / 2.0e0;  //divide the selected box equally in y direction
				
				//create 4 boxes inside the selected boxes
				for (int ib1 = 1; ib1 <= 2; ib1++)
				{
					for (int ib2 = 1; ib2 <= 2; ib2++)
					{
						xb_mint = xb_min[ib][k] + (ib1 - 1) * dbx;
						xb_maxt = xb_min[ib][k] + ib1 * dbx;
						yb_mint = yb_min[ib][k] + (ib2 - 1) * dby;
						yb_maxt = yb_min[ib][k] + ib2 * dby;

						//calculate the particles present inside the box
						//amount_inbox(n0, n1, npi, xi, yi, n2, n3, npj, xj, yj, xb_mint, xb_maxt,yb_mint, yb_maxt, npr_in, n_inter);
						amount_inbox_new(n0, n1, npi, xi, yi, n2, n3, npj, xj, yj, xb_mint, xb_maxt,
									 yb_mint, yb_maxt, npr_in, n_inter, particleinbox, k, ib, k1, nb[k1]+1 );
						
						//If there are particles present in the newly created box (new layer)
						if (npr_in > 0)
						{
							nb[k1] = nb[k1] + 1;			          			//update the number of boxes in the new layer
							nchild[ib][k] = nchild[ib][k] + 1;					//update the number of child in the previous layer for the selected box
							ichild[ib][k][nchild[ib][k]] = nb[k1];				//update the child for the selected box in the previous layer
							nprin[nb[k1]][k1] = npr_in;							//update the number of particles in the new box
							xb_min[nb[k1]][k1] = xb_mint;						//save the box minimum x coordinate
							xb_max[nb[k1]][k1] = xb_maxt;						//save the box maximum x coordinate
							yb_min[nb[k1]][k1] = yb_mint;						//save the box minimum y coordinate
							yb_max[nb[k1]][k1] = yb_maxt;						//save the box maximum y coordinate
							xb_cen[nb[k1]][k1] = (xb_maxt + xb_mint) / 2.0e0;   //center of box x-direction
							yb_cen[nb[k1]][k1] = (yb_maxt + yb_mint) / 2.0e0;   //center of box y-direction
							iparent[nb[k1]][k1] = ib;							//set the parent for the new generated box
						}
					}
				}
			}
			
		}
		lev = k1;
		if (n_stop == 0)
		{
			lev = k;
			goto label_1; // exit in fortran90
		}
	}

label_1:
	for (int i = 1; i <= nb[lev]; i++)
	{
		nchild[i][lev] = 0;
	} // nchild(span(1,nb(lev)),lev) = 0;

	// -- delete internal variables
	for (size_t i = 0; i < Pars::nbmrl + 1; i++)
	{
		delete[] nprin[i];
	}
	delete[] nprin;
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::list_one(int ib, int k, const std::vector<int> &nb, int lev,
							const std::vector<std::vector<double>> &xb_min, const std::vector<std::vector<double>> &yb_min,
							const std::vector<std::vector<double>> &xb_max, const std::vector<std::vector<double>> &yb_max,
							const std::vector<std::vector<int>> &nchild,
							int &nls1, std::vector<int> &ils1, std::vector<int> &kls1)
{
	// xb_min.resize(nbmrl+1,lmax+1);	xb_max.resize(nbmrl+1,lmax+1);
	// yb_min.resize(nbmrl+1,lmax+1);	yb_max.resize(nbmrl+1,lmax+1);
	// nb.resize(lmax+1); 		nchild.resize(nbmrl+1,lmax+1);
	// ils1.resize(nbl1+1); 	kls1.resize(nbl1+1);
	int istop;
	double dx_xn, dx_nx, dy_xn, dy_nx;

	if (nchild[ib][k] == 0)
	{
		nls1 = 1;
		ils1[nls1] = ib;
		kls1[nls1] = k;
		for (int k2 = 1; k2 <= lev; k2++)
		{
			for (int ib2 = 1; ib2 <= nb[k2]; ib2++)
			{
				istop = 0;
				if ((k2 == k) && (ib2 == ib))
				{
					istop = 1;
				}
				if ((nchild[ib2][k2] == 0) && (istop == 0))
				{
					dx_xn = std::abs(xb_max[ib2][k2] - xb_min[ib][k]);
					dx_nx = std::abs(xb_min[ib2][k2] - xb_max[ib][k]);
					dy_xn = std::abs(yb_max[ib2][k2] - yb_min[ib][k]);
					dy_nx = std::abs(yb_min[ib2][k2] - yb_max[ib][k]);

					if ((((dy_xn <= Pars::tol2) || (dy_nx <= Pars::tol2)) &&
						 (xb_min[ib2][k2] <= xb_max[ib][k] + Pars::tol2) &&
						 (xb_max[ib2][k2] >= xb_min[ib][k] - Pars::tol2)) ||
						(((dx_xn <= Pars::tol2) || (dx_nx <= Pars::tol2)) &&
						 (yb_min[ib2][k2] <= yb_max[ib][k] + Pars::tol2) &&
						 (yb_max[ib2][k2] >= yb_min[ib][k] - Pars::tol2)))
					{
						nls1 = nls1 + 1;
						ils1[nls1] = ib2;
						kls1[nls1] = k2;
					}
				}
			}
		}
	}
	else if (nchild[ib][k] > 0)
	{
		nls1 = 0;
	}
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::list_two(int ib, int k, const std::vector<int> &nb,
							const std::vector<std::vector<double>> &xb_min, const std::vector<std::vector<double>> &yb_min,
							const std::vector<std::vector<double>> &xb_max, const std::vector<std::vector<double>> &yb_max,
							const std::vector<std::vector<int>> &iparent, const std::vector<std::vector<int>> &nchild,
							const std::vector<std::vector<std::vector<int>>> &ichild,
							int &nls2, std::vector<int> &ils2, std::vector<int> &kls2)
{
	// xb_min.resize(nbmrl+1,lmax+1); xb_max.resize(nbmrl+1,lmax+1);
	// yb_min.resize(nbmrl+1,lmax+1); yb_max.resize(nbmrl+1,lmax+1);
	// nb.resize(lmax+1);
	// iparent.resize(nbmrl+1,lmax+1); nchild.resize(nbmrl+1,lmax+1);
	// ichild.resize(nbmrl+1,lmax+1,4+1);
	// ils2.resize(nbl2+1); kls2.resize(nbl2+1);

	int ipt, kpt, ic;
	double dx_xn, dx_nx, dy_xn, dy_nx;

	if (k == 1) // if the highest level (4 boxes) (THE BOX at LEVEL 1)
	{
		nls2 = 0;
		for (int ib1 = 1; ib1 <= nb[1]; ib1++) //check all boxes
		{
			if (ib1 != ib)
			{
				if ((xb_min[ib1][1] > xb_max[ib][1] + Pars::tol2) ||
					(xb_max[ib1][1] < xb_min[ib][1] - Pars::tol2) ||
					(yb_min[ib1][1] > yb_max[ib][1] + Pars::tol2) ||
					(yb_max[ib1][1] < yb_min[ib][1] - Pars::tol2))
				{
					nls2 = nls2 + 1;
					ils2[nls2] = ib1;
					kls2[nls2] = 1;
				}
			}
		}
	}
	else if (k > 1) //if not the highest level
	{
		nls2 = 0;
		ipt = iparent[ib][k];
		kpt = k - 1;
		for (int ib2 = 1; ib2 <= nb[kpt]; ib2++)
		{
			if ((nchild[ib2][kpt] > 0) && (ib2 != ipt))
			{
				// Check wheter it is the adjacent parent box
				dx_xn = std::abs(xb_max[ib2][kpt] - xb_min[ipt][kpt]);
				dx_nx = std::abs(xb_min[ib2][kpt] - xb_max[ipt][kpt]);
				dy_xn = std::abs(yb_max[ib2][kpt] - yb_min[ipt][kpt]);
				dy_nx = std::abs(yb_min[ib2][kpt] - yb_max[ipt][kpt]);

				if ((((dy_xn <= Pars::tol2) || (dy_nx <= Pars::tol2)) &&
					 (xb_min[ib2][kpt] <= xb_max[ipt][kpt] + Pars::tol2) &&
					 (xb_max[ib2][kpt] >= xb_min[ipt][kpt] - Pars::tol2)) ||
					(((dx_xn <= Pars::tol2) || (dx_nx <= Pars::tol2)) &&
					 (yb_min[ib2][kpt] <= yb_max[ipt][kpt] + Pars::tol2) &&
					 (yb_max[ib2][kpt] >= yb_min[ipt][kpt] - Pars::tol2)))
				{
					for (int ic2 = 1; ic2 <= nchild[ib2][kpt]; ic2++)
					{
						ic = ichild[ib2][kpt][ic2];
						if (ic != ib)
						{
							if ((xb_min[ic][k] > xb_max[ib][k] + Pars::tol2) ||
								(xb_max[ic][k] < xb_min[ib][k] - Pars::tol2) ||
								(yb_min[ic][k] > yb_max[ib][k] + Pars::tol2) ||
								(yb_max[ic][k] < yb_min[ib][k] - Pars::tol2))
							{
								nls2 = nls2 + 1;
								ils2[nls2] = ic;
								kls2[nls2] = k;
							}
						}
					}
				}
			}
		}
	}
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::list_three(int ib, int k, const std::vector<int> &nb, int lev,
							  const std::vector<std::vector<double>> &xb_min, const std::vector<std::vector<double>> &yb_min,
							  const std::vector<std::vector<double>> &xb_max, const std::vector<std::vector<double>> &yb_max,
							  const std::vector<std::vector<int>> &nchild, const std::vector<std::vector<std::vector<int>>> &ichild,
							  int &nls3, std::vector<int> &ils3, std::vector<int> &kls3)
{
	// xb_min.resize(nbmrl+1,lmax+1);	xb_max.resize(nbmrl+1,lmax+1);
	// yb_min.resize(nbmrl+1,lmax+1);	yb_max.resize(nbmrl+1,lmax+1);
	// nb.resize(lmax+1); 		nchild.resize(nbmrl+1,lmax+1);
	// ils3.resize(nbl3+1); 	kls3.resize(nbl3+1);
	// ichild.resize(nbmrl+1,lmax+1,4+1);

	int istop, ic, kc;
	double dx_xn, dx_nx, dy_xn, dy_nx;

	if (nchild[ib][k] == 0)
	{
		nls3 = 0;
		for (int k2 = k; k2 <= lev; k2++)
		{
			for (int ib2 = 1; ib2 <= nb[k2]; ib2++)
			{
				istop = 0;
				if ((k2 == k) && (ib2 == ib))
				{
					istop = 1;
				}

				if ((nchild[ib2][k2] > 0) && (istop == 0))
				{
					dx_xn = std::abs(xb_max[ib2][k2] - xb_min[ib][k]);
					dx_nx = std::abs(xb_min[ib2][k2] - xb_max[ib][k]);
					dy_xn = std::abs(yb_max[ib2][k2] - yb_min[ib][k]);
					dy_nx = std::abs(yb_min[ib2][k2] - yb_max[ib][k]);

					if ((((dy_xn <= Pars::tol2) || (dy_nx <= Pars::tol2)) &&
						 (xb_min[ib2][k2] <= xb_max[ib][k] + Pars::tol2) &&
						 (xb_max[ib2][k2] >= xb_min[ib][k] - Pars::tol2)) ||
						(((dx_xn <= Pars::tol2) || (dx_nx <= Pars::tol2)) &&
						 (yb_min[ib2][k2] <= yb_max[ib][k] + Pars::tol2) &&
						 (yb_max[ib2][k2] >= yb_min[ib][k] - Pars::tol2)))
					{
						for (int ic2 = 1; ic2 <= nchild[ib2][k2]; ic2++)
						{
							ic = ichild[ib2][k2][ic2];
							kc = k2 + 1;
							if ((xb_min[ic][kc] > xb_max[ib][k] + Pars::tol2) ||
								(xb_max[ic][kc] < xb_min[ib][k] - Pars::tol2) ||
								(yb_min[ic][kc] > yb_max[ib][k] + Pars::tol2) ||
								(yb_max[ic][kc] < yb_min[ib][k] - Pars::tol2))
							{
								nls3 = nls3 + 1;
								ils3[nls3] = ic;
								kls3[nls3] = kc;
							}
						}
					}
				}
			}
		}
	}
	else if (nchild[ib][k] > 0)
	{
		nls3 = 0;
	}
}

// !!!!============================================================
// !!!!============================================================
void velocity_biot_savart::list_four(int ib, int k, const std::vector<int> &nb,
							 const std::vector<std::vector<double>> &xb_min, const std::vector<std::vector<double>> &yb_min,
							 const std::vector<std::vector<double>> &xb_max, const std::vector<std::vector<double>> &yb_max,
							 const std::vector<std::vector<int>> &iparent, const std::vector<std::vector<int>> &nchild,
							 int &nls4, std::vector<int> &ils4, std::vector<int> &kls4)
{
	// xb_min.resize(nbmrl+1,lmax+1);	xb_max.resize(nbmrl+1,lmax+1);
	// yb_min.resize(nbmrl+1,lmax+1);	yb_max.resize(nbmrl+1,lmax+1);
	// nb.resize(nbmrl+1);		nchild.resize(nbmrl+1,lmax+1);
	// ils4.resize(nbl4+1); 	kls4.resize(nbl4+1);
	// iparent.resize(nbmrl+1,lmax+1);

	int ipt, kpt;
	double dx_mxmn, dx_mnmx, dy_mxmn, dy_mnmx;

	ipt = iparent[ib][k];
	kpt = k - 1;

	if (k == 1)
	{
		nls4 = 0;
	}
	else if (k > 1)
	{
		nls4 = 0;
		for (int k2 = 1; k2 <= k - 1; k2++)
		{
			for (int ib2 = 1; ib2 <= nb[k2]; ib2++)
			{
				if (nchild[ib2][k2] == 0)
				{
					dx_mxmn = std::abs(xb_max[ib2][k2] - xb_min[ipt][kpt]);
					dx_mnmx = std::abs(xb_min[ib2][k2] - xb_max[ipt][kpt]);
					dy_mxmn = std::abs(yb_max[ib2][k2] - yb_min[ipt][kpt]);
					dy_mnmx = std::abs(yb_min[ib2][k2] - yb_max[ipt][kpt]);

					if ((((dy_mxmn <= Pars::tol2) || (dy_mnmx <= Pars::tol2)) &&
						 (xb_min[ib2][k2] <= xb_max[ipt][kpt] + Pars::tol2) &&
						 (xb_max[ib2][k2] >= xb_min[ipt][kpt] - Pars::tol2)) ||
						(((dx_mxmn <= Pars::tol2) || (dx_mnmx <= Pars::tol2)) &&
						 (yb_min[ib2][k2] <= yb_max[ipt][kpt] + Pars::tol2) &&
						 (yb_max[ib2][k2] >= yb_min[ipt][kpt] - Pars::tol2)))
					{
						if ((xb_min[ib2][k2] > xb_max[ib][k] + Pars::tol2) ||
							(xb_max[ib2][k2] < xb_min[ib][k] - Pars::tol2) ||
							(yb_min[ib2][k2] > yb_max[ib][k] + Pars::tol2) ||
							(yb_max[ib2][k2] < yb_min[ib][k] - Pars::tol2))
						{
							nls4 = nls4 + 1;
							ils4[nls4] = ib2;
							kls4[nls4] = k2;
						}
					}
				}
			}
		}
	}
}

// !!!!==============================================================
// !!!!==============================================================
void velocity_biot_savart::bico(int n, int &k, double &c)     // Combinatoric nCk = n!/(k!(n-k)!)
{
	int n_k;

	if ((k < 0) || (k > n) || (n < 0))
	{
		c = 0.0e0;
		std::cout << "\nn = " << n << " " << "k = " << k;
		std::cout << "\nerror: binomial input are out of range\n";
	}
	else if ((k == 0) || (k == n) || (n == 0))
	{
		c = 1.0e0;
	}
	else if ((k == 1) || (k == n - 1))
	{
		c = (double)(n);
	}
	else
	{
		if (k > n / 2)
		{
			k = n - k;
		}
		n_k = n - k;
		c = 1.0e0;
		for (int i = 1; i <= k; i++)
		{
			c = c * ((double)(n_k + i) / (double)(i));
		}
	}
}

// !!!!==============================================================
// !!!!==============================================================
void velocity_biot_savart::direct_sum(int nip, int npi, const std::vector<int> &ipp, const std::vector<double> &xi,
							  const std::vector<double> &yi, const std::vector<double> &si, std::vector<double> &ui, std::vector<double> &vi,
							  int njp, int npj, const std::vector<int> &jpp, const std::vector<double> &xj, const std::vector<double> &yj,
							  const std::vector<double> &sj, const std::vector<double> &gj, int icutoff)
{
	//int i, j;
	
	//#pragma omp parallel for
	for (int i2 = 1; i2 <= nip; i2++)
	{
		int i;
		i = ipp[i2-1]; //check again
		for (int j2 = 1; j2 <= njp; j2++)
		{
			int j;
			double q, dxij, dyij, rij2, sij2, rij, sij;
			j = jpp[j2-1]; //check again

			dxij = xi[i] - xj[j];
			dyij = yi[i] - yj[j];
			rij2 = std::pow(dxij, 2) + std::pow(dyij, 2);
			sij2 = 0.25e0 * (std::pow(si[i], 2) + std::pow(sj[j], 2)); // ! previously incorrect

			rij = std::sqrt(rij2);
			sij = std::sqrt(sij2);

			if (rij2 > 0.0e0)
			{
				regul_func_2d(rij, sij, q);
				ui[i] = ui[i] - q * dyij * gj[j] / rij2;
				vi[i] = vi[i] + q * dxij * gj[j] / rij2;
			}
		}
	}
}

// allocate-deallocate internal variables
// !!!! =============end of subroutines ============================
// !!!! =============end of subroutines ============================
// !!!! =============end of subroutines ============================
void velocity_biot_savart::amount_inbox_new(int n0, int n1, int npi, const std::vector<double> &xi, const std::vector<double> &yi,
								int n2, int n3, int npj, const std::vector<double> &xj, const std::vector<double> &yj, double xb_mini, double xb_maxi,
								double yb_mini, double yb_maxi, int &npr_in, int n_inter, std::vector<std::vector<std::vector<int>>> &numberofparticle, int parentx, int parenty, int x, int y)
{
	// internal variables
	int isame;

	if (n_inter == 1 && (parentx == -1 || parenty == -1) )
	{
		npr_in = 0;
		//printf("masuk A\n");
		//#pragma omp parallel for reduction(+:npr_in)
		//printf("%d\n", numberofparticle[x][y].size());
		for (int i = n0; i <= n1; i++)
		{
			if ((xi[i] >= xb_mini) && (xi[i] < xb_maxi) &&
				(yi[i] >= yb_mini) && (yi[i] < yb_maxi))
			{
				npr_in = npr_in + 1;
				numberofparticle[x][y].push_back(i);
			}
		}
	}

	if (n_inter == 1 && parentx != -1 )
	{
		npr_in = 0;
		//printf("masuk b\n");
		//printf("banyak partikel di cell ini %d ", numberofparticle[parentx][parenty].size());
		for (int i = 0; i < numberofparticle[parentx][parenty].size(); i++)
		{
			int ii = numberofparticle[parentx][parenty][i];
			//printf("%d ", ii);
			if ((xi[ii] >= xb_mini) && (xi[ii] < xb_maxi) &&
				(yi[ii] >= yb_mini) && (yi[ii] < yb_maxi))
			{

				npr_in = npr_in + 1;
				numberofparticle[x][y].push_back(ii);
			}
		}
		//printf("di cell %d, %d punya %d\n", x, y, numberofparticle[x][y].size());
	}
}
// !!!!============================================================