#include "data_saving.hpp"

void save_data::grid_data(int np, vector<double> &xp, vector<double> &yp, vector<double> &sp, vector<double> &gpz,
	vector<double> &up, vector<double> &vp, int &nvertn, vector<vector<double>> &vertex)
{
	// ===accessign struct data
	int ngrd, ngrx, ngry;
	double grid_space;
	vector<double> xmax, xmin, xg, yg, ug, vg, ggz;

	// // internal variables
	// vector<int> ncgrid(2);
	// vector<double> sg; // core size
	// int ix,iy,ip_grid, kernel_type, n_inter;

	// // === Calculate size of Grid for Allocation
	// // =====================================
	// grid_space = Pars::sigma; // New grids for redistribute_particles, spacing = core size
	//                   // =sigma_max  <=> grid spacing=core size then we can remesh VOTEX STRENGTH, 
	//                   // if not we have to /core-size^3=omega --> remesh VORTICITY omega
	// // =====================================
	// // === grid limitation -Preliminary ===
	// //  We can change size of box in here
	// xmin[0] = *min_element(xp.begin(), xp.begin()+np) - 0.8e0;
	// xmin[1] = *min_element(yp.begin(), yp.begin()+np) - 0.3e0;
	// xmax[0] = *max_element(xp.begin(), xp.begin()+np) + 0.5e0;
	// xmax[1] = *max_element(yp.begin(), yp.begin()+np) + 0.3e0;
	// //  Calculate grid size first for saving initial generating memory:
	// ncgrid[0] = std::ceil( ( xmax[0] - xmin[0] )/grid_space ) +1+2; // +2 to makesure, sometime only +1 not enough
	// ncgrid[1] = std::ceil( ( xmax[1] - xmin[1] )/grid_space ) +1+2; // +2 to makesure, sometime only +1 not enough
	// // === allocate redistribution grids ===
	// vector<double> xgrid(ncgrid[0]);
	// vector<double> ygrid(ncgrid[1]);
	// // === create grids for redistribution ===
	// d_base_grid.create_grid(0.0e0, grid_space, xmin[0], xmax[0], ncgrid[0], ngrx, xgrid);
	// d_base_grid.create_grid(0.0e0, grid_space, xmin[1], xmax[1], ncgrid[1], ngry, ygrid);
	// //  ============ Generating Grid ==================================
	// // - make Grid from Grid spacing we want, and Preliminary max min size of grid 
	// printf("Extruding Box - ngrx x ngry : %d x %d\n", ngrx, ngry);
	// // --------------------------------------------- 
	// ngrd = ngrx*ngry;
	// // -- resize grid data to be saved
	// xg.resize(ngrd); yg.resize(ngrd); ggz.resize(ngrd); ug.resize(ngrd); vg.resize(ngrd);
	// sg.resize(ngrd);

	// vector<vector<double>> ggz_grid(ngrx, vector<double>(ngry));
	// vector<vector<double>> upx_grid(ngrx, vector<double>(ngry));
	// vector<vector<double>> upy_grid(ngrx, vector<double>(ngry));
	// vector<double> ugp(ngrd), vgp(ngrd);

	// // ============ Coordinates =================
	// // for matlab reshape(f,ngry,ngrx) to get back the array data.
	// ip_grid = 0;
	// for (int ix = 0; ix < ngrx; ix++)
	// {
	// 	for (int iy = 0; iy < ngry; iy++)
	// 	{
	// 		xg[ip_grid]	= xgrid[ix];
	// 		yg[ip_grid]	= ygrid[iy];
	// 		ip_grid++;
	// 	}
	// }

	// // ============ Calculating VORTICITY on Grid-nodes ==================================================
	// //  Using iterpolation M4 to interpolate Strength from particle to grid: 
	// //   gpz --> ggz, 
	// kernel_type = Pars::opt_remesh_W;
	// d_base_remeshing.P2G(np, xp, yp, gpz, ngrx, ngry, xgrid, ygrid, grid_space, kernel_type, ggz_grid);
	// // =============Vorticity strengths =========
	// // saving 1D-array data
	// ip_grid = 0;
	// for (int ix = 0; ix < ngrx; ix++)
	// {
	// 	for (int iy = 0; iy < ngry; iy++)
	// 	{
	// 		ggz[ip_grid] = ggz_grid[ix][iy];
	// 		sg[ip_grid]	 = sp[0]; // In case of sigma not changed
	// 		ip_grid++;
	// 	}
	// }

	// // =========== Calculating VELOCITY on GRID-nodes ==============
	// if (Pars::opt_gdata_vel == 1)
	// {
	// 	// Using Biot-Savart to calculate grid velocity from grid Strength
	// 	// ggz and xg,yg --> ugp,vgp,wgp
	// 	std::fill(sp.begin(), sp.begin()+np, grid_space);
	// 	if (Pars::ioptfmm == 0) // still wrong, please do create poisson_2 that be able to take care two different distribution
	// 	{
	// 		d_base_poisson.biotsavart_direct_2d(1,ngrd,ngrd,xg,yg,sg,ugp,vgp,
	// 		1,np,np,xp,yp,sp,gpz,Pars::icutoff);
	// 	}
	// 	else if (Pars::ioptfmm == 1)
	// 	{
	// 		n_inter = 1;
	// 		d_base_poisson.biotsavart_fmm_2d(1,ngrd,ngrd,xg,yg,sg,ugp,vgp,
	// 		1,np,np,xp,yp,sp,gpz,Pars::icutoff,Pars::n_s,n_inter,Pars::ndp);
	// 	}
	// 	for (int i = 0; i < ngrd; i++)
	// 	{
	// 		ug[i] = ugp[i] + Pars::uin; // Biot-Savart calculate rotational velocity part only
	// 		vg[i] = vgp[i] + Pars::vin; // must be added irrotational velocity part
	// 	}
	// }

	// // --> calculate ug,vg,wg velocity grid nodes 
	// else if (Pars::opt_gdata_vel == 2)
	// {
	// 	// Using interpolation for Calculating Velocity of Grid-nodes from particles
	// 	d_base_remeshing.P2G(np,xp,yp,up, ngrx,ngry, xgrid,ygrid, grid_space,kernel_type, upx_grid); // x-dir
	// 	d_base_remeshing.P2G(np,xp,yp,vp, ngrx,ngry, xgrid,ygrid, grid_space,kernel_type, upy_grid); // y-dir
	// 	// ========== Velocity vectors =============
	// 	// saving 1D-array data
	// 	ip_grid = 0;
	// 	for (int ix = 0; ix < ngrx; ix++)
	// 	{
	// 		for (int iy = 0; iy < ngry; iy++)
	// 		{
	// 			ug[ip_grid]	= upx_grid[ix][iy];
	// 			vg[ip_grid]	= upy_grid[ix][iy];
	// 			ip_grid++;
	// 		}
	// 	}
	// }
	// // correcting near body
	// Body _b;
	// _b.uT.resize(1);
	// _b.vT.resize(1);
	// d_penalization.Penalization(0/*it*/, _b, ngrd, xg, yg, ug, vg, sg, ggz);

	// // deallocating memory
	// xgrid.clear(); ygrid.clear();
	// ggz_grid.clear();
	// upx_grid.clear(); upy_grid.clear();
	// ugp.clear(); vgp.clear();
	// sg.clear();
} // end of function
