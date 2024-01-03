#include "velocity_biot_savart.hpp"
#include <complex>

void velocity_biot_savart::biotsavart_fmm_2d(Particle &pi, Particle &pj, const int icutoff, const int n_s, const int n_inter, const int ndp)
{

	// #pragma omp declare reduction(vec_float_plus : std::vector<double> : \
    //               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    //                 initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	
	// #pragma omp declare reduction(vec_float_minus : std::vector<double> : \
    //               std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::minus<double>())) \
    //                 initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

	// Shifting data in order to start indexing from 1
	int n0 = 1;
	int &n1 = pi.num;
	int &npi = pi.num;
	int n2 = 1;
	int &n3 = pj.num;
	int &npj = pj.num;
	std::vector<double> xi(npi + 1), yi(npi + 1), si(npi + 1), ui(npi + 1), vi(npi + 1);
	std::vector<double> xj(npj + 1), yj(npj + 1), sj(npj + 1), gj(npj + 1);
	double start,finish,total;
	
	//#pragma omp parallel for
	for (int i = 1; i <= npi; i++)
	{
		xi[i] = pi.x[i - 1];
		yi[i] = pi.y[i - 1];
		si[i] = pi.s[i - 1];
		ui[i] = 0.0e0;
		vi[i] = 0.0e0;
	}

	//#pragma omp parallel for 
	for (int i = 1; i <= npj; i++)
	{
		xj[i] = pj.x[i - 1];
		yj[i] = pj.y[i - 1];
		sj[i] = pj.s[i - 1];
		gj[i] = pj.gz[i - 1];
	}
	

	// internal variables
	std::complex<double> cmw, cmw2, ci, z, zi, zj, zp, zc, zo, zg, zm;
	int lev, ibtt, njp, nls1, nls2, nls3, nls4, nip;
	double xmin, ymin, xmax, ymax, cbi;
	std::vector<std::vector<std::vector<int>>> particleinbox;
	particleinbox.resize(Pars::lmax + 1, std::vector<std::vector<int>>(Pars::nbmrl+1));
	std::vector<std::vector<double>> xb_min(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<std::vector<double>> xb_max(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<std::vector<double>> yb_min(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<std::vector<double>> yb_max(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<std::vector<double>> xb_cen(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<std::vector<double>> yb_cen(Pars::nbmrl + 1, std::vector<double>(Pars::lmax + 1, 0.0e0));
	std::vector<int> nb(Pars::lmax + 1, 0);
	std::vector<std::vector<int>> iparent(Pars::nbmrl + 1, std::vector<int>(Pars::lmax + 1, 0));
	std::vector<std::vector<int>> nchild(Pars::nbmrl + 1, std::vector<int>(Pars::lmax + 1, 0));
	std::vector<std::vector<std::vector<int>>> ichild(Pars::nbmrl + 1, std::vector<std::vector<int>>(Pars::lmax + 1, std::vector<int>(4 + 1, 0))); 
	std::vector<int> ipp(npi + 1, 0);
	std::vector<int> jpp(npj + 1, 0);
	std::vector<int> ils1(Pars::nbl1 + 1, 0);
	std::vector<int> kls1(Pars::nbl1 + 1, 0);
	std::vector<int> ils2(Pars::nbl2 + 1, 0);
	std::vector<int> kls2(Pars::nbl2 + 1, 0);
	std::vector<int> ils3(Pars::nbl3 + 1, 0);
	std::vector<int> kls3(Pars::nbl3 + 1, 0);
	std::vector<int> ils4(Pars::nbl4 + 1, 0);
	std::vector<int> kls4(Pars::nbl4 + 1, 0);
	std::vector<std::vector<int>> ibt(Pars::nbmrl + 1, std::vector<int>(Pars::lmax + 1, 0));

	std::vector<std::vector<std::complex<double>>> ak(Pars::nbmax + 1, std::vector<std::complex<double>>(Pars::npcm + 1, std::complex<double>(0.0e0, 0.0e0)));
	std::vector<std::vector<std::complex<double>>> bk(Pars::nbmax + 1, std::vector<std::complex<double>>(Pars::npcm + 1, std::complex<double>(0.0e0, 0.0e0)));

	ci = std::complex<double>(0.0, 1.0);

	// !!!!=== step 1 === //Build a Tree (bisa dioptimize dengan nyimpan treenya di private variable (sepertinya))
	start = omp_get_wtime();
	hierarchy_mesh(n0, n1, npi, xi, yi, n2, n3, npj, xj, yj, n_s, n_inter, xmin, ymin, xmax, ymax,
				   nb, xb_min, yb_min, xb_max, yb_max, xb_cen, yb_cen, lev, nchild, ichild, iparent, particleinbox);
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> Tree finished in : %f s\n", finish-start);

	
	int nTotal = lev * 6 - 2;
	//ProgressBar progressBar(nTotal, 40, '|', ' ');
	// !!!!===============================================
	
	start = omp_get_wtime();
	ibtt = 0;
	//ndp = number of multipole expansion
	for (int k = 1; k <= lev; k++) //check every level
	{
		//++progressBar;		   // record the tick
		//progressBar.display(); // display the bar only at certain steps
		
		for (int ib = 1; ib <= nb[k]; ib++) //check every boxes in the level
		{
			ibtt = ibtt + 1;
			ibt[ib][k] = ibtt;
			for (int kp = 1; kp <= (ndp + 1); kp++)
			{
				if (kp <= ndp)
				{
					ak[ibt[ib][k]][kp] = std::complex<double>(0.0e0, 0.0e0);
				}
				bk[ibt[ib][k]][kp] = std::complex<double>(0.0e0, 0.0e0);
			}
		}
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 1 : %f s\n", finish-start);
	//============= finished initialization =============


	//====================== Multipole construction==============================
	// !!!!=== step 2.1 ====================== Multipole expansion
	start = omp_get_wtime();
	for (int k = 1; k <= lev; k++) //check every level
	{
		//++progressBar;		   // record the tick
		//progressBar.display(); // display the bar only at certain steps
		
		for (int ib = 1; ib <= nb[k]; ib++) //check every boxes in the level
		{
			if (nchild[ib][k] == 0)	//if the box is childless
			{
				//par_loc(n2, n3, npj, xj, yj, xb_min[ib][k], xb_max[ib][k],
				//		yb_min[ib][k], yb_max[ib][k], njp, jpp);
				njp = particleinbox[k][ib].size(); 							//number of particles inside the box
				jpp = particleinbox[k][ib]; 								//it starts from 0 - (n-1)
			
				zm = std::complex<double>(xb_cen[ib][k], yb_cen[ib][k]);    //childless box center

				if (njp > 0)     //if there is particles inside the box
				{
					for (int kp = 1; kp <= ndp; kp++)  // There are 10 expansion of multipoles
					{
						for (int jp2 = 1; jp2 <= njp; jp2++) //for all particles inside the box
						{
							int jp = jpp[jp2-1]; 															//the selected particle index
							zj = std::complex<double>(xj[jp], yj[jp]);										//the location of the selected particle
							ak[ibt[ib][k]][kp] = ak[ibt[ib][k]][kp] + gj[jp] * pow((zj - zm), (kp - 1));	//multipole expansion for the selected childless box, <ibt[ib][k], is the index of the box and kp is its multipole>
						}
					}
				}
			}
		}
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 2 : %f s\n", finish-start);
	//printf("checkpoint 2\n");

	
	// !!!!=== step 2.2 === Multipole shifting bottom up 
	start = omp_get_wtime();
	for (int k = lev - 1; k >= 1; k--) //check every level from bottom to up
	{
		//++progressBar;		   // record the tick
		//progressBar.display(); // display the bar only at certain steps
		
		for (int ib = 1; ib <= nb[k]; ib++)  //check all the boxes at the current level
		{
			if (nchild[ib][k] > 0)			 //if the box has children, shift the multipole from children to the parent
			{
				zp = std::complex<double>(xb_cen[ib][k], yb_cen[ib][k]);   //the parent box center coordinate
				for (int ic2 = 1; ic2 <= nchild[ib][k]; ic2++)			   //check the child box of the selected parent 
				{
					int ic = ichild[ib][k][ic2];										//the child box index
					zc = std::complex<double>(xb_cen[ic][k + 1], yb_cen[ic][k + 1]);	//the child box center coordinate
					for (int kl = 1; kl <= ndp; kl++) 	//From all of the multipole
					{
						for (int kp = 1; kp <= kl; kp++) 
						{
							int kpm1 = kp - 1;
							bico(kl - 1, kpm1, cbi);
							ak[ibt[ib][k]][kl] = ak[ibt[ib][k]][kl] + cbi * pow((zc - zp), (kl - kp)) * ak[ibt[ic][k + 1]][kp];
						}
					}
				}
			}
		}
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 3 : %f s\n", finish-start);
	//=============================== Finished Multipole expansion =============================

	double total1 = 0, total2 = 0, total3 = 0, start1 = 0;
	double total4 = 0, total5 = 0;
 
	// !!!!=== step 4 & 6 === Local expansion ??
	start = omp_get_wtime(); 
	for (int k = 1; k <= lev; k++){//check every level
		//++progressBar;		   // record the tick
		//progressBar.display(); // display the bar only at certain steps

		for (int ib = 1; ib <= nb[k]; ib++) //check all of the boxes at the level
		{
			zg = std::complex<double>(xb_cen[ib][k], yb_cen[ib][k]); //the center of the selected box

			start1 = omp_get_wtime();
			// !!!! ==== step 4 ===
			list_two(ib, k, nb, xb_min, yb_min, xb_max, yb_max,
					 iparent, nchild, ichild, nls2, ils2, kls2);
			
			total1 += omp_get_wtime() - start1;
			start1 = omp_get_wtime();

			if (nls2 > 0)
			{
				for (int ib1 = 1; ib1 <= nls2; ib1++)
				{
					int ib2 = ils2[ib1];
					int k2 = kls2[ib1];
					zo = std::complex<double>(xb_cen[ib2][k2], yb_cen[ib2][k2]);
					for (int kl = 1; kl <= ndp + 1; kl++)
					{
						int ll = kl - 1;
						for (int kp = 1; kp <= ndp; kp++)
						{
							int kpm1 = kp - 1;
							bico(ll + kp - 1, kpm1, cbi);
							bk[ibt[ib][k]][kl] = bk[ibt[ib][k]][kl] + pow((-1.0e0), kp) * cbi *
																		  ak[ibt[ib2][k2]][kp] / (pow((zo - zg), (kp + ll)));
						}
					}
				}
			}
			total2 += omp_get_wtime() - start1;
			start1 = omp_get_wtime();
			// !!!! ==== step 6 ===
			list_four(ib, k, nb, xb_min, yb_min, xb_max, yb_max,
					  iparent, nchild, nls4, ils4, kls4);
			
			total3 += omp_get_wtime() - start1;
			start1 = omp_get_wtime();
			if (nls4 > 0)
			{
				for (int ib1 = 1; ib1 <= nls4; ib1++)
				{
					int ib2 = ils4[ib1];
					int k2 = kls4[ib1];
					//par_loc(n2, n3, npj, xj, yj, xb_min[ib2][k2], xb_max[ib2][k2],
					//		yb_min[ib2][k2], yb_max[ib2][k2], njp, jpp);
					njp = particleinbox[k2][ib2].size();
					jpp = particleinbox[k2][ib2]; 
					
					
					for (int jp2 = 1; jp2 <= njp; jp2++)
					{
						int jp = jpp[jp2-1];
						zo = std::complex<double>(xj[jp], yj[jp]);
						for (int kl = 1; kl <= ndp + 1; kl++)
						{
							int ll = kl - 1;
							bk[ibt[ib][k]][kl] = bk[ibt[ib][k]][kl] - gj[jp] / (pow((zo - zg), (1 + ll)));
						}
					}
				}
			}
			
			total4 += omp_get_wtime() - start1;
			start1 = omp_get_wtime();
		}
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 4 : %f s\n", finish-start);
	printf("List 2  = %f s\n", total1);
	printf("M2L     = %f s\n", total2);
	printf("List 4  = %f s\n", total3);
	printf("Outsider= %f s\n", total4);
	//printf("checkpoint 4\n");
	
	// !!!!=== step 7 ===
	start = omp_get_wtime();
	for (int k = 1; k <= lev - 1; k++)
	{
		//++progressBar; // record the tick
		// usleep(200); // simulate work
		// if (i % 10 == 0)
		//progressBar.display(); // display the bar only at certain steps
		//#pragma omp parallel for shared(nchild, bk, xb_cen, yb_cen, ibt) private(zp, zc, cbi)
		for (int ib = 1; ib <= nb[k]; ib++)
		{
			if (nchild[ib][k] > 0)
			{
				zp = std::complex<double>(xb_cen[ib][k], yb_cen[ib][k]);
				for (int ic = 1; ic <= nchild[ib][k]; ic++)
				{
					int ib2 = ichild[ib][k][ic];
					zc = std::complex<double>(xb_cen[ib2][k + 1], yb_cen[ib2][k + 1]);
					for (int kl = 1; kl <= ndp + 1; kl++)
					{
						int ll = kl - 1;
						for (int kp2 = kl; kp2 <= ndp + 1; kp2++)
						{
							int kp = kp2 - 1;
							bico(kp, ll, cbi);
							bk[ibt[ib2][k + 1]][kl] = bk[ibt[ib2][k + 1]][kl] + cbi * pow((zc - zp), (kp - ll)) * bk[ibt[ib][k]][kp2];
						}
					}
				}
			}
		}
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 5 : %f s\n", finish-start);
	//printf("checkpoint 5\n");
	total1 = 0, total2 = 0, total3 = 0, start1 = 0;
	total4 = 0, total5 = 0;
	// !!!!=== step 3 & 5 & 8 ===
	start = omp_get_wtime();
	for (int k = 1; k <= lev; k++)
	{
		for (int ib = 1; ib <= nb[k]; ib++)
		{
			if (nchild[ib][k] == 0)
			{
				start1 = omp_get_wtime();
				//par_loc(n0, n1, npi, xi, yi, xb_min[ib][k], xb_max[ib][k],
				//		yb_min[ib][k], yb_max[ib][k], nip, ipp);
				nip = particleinbox[k][ib].size();
				ipp = particleinbox[k][ib];
				
				zg = std::complex<double>(xb_cen[ib][k], yb_cen[ib][k]);

				// !!!! ==== step 3 ===
				list_one(ib, k, nb, lev, xb_min, yb_min, xb_max, yb_max,
						 nchild, nls1, ils1, kls1);
				
				total1 += omp_get_wtime() - start1;
				start1 = omp_get_wtime();

				if ((nip > 0) && (nls1 > 0))
				{
					for (int il = 1; il <= nls1; il++)
					{
						int ib2 = ils1[il];
						int k2 = kls1[il];
						//par_loc(n2, n3, npj, xj, yj, xb_min[ib2][k2], xb_max[ib2][k2],
						//		yb_min[ib2][k2], yb_max[ib2][k2], njp, jpp);
						njp = particleinbox[k2][ib2].size();
						jpp = particleinbox[k2][ib2];
						
						if (njp > 0)
						{
							direct_sum(nip, npi, ipp, xi, yi, si, ui, vi, njp, npj, jpp, xj, yj, sj, gj, Pars::icutoff);
						}
					}
				}
				total2 += omp_get_wtime() - start1;
				start1 = omp_get_wtime();
				// !!!! ==== step 5 ===
				list_three(ib, k, nb, lev, xb_min, yb_min, xb_max, yb_max,
						   nchild, ichild, nls3, ils3, kls3);
				
				total3 += omp_get_wtime() - start1;
				start1 = omp_get_wtime();

				if ((nip > 0) && (nls3 > 0))
				{
					//printf("%d, %d\n", nip, nls3);
					//#pragma omp parallel for reduction(vec_float_plus: ui) reduction(vec_float_minus: vi) schedule(dynamic)
					for (int ip2 = 1; ip2 <= nip; ip2++)
					{
						int ip = ipp[ip2-1]; //check each particle inside the evaluated box
						cmw = std::complex<double>(0.0e0, 0.0e0);
						z = std::complex<double>(xi[ip], yi[ip]);
						
						for (int il = 1; il <= nls3; il++)
						{
							//int ic = ils3[il];
							//int ikc = kls3[il];
							zm = std::complex<double>(xb_cen[ils3[il]][kls3[il]], yb_cen[ils3[il]][kls3[il]]);
							std::complex<double> zmm = std::complex<double>(xb_cen[ils3[il]][kls3[il]], yb_cen[ils3[il]][kls3[il]]);
							for (int kp = 1; kp <= ndp; kp++) 
							{
								cmw = cmw + ak[ibt[ils3[il]][kls3[il]]][kp] / (pow((z - zm), kp));
							}
						}

						cmw2 = -ci * cmw / (2.0e0 * M_PI);

						ui[ip] = ui[ip] + std::real(cmw2);
						vi[ip] = vi[ip] - std::imag(cmw2);						
					}
				}
				total4 += omp_get_wtime() - start1;
				start1 = omp_get_wtime();
				// !            !!!! ==== step 8 ===
				if (nip > 0)
				{
					//#pragma omp parallel for reduction(vec_float_plus: ui) reduction(vec_float_minus: vi) schedule(dynamic)
					for (int ip2 = 1; ip2 <= nip; ip2++) //check
					{
						//int ip = ipp[ip2-1];
						zi = std::complex<double>(xi[ipp[ip2-1]], yi[ipp[ip2-1]]);
						cmw = std::complex<double>(0.0e0, 0.0e0);

						//#pragma omp paralel for 
						for (int kl = 1; kl <= ndp + 1; kl++)
						{
							cmw = cmw + bk[ibt[ib][k]][kl] * pow((zi - zg), kl-1);
						}

						cmw2 = -ci * cmw / (2.0e0 * M_PI);
						//printf("%f, %f\n", std::real(cmww), std::imag(cmww));
						ui[ipp[ip2-1]] = ui[ipp[ip2-1]] + std::real(cmw2);
						vi[ipp[ip2-1]] = vi[ipp[ip2-1]] - std::imag(cmw2);
					}
				}
				total5 += omp_get_wtime() - start1;
			}
		}
		//++progressBar; // record the tick
		// usleep(200); // simulate work
		// if (i % 10 == 0)
		//progressBar.display(); // display the bar only at certain steps
	}
	finish = omp_get_wtime();
	total += (finish-start);
	printf("<+> checkpoint 6 : %f s\n", finish-start);
	printf("List 1       = %f s\n", total1);
	printf("Direct Calc  = %f s\n", total2);
	printf("List 3       = %f s\n", total3);
	printf("Far outsider = %f s\n", total4);
	printf("Farfield     = %f s\n", total5);
	printf("<+> Total Time : %f s\n", total);
	
	//progressBar.done(); // terminating progression display
	//printf("checkpoint 6\n");

	#pragma omp parallel for schedule(static)
	for (int i = 1; i <= npi; i++)
	{
		pi.u[i - 1] = ui[i];
		pi.v[i - 1] = vi[i];
	}

	//save quadtree data//
	
	// std::ofstream ofs;
	// ofs.open("treedata.csv");
	// ofs << "level" << "," << "Cell number" << "," << "number of child" << "," << "Parent" <<"," << "number of particles\n";
	// for (int i = 1; i <= Pars::lmax; i++){
	// 	for (int jj = 1; jj <= nb[i]; jj++){
	// 		//printf("%d\n", particleinbox[i][jj].size());
	// 		ofs<< i << "," << jj << "," << nchild[jj][i] << "," <<iparent[jj][i] << "," << particleinbox[i][jj].size()<<"\n";
	// 	}
	// } 
	// ofs.close(); 

	// -- deallocating memory
	xb_min.clear();
	xb_max.clear();
	yb_min.clear();
	yb_max.clear();
	xb_cen.clear();
	yb_cen.clear();
	nb.clear();
	
	iparent.clear();
	nchild.clear();
	ichild.clear();
	ipp.clear();
	jpp.clear();
	ils1.clear();
	kls1.clear();
	ils2.clear();
	kls2.clear();

	ils3.clear();
	kls3.clear();
	ils4.clear();
	
	kls4.clear();
	
	ibt.clear();
	//printf("checkpoint 8\n");
	ak.clear();
	bk.clear();


	xi.clear();
	yi.clear();
	si.clear();
	ui.clear();
	vi.clear();
	xj.clear();
	yj.clear();
	sj.clear();
	gj.clear();
}
