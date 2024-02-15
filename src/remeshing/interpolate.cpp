#include "remeshing.hpp"
#include "../grid_block/generateGrid.hpp"
#include "../grid_block/gridNodeNgh.hpp"
#include "../LSMPS/LSMPSb.hpp"

/**
 *  @brief  A particle interpolation subroutine. Interpolate the property data value into new 
 *  distribution from the source distribution.
 *         
 *  @param  _targetValue  [OUTPUT] Target property container.
 *  @param  _targetPar  Target particle data container.
 *  @param  _sourcePar  Source particle data container.
 *  @param  _sourceValue  Source property data.
*/
void remeshing::interpolate(
    std::vector<double> &_trgVal,
	const Particle &_trgPar,
	const Particle &_srcPar,
	const std::vector<double> &_srcVal)
{
    // Perform the LSMPS type b interpolation
    LSMPSb lsmpsCalc;
    #if (DIM == 2)
        lsmpsCalc.set_LSMPS_2D(_trgPar.x, _trgPar.y, _trgPar.s,
                            _srcPar.x, _srcPar.y, _srcPar.s, _srcVal, 
                            _trgPar.neighbor);
        _trgVal = lsmpsCalc.get_d0();
    #elif (DIM == 3)
        lsmpsCalc.set_LSMPS_3D(_trgPar.x, _trgPar.y, _trgPar.z, _trgPar.s, _trgVal,
                               _srcPar.x, _srcPar.y, _srcPar.z, _srcPar.s, _srcVal, 
                               _trgPar.neighbor, _trgPar.neighbor);
        _trgVal = lsmpsCalc.get_d0();
    #endif

    return;
}

/**
 *  A type of interpolation data for 3D simulation
 *    0:= Only interpolate velocity
 *    1:= Interpolate velocity and vorticity data
 *    2:= Addition of Q criterion data
*/
#define INTERPOLATE_3D_TYPE	0

/**
 *  @brief  A particle interpolation from source distribution data into a new
 *  distribution data. The interpolation is carried out using the LSMPS standard B.
 *         
 *  @param  _targetPar  [OUTPUT] Target particle data container.
 *  @param  _sourcePar  Source particle data container.
 *  @param  _baseGrid  	Grid node data of source particle.
*/
void remeshing::re_arrange_distribution(Particle &trgPar, const Particle &srcPar, const GridNode &nGrd){
	/** PROCEDURE:
	 * 	1. Generate the target particle in single resolution based on the user defined properties
	 *  2. Assign the leaf node that contains the particle
	 *  3. Evaluate the neighbor toward the source particle.
	 *  4. Interpolate the necessary properties.
	*/

	// Print prompt
	printf("%s\nRe-arrange distribution ... %s\n", FONT_TORQUOISE, FONT_RESET);
	
	// PROCEDURE 0:
	// ***********
	// Define the target particle distribution parameter
	
    // ========== USER DEFINED SECTION ========== 
    // ******************************************
    const double sigma = Pars::sigmaInt;	// The basic particle size (distribution spacing)
	const double _xdom = Pars::xdomInt;		// The x upstream gap distance
	std::vector<double> domLen = 			// The domain length
		{Pars::lxdomInt, Pars::lydomInt, Pars::lzdomInt};	
		// {5.0e0, 2.0e0, 2.0e0};
    // ************** END SECTION ***************

	int parCnt[DIM];		// The count of particle in each direction
	int parNum = 1;			// The total number of particle
	
    // Calculate the particle count
    basis_loop(d){
        // Update the particle count
        parCnt[d] = domLen[d] / sigma;
        // Update the domain length
        domLen[d] = parCnt[d] * sigma;
    }
    
    // Calculate the total particle number
	basis_loop(d) parNum *= parCnt[d];


	// PROCEDURE 1:
	// ***********
	// Set the basis parameter for particle generation
	// Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif

	double globalMinCoord[DIM];		// The global minimum coordinate
	globalMinCoord[0] = _xdom;		// Set the x coordinate basis
	for (int d = 1; d < DIM; d++){	// Set for the other basis
		globalMinCoord[d] = -domLen[d] / 2.0;
	}
	
	// Update particle data size
	trgPar.num = parNum;
	
	// Generate the target particle data
	#if (DIM == 2)
		// Resize the particle data
		trgPar.x = std::vector<double>(parNum, 0.0e0);
		trgPar.y = std::vector<double>(parNum, 0.0e0);
		
        int ID = 0;
		for (int j = 0; j < parCnt[1]; j++){
			for (int i = 0; i < parCnt[0]; i++){
				trgPar.x[ID] = globalMinCoord[0] + (i+0.5) * sigma;
				trgPar.y[ID] = globalMinCoord[1] + (j+0.5) * sigma;
				ID++;
			}
		}
	#elif (DIM == 3)
		// Resize the particle data
		trgPar.x = std::vector<double>(parNum, 0.0e0);
		trgPar.y = std::vector<double>(parNum, 0.0e0);
		trgPar.z = std::vector<double>(parNum, 0.0e0);

		int ID = 0;
		for (int k = 0; k < parCnt[2]; k++){
			for (int j = 0; j < parCnt[1]; j++){
				for (int i = 0; i < parCnt[0]; i++){
					trgPar.x[ID] = globalMinCoord[0] + (i+0.5) * sigma;
					trgPar.y[ID] = globalMinCoord[1] + (j+0.5) * sigma;
					trgPar.z[ID] = globalMinCoord[2] + (k+0.5) * sigma;
					ID++;
				}
			}
		}
	#endif
	
	// [DEBUG] Update particle data size
	if (parNum != ID)
	ERROR_LOG << "Something wrong with the ID calculation, not consitent!\n";

	// Particle generation summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Init. Generate new particle        [%f s]\n", _time);
	
	
	// PROCEDURE 2:
	// ***********
	// Assign the particle corresponding node
	// Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

	// Carry out the procedure below !
    generateGrid grid_tool;
    std::unordered_map<int, std::vector<int>> nodeParMap;   // A container of target particle
    grid_tool.assignNodeID(nGrd, nodeParMap, trgPar);

	// Particle generation summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Step 1: Node grouping              [%f s]\n", _time);


    // PROCEDURE 3:
	// ***********
	// Evaluate the neighbor toward the source particle
	// Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

	// Carry out the procedure below !
    GridNodeNgh ngh_tool;
    ngh_tool.find_inter_neighbor(trgPar.neighbor, trgPar, nodeParMap, nGrd, srcPar);

	// Particle generation summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Step 2: Evaluate neighbor          [%f s]\n", _time);


    // PROCEDURE 4:
	// ***********
	// Interpolate the necessary properties
	// Computational time accumulation
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

	// Carry out the procedure below !
    #if (DIM == 2)
        this->interpolate(trgPar.vorticity, trgPar, srcPar, srcPar.vorticity);  // Vorticity
        this->interpolate(trgPar.u, trgPar, srcPar, srcPar.u);                  // Velocity x direction
        this->interpolate(trgPar.v, trgPar, srcPar, srcPar.v);                  // Velocity y direction
    // #elif (DIM == 3)
        // Only interpolate the velocity
        if (INTERPOLATE_3D_TYPE == 0){
            this->interpolate(trgPar.u, trgPar, srcPar, srcPar.u);      // Velocity x direction
            this->interpolate(trgPar.v, trgPar, srcPar, srcPar.v);      // Velocity y direction
            this->interpolate(trgPar.w, trgPar, srcPar, srcPar.w);      // Velocity z direction
            
            // this->interpolate(trgPar.vortx, trgPar, srcPar, srcPar.vortx);   // Vorticity x direction
            // this->interpolate(trgPar.vorty, trgPar, srcPar, srcPar.vorty);   // Vorticity y direction
            // this->interpolate(trgPar.vortz, trgPar, srcPar, srcPar.vortz);   // Vorticity z direction
        }
        // Only interpolate the velocity and vorticity
        else if(INTERPOLATE_3D_TYPE == 1 || INTERPOLATE_3D_TYPE == 2){
            // LSMPSb lsmpsdu, lsmpsdv, lsmpsdw;
            LSMPSb lsmpsCalc;

            // Container of each velocity differential
            std::vector<double> dudx, dudy, dudz;
            std::vector<double> dvdx, dvdy, dvdz;
            std::vector<double> dwdx, dwdy, dwdz;

            // [Sec] Interpolation of differential of x velocity
            lsmpsCalc.set_LSMPS_3D(trgPar.x, trgPar.y, trgPar.z, trgPar.s, trgPar.u,
                                srcPar.x, srcPar.y, srcPar.z, srcPar.s, srcPar.u, 
                                trgPar.neighbor, trgPar.neighbor);
            // Interpolation
            trgPar.u = lsmpsCalc.get_d0();
            // Differential
            dudx = lsmpsCalc.get_ddx();
            dudy = lsmpsCalc.get_ddy();
            dudz = lsmpsCalc.get_ddz();


            // [Sec] Interpolation of differential of y velocity
            lsmpsCalc.set_LSMPS_3D(trgPar.x, trgPar.y, trgPar.z, trgPar.s, trgPar.v,
                                srcPar.x, srcPar.y, srcPar.z, srcPar.s, srcPar.v, 
                                trgPar.neighbor, trgPar.neighbor);
            // Interpolation
            trgPar.v = lsmpsCalc.get_d0();
            // Differential
            dvdx = lsmpsCalc.get_ddx();
            dvdy = lsmpsCalc.get_ddy();
            dvdz = lsmpsCalc.get_ddz();

            
            // [Sec] Interpolation of differential of z velocity
            lsmpsCalc.set_LSMPS_3D(trgPar.x, trgPar.y, trgPar.z, trgPar.s, trgPar.w,
                                srcPar.x, srcPar.y, srcPar.z, srcPar.s, srcPar.w, 
                                trgPar.neighbor, trgPar.neighbor);
            // Interpolation
            trgPar.w = lsmpsCalc.get_d0();
            // Differential
            dwdx = lsmpsCalc.get_ddx();
            dwdy = lsmpsCalc.get_ddy();
            dwdz = lsmpsCalc.get_ddz();


            // [Sec] Calculation of vorticity
            trgPar.vortx = std::vector<double>(trgPar.num);
            trgPar.vorty = std::vector<double>(trgPar.num);
            trgPar.vortz = std::vector<double>(trgPar.num);
            trgPar.vorticity = std::vector<double>(trgPar.num);
            #pragma omp parallel for
            for (int i = 0; i < trgPar.num; i++){
                // Calculation of each vorticity
                double vortX = dwdy[i] - dvdz[i];
                double vortY = dudz[i] - dwdx[i];
                double vortZ = dvdx[i] - dudy[i];

                // Assign the data into container
                trgPar.vortx[i] = vortX;
                trgPar.vorty[i] = vortY;
                trgPar.vortz[i] = vortZ;
                trgPar.vorticity[i] = std::sqrt(vortX*vortX + vortY*vortY + vortZ*vortZ);
            }


            // Additional data calculation of Q criterion (2nd Invariant of velocity gradient)
            if(INTERPOLATE_3D_TYPE == 2){
                // [Sec] Calculation of Q Criterion
                trgPar.P = std::vector<double>(trgPar.num);
                #pragma omp parallel for
                for (int i = 0; i < trgPar.num; i++){
                    // Calculation of each vorticity
                    double Q = -0.5 * (dudx[i]*dudx[i] + dvdy[i]*dvdy[i] + dwdz[i]*dwdz[i])
                                    - (dudy[i]*dvdx[i] + dudz[i]*dwdx[i] + dvdz[i]*dwdy[i]);
                    
                    // Assign the data into container
                    trgPar.P[i] = Q;
                }
            }
        }

    #endif

	// Particle generation summary time display
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
	printf("<-> Step 3: Interpolate property       [%f s]\n", _time);
	
	return;
}
