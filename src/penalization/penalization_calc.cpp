#include "penalization.hpp"
#include "../LSMPS/LSMPSa.hpp"
#include "../velocity_calculation/velocity_calc.hpp"

/**
 *  @brief Single penalization calculation. Update the velocity and vorticity 
 *  based on penalization calculation.
 *  
 *  @param	_tempParticle  The temporary particle data for calculation.
 *  @param	_baseParticle  The original particle data in calculation.
 *  @param  _bodyList  The body data list.
 *  @param  _step  The current simulation iteration step.
 */
void penalization::no_slip(Particle &_p, 
						   Particle &_basePar, 
						   const std::vector<Body> &_bodyList, 
						   int step) const
{
	/* Procedure:
        1. Define the body velocity
        2. Calculate velocity penalization
        3. Update vorticity by penalization
        4. Update the velocity and vorticity into the origin particle
	*/

	// Declaring internal variable
	// Velocity differential
	std::vector<double> du(_p.num, 0.0e0);
	std::vector<double> dv(_p.num, 0.0e0);
	// Penalized velocity
	std::vector<double> u_pen(_p.num, 0.0e0);
	std::vector<double> v_pen(_p.num, 0.0e0);
	// Solid body velocity
	std::vector<double> uS(_p.num, 0.0e0);
	std::vector<double> vS(_p.num, 0.0e0);


	// PROCEDURE 1: Define the body velocity
    // ********************************************************************
	// Calculate the solid velocity
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++)
	{
		// Rotation still not modeled
		double uR = -(_p.y[i] - Pars::ycenter) * Pars::omega;      // Body node x-velocity influenced by body rotation
		double vR = (_p.x[i] - Pars::xcenter) * Pars::omega;      // Body node y-velocity influenced by body rotation

		// Tranlation take from the body variable
		double uT = _bodyList[_p.bodyPart[i]].uT[step];    // Body node x-velocity influenced by body translation at each iteration
		double vT = _bodyList[_p.bodyPart[i]].vT[step];    // Body node y-velocity influenced by body translation at each iteration
		
		// // Still not consider the U_deformation yet!
		// double uDEF = 0.0e0;    // Body node x-velocity influenced by body deformation
		// double vDEF = 0.0e0;    // Body node y-velocity influenced by body deformation
		
		// Node resultant velocity
		uS[i] = uR + uT/* + uDEF*/;
		vS[i] = vR + vT/* + vDEF*/;
	}


    // PROCEDURE 2: Calculate velocity penalization
    // ********************************************************************
	/* Note:
		> There are implicit, semi-implicit, and explicit schemes
		> Velocity penalization calculation 1 for: Implicit and semi implicit
		> Velocity penalization calculation 2 for: Explicit
	*/
	// Velocity penalization calculation 1 (IMPLICIT + SEMI)
	if (Pars::opt_pen == 1 || Pars::opt_pen == 2)
	{
		// if    (Pars::opt_pen == 1)  {std::cout << "<+> Penalization type               :      IMPLICIT\n";}
		// else/*(Pars::opt_pen == 2)*/{std::cout << "<+> Penalization type               : SEMI-IMPLICIT\n";}
		
		#pragma omp parallel for
		for (int i = 0; i < _p.num; i++)
		{
			// Consider all velocity
			u_pen[i] = (_p.u[i] + Pars::lambda*Pars::dt* _p.chi[i] * uS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			v_pen[i] = (_p.v[i] + Pars::lambda*Pars::dt* _p.chi[i] * vS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));

			// // Only consider the rotational velocity
			// u_pen[i] = (_p.u[i] + Pars::lambda * Pars::dt * _p.chi[i] * (uS[i] - Pars::u_inf)) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			// v_pen[i] = (_p.v[i] + Pars::lambda * Pars::dt * _p.chi[i] * (vS[i] - Pars::v_inf)) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
		}
	}
	// Velocity penalization calculation 2 (EXPLICIT)
	else if (Pars::opt_pen == 3)
	{
		// std::cout << "<+> Penalization type               :      EXPLICIT\n";
		#pragma omp parallel for
		for (int i = 0; i < _p.num; i++)
		{
			// Consider all velocity
			u_pen[i] = (1 - _p.chi[i]) * _p.u[i] + _p.chi[i] * (uS[i]);
			v_pen[i] = (1 - _p.chi[i]) * _p.v[i] + _p.chi[i] * (vS[i]);
			
			// // Only consider the rotational velocity
			// u_pen[i] = (1 - _p.chi[i]) * _p.u[i] + _p.chi[i] * (uS[i] - Pars::u_inf);
			// v_pen[i] = (1 - _p.chi[i]) * _p.v[i] + _p.chi[i] * (vS[i] - Pars::v_inf);
		}
	}
	
	// Calculate the velocity difference
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		du[i] = u_pen[i] - _p.u[i];
		dv[i] = v_pen[i] - _p.v[i];
	}
	
	// Update temporary particle velocity
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		_p.u[i] = u_pen[i];
		_p.v[i] = v_pen[i];
	}


	// PROCEDURE 3: Update vorticity by penalization
    // ********************************************************************
	// Updating vorticity from calculation
	// Performing LSMPS calculation
	LSMPSa lsmpsa_du;
	LSMPSa lsmpsa_dv;
	lsmpsa_du.set_LSMPS(_p.x, _p.y, _p.s, du, _p.neighbor);
	lsmpsa_dv.set_LSMPS(_p.x, _p.y, _p.s, dv, _p.neighbor);
	std::vector<double> _ddvdx = lsmpsa_dv.get_ddx();
	std::vector<double> _ddudy = lsmpsa_du.get_ddy();
	
	// Evaluating added penalization vorticity field [dw = grad (cross) du]
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		// Calculate the vorticity and vortex strength
		double _dvor = _ddvdx[i] - _ddudy[i];         // w = del x u
		double _dgpz = _dvor * std::pow(_p.s[i], 2);  // γ = w * Area

		// Update the vorticity and vortex strength
		_p.vorticity[i] += _dvor;
		_p.gz[i] += _dgpz;
	}


	// PROCEDURE 4: Update the penalization result into the original particle
    // ********************************************************************
    // Update the penalization value to original particle variable
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++){
        // Aliasing the original ID
        const int &ori_ID = this->baseID[i];

        // Assign each new data to original data
        _basePar.u[ori_ID]  = _p.u[i];
        _basePar.v[ori_ID]  = _p.v[i];
        _basePar.gz[ori_ID] = _p.gz[i];
        _basePar.vorticity[ori_ID] = _p.vorticity[i];
    }
	
	return;
}

/**
 *  @brief Single penalization calculation for 3D simulation. Update the velocity
 *  and vorticity based on penalization calculation.
 *  
 *  @param	_tempParticle  The temporary particle data for calculation.
 *  @param	_baseParticle  The original particle data in calculation.
 *  @param  _bodyList  The body data list.
 *  @param  _step  The current simulation iteration step.
 */
void penalization::no_slip_3d(Particle &_p, 
						      Particle &_basePar, 
						      const std::vector<Body> &_bodyList, 
						      int step) const
{
	/* Procedure:
        1. Define the body velocity
        2. Calculate velocity penalization
        3. Update vorticity by penalization
        4. Update the velocity and vorticity into the origin particle
	*/

	// Declaring internal variable
	// Velocity differential
	std::vector<double> du(_p.num, 0.0e0);
	std::vector<double> dv(_p.num, 0.0e0);
	std::vector<double> dw(_p.num, 0.0e0);
	// Penalized velocity
	std::vector<double> u_pen(_p.num, 0.0e0);
	std::vector<double> v_pen(_p.num, 0.0e0);
	std::vector<double> w_pen(_p.num, 0.0e0);
	// Solid body velocity
	std::vector<double> uS(_p.num, 0.0e0);
	std::vector<double> vS(_p.num, 0.0e0);
	std::vector<double> wS(_p.num, 0.0e0);


	// PROCEDURE 1: Define the body velocity
    // ********************************************************************
	// Calculate the solid velocity
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++)
	{
		// // Rotation still not modeled
		// double uR = -(_p.y[i] - Pars::ycenter) * Pars::omega;      // Body node x-velocity influenced by body rotation
		// double vR = (_p.x[i] - Pars::xcenter) * Pars::omega;      // Body node y-velocity influenced by body rotation
		// double wR = (_p.x[i] - Pars::xcenter) * Pars::omega;      // Body node z-velocity influenced by body rotation

		// Tranlation take from the body variable
		double uT = _bodyList[_p.bodyPart[i]].uT[step];    // Body node x-velocity influenced by body translation at each iteration
		double vT = _bodyList[_p.bodyPart[i]].vT[step];    // Body node y-velocity influenced by body translation at each iteration
		double wT = _bodyList[_p.bodyPart[i]].wT[step];    // Body node z-velocity influenced by body translation at each iteration
		
		// // Still not consider the U_deformation yet!
		// double uDEF = 0.0e0;    // Body node x-velocity influenced by body deformation
		// double vDEF = 0.0e0;    // Body node y-velocity influenced by body deformation
		// double wDEF = 0.0e0;    // Body node z-velocity influenced by body deformation
		
		// Node resultant velocity
		uS[i] = /*uR + */uT/* + uDEF*/;
		vS[i] = /*vR + */vT/* + vDEF*/;
		wS[i] = /*wR + */wT/* + wDEF*/;
	}


    // PROCEDURE 2: Calculate velocity penalization
    // ********************************************************************
	/* Note:
		> There are implicit, semi-implicit, and explicit schemes
		> Velocity penalization calculation 1 for: Implicit and semi implicit
		> Velocity penalization calculation 2 for: Explicit
	*/
	// Velocity penalization calculation 1 (IMPLICIT and SEMI)
	if (Pars::opt_pen == 1 || Pars::opt_pen == 2)
	{
		#pragma omp parallel for
		for (int i = 0; i < _p.num; i++)
		{
			// Consider velocity resultant
			u_pen[i] = (_p.u[i] + Pars::lambda*Pars::dt* _p.chi[i] * uS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			v_pen[i] = (_p.v[i] + Pars::lambda*Pars::dt* _p.chi[i] * vS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
			w_pen[i] = (_p.w[i] + Pars::lambda*Pars::dt* _p.chi[i] * wS[i]) / (1.0e0 + (Pars::lambda * Pars::dt * _p.chi[i]));
		}
	}
	// Velocity penalization calculation 2 (EXPLICIT)
	else if (Pars::opt_pen == 3)
	{
		#pragma omp parallel for
		for (int i = 0; i < _p.num; i++)
		{
			// Consider velocity resultant
			u_pen[i] = (1 - _p.chi[i]) * _p.u[i] + _p.chi[i] * (uS[i]);
			v_pen[i] = (1 - _p.chi[i]) * _p.v[i] + _p.chi[i] * (vS[i]);
			w_pen[i] = (1 - _p.chi[i]) * _p.w[i] + _p.chi[i] * (wS[i]);
		}
	}
	
	// Calculate the velocity difference
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		du[i] = u_pen[i] - _p.u[i];
		dv[i] = v_pen[i] - _p.v[i];
		dw[i] = w_pen[i] - _p.w[i];
	}
	
	// Update temporary particle velocity
	#pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		_p.u[i] = u_pen[i];
		_p.v[i] = v_pen[i];
		_p.w[i] = w_pen[i];
	}


	// PROCEDURE 3: Update vorticity by penalization
    // ********************************************************************
	// Updating vorticity from calculation
	// Performing LSMPS calculation
	LSMPSa lsmpsa;
	lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, du, 
						_p.x, _p.y, _p.z, _p.s, du, _p.neighbor);
	std::vector<double> _ddudy = lsmpsa.get_ddy();	// [d(du)/dz]
	std::vector<double> _ddudz = lsmpsa.get_ddz();	// [d(du)/dy]

	lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, dv, 
						_p.x, _p.y, _p.z, _p.s, dv, _p.neighbor);
	std::vector<double> _ddvdx = lsmpsa.get_ddx();	// [d(dv)/dx]
	std::vector<double> _ddvdz = lsmpsa.get_ddz();	// [d(dv)/dz]

	lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, dw, 
						_p.x, _p.y, _p.z, _p.s, dw, _p.neighbor);
	std::vector<double> _ddwdx = lsmpsa.get_ddx();	// [d(dw)/dx]
	std::vector<double> _ddwdy = lsmpsa.get_ddy();	// [d(dw)/dy]


	
	// Evaluating added penalization vorticity field [dw = grad (cross) du]
	// #pragma omp parallel for
	for (int i = 0; i < _p.num; i++){
		// Calculate the vorticity
		double _dvortx = _ddwdy[i] - _ddvdz[i];         // w = del x u
		double _dvorty = _ddudz[i] - _ddwdx[i];         // w = del x u
		double _dvortz = _ddvdx[i] - _ddudy[i];         // w = del x u

		// Update the vorticity
		_p.vortx[i] += _dvortx;
		_p.vorty[i] += _dvorty;
		_p.vortz[i] += _dvortz;
	}


	// PROCEDURE 4: Update the penalization result into the original particle
    // ********************************************************************
    // Update the penalization value to original particle variable
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++){
        // Aliasing the original ID
        const int &ori_ID = this->baseID[i];

        // Assign each new data to original data
        _basePar.u[ori_ID]   = _p.u[i];
        _basePar.v[ori_ID]   = _p.v[i];
		_basePar.w[ori_ID]   = _p.w[i];
        _basePar.vortx[ori_ID]  = _p.vortx[i];
		_basePar.vorty[ori_ID]  = _p.vorty[i];
		_basePar.vortz[ori_ID]  = _p.vortz[i];
    }

	// Update the vorticity
	_basePar.vorticity.resize(_basePar.num,0.0);
	#pragma omp parallel for
    for (int i = 0; i < _p.num; i++){
		_basePar.vorticity[i] = std::sqrt(
			_basePar.vortx[i]*_basePar.vortx[i] +
			_basePar.vorty[i]*_basePar.vorty[i] +
			_basePar.vortz[i]*_basePar.vortz[i]
		);
	}
	
	return;
}

/**
 *  @brief Iterative penalization calculation. Update the velocity and vorticity 
 *  based on penalization calculation.
 *  
 *  @param	_tempParticle  The temporary particle data for calculation.
 *  @param	_baseParticle  The original particle data in calculation.
 *  @param  _bodyList  The body data list.
 *  @param  _step  The current simulation iteration step.
 */
void penalization::no_slip_iterative(Particle &_p, 
									 Particle &_basePar, 
									 const std::vector<Body> &_bodyList, 
									 int step) const
{
	/* Procedure Iteration:
        1. Apply penalization calculation
        2. Calculate velocity after penalization
        3. Apply the temporary data
	*/

	// Create the velocity calculation manager
	VelocityCalc velocity_tool;

	// Calculate the single no slip boudanry condition penalization
	for (int _it = 0; _it < Pars::pen_iter - 1; _it++){
		// // <!> Procedure 1 [Move to the bottom of the loop]
		// // Calculate penalization
		// if (DIM == 2){
		// 	this->no_slip(_p, _basePar, _bodyList, step);
		// }
        // else if (DIM == 3){
		// 	this->no_slip_3d(_p, _basePar, _bodyList, step);
		// }

		// <!> Procedure 1.2
		// Update the vorticity inside the body to zero (cancel the vorticity residual by error LSMPS calculation)
		#if (DIM == 2)
			#pragma omp parallel for
			for(int ID = 0; ID < _basePar.num; ID++){
				if (_basePar.insideBody[ID] == true){
					// Set the vorticity inside body to zero (0) to counter truncation error
					_basePar.vorticity[ID] = 0;
					_basePar.gz[ID] = 0;
				}
			}
		#elif (DIM == 3)
			#pragma omp parallel for
			for(int ID = 0; ID < _basePar.num; ID++){
				if (_basePar.insideBody[ID] == true){
					// Set the vorticity inside body to zero (0) to counter truncation error
					_basePar.vortx[ID] = 0;
					_basePar.vorty[ID] = 0;
					_basePar.vortz[ID] = 0;
				}
			}
		#endif


		// <!> Procedure 2
		// Calculate velocity
		velocity_tool.get_velocity(_basePar, step);


		// <!> Procedure 3
		// Update the vorticity and velocity data into temporary particle (_p)
		#if (DIM == 2)
			#pragma omp parallel for
			for (int i = 0; i < _p.num; i++){
				// Aliasing the original ID
				const int &ori_ID = this->baseID[i];

				// Assign each new data to original data
				_p.u[i]  = _basePar.u[ori_ID];
				_p.v[i]  = _basePar.v[ori_ID];
				_p.gz[i] = _basePar.gz[ori_ID];
				_p.vorticity[i] = _basePar.vorticity[ori_ID];
			}
		#elif (DIM == 3)
			#pragma omp parallel for
			for (int i = 0; i < _p.num; i++){
				// Aliasing the original ID
				const int &ori_ID = this->baseID[i];

				// Assign each new data to original data
				_p.u[i] = _basePar.u[ori_ID];
				_p.v[i] = _basePar.v[ori_ID];
				_p.w[i] = _basePar.w[ori_ID];
				_p.vortx[i] = _basePar.vortx[ori_ID];
				_p.vorty[i] = _basePar.vorty[ori_ID];
				_p.vortz[i] = _basePar.vortz[ori_ID];
			}
		#endif

		// <!> Procedure 1
		// Calculate penalization
		if (DIM == 2){
			this->no_slip(_p, _basePar, _bodyList, step);
		}
        else if (DIM == 3){
			this->no_slip_3d(_p, _basePar, _bodyList, step);
		}
	}
	
	return;
}