#include "stretching.hpp"

#ifndef INCLUDED_LSMPSa
#include "../LSMPS/LSMPSa.hpp"
#endif

/**
 *  @brief A particle stretching calculation.
 * 
 *  @param _particle The particle data for stretching calculation.
 *
 *  @headerfile stretching.hpp
 */
void stretching::main_stretching(Particle &p)
{
    // Internal method
    LSMPSa lsmpsa;		          // To calculate velocity differential

    // Internal variables
    Particle _particle;           // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;      // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _eval;      // The flag list of evaluated particle
    _eval = std::vector<bool>(p.num, false);    // At initial still no evaluated particle list
    
    // Stretching computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
    
    // Stretching prompt
    printf("\nCalculating stretching ...\n");
    
    
    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data

    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _eval[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _eval[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_eval[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _eval[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _particle
    _particle.num = _index.size();
    _particle.x.resize(_particle.num);
    _particle.y.resize(_particle.num);
    _particle.z.resize(_particle.num);
    _particle.u.resize(_particle.num);
    _particle.v.resize(_particle.num);
    _particle.w.resize(_particle.num);
    _particle.s.resize(_particle.num);
    _particle.vortx.resize(_particle.num);
    _particle.vorty.resize(_particle.num);
    _particle.vortz.resize(_particle.num);
    _particle.neighbor.resize(_particle.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++)
    {
        const int &_ID = _index[i];
        _particle.x[i] = (p.x[_ID]);
        _particle.y[i] = (p.y[_ID]);
        _particle.z[i] = (p.z[_ID]);
        _particle.u[i] = (p.u[_ID]);
        _particle.v[i] = (p.v[_ID]);
        _particle.w[i] = (p.w[_ID]);
        _particle.s[i] = (p.s[_ID]);
        _particle.vortx[i] = (p.vortx[_ID]);
        _particle.vorty[i] = (p.vorty[_ID]);
        _particle.vortz[i] = (p.vortz[_ID]);
        _particle.neighbor[i] = (p.neighbor[_ID]);
    }

    
    // PROCEDURE 2:
    // ************
    // Calculating the stretching

    // Calculate the first order differential of velocity
    Particle &_p = _particle;       // Aliasing
        
    // Calculate the velocity x (u) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.u,
                        p.x, p.y, p.z, p.s, p.u, _p.neighbor);
    std::vector<double> _dudx = lsmpsa.get_ddx();
    std::vector<double> _dudy = lsmpsa.get_ddy();
    std::vector<double> _dudz = lsmpsa.get_ddz();

    // Calculate the velocity y (v) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.v,
                        p.x, p.y, p.z, p.s, p.v, _p.neighbor);
    std::vector<double> _dvdx = lsmpsa.get_ddx();
    std::vector<double> _dvdy = lsmpsa.get_ddy();
    std::vector<double> _dvdz = lsmpsa.get_ddz();

    // Calculate the velocity z (w) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.w,
                        p.x, p.y, p.z, p.s, p.w, _p.neighbor);
    std::vector<double> _dwdx = lsmpsa.get_ddx();
    std::vector<double> _dwdy = lsmpsa.get_ddy();
    std::vector<double> _dwdz = lsmpsa.get_ddz();

    // Calculate the vortex stretching (dω/dt = ω.∇u) := omega(dot)(nabla)u
    #pragma omp parallel for
    for (int i = 0; i < _particle.num; i++){
        const int &_ID = _index[i];
        p.vortx[_ID] += Pars::dt * ((p.vortx[_ID] * _dudx[i]) + (p.vorty[_ID] * _dudy[i]) + (p.vortz[_ID] * _dudz[i]));
        p.vorty[_ID] += Pars::dt * ((p.vortx[_ID] * _dvdx[i]) + (p.vorty[_ID] * _dvdy[i]) + (p.vortz[_ID] * _dvdz[i]));
        p.vortz[_ID] += Pars::dt * ((p.vortx[_ID] * _dwdx[i]) + (p.vorty[_ID] * _dwdy[i]) + (p.vortz[_ID] * _dwdz[i]));
    }

    // *Note: At this point dont need to calculate the vorticity absolute value,
    //         because it will be calculated in the next section (particle redistribution)

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Stretching total computation time: [%f s]\n", _time);
    
    return;
}


/**
 *  @brief Combination of diffusion and stretching only for 3D simulation.
 *  NOTE: Supposed to be more efficient, than calculating each separately.
 * 
 *  @param _particle The particle data for diffusion and stretching calculation.
 *
 *  @headerfile stretching.hpp
 */
void stretching::calc_diff_stretch(Particle &p){
    // Violation of program
    if (DIM != 3){
        ERROR_LOG << "Abuse of program, can't do stretching for 2D simulation!";
        throw std::runtime_error("[ERROR] An attempt to calculate stretching for 2D simulation");
    }
    
    // Internal method
    LSMPSa lsmpsa;		          // To calculate velocity differential

    // Internal variables
    Particle _p;           // The data of particle inside 'active particle + buffer zone'
    std::vector<int> _index;      // The index list of particle inside 'active particle + buffer zone'
    std::vector<bool> _eval;      // The flag list of evaluated particle
    _eval = std::vector<bool>(p.num, false);    // At initial still no evaluated particle list
    
    // Diffusion and stretching computational time manager
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::_V2::system_clock::time_point tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        double _time = omp_get_wtime();
    #endif
    
    
    // PROCEDURE 1:
    // ************
    // Collecting the only active particle data
    printf("\nCollect data for diffusion and stretching ...\n");

    // Assign the index of active particles
    for (int ID = 0; ID < p.num; ID++)
    {
        // Check if the particle ID is active
        if (p.isActive[ID] == true)
        {
            // [1] The current ID particle
            if(!Pars::flag_ngh_include_self && _eval[ID] == false)
            {
                // Assign ID into index list if not in the list (if flag == false)
                _index.push_back(ID);    // Assign ID into the index list
                _eval[ID] = true;        // Update the evaluation flag of particle ID
            }

            // [2] The neighbor particles of current ID particle
            for (auto _ngh_ID : p.neighbor[ID])
            {
                if(_eval[_ngh_ID] == false)
                {
                    // Assign into index list if not in the list (if flag == false)
                    _index.push_back(_ngh_ID);  // Assign _ngh_ID particle into the index list
                    _eval[_ngh_ID] = true;      // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Resize the _p
    _p.num = _index.size();
    _p.x.resize(_p.num);
    _p.y.resize(_p.num);
    _p.z.resize(_p.num);
    _p.u.resize(_p.num);
    _p.v.resize(_p.num);
    _p.w.resize(_p.num);
    _p.s.resize(_p.num);
    _p.vortx.resize(_p.num);
    _p.vorty.resize(_p.num);
    _p.vortz.resize(_p.num);
    _p.neighbor.resize(_p.num);

    // Store the particle data of each particle inside _index list
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++)
    {
        const int &_ID = _index[i];
        _p.x[i] = (p.x[_ID]);
        _p.y[i] = (p.y[_ID]);
        _p.z[i] = (p.z[_ID]);
        _p.u[i] = (p.u[_ID]);
        _p.v[i] = (p.v[_ID]);
        _p.w[i] = (p.w[_ID]);
        _p.s[i] = (p.s[_ID]);
        _p.vortx[i] = (p.vortx[_ID]);
        _p.vorty[i] = (p.vorty[_ID]);
        _p.vortz[i] = (p.vortz[_ID]);
        _p.neighbor[i] = (p.neighbor[_ID]);
    }
    
    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        std::chrono::duration<double> span = std::chrono::system_clock::now() - tick;
        double _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Collect data computing time:       [%f s]\n", _time);

    

    // PROCEDURE 2:
    // ************
    // Calculating the diffusion
    printf("\nCalculating diffusion ...\n");
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif
    
    // Calculate the second order differential of vorticity
    // Calculate the x directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortx, 
                        p.x, p.y, p.z, p.s, p.vortx, _p.neighbor);
    std::vector<double> _d2vortxd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortxd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortxd2z = lsmpsa.get_d2d2z();

    // Calculate the y directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vorty, 
                        p.x, p.y, p.z, p.s, p.vorty, _p.neighbor);
    std::vector<double> _d2vortyd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortyd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortyd2z = lsmpsa.get_d2d2z();

    // Calculate the z directing vorticity
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.vortz, 
                        p.x, p.y, p.z, p.s, p.vortz, _p.neighbor);
    std::vector<double> _d2vortzd2x = lsmpsa.get_d2d2x();
    std::vector<double> _d2vortzd2y = lsmpsa.get_d2d2y();
    std::vector<double> _d2vortzd2z = lsmpsa.get_d2d2z();
    
    // Calculate the vorticity laplacian to get the vorticity time differential
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++)
    {
        const int &_ID = _index[i];     // ID alliasing

        p.vortx[_ID] += Pars::dt * Pars::NU * ((_d2vortxd2x[i] + _d2vortxd2y[i] + _d2vortxd2z[i]));
        p.vorty[_ID] += Pars::dt * Pars::NU * ((_d2vortyd2x[i] + _d2vortyd2y[i] + _d2vortyd2z[i]));
        p.vortz[_ID] += Pars::dt * Pars::NU * ((_d2vortzd2x[i] + _d2vortzd2y[i] + _d2vortzd2z[i]));
    }

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Diffusion total computation time:  [%f s]\n", _time);



    // PROCEDURE 3:
    // ************
    // Calculating the stretching
    printf("\nCalculating stretching ...\n");
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        tick = std::chrono::system_clock::now();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime();
    #endif

    // Calculate the first order differential of velocity
    // Calculate the velocity x (u) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.u,
                        p.x, p.y, p.z, p.s, p.u, _p.neighbor);
    std::vector<double> _dudx = lsmpsa.get_ddx();
    std::vector<double> _dudy = lsmpsa.get_ddy();
    std::vector<double> _dudz = lsmpsa.get_ddz();

    // Calculate the velocity y (v) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.v,
                        p.x, p.y, p.z, p.s, p.v, _p.neighbor);
    std::vector<double> _dvdx = lsmpsa.get_ddx();
    std::vector<double> _dvdy = lsmpsa.get_ddy();
    std::vector<double> _dvdz = lsmpsa.get_ddz();

    // Calculate the velocity z (w) differential
    lsmpsa.set_LSMPS_3D(_p.x, _p.y, _p.z, _p.s, _p.w,
                        p.x, p.y, p.z, p.s, p.w, _p.neighbor);
    std::vector<double> _dwdx = lsmpsa.get_ddx();
    std::vector<double> _dwdy = lsmpsa.get_ddy();
    std::vector<double> _dwdz = lsmpsa.get_ddz();

    // Calculate the vortex stretching
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++)
    {
        const int &_ID = _index[i];     // ID alliasing

        p.vortx[_ID] += Pars::dt * ((p.vortx[_ID] * _dudx[i]) + (p.vorty[_ID] * _dudy[i]) + (p.vortz[_ID] * _dudz[i]));
        p.vorty[_ID] += Pars::dt * ((p.vortx[_ID] * _dvdx[i]) + (p.vorty[_ID] * _dvdy[i]) + (p.vortz[_ID] * _dvdz[i]));
        p.vortz[_ID] += Pars::dt * ((p.vortx[_ID] * _dwdx[i]) + (p.vorty[_ID] * _dwdy[i]) + (p.vortz[_ID] * _dwdz[i]));
    }

    // *Note: At this point dont need to calculate the vorticity absolute value,
    //         because it will be calculated in the next section (particle redistribution)

    // Display computational time
    #if (TIMER_PAR == 0)
        // Timer using super clock (chrono)
        span = std::chrono::system_clock::now() - tick;
        _time = span.count();
    #elif (TIMER_PAR == 1)
        // Timer using paralel package
        _time = omp_get_wtime() - _time;
    #endif
    printf("<-> Stretching total computation time: [%f s]\n", _time);
    
    return;

}