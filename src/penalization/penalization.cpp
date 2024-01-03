#include "penalization.hpp"
#include <unordered_map>

/**
 *  @brief Penalization calculation manager. Update the velocity and vorticity 
 *  based on penalization calculation.
 *  NOTE: Works only for 2D simulation.
 *  
 *  @param	_particle  The particle data for calculation.
 *  @param  _bodyList  The body data list.
 *  @param  _step  The current simulation iteration step.
 */
void penalization::get_penalization(Particle &_currPar, const std::vector<Body> &bL, int step)
{
    /* Procedure:
       1. Collect all particle near the body (offset by Pars::numpen)
       2. Calculate the chi parameter
       3. Perform the vorticity penalization
       4. **Update the penalization result into the original particle [Put inside the no_slip function]
    */
    
    // Penalization time manager
    double _time = omp_get_wtime();
    printf("\nCalculating penalization ... \n");

    // PROCEDURE 1: Collect near body particle
    // ********************************************************************
    // Internal variables
    std::unordered_map<int,int> tempIdx;  	// Index convertion from original to new index (for obtaining neighbor list)
    Particle _tempPar;          // Temporary particle for evaluation (Near body particle only)
    _tempPar.num = 0;           // Initial particle number value
    

    // Take only the particle inside the penalization evaluation domain
    this->baseID.clear();       // Reserve the data first

    // Take the index converter
    for (int ID = 0; ID < _currPar.num; ID++){
        if (_currPar.bodyPart[ID] != -1){
            // Put the ID conversion between original and temporary
            this->baseID.push_back(ID);
            tempIdx.insert({ID, _tempPar.num});
            // Update the number
            _tempPar.num++;
        }
    }
    
    // Update the other data
    if (DIM == 2){
        _tempPar.x.resize(_tempPar.num);
        _tempPar.y.resize(_tempPar.num);
        _tempPar.u.resize(_tempPar.num);
        _tempPar.v.resize(_tempPar.num);
        _tempPar.s.resize(_tempPar.num);
        _tempPar.gz.resize(_tempPar.num);
        _tempPar.vorticity.resize(_tempPar.num);
        _tempPar.bodyPart.resize(_tempPar.num);
        _tempPar.chi.resize(_tempPar.num);

        #pragma omp parallel for
        for (int i = 0; i < _tempPar.num; i++){
            // Aliasing to the original index
            int &ori_ID = this->baseID[i];

            // Put the 2D data
            _tempPar.x[i] = _currPar.x[ori_ID];
            _tempPar.y[i] = _currPar.y[ori_ID];
            _tempPar.u[i] = _currPar.u[ori_ID];
            _tempPar.v[i] = _currPar.v[ori_ID];
            _tempPar.s[i] = _currPar.s[ori_ID];
            _tempPar.gz[i] = _currPar.gz[ori_ID];
            _tempPar.vorticity[i] = _currPar.vorticity[ori_ID];
            _tempPar.bodyPart[i] = _currPar.bodyPart[ori_ID];
            _tempPar.chi[i] = _currPar.chi[ori_ID];
        }
    }
    else if (DIM == 3){
        _tempPar.x.resize(_tempPar.num);
        _tempPar.y.resize(_tempPar.num);
        _tempPar.z.resize(_tempPar.num);
        _tempPar.u.resize(_tempPar.num);
        _tempPar.v.resize(_tempPar.num);
        _tempPar.w.resize(_tempPar.num);
        _tempPar.s.resize(_tempPar.num);
        _tempPar.vortx.resize(_tempPar.num);
        _tempPar.vorty.resize(_tempPar.num);
        _tempPar.vortz.resize(_tempPar.num);
        _tempPar.bodyPart.resize(_tempPar.num);
        _tempPar.chi.resize(_tempPar.num);

        #pragma omp parallel for
        for (int i = 0; i < _tempPar.num; i++){
            // Aliasing to the original index
            int &ori_ID = this->baseID[i];

            // Put the 3D data
            _tempPar.x[i] = _currPar.x[ori_ID];
            _tempPar.y[i] = _currPar.y[ori_ID];
            _tempPar.z[i] = _currPar.z[ori_ID];
            _tempPar.u[i] = _currPar.u[ori_ID];
            _tempPar.v[i] = _currPar.v[ori_ID];
            _tempPar.w[i] = _currPar.w[ori_ID];
            _tempPar.s[i] = _currPar.s[ori_ID];
            _tempPar.vortx[i] = _currPar.vortx[ori_ID];
            _tempPar.vorty[i] = _currPar.vorty[ori_ID];
            _tempPar.vortz[i] = _currPar.vortz[ori_ID];
            _tempPar.bodyPart[i] = _currPar.bodyPart[ori_ID];
            _tempPar.chi[i] = _currPar.chi[ori_ID];
        }
    }


    // Get the local neighbor ID
    _tempPar.neighbor.resize(_tempPar.num);         // Take from original neighbor list, transform the index to new index
    #pragma omp parallel for
    for (int i = 0; i < _tempPar.num; i++){
        // Aliasing to the original index
        int &original_ID = this->baseID[i];

        // Iterate through the neighbor 
        for (const auto& oldNghID : _currPar.neighbor[original_ID]){
            // Check whether the current particle is included to temporary particle data
            if (tempIdx.count(oldNghID) != 0){
                _tempPar.neighbor[i].push_back(tempIdx.at(oldNghID));
            }
        }
    }
    

    // PROCEDURE 2: Set up initial penalization parameter
    // ********************************************************************
    // Calculate the chi (Only if the particle is adapted) [*] Already calculated in redistribution step
    // if (_currPar.isAdapt){
    //     // Calculate the chi data
    //     this->calculate_chi(_tempPar, bL);

    //     // Assign the penalization mark data
    //     _currPar.chi.clear(); _currPar.chi.resize(_currPar.num,0.0);
    //     #pragma omp parallel for
    //     for (int i = 0; i < _tempPar.num; i++){
    //         _currPar.chi[this->baseID[i]] = _tempPar.chi[i];       
    //     }
    // }

    // PROCEDURE 3: Perform the vorticity penalization
    // ********************************************************************
    // Perform the penalization
    if (Pars::opt_pen_iter == 1){
        if (DIM == 2)
        this->no_slip(_tempPar, _currPar, bL, step);
        else if (DIM == 3)
        this->no_slip_3d(_tempPar, _currPar, bL, step);
    }
    else if (Pars::opt_pen_iter == 2){
        this->no_slip_iterative(_tempPar, _currPar, bL, step);
    }

    // PROCEDURE 4: [ADDITIONAL] Set vorticity inside body to zero
    // ********************************************************************
    // Update the vorticity inside the body to zero (cancel the vorticity residual by error LSMPS calculation)
    if (DIM == 2){
        #pragma omp parallel for
        for(int ID = 0; ID < _currPar.num; ID++){
            if (_currPar.insideBody[ID] == true){
                // Set the vorticity inside body to zero (0) to counter truncation error
                _currPar.vorticity[ID] = 0;
                _currPar.gz[ID] = 0;
            }
        }
    }
    else if (DIM == 3){
        #pragma omp parallel for
        for(int ID = 0; ID < _currPar.num; ID++){
            if (_currPar.insideBody[ID] == true){
                // Set the vorticity inside body to zero (0) to counter truncation error
                _currPar.vortx[ID] = 0;
                _currPar.vorty[ID] = 0;
                _currPar.vortz[ID] = 0;
            }
        }
    }

    // Display penalization computational time
    _time = omp_get_wtime() - _time;
    printf("<-> Penalization computational time:   [%f s]\n", _time);
}
