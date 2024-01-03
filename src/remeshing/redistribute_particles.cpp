#include "remeshing.hpp"
#include "../LSMPS/LSMPSb.hpp"

// A debugger to check the neighbor evaluation after adaptation and redistribution
void save_ngh_all(const Particle &_myPar, int ID);      // For DEBUGGING
static int value = 0;                                   // For DEBUGGING

/**
 *  @brief  Redistribution all particle using LSMPS B interpolation on vorticity.
 *  Interpolate the vorticuty from most updated particle into the base particle.
 *  Assign the interpolated data into the particle.
 *         
 *  @param  _particle  Particle most updated data container.
*/
void remeshing::redistribute_particles(Particle &parEval)
{   
    // ============================================= //
	// ========== Particle Redistribution ========== //
	// ============================================= //
    /** To Do:
     *  > Interpolate particle [vorticity] to base particle\
     *  > Update the other particle properties
     *     - Coordinate, num, : Take the value from base particle value
     *        level, ngh, node relation
     *     - Particle size    : Calculate from level
     *     - Velocity, P, chi : Resize with 0.0
     *     - Vortex strength  : Calculate from interpolated vorticity
     *     - Active Flag      : Re-evaluate from vorticity value
     *  > In this set up -> Update the base Particle -> assign the updated data to current particle
    */

    // Set up the neighbor for data interpolation
    std::vector<std::vector<int>> sourceNghList; // The neighbor of base particle using the ID of new particle (inter neighbor)
    std::vector<std::vector<int>> selfNghList;   // The neighbor of base particle relative to itself (self neihgbor)

    // Create an alias for each data
    std::vector<double> &_tempX = this->particleBase->x;
    std::vector<double> &_tempY = this->particleBase->y;
    std::vector<double> &_tempZ = this->particleBase->z;

    std::vector<double> _tempSize = this->particleBase->s;      // Not a reference, because is altered when adaptation flag is on

    // For 2D simulation
    std::vector<double> &_tempVor = this->particleBase->vorticity;
    
    // For 3D simulation
    std::vector<double> &_tempVortx = this->particleBase->vortx;
    std::vector<double> &_tempVorty = this->particleBase->vorty;
    std::vector<double> &_tempVortz = this->particleBase->vortz;

    // NOTE:
    //  > Both [self neighbor] and [inter neighbor] hold the same value if not undergo adaptation
    //  > If adaptation is done "Base Particle" hold the [inter neighbor] data, while 
    //     [self neighbor] need to be calculated next
    //  > The [inter neighbor] is different if adapted and not adapted
    //     If adapted it already obtain all

    // Computational time manager
    double _time;


    // PROCEDURE 1: Neighbor Alliasing
    // ************
    // **Evaluate the _par2grd neighbor list [inter neighbor]
    sourceNghList = this->particleBase->neighbor;
    // Add the particle ID itself into the neighborID
    if ((this->adap_flag == false) && (Pars::flag_ngh_include_self == false)){
        // ERROR_LOG << "GETTING INTO FIRST PAGE\n";
        for (int i = 0; i < this->particleBase->num; i++){
            sourceNghList[i].push_back(i);
        }
    }
    

    // **Evaluate the _grd2grd neighbor list  [self neighbor]
    if (this->adap_flag == true){
        // ERROR_LOG << "GETTING INTO SECOND PAGE\n";
        // Update the particle size
        for (int i = 0; i < this->particleBase->num; i++)
        this->particleBase->s[i] = this->tempGrid->gridSize/(Pars::intPow(2,this->particleBase->level[i]) * this->tempGrid->baseParNum);
        
        // Calculate the neighbor (At the same time update the neighbor for next data)
        neighbor_step.neigbor_search(*(this->particleBase), *(this->tempGrid));
    }

    // Assign the current stored neighbor on Base Particle
    selfNghList = this->particleBase->neighbor;
    // Add the particle ID itself into the neighborID
    if (!Pars::flag_ngh_include_self){
        // ERROR_LOG << "GETTING INTO THIRD PAGE\n";
        for (int i = 0; i < this->particleBase->num; i++){
            selfNghList[i].push_back(i);
        }
    }


    // // [DEBUG]
    // ERROR_LOG << "CHECK BEFOFE INTERPOLATED\n";
    // printf(" > Size of target particle (x,y,s,f) = (%ld,%ld,%ld, %ld)\n",_tempX.size(), _tempY.size(), _tempSize.size(), _tempVor.size());
    // printf(" > Size of source particle (x,y,s,f) = (%ld,%ld,%ld, %ld)\n",parEval.x.size(), parEval.y.size(), parEval.s.size(), parEval.vorticity.size());
    // printf(" > Size of neighbor source = (%ld)\n",sourceNghList.size());
    // printf(" > Size of neighbor self = (%ld)\n",selfNghList.size());

    // int NN = 20;
    // int timer[NN];
    // for (int i = 0; i < NN; i++){
    //     timer[i] = (((((i+123) * 13) %143) + ((i+97) * 17))) * i;
    // }

    // WARNING_LOG << "Data Number : " << this->particleBase->num << "\n";
    // clock_t time = clock();
    // for (const auto &mul : timer){
    //     int ID = (static_cast<int>(mul * (static_cast<double>(clock()-time)))) % this->particleBase->num;
    //     // MESSAGE_LOG << "List of random ID " << ID << "\n";
    //     save_ngh_all(*this->particleBase,std::abs(ID));
    // }
    // value++;

    // Problem:
    //  1. Neighbor missmatch       [DONE wrong definition at adaptation step]

    

    // PROCEDURE 2: Interpolation of Vorticity
    // ************
    _time = omp_get_wtime();

    // Perform LSMPS interpolation of vorticity
    LSMPSb lsmpsb;
    if (DIM == 2){
        lsmpsb.set_LSMPS(_tempX, _tempY, _tempSize, _tempVor,
                         parEval.x, parEval.y, parEval.s, parEval.vorticity,
                         sourceNghList, selfNghList);
        // Update the particle base data vorticity and vortex strength
        this->particleBase->vorticity = lsmpsb.get_d0();
    }
    else if (DIM == 3){
        // Interpolate the vorticity in x direction
        lsmpsb.set_LSMPS_3D(_tempX, _tempY, _tempZ, _tempSize, _tempVortx,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vortx,
                            sourceNghList, selfNghList);
        this->particleBase->vortx = lsmpsb.get_d0();

        // Interpolate the vorticity in y direction
        lsmpsb.set_LSMPS_3D(_tempX, _tempY, _tempZ, _tempSize, _tempVorty,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vorty,
                            sourceNghList, selfNghList);
        this->particleBase->vorty = lsmpsb.get_d0();

        // Interpolate the vorticity in y direction
        lsmpsb.set_LSMPS_3D(_tempX, _tempY, _tempZ, _tempSize, _tempVortz,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vortz,
                            sourceNghList, selfNghList);
        this->particleBase->vortz = lsmpsb.get_d0();
    }
    
    // Time manager console log
    _time = omp_get_wtime() - _time;
    printf("<-> Calculating LSMPS interpolation:   [%f s]\n", _time);
    
    
    if (DIM == 2){
        // Update the vortex strength <?> Should be changed next
        this->particleBase->gz.clear();this->particleBase->gz.resize(this->particleBase->num,0.0);
        #pragma omp parallel for
        for (int _i = 0; _i < this->particleBase->num; _i++){
            this->particleBase->gz[_i] = this->particleBase->vorticity[_i] * std::pow(this->particleBase->s[_i],2.0);
        }
    }
    else if (DIM == 3){
        // Update the vorticity
        // Create base particle alias
        Particle &_p = *this->particleBase;
        
        // Reset the size of vorticity
        _p.vorticity.clear();_p.vorticity.resize(_p.num,0.0);
        
        // Calculate the vorticity
        #pragma omp parallel for
        for (int _i = 0; _i < _p.num; _i++){
            _p.vorticity[_i] = std::sqrt(_p.vortx[_i]*_p.vortx[_i] + 
                                         _p.vorty[_i]*_p.vorty[_i] + 
                                         _p.vortz[_i]*_p.vortz[_i]);
        }
    }

    
    // PROCEDURE 3: Update Active Particle
    // ************    
    _time = omp_get_wtime();

    // Update the redistribute data into the particle variable
    this->set_particle_sign();

    // Update the redistribute data into the particle variable
    parEval = *(this->particleBase);
    parEval.isAdapt = this->adap_flag;

    _time = omp_get_wtime() - _time;
    printf("<-> Re-assigning the particle data:    [%f s]\n", _time);


    // // [DEBUG]
    // time = clock();
    // WARNING_LOG << "Data Number : " << this->particleBase->num << "\n";
    // for (const auto &mul : timer){
    //     int ID = (static_cast<int>(mul * (static_cast<double>(clock()-time)))) % this->particleBase->num;
    //     // MESSAGE_LOG << "List of random ID " << ID << "\n";
    //     save_ngh_all(*this->particleBase,std::abs(ID));
    // }
    // value++;
}


/**
 *  @brief  Redistribution active particle only using LSMPS B interpolation on vorticity.
 *  Interpolate the vorticuty from most updated particle into the base particle.
 *  Assign the interpolated data into the particle.
 *         
 *  @param  _particle  Particle most updated data container.
*/
void remeshing::redistribute_active_particles(Particle &parEval)
{   
    // ============================================= //
	// ========== Particle Redistribution ========== //
	// ============================================= //
    /** To Do:
     *  > Interpolate particle [vorticity] to base particle\
     *  > Update the other particle properties
     *     - Coordinate, num, : Take the value from base particle value
     *        level, ngh, node relation
     *     - Particle size    : Calculate from level
     *     - Velocity, P, chi : Resize with 0.0
     *     - Vortex strength  : Calculate from interpolated vorticity
     *     - Active Flag      : Re-evaluate from vorticity value
     *  > In this set up -> Update the base Particle -> assign the updated data to current particle
    */

    // Set up the neighbor for data interpolation
    std::vector<std::vector<int>> sourceNghList; // The neighbor of base particle using the ID of new particle (inter neighbor)
    std::vector<std::vector<int>> selfNghList;   // The neighbor of base particle relative to itself (self neihgbor)

    // Create a duplicate of size
    std::vector<double> _tempSize = this->particleBase->s;      // Not a reference, because is altered when adaptation flag is on

    // NOTE:
    //  > Both [self neighbor] and [inter neighbor] hold the same value if not undergo adaptation
    //  > If adaptation is done "Base Particle" hold the [inter neighbor] data, while 
    //     [self neighbor] need to be calculated next
    //  > The [inter neighbor] is different if adapted and not adapted
    //     If adapted it already obtain all

    // Computational time manager
    double _time;


    // PROCEDURE 1: Neighbor Alliasing
    // ************
    // **Evaluate the _par2grd neighbor list [inter neighbor]
    sourceNghList = this->particleBase->neighbor;
    // Add the particle ID itself into the neighborID
    if ((this->adap_flag == false) && (Pars::flag_ngh_include_self == false)){
        // ERROR_LOG << "GETTING INTO FIRST PAGE\n";
        for (int i = 0; i < this->particleBase->num; i++){
            sourceNghList[i].push_back(i);
        }
    }
    

    // **Evaluate the _grd2grd neighbor list  [self neighbor]
    if (this->adap_flag == true){
        // ERROR_LOG << "GETTING INTO SECOND PAGE\n";
        // Update the particle size
        for (int i = 0; i < this->particleBase->num; i++)
        this->particleBase->s[i] = this->tempGrid->gridSize/(Pars::intPow(2,this->particleBase->level[i]) * this->tempGrid->baseParNum);
        
        // Calculate the neighbor (At the same time update the neighbor for next data)
        neighbor_step.neigbor_search(*(this->particleBase), *(this->tempGrid));
    }

    // Assign the current stored neighbor on Base Particle
    selfNghList = this->particleBase->neighbor;
    // Add the particle ID itself into the neighborID
    if (!Pars::flag_ngh_include_self){
        // ERROR_LOG << "GETTING INTO THIRD PAGE\n";
        for (int i = 0; i < this->particleBase->num; i++){
            selfNghList[i].push_back(i);
        }
    }



    // PROCEDURE 2: Colecting active particle
    // ************
    _time = omp_get_wtime();

    std::vector<int> activeID;          // Point to the original particle ID
    std::vector<int> tempID(this->particleBase->num, -1);    // Point the ID to active index
    std::vector<bool> activeCheckFlag(this->particleBase->num, false);   // A marker for already checked particle

    // Assign the index of active particle only
    for (int ID = 0; ID < this->particleBase->num; ID++)
    {
        // Check if the particle ID is active
        if (this->particleBase->isActive[ID] == true)       // Need further check **
        {
            // [1] Put the current ID particle to the active ID list
            if(!Pars::flag_ngh_include_self && activeCheckFlag[ID] == false){
                // Assign ID into index list if not in the list (if flag == false)
                activeID.push_back(ID);         // Assign ID into the index list
                activeCheckFlag[ID] = true;     // Update the evaluation flag of particle ID
            }

            // [2] Put the neighbor particles of current ID particle
            for (auto _ngh_ID : this->particleBase->neighbor[ID]){
                if(activeCheckFlag[_ngh_ID] == false){
                    // Assign into index list if not in the list (if flag == false)
                    activeID.push_back(_ngh_ID);        // Assign _ngh_ID particle into the index list
                    activeCheckFlag[_ngh_ID] = true;    // Update the evaluation flag of particle _ngh_ID
                }
            }
        }
    }

    // Create the vector array of interpolated particle
    int actNum = activeID.size();       // Active particle number
    std::vector<double> _actX, _actY, _actZ, _actSize;

    // New neighbor data
    std::vector<std::vector<int>> actSrcNghList(actNum);
    std::vector<std::vector<int>> actSelfNghList(actNum);

    // For 2D simulation
    std::vector<double> _actVor;
    
    // For 3D simulation
    std::vector<double> _actVortx, _actVorty, _actVortz;

    // Update the temporary ID
    #pragma omp parallel for
    for (int ID = 0; ID < actNum; ID++){
        // Alias to the original ID
        const int &_oriID = activeID[ID];
        // Put the temporary ID
        tempID[_oriID] = ID;
    }

    // Store the particle data of each particle inside activeID list
    // Resize the container
    _actSize.resize(actNum, 0.0);
    _actX.resize(actNum, 0.0);
    _actY.resize(actNum, 0.0);
    if (DIM == 3) _actZ.resize(actNum, 0.0);
    
    #pragma omp parallel for
    for (int i = 0; i < actNum; i++){
        // Alias to the original ID
        const int &_ID = activeID[i];
        
        // Store the coordinate data
        _actSize[i] = _tempSize[_ID];
        _actX[i] = this->particleBase->x[_ID];
        _actY[i] = this->particleBase->y[_ID];
        if (DIM == 3) _actZ[i] = this->particleBase->z[_ID];
        
        // Store the neighbor data
        actSrcNghList[i] = sourceNghList[_ID];      // Directly take from the temporary data
        actSelfNghList[i].clear();                  // Reserve the current data (redundant but for making sure)
        for (const auto &ori_ngh_ID : selfNghList[_ID]){
            // Alias to the temporary ID
            const int &_tmpID = tempID[ori_ngh_ID];
            // Skip if the corresponding ID is not existed in the active data
            if (_tmpID == -1) continue;
            actSelfNghList[i].push_back(_tmpID);
        }
    }

    // Time manager console log
    _time = omp_get_wtime() - _time;
    printf("<-> Collecting the active particle:    [%f s]\n", _time);


    // // [DEBUG]
    // ERROR_LOG << "CHECK BEFOFE INTERPOLATED\n";
    // printf(" > Size of target particle (x,y,s,f) = (%ld,%ld,%ld, %ld)\n",_tempX.size(), _tempY.size(), _tempSize.size(), _tempVor.size());
    // printf(" > Size of source particle (x,y,s,f) = (%ld,%ld,%ld, %ld)\n",parEval.x.size(), parEval.y.size(), parEval.s.size(), parEval.vorticity.size());
    // printf(" > Size of neighbor source = (%ld)\n",sourceNghList.size());
    // printf(" > Size of neighbor self = (%ld)\n",selfNghList.size());

    // int NN = 20;
    // int timer[NN];
    // for (int i = 0; i < NN; i++){
    //     timer[i] = (((((i+123) * 13) %143) + ((i+97) * 17))) * i;
    // }

    // WARNING_LOG << "Data Number : " << this->particleBase->num << "\n";
    // clock_t time = clock();
    // for (const auto &mul : timer){
    //     int ID = (static_cast<int>(mul * (static_cast<double>(clock()-time)))) % this->particleBase->num;
    //     // MESSAGE_LOG << "List of random ID " << ID << "\n";
    //     save_ngh_all(*this->particleBase,std::abs(ID));
    // }
    // value++;

    // Problem:
    //  1. Neighbor missmatch       [DONE wrong definition at adaptation step]

    

    // PROCEDURE 3: Interpolation of Vorticity
    // ************
    _time = omp_get_wtime();

    // Perform LSMPS interpolation of vorticity
    LSMPSb lsmpsb;
    if (DIM == 2){
        lsmpsb.set_LSMPS(_actX, _actY, _actSize, _actVor,
                         parEval.x, parEval.y, parEval.s, parEval.vorticity,
                         actSrcNghList, actSelfNghList);
        // Update the particle base data vorticity and vortex strength
        // this->particleBase->vorticity = lsmpsb.get_d0();
        _actVor = lsmpsb.get_d0();
    }
    else if (DIM == 3){
        // Interpolate the vorticity in x direction
        lsmpsb.set_LSMPS_3D(_actX, _actY, _actZ, _actSize, _actVortx,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vortx,
                            actSrcNghList, actSelfNghList);
        // this->particleBase->vortx = lsmpsb.get_d0();
        _actVortx = lsmpsb.get_d0();

        // Interpolate the vorticity in y direction
        lsmpsb.set_LSMPS_3D(_actX, _actY, _actZ, _actSize, _actVorty,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vorty,
                            actSrcNghList, actSelfNghList);
        // this->particleBase->vorty = lsmpsb.get_d0();
        _actVorty = lsmpsb.get_d0();

        // Interpolate the vorticity in y direction
        lsmpsb.set_LSMPS_3D(_actX, _actY, _actZ, _actSize, _actVortz,
                            parEval.x, parEval.y, parEval.z, parEval.s, parEval.vortz,
                            actSrcNghList, actSelfNghList);
        // this->particleBase->vortz = lsmpsb.get_d0();
        _actVortz = lsmpsb.get_d0();
    }
    
    // Time manager console log
    _time = omp_get_wtime() - _time;
    printf("<-> Calculating LSMPS interpolation:   [%f s]\n", _time);
    
    
    if (DIM == 2){
        // Update the vorticity and vortex strength <?> Should be changed next 
        // <!> no need: 2D use gz and vorticity; 3D use vortx, vorty, vortz, and vorticity
        Particle &_p = *this->particleBase; // Create base particle alias

        // Reset the size of vorticity and vortex strength
        _p.gz.clear();          _p.gz.resize(_p.num,0.0);
        _p.vorticity.clear();   _p.vorticity.resize(_p.num,0.0);

        // Update the vortex strength 
        #pragma omp parallel for
        for (int _i = 0; _i < actNum; _i++){
            // Alias to the original ID
            const int &_oriID = activeID[_i];
            
            // Alias to the vorticity
            const double &__vor = _actVor[_i];

            // Assign the vorticity data into the container
            _p.gz[_oriID] = __vor * std::pow(_p.s[_oriID],2.0);
            _p.vorticity[_oriID] = __vor;
        }
    }
    else if (DIM == 3){
        // Update the vorticity
        // Create base particle alias
        Particle &_p = *this->particleBase;
        
        // Reset the size of vorticity
        _p.vortx.clear();       _p.vortx.resize(_p.num,0.0);
        _p.vorty.clear();       _p.vorty.resize(_p.num,0.0);
        _p.vortz.clear();       _p.vortz.resize(_p.num,0.0);
        _p.vorticity.clear();   _p.vorticity.resize(_p.num,0.0);
        
        // Calculate the vorticity
        #pragma omp parallel for
        for (int _i = 0; _i < actNum; _i++){
            // Alias to the original ID
            const int &_oriID = activeID[_i];

            // Alias to the vorticity
            const double &__vortX = _actVortx[_i];
            const double &__vortY = _actVorty[_i];
            const double &__vortZ = _actVortz[_i];

            // Assign the vorticity data into the container
            _p.vortx[_oriID] = __vortX;
            _p.vorty[_oriID] = __vortY;
            _p.vortz[_oriID] = __vortZ;
            _p.vorticity[_oriID] = std::sqrt(__vortX*__vortX + __vortY*__vortY + __vortZ*__vortZ);
        }
    }

    
    // PROCEDURE 4: Update Active Particle
    // ************    
    _time = omp_get_wtime();

    // Update the redistribute data into the particle variable
    this->set_particle_sign();

    // Update the redistribute data into the particle variable
    parEval = *(this->particleBase);
    parEval.isAdapt = this->adap_flag;

    _time = omp_get_wtime() - _time;
    printf("<-> Re-assigning the particle data:    [%f s]\n", _time);
    return;
}


/**
 *  @brief  Update the particle active sign by considering the value of particle
 *  vorticity. Set the vorticity property to 0 for non active particle.
 *  NOTE: The updated particle is the base particle.
*/
void remeshing::set_particle_sign(){
    // **Find the maximum vorticity
    double vor_max = 0.0e0;
    for (int i = 0; i < this->particleBase->num; i++){
        if (std::abs(this->particleBase->vorticity[i]) > vor_max)
        vor_max = std::abs(this->particleBase->vorticity[i]);
    }

    // **Evaluate the active flag particle
    this->particleBase->isActive.clear();this->particleBase->isActive.resize(this->particleBase->num,false);
    double SIGNIFICANCE_LIMIT = Pars::active_sig*vor_max;
    #pragma omp parallel for
    for (int i = 0; i < this->particleBase->num; i++){
        double vor = std::abs(this->particleBase->vorticity[i]);
        this->particleBase->isActive[i] = (vor >= SIGNIFICANCE_LIMIT) ? true : false;
    }
    
    
    // **Set the vorticity to zero for non active particle
    if (DIM == 2){
        #pragma omp parallel for
        for (int i = 0; i < this->particleBase->num; i++){
            if (this->particleBase->isActive[i] == false){
                // To prevent truncated error, set the value of vorticity to 0.0
                this->particleBase->vorticity[i] = 0.0e0;
                this->particleBase->gz[i] = 0.0e0;
            }
        }
        // // [OLD SECTION] not use but for reference
        // // Find the maximum absolute value of vorticity
        // double vor_max = 0.0e0;
        // for (int i = 0; i < this->particleBase->num; i++){
        //     if (std::abs(this->particleBase->vorticity[i]) > vor_max)
        //     vor_max = std::abs(this->particleBase->vorticity[i]);
        // }

        // // Evaluate the active flag particle
        // this->particleBase->isActive.clear();this->particleBase->isActive.resize(this->particleBase->num,false);
        // #pragma omp parallel for
        // for (int i = 0; i < this->particleBase->num; i++){
        //     double vor = std::abs(this->particleBase->vorticity[i]);
        //     if (vor >= Pars::active_sig*vor_max /* || this->particleBase->isActive[i] == true*/){
        //         this->particleBase->isActive[i] = true;
        //     }
        //     else{
        //         this->particleBase->isActive[i] = false;
        //         // To prevent truncated error, set the value of vorticity to 0.0
        //         this->particleBase->vorticity[i] = 0.0e0;
        //         this->particleBase->gz[i] = 0.0e0;  
        //     }
        // }
    }
    else if (DIM == 3){
        #pragma omp parallel for
        for (int i = 0; i < this->particleBase->num; i++){
            if (this->particleBase->isActive[i] == false){
                // To prevent truncated error, set the value of vorticity to 0.0
                this->particleBase->vorticity[i] = 0.0e0;
                this->particleBase->vortx[i] = 0.0e0;
                this->particleBase->vorty[i] = 0.0e0;
                this->particleBase->vortz[i] = 0.0e0;
            }
        }
    }
    return;
}

// [DEBUG TOOLS] Saving data each neighbor
void save_ngh_all(const Particle &_myPar, int ID){
    std::ofstream write;
    std::string name = "data_result/nghData_" + std::to_string(value) + "_" + std::to_string(ID) + ".csv";
    write.open(name.c_str(), std::ofstream::out);

    // Write the table header
    write << "x"
          << "," << "y";
    if (DIM > 2)
    write << "," << "z";
    write << "," << "s"
          << "\n";

    // Write the table data
    for (size_t i = 0; i < _myPar.neighbor.at(ID).size(); i++){
    const int &_nghID = _myPar.neighbor.at(ID)[i];
    write << _myPar.x[_nghID]
          << "," << _myPar.y[_nghID];
    if (DIM > 2)
    write << "," << _myPar.z[_nghID];
    write << "," << _myPar.s[_nghID]
          << "\n";    
    }

    // Write the circular cylinder of the support radius
    int N = 70;
    double dAlp = 2*M_PI/N;
    double R = _myPar.s[ID] * Pars::r_sup * Pars::r_buff;
    
    for (int i = 0; i < N; i++){
    write << (_myPar.x[ID] + R * cos(i*dAlp))
          << "," << (_myPar.y[ID] + R * sin(i*dAlp));
    if (DIM > 2)
    write << "," << (_myPar.z[ID]);
    write << "," << (Pars::sigma/2)
          << "\n";    
    }

    write.close();
    return;
}