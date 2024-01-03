#include "neighbor.hpp"

/**
 *  @brief A grid generator, used as a temporary data grouping for neighbor evaluation.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_parEval  The particle data for grid construction.
 *  @param	_tempGrid  [OUTPUT] The temporary grid to be generated.
 */
void neighbor::create_temp_grid(const Particle& _par, NghBaseGrid &_grid){
    // Evaluate extreme coordinate
	double maxCoor[DIM];
	double minCoor[DIM];
	// Evaluate extreme in x direction
		minCoor[0] = *(std::min_element(_par.x.begin(),_par.x.end()));
		maxCoor[0] = *(std::max_element(_par.x.begin(),_par.x.end()));
	// Evaluate extreme in y direction
		minCoor[1] = *(std::min_element(_par.y.begin(),_par.y.end()));
		maxCoor[1] = *(std::max_element(_par.y.begin(),_par.y.end()));
	// Evaluate extreme in z direction
	if (DIM > 2){
		minCoor[2] = *(std::min_element(_par.z.begin(),_par.z.end()));
		maxCoor[2] = *(std::max_element(_par.z.begin(),_par.z.end()));
	}
    // Perform grid initialization
    _grid.initialize_grid(maxCoor, minCoor, _par.s);

    return;
}

/**
 *  @brief Neighbor search function manager for neighbor search toward 
 *  the particle data itself.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_evalPar  The particle data for neighbor evaluation.
 *  @param  _baseGrid  The grid node data for Grid Node type neighbor search only.
 */
void neighbor::neigbor_search(Particle &parEval, const GridNode &_baseGrid){
    // Print log to console
    printf("Evaluating neighbor ...\n");
    double _time = omp_get_wtime();

    // Direct neighbor search evaluation
    if (Pars::opt_neighbor == 0){
        printf("%s<+> Type 0: Direct neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Evaluate the neighbor
        directFindNgh _directFind;
        _directFind.find_neighbor(parEval.neighbor, parEval.s, 
                                 parEval.x, parEval.y, parEval.z);
    }
    
    // Link list neighbor search evaluation
    else if (Pars::opt_neighbor == 1){
        printf("%s<+> Type 1: Link list neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Create temporary grid
        NghBaseGrid _tempGrid;
        this->create_temp_grid(parEval, _tempGrid);

        // Evaluate the neighbor
        LinkListNgh _linkList;
        _linkList.find_neighbor(parEval.neighbor, _tempGrid, parEval.s,
                                parEval.x, parEval.y, parEval.z);
    }
    
    // [NEED FURTHER REFACTOR]
    // Included neighbor search package from Cell List
    else if (Pars::opt_neighbor == 2){
        printf("%s<+> Type 2: Cell list neighbor search %s\n", FONT_CYAN, FONT_RESET);
        this->cellListData.findNeighbor(parEval);
        // this->cellListData.checkNGH(parEval);
    }

    // Spatial hash neighbor search evaluation
    else if (Pars::opt_neighbor == 3){
        printf("%s<+> Type 3: Spatial hash neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Create temporary grid
        NghBaseGrid _tempGrid;
        this->create_temp_grid(parEval, _tempGrid);

        // Evaluate the neighbor
        SpatialHashNgh _spatialHash;
        _spatialHash.find_neighbor(parEval.neighbor, _tempGrid, parEval.s,
                                   parEval.x, parEval.y, parEval.z);
    }

    // Grid node neighbor search evaluation
    else if (Pars::opt_neighbor == 4){
        printf("%s<+> Type 4: Grid node neighbor search %s\n", FONT_CYAN, FONT_RESET);
        
        // Evaluate the neighbor
        GridNodeNgh _gridNode;
        _gridNode.find_neighbor(parEval.neighbor, _baseGrid, parEval);
    }
    
    // Neighbor search summary time display
    _time = omp_get_wtime() - _time;
    printf("<-> Neighbor search computation time:  [%f s]\n", _time);
}

/**
 *  @brief Neighbor search function manager for two different particle data set.
 *  Evaluation is done in one direction to evaluate the source particle data as
 *  the neighbor of the target particle data.
 *  NOTE: Works well both for 2D and 3D simulation.
 *  
 *  @param	_evalPar  The particle target particle data for neighbor evaluation.
 *  @param	_srcPar  The source particle data, work as a data format.
 */
void neighbor::neigbor_search(Particle& _evalPar, Particle &_srcPar){
    // Print log to console
    printf("Evaluating inter search neighbor ...\n");
    double _time = omp_get_wtime();

    // Evaluate the neighbor
    InterSearchNgh _interSearch;
    _interSearch.find_neighbor(_srcPar, _evalPar, _evalPar.neighbor);
    
    // Neighbor search summary time display
    _time = omp_get_wtime() - _time;
    printf("<-> Neighbor search computation time:  [%f s]\n", _time);
}


// ************************************************************************************
// ====================================================================================

// [OLD PACKAGE] The particle neighbor list
void neighbor::cell_list_init(Particle& parEval){
    // Generate the Cell List
    printf("Generating cell list ...\n");
    clock_t t = clock();

    // Initialize the Cell List
    cellListData.initCellList(parEval);
    
    // Evaluate the particle cell ID
    cellListData.createCellList(parEval);

    // // Show particle neighbor
    // cellListData.showBasisNgh();
    
    // particle generation initialization summary time display
    t = clock() - t;
	printf("<-> Cell list initialization\n");
    printf("    comp. time:                        [%f s]\n", (double)t/CLOCKS_PER_SEC);

    // Evaluate the particle cell ID
    if (Pars::flag_save_cell){cellListData.saveCellData();}
    
}

// ************************************************************************************
// ====================================================================================

// [OLD PACKAGE] Particle adaptation 
bool neighbor::particle_adaptation(const Particle& parEval, Particle& baseParticle, std::vector<double>& PARsize){
    // Note: At initial the distribution of parEval and baseParticle is the same but not the properties

    // computational time accumulation
    clock_t t;
    t = clock();
    /* Procedure of particle adaptation:
       [1] Evaluate each particle, determine if need adaptation of not
       [2] Perform the particle adaptation if necesarry (based on the procedure 1)
    */

    // PROCEDURE 1 ! : Evaluating particle adaptation
    // *************
    // The marker whether the adaptation need to performed or not
    bool _adaptation = false;

    // TODO: Find maximum vorticity
    double vor_max = 0.0e0;
    for (int i = 0; i < parEval.num; i++)
    {
        vor_max = vor_max > std::abs(parEval.vorticity[i]) ? 
                  vor_max : std::abs(parEval.vorticity[i]);
    }

    // Initialize the Cell Level Up setting (Setting for particle adaptation)
    this->cellListData.setLevelUp(0, 0, 0);

    // Assign the cell for Level Up (evaluate each particle)
    for (int i = 0; i < parEval.num; i++)
    {   
        // ... LEVEL 1 CHECK ...
        // Adaptation evaluation still based on vorticity value of the evaluated particle
        if (std::abs(parEval.vorticity[i]) >= 1.0e-5 * vor_max)  // Threshold set to be 1.0e-5 of the max value,
                                                                 //   set lower than the particle redistribution
        {
            // ... LEVEL 2 CHECK ...
            // Do adaptation if the particle level != maxlevel (or lower) and the
            if (parEval.level[i] < Pars::max_level){
                // Assign the cell for Level Up
                this->cellListData.setLevelUp(parEval.basis_label[i], parEval.cell_label[i], 1);

                // Do the adaptation
                _adaptation = true;
            }
        }
    }
    
    // Adaptation procedure:
    // -> List all the cell to be devided (save the basis-cell-ID pair)
    // -> Find the neighboring cell to be divided (save the basis-cell-ID pair)
    // -> At each cell: divide the particle into the target level (only if par level < cell level)
    // -> Put the particle into the corresponding cell
    // -> Devide the cell
    
    // PROCEDURE 2 ! : Performing particle adaptation
    // *************
    if (_adaptation == true)
    {
        // Perform the adaptive algorithm
        this->cellListData.performAdaptive(baseParticle, 0, Pars::max_level, PARsize);

        // Particle adaptation summary time display
        t = clock() - t;
        printf("<+> Particle adaptation done ...\n");
        printf("<+> Number of particle after adaptation : %8d\n", baseParticle.num);
        printf("<-> Particle adaptation calculation \n");
        printf("    comp. time:                        [%f s]\n\n", (double)t/CLOCKS_PER_SEC);
        
        // Evaluate the verlet list of new particle distribution
        // this->neigbor_search(baseParticle);		// Already taken on the particle adaptation
    }else{
        printf("<+> No particle adaptation\n");
    }

    return _adaptation;
}

