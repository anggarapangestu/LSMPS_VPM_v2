#include "gridNodeNgh.hpp"

/**
 *  @brief  Evaluate the neighboring interaction for each particle. 
 *  Interaction pairs are determined by grid node relation.
 *
 *  @param  _nghIDList	[OUTPUT] Neighbor ID list.
 *  @param  _baseGrid 	The node container as neighbor evaluation tools.
 *  @param  _evalPar  	The particle to evaluate neighbor.
*/
void GridNodeNgh::find_neighbor(std::vector<std::vector<int>> &_nghIDList, 
								const GridNode &_baseGrid, const Particle &_par)
{
    // **Evaluate neighbor interation (spatial hash)
    _nghIDList.clear();
    _nghIDList.resize(_par.num);
    std::vector<bool> evalFlag(_par.num, false);

	// Collect all leaf node only (for the sake of parallel programming)
	std::vector<const Node*> leafNodeList;
	for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		if(!(_currNode->isLeaf)) continue;
		leafNodeList.push_back(_currNode);
	}
	
	// Iterate through all node in the grid
	#pragma omp parallel for			// [-> Create a problem in the server computer [RESOLVED by create a container 'leafNodeList' first]]
	for (size_t n = 0; n < leafNodeList.size(); n++){
	// for (const auto&[_nodeID,_currNode] : _baseGrid.nodeMap){
		// Alias to the current node
		const Node *&_currNode = leafNodeList[n];

		// // Only evaluate the leaf node
        // if(!(_currNode->isLeaf)) continue;
        
        // Get the list of neighbor node
        std::vector<Node*> _nghNodeList;
        _baseGrid.findNghAll(_nghNodeList, _currNode);

		// Iterate through all point inside the current node (it will share the same neighbor node)
		for (size_t i = 0; i < _currNode->parList.size(); i++){
			// Aliasing the current particle ID
			const int &ID_i = _currNode->parList[i];
			
            // Evaluate to each point inside the grid neighbor
			for (const auto &nghNode : _nghNodeList){
				for (const auto &ID_j : nghNode->parList){		
					// An exception for include itself on neighbor list
					if (!Pars::flag_ngh_include_self && (ID_j == ID_i)) continue;

                    // Calculate the distance square between two points
					double _dr2 = 0;
					// Calculate in x direction
						double _dx = _par.x[ID_i] - _par.x[ID_j];
						_dr2 += (_dx*_dx);
					// Calculate in y direction
						double _dy = _par.y[ID_i] - _par.y[ID_j];
						_dr2 += (_dy*_dy);
					// Calculate in z direction
					if (DIM > 2){
						double _dz = _par.z[ID_i] - _par.z[ID_j];
						_dr2 += (_dz*_dz);
					}

					// Evaluate distance
					if (Pars::opt_ngh_interact == 1){
						// Check the particle i
						double _rSup = _par.s[ID_i] * Pars::r_sup * Pars::r_buff;   // Support radius of particle i
						if (_dr2 < (_rSup*_rSup)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
					else if (Pars::opt_ngh_interact == 2){
						// Support radius of average size
						double _rSupAve = ((_par.s[ID_i] + _par.s[ID_j]) / 2.0e0) * Pars::r_sup * Pars::r_buff;
						if (_dr2 < (_rSupAve*_rSupAve)){
							_nghIDList[ID_i].push_back(ID_j);
						}
					}
				}
			}
		}
	}

    return;
}