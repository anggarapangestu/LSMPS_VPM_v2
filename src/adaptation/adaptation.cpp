#include "adaptation.hpp"

/**
 *  @brief  Particle adaptation manager distribution manager. Generate the particle data
 *  using the selected type of distribution by user in global.cpp file.
 *         
 *  @param  _particle  Particle data container.
 *  @param  _bodyList  The list of body data container used for particle generation.
 *  @param  _gridNode  Grid node data container used for particle generation.
*/
bool adaptation::get_adaptation(const Particle &_parEval, Particle *&_parBase, GridNode &_baseGrid){
    // Perform particle adaptation
    // ***************************
    // Adaptation using gridNode
    GridNodeAdapt _gridNode;
    bool adapt_check;
    adapt_check = _gridNode.get_adaptation(_baseGrid, *_parBase, _parEval);

    // Check if there is no adaptation
    if (adapt_check == false){
        return adapt_check;
    }
    
    // Retrieve the data to base particle
    delete _parBase;        // Avoid memory leaking (The old data template particle is not necessary)
    _gridNode.take_particle_pointer(_parBase);

    // Update the particle data
    // ************************
    int &num = _parBase->num;
    // Resize the other particle data
    // [NOTE] Modify if there is any properties update
    _parBase->u.resize(num);
    _parBase->v.resize(num);
    if (DIM == 3){
        _parBase->w.resize(num);
        _parBase->gx.resize(num);
        _parBase->gy.resize(num);
    }
    _parBase->gz.resize(num);
    _parBase->vorticity.resize(num);

    return adapt_check;
}