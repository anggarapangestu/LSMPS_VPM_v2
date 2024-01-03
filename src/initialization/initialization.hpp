#ifndef INCLUDED_INITIALIZATION
#define INCLUDED_INITIALIZATION

#include "../../Utils.hpp"
#include "../grid_block/gridNode.hpp"
#include "../geometry/geometry.hpp"

/**
 *  @brief  Initialization class consisted of method to generate 
 *  initial particle distribution.
 * 
 *  @headerfile initialization.hpp
*/
class initialization
{
private:
    // Private method
    void generate_particle(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &_baseGrid);
    void read_2d_particle(Particle &_particle, int _iteration);
    void read_3d_particle(Particle &_particle, int _iteration);

    // 2D Particles initialization
    void init_2d_single_res(Particle &_particle);                                                   // [DONE]
    void init_2d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_2d_multi_res_multi_block(Particle &_particle, const std::vector<Body> &_bodyList);     // Ongoing
    void init_2d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]

    // 3D Particles initialization
    void init_3d_single_res(Particle &_particle);                                                    // [DONE]
    void init_3d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_3d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]

    // Particles initialization by grid block method
    void init_par_grid_block(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &g);

public:
    // Initialization global procedure
    void initialize_particle(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &_baseGrid);
};
#endif