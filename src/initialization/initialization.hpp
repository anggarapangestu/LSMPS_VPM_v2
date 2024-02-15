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

    // Testing Initialization
    void init_2d_test_1(Particle &_particle);
    void init_2d_test_2(Particle &_particle);

    // 2D Particles initialization
    void init_2d_single_res(Particle &_particle);                                                   // [DONE]
    void init_2d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_2d_multi_res_multi_block(Particle &_particle, const std::vector<Body> &_bodyList);     // Ongoing
    void init_2d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]

    // 3D Particles initialization
    void init_3d_single_res(Particle &_particle);                                                    // [DONE]
    void init_3d_multi_res_single_block(Particle &_particle, const std::vector<Body> &_bodyList);    // [DONE]
    void init_3d_multi_res_body_adjusted(Particle &_particle, const std::vector<Body> &_bodyList);   // [DONE]

    // Vorticity Initialization type
    void perlman_vorticity(Particle &_particle);

    // Particles initialization by grid block method
    void init_par_grid_block(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &g);

public:
    // Particle reader
    void read_2d_particle(Particle &_particle, int _iteration);
    void read_3d_particle(Particle &_particle, int _iteration);

    // Solution to initial vorticity
    void perlman_velocity_solution(Particle &_particle);

    // Update boundary
    void update_domain_boundary(Particle &_particle);

    // Initialization global procedure
    void initialize_particle(Particle &_particle, const std::vector<Body> &_bodyList, GridNode &_baseGrid);
    void initialize_vorticity(Particle &_particle);
};
#endif