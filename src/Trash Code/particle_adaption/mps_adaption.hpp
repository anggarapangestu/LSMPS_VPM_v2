#ifndef INCLUDED_MPS
#define INCLUDED_MPS

#ifndef INCLUDED_GLOBAL
#include "../../../global.hpp"
#endif

#ifndef INCLUDED_UTILS
#include "../../../Utils.hpp"
#endif

// #ifndef INCLUDED_PARTICLE_ADAPTION
// #include "particle_adaption.hpp"
// #endif

#ifndef INCLUDED_BASE_REMESHING
#include "../base_remeshing.hpp"
#endif

#ifndef INCLUDED_NEIGHBOR
#include "../../neighbor/neighbor.hpp"
#endif

#ifndef INCLUDED_LSMPS
#include "../../LSMPS/LSMPSa.hpp"
#endif

// class particle_adaption;
class base_remeshing;
class neighbor;
class mps_adaption
{
private:
    // * declaring instance
    // particle_adaption _particle_adaption;
    base_remeshing _base_remeshing;
    neighbor d_neighbor;

    std::vector<double> _Dp;     // target diameter
    std::vector<bool> _isMerged; // merging condition

    // TODO: generate target diameter from spatial function
    void generate_target_diameter(Particle &particle, const Body &body);
    // TODO: particle merging algorithms
    void particle_merging(Particle &particle, const Body &body);
    // TODO: weighting function
    double weighting_function(const double &rij, const double &Rij);

public:
    // TODO: perform MPS adaption 
    void get_mps_adaption(Particle &particle, const Body &body);
    // TODO: calculate packing ratio of particles
    std::vector<double> get_packing_ratio(const Particle &particle);
    // TODO; move particle algoritm
    void particle_moved(Particle &particle);
};
#endif
