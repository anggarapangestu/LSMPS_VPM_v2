#ifndef INCLUDED_MULTIBLOCK
#define INCLUDED_MULTIBLOCK

#ifndef INCLUDED_GLOBAL
#include "../../../global.hpp"
#endif

#ifndef INCLUDED_UTILS
#include "../../../Utils.hpp"
#endif

class multiblock_adaption
{
private:
    /* data */
    double _rSigmaIn;
    double _rSigmaOut;
    double _sigmaMin;
    double _sigmaMax;
    double _maxBodyLength;           // maximum length of the obstacle
    std::vector<double> _bodycenter; // center coordinate of obstacle
    std::vector<double> _bodyLength; // length of the obstacle

    std::vector<double> _minBodyCoordinate; // minimum limit coordinate of obstacle
    std::vector<double> _maxBodyCoordinate; // maximum limit coordinate of obstacle

    std::vector<double> _centerPanelX;  // Center x coordinate of each obstacle panel
    std::vector<double> _centerPanelY;  // Center y coordinate of each obstacle panel
    std::vector<double> _normalPanelX;  // Unit vector length of each obstacle panel in x-direction
    std::vector<double> _normalPanelY;  // Unit vector length of each obstacle panel in y-direction

    std::vector<double> _xB; // x coordinate of generated block
    std::vector<double> _yB; // y coordinate of generated block
    std::vector<double> _sB; // equivalent core size of generated block

    double target_resolution_function(const int idx, const Particle &particle, const Body &body);
    // ! 1st step
    void generate_multiblock(const Particle &particle, const Body &body);
    // ! 2nd step [optional]
    void generate_multiblock_wake(const Particle &particle, const Body &body);

public:
    void get_multiblock_adaption(Particle &particle, const Body &body);
};
#endif