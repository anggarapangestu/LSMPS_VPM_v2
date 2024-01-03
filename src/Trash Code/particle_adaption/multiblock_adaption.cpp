#include "multiblock_adaption.hpp"

#pragma region public methods
void multiblock_adaption::get_multiblock_adaption(Particle &particle, const Body &body)
{
    this->_bodycenter.resize(2);
    this->_bodyLength.resize(2);
    this->_minBodyCoordinate.resize(2);
    this->_maxBodyCoordinate.resize(2);
    //0 = arah sumbu x di cartesian
    //1 = arah sumbu y di cartesian
    //menentukan remeshing yang dekat sama object
    this->_minBodyCoordinate[0] = *std::min_element(body.x.begin(), body.x.end()) - 30*Pars::sigma;//- 1.25; // Changed 13:13 March 23 by Angga, initial state just uncomment to change the constant
    this->_minBodyCoordinate[1] = *std::min_element(body.y.begin(), body.y.end()) - 40*Pars::sigma;//- 2;
    this->_maxBodyCoordinate[0] = *std::max_element(body.x.begin(), body.x.end()) + 30*Pars::sigma;//+ 1.25;
    this->_maxBodyCoordinate[1] = *std::max_element(body.y.begin(), body.y.end()) + 40*Pars::sigma;//+ 2;// +2.5 nxt time

    this->_bodycenter[0] = (this->_maxBodyCoordinate[0] + this->_minBodyCoordinate[0]) * 0.5;
    this->_bodycenter[1] = (this->_maxBodyCoordinate[1] + this->_minBodyCoordinate[1]) * 0.5;

    this->_bodyLength[0] = (this->_maxBodyCoordinate[0] - this->_minBodyCoordinate[0]);
    this->_bodyLength[1] = (this->_maxBodyCoordinate[1] - this->_minBodyCoordinate[1]);

    this->_maxBodyLength = _bodyLength[0] > _bodyLength[1] ? _bodyLength[0] : _bodyLength[1];

    double _nPanel = body.x.size() - 1;
    this->_centerPanelX.resize(_nPanel);
    this->_centerPanelY.resize(_nPanel);
    this->_normalPanelX.resize(_nPanel);
    this->_normalPanelY.resize(_nPanel);

    //buat center panel dan normal panelnya per panel 
    for (int i = 0; i < _nPanel; i++)
    {
        // determine midpoint of panel
        this->_centerPanelX[i] = body.x[i] + ((body.x[i + 1] - body.x[i]) * 0.5e0);
        this->_centerPanelY[i] = body.y[i] + ((body.y[i + 1] - body.y[i]) * 0.5e0);
        // determine normal panel
        double _dx = body.x[i + 1] - body.x[i];
        double _dy = body.y[i + 1] - body.y[i];
        double _theta = atan2(_dy, _dx);
        this->_normalPanelX[i] = sin(_theta);  // unit normal vector in x-direction
        this->_normalPanelY[i] = -cos(_theta); // unit normal vector in y-direction
    }
    
    // --------------------------------------------------------------------------------- //
    this->_rSigmaIn = this->_maxBodyLength * 0.5 + 10 * Pars::sigma; // ! still CHANGEABLE
    this->_rSigmaOut = this->_maxBodyLength * 2;                     // ! still CHANGEABLE
    this->_sigmaMin = Pars::sigma * 1;                               // ! still CHANGEABLE
    this->_sigmaMax = Pars::sigma * 2;                               // ! still CHANGEABLE, should be changed to use RATIO
    // --------------------------------------------------------------------------------- //

    // TODO: generate multiblock mesh
    this->generate_multiblock(particle, body);
    // this->generate_multiblock_wake(particle, body);// !!!

    particle.x = this->_xB;
    particle.y = this->_yB;
    particle.s = this->_sB;
    particle.num = this->_xB.size();
    particle.u.resize(particle.num, 0.0e0);
    particle.v.resize(particle.num, 0.0e0);
    particle.gz.resize(particle.num, 0.0e0);
    particle.neighbor.resize(particle.num, std::vector<int>());
    particle.isActive.resize(particle.num, false);
    particle.chi.resize(particle.num, 0.0);

}
#pragma endregion

#pragma region private methods
void multiblock_adaption::generate_multiblock(const Particle &particle, const Body &body)
{
    // copy to local variables
    const std::vector<double> &_xp = particle.x;
    const std::vector<double> &_yp = particle.y;
    const std::vector<double> &_sp = particle.s;
    // internal variables
    int *npar = new int[2];
    double *maxpar = new double[2];
    double *minpar = new double[2];
    
    maxpar[0] = *std::max_element(_xp.begin(), _xp.end()) + 10 * Pars::sigma; // !
    maxpar[1] = *std::max_element(_yp.begin(), _yp.end()) + 10 * Pars::sigma; // !
    minpar[0] = *std::min_element(_xp.begin(), _xp.end()) - 10 * Pars::sigma; // !
    minpar[1] = *std::min_element(_yp.begin(), _yp.end()) - 10 * Pars::sigma; // !
    
    for (size_t i = 0; i < 2; i++)
    {
        npar[i] = std::ceil(std::abs(maxpar[i] - minpar[i]) / Pars::D_0);
    }

    std::vector<double> _xbox0;
    std::vector<double> _ybox0;
    
    //bikin titik titiknya di dummy
    for (size_t i = 0; i < npar[0]; i++)
    {
        double _xBoxPos = minpar[0] + static_cast<double>(i) * Pars::D_0;
        _xbox0.push_back(_xBoxPos);
    }
    
    for (size_t i = 0; i < npar[1]; i++)
    {
        double _yBoxPos = minpar[1] + static_cast<double>(i) * Pars::D_0;
        _ybox0.push_back(_yBoxPos);
    }

    this->_xB.clear();
    this->_yB.clear();
    this->_sB.clear();
    
    for (size_t i = 0; i < npar[0]; i++)
    {
        for (size_t j = 0; j < npar[1]; j++)
        {
            this->_xB.push_back(_xbox0[i]);
            this->_yB.push_back(_ybox0[j]);
            this->_sB.push_back(Pars::D_0);
        }
    }

    bool _isFinish = false;
    do
    {
        std::vector<double> _xBtemp;
        std::vector<double> _yBtemp;
        std::vector<double> _sBtemp;

        Particle _particle;
        _particle.x = _xB;
        _particle.y = _yB;
        _particle.s = _sB;

        for (size_t i = 0; i < _xB.size(); i++)
        {
            double _sTarget = target_resolution_function(i, _particle, body);
            if (_sB[i] > _sTarget)
            {
                // =========== Devide And Conquer ================
                _xBtemp.push_back(_xB[i] + 0.25 * _sB[i]);
                _xBtemp.push_back(_xB[i] + 0.25 * _sB[i]);
                _xBtemp.push_back(_xB[i] - 0.25 * _sB[i]);
                _xBtemp.push_back(_xB[i] - 0.25 * _sB[i]);
                
                _yBtemp.push_back(_yB[i] + 0.25 * _sB[i]);
                _yBtemp.push_back(_yB[i] - 0.25 * _sB[i]);
                _yBtemp.push_back(_yB[i] + 0.25 * _sB[i]);
                _yBtemp.push_back(_yB[i] - 0.25 * _sB[i]);
                for (int nS = 0; nS < 4; nS++)
                    _sBtemp.push_back(0.5 * _sB[i]);
            }
            else
            {
                _xBtemp.push_back(_xB[i]);
                _yBtemp.push_back(_yB[i]);
                _sBtemp.push_back(_sB[i]);
            }
        }

        this->_xB = _xBtemp;
        this->_yB = _yBtemp;
        this->_sB = _sBtemp;

        _xBtemp.clear();
        _yBtemp.clear();
        _sBtemp.clear();

        double minSb = *std::min_element(_sB.begin(), _sB.end());
        // cout << minSb << endl;
        if (minSb <= Pars::sigma * 1) // !!!!!!!!!!!!!
        {
            _isFinish = true;
        }
    } while (_isFinish == false);
}

double multiblock_adaption::target_resolution_function(const int idx, const Particle &particle, const Body &body)
{
    // const std::vector<double> &xp = particle.x;
    // const std::vector<double> &yp = particle.y;
    // const std::vector<double> &sp = particle.s;
    double _xp = particle.x[idx];
    double _yp = particle.y[idx];
    double _sp = particle.s[idx];
    double _sTarget;

    // double r2bodycenter = std::sqrt(std::pow(xp[idx] - this->_bodycenter[0], 2) + std::pow(yp[idx] - this->_bodycenter[1], 2));
    // double _dx = xp[idx] - this->_bodycenter[0];
    // double _dy = yp[idx] - this->_bodycenter[1];
    // if (r2bodycenter < _rSigmaIn)
    // {
    //     _sTarget = _sigmaMin;
    // }
    // else if (r2bodycenter >= _rSigmaIn && r2bodycenter < _rSigmaOut)
    // {
    //     double _r = r2bodycenter - _rSigmaIn;
    //     double _L = _rSigmaOut - _rSigmaIn;
    //     _sTarget = (_r / _L) * (_sigmaMax - _sigmaMin) + _sigmaMin;
    // }
    // else
    // {
    //     _sTarget = _sigmaMax;
    // }
    if (_xp >= this->_minBodyCoordinate[0] - 10 * Pars::sigma && _xp <= this->_maxBodyCoordinate[0] + 10 * Pars::sigma &&
        _yp >= this->_minBodyCoordinate[1] - 10 * Pars::sigma && _yp <= this->_maxBodyCoordinate[1] + 10 * Pars::sigma)
    {
        // initial guess value
        // int _k1 = 0;
        // double _rmins = 100.0e0;
        // double _signNormal = 0.0e0;
        // double _testPin = 0.0e0;

        // for (int k = 0; k < body.x.size() - 1; k++)
        // {
        //     // Distance between center panel and the point
        //     double _r2panel = std::sqrt(std::pow((_xp - this->_centerPanelX[k]), 2) +
        //                                 std::pow((_yp - this->_centerPanelY[k]), 2));
        //     if (std::abs(_r2panel) <= _rmins)
        //     {
        //         _rmins = _r2panel; // find the center panel closest to the point
        //         _k1 = k;           // Name of closest panel
        //     }
        // }

        // /* normal distance from considered point to panel
        // * = 0 when grid node lie on the panel, or extended-panel
        // */
        // double _rNormal = (((_xp - this->_centerPanelX[_k1]) * (this->_normalPanelX[_k1])) +
        //                    ((_yp - this->_centerPanelY[_k1]) * (this->_normalPanelY[_k1])));

        // if (_rNormal == 0.0e0)
        //     _signNormal = 0.0e0;
        // else
        //     _signNormal = _rNormal / (std::abs(_rNormal));

        // ==== Midpoint Distance: using distance from considered point to midpoint of triangle face
        // read more .../Report_progress/April2017_Airfoil_Penalization_inOutside_problem.odp, part B
        // _testPin = _rmins * _signNormal;
        // ==== Normal Distance: Using normal distance from considered point to triangle face
        // testpin(i) = r_normal;

        // if (std::abs(_testPin) >= 10 * Pars::sigma)
        //     _sTarget = _sigmaMax;
        // else
            _sTarget = _sigmaMin;
    }
    else
    {
        _sTarget = _sigmaMax;
    }

    return _sTarget;
}

void multiblock_adaption::generate_multiblock_wake(const Particle &particle, const Body &body)
{
    // copy to local variables
    const std::vector<double> &_xp = particle.x; // copy x coordinate into local variables
    const std::vector<double> &_yp = particle.y; // copy y coordinate into local variables
    const std::vector<double> &_sp = particle.s; // copy blob size into local variables
    // internal variables
    double *_maxpar = new double[2]; // maximum coordinate of particles
    double *_minpar = new double[2]; // minimum coordinate of particles

    _maxpar[0] = *std::max_element(_xp.begin(), _xp.end()) + 10 * Pars::sigma; // !
    _maxpar[1] = *std::max_element(_yp.begin(), _yp.end()) + 10 * Pars::sigma; // !
    _minpar[0] = *std::min_element(_xp.begin(), _xp.end()) - 10 * Pars::sigma; // !
    _minpar[1] = *std::min_element(_yp.begin(), _yp.end()) - 10 * Pars::sigma; // !

    const double _r0 = this->_rSigmaIn;            // inner circle radius
    const double _r1 = this->_maxBodyLength * 1.5; // outer circle radius

    // TODO: 1. particle inside main circle
    std::vector<double> _xInner;
    std::vector<double> _yInner;
    std::vector<double> _sInner;
    int _nx = 2 * std::ceil(this->_rSigmaIn / Pars::sigma);
    int _ny = 2 * std::ceil(this->_rSigmaIn / Pars::sigma);
    for (int i = 0; i < _nx; i++)
    {
        double _x = -this->_rSigmaIn + static_cast<double>(i) * Pars::sigma;
        double _s = Pars::sigma;
        for (int j = 0; j < _ny; j++)
        {
            double _y = -this->_rSigmaIn + static_cast<double>(j) * Pars::sigma;
            _xInner.push_back(_x + this->_bodycenter[0]);
            _yInner.push_back(_y + this->_bodycenter[1]);
            _sInner.push_back(_s);
        }
    }
    std::vector<bool> _isDeleted(_xInner.size(), false);

    // TODO: 2. particle outside main circle
    std::vector<double> _xOuter;
    std::vector<double> _yOuter;
    std::vector<double> _sOuter;
    std::vector<double> _rDistribution;                                    // to store r for outer domain
    std::vector<double> _sDistribution;                                    // to store variable blob size for outer domain
    double _d_r_rho = std::log(Pars::sigma * (std::exp(1) - 1) / _r1 + 1); // distance in mapped radius for outer domain

    double _ri, _i = 0;
    do
    {
        double _rho = static_cast<double>(_i * _d_r_rho);           // assign rho (mapped radius)
        double _riRatio = (std::exp(_rho) - 1) / (std::exp(1) - 1); // logaritmic distribution
        _ri = _r0 + (_riRatio * _r1);                               // map to cartesian coordinate
        _rDistribution.push_back(_ri);                              // assign r
        _i++;
    } while (_ri <= _maxpar[0]);

    _sDistribution.push_back(Pars::sigma); // initial value
    for (size_t i = 1; i < _rDistribution.size(); i++)
    {
        _sDistribution.push_back(_rDistribution[i] - _rDistribution[i - 1]);
    }
    _sDistribution.pop_back(); // remove last element

    // TODO: generate mesh in wake region
    for (size_t i = 0; i < _rDistribution.size(); i++)
    {
        double _perimeter = 2 * Pars::pi * _rDistribution[i];
        int _ny = std::ceil(_perimeter / _sDistribution[i]);
        double _dteta = (Pars::pi * 2) / _ny;

        for (size_t j = 0; j < _ny; j++)
        {
            double _teta = static_cast<double>(j) * _dteta;
            double _x = _rDistribution[i] * std::cos(_teta);
            double _y = _rDistribution[i] * std::sin(_teta);
            double _s = _sDistribution[i];

            _xOuter.push_back(_x + this->_bodycenter[0]);
            _yOuter.push_back(_y + this->_bodycenter[1]);
            _sOuter.push_back(_s);
        }
    }

    // labelling deleted particles from inner circle.
    // thus it's not intersecting with outer particles
    for (size_t i = 0; i < _xInner.size(); i++)
    {
        double _r = std::sqrt(std::pow(_xInner[i] - this->_bodycenter[0], 2) + std::pow(_yInner[i] - this->_bodycenter[1], 2));
        if (_r > (this->_rSigmaIn - (Pars::sigma * 0.4))) // !
        {
            _isDeleted[i] = true;
        }
    }

    // TODO: assign full mesh
    std::vector<double> _xBtemp;
    std::vector<double> _yBtemp;
    std::vector<double> _sBtemp;
    for (size_t i = 0; i < _xInner.size(); i++)
    {
        if (_isDeleted[i] == false)
        {
            _xBtemp.push_back(_xInner[i]);
            _yBtemp.push_back(_yInner[i]);
            _sBtemp.push_back(_sInner[i]);
        }
    }
    for (size_t i = 0; i < _xOuter.size(); i++)
    {
        _xBtemp.push_back(_xOuter[i]);
        _yBtemp.push_back(_yOuter[i]);
        _sBtemp.push_back(_sOuter[i]);
    }
    this->_xB = _xBtemp;
    this->_yB = _yBtemp;
    this->_sB = _sBtemp;
    std::cout << "\ncurrent particle number : " << _xBtemp.size() << std::endl;

    // delete local variables
    delete[] _maxpar;
    delete[] _minpar;
}
#pragma endregion
