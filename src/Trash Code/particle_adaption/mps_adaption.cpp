#include "mps_adaption.hpp"

#pragma region public_methods
void mps_adaption::get_mps_adaption(Particle &particle, const Body &body)
{
    this->generate_target_diameter(particle, body);
    this->particle_merging(particle, body);
}
#pragma endregion

#pragma region private_methods_merging
void mps_adaption::generate_target_diameter(Particle &particle, const Body &body)
{
    vector<double> &xp = particle.x;
    vector<double> &yp = particle.y;
    vector<double> &sp = particle.s;

    uint ntotal = xp.size();

    double *bodycenter = new double[2];
    bodycenter[0] = (*std::max_element(body.x.begin(), body.x.end()) + *std::min_element(body.x.begin(), body.x.end())) * 0.5;
    bodycenter[1] = (*std::max_element(body.y.begin(), body.y.end()) + *std::min_element(body.y.begin(), body.y.end())) * 0.5;

    // --------------------------------------------------------------------------------- //
    double _rSigmaIn = Pars::ly * 0.5 + Pars::sigma; // EDITABLE, still unstable
    double _rSigmaOut = Pars::ly * 2;                // EDITABLE, still unstable
    // --------------------------------------------------------------------------------- //
    double _sigmaMin = Pars::sigma * 1; // ! use from global variable
    double _sigmaMax = Pars::D_0 * 1;   // ! use from global variable

    // == reset variables
    _Dp.clear();
    // == assign new size of variables
    _Dp.resize(ntotal);

    for (int i = 0; i < ntotal; i++)
    {
        double r2bodycenter = std::sqrt(std::pow(xp[i] - bodycenter[0], 2) + std::pow(yp[i] - bodycenter[1], 2));
        if (r2bodycenter < _rSigmaIn)
        {
            _Dp[i] = _sigmaMin;
        }
        else if (r2bodycenter >= _rSigmaIn && r2bodycenter < _rSigmaOut)
        {
            double _r = r2bodycenter - _rSigmaIn;
            double _L = _rSigmaOut - _rSigmaIn;
            _Dp[i] = (_r / _L) * (_sigmaMax - _sigmaMin) + _sigmaMin;
        }
        else
        {
            _Dp[i] = _sigmaMax;
        }
    }
    delete[] bodycenter;
}

void mps_adaption::particle_merging(Particle &particle, const Body &body)
{
    // vector<double> &xp = particle.x;
    // vector<double> &yp = particle.y;
    // vector<double> &sp = particle.s;
    vector<double> xp{};
    vector<double> yp{};
    vector<double> sp{};

    vector<double> _xp = particle.x;
    vector<double> _yp = particle.y;
    vector<double> _sp = particle.s;

    // == reset variables
    _isMerged.clear();
    // == assign new size of variables
    _isMerged.resize(_sp.size(), false);

    // generating neighorhood
    std::vector<std::vector<int>> _par2gridNeighbor = _base_remeshing.inter_search(particle, particle); // neighbor : particle --> grid

    // internal variables
    std::vector<int> _pair;
    std::vector<double> _xpMerge;
    std::vector<double> _ypMerge;
    std::vector<double> _spMerge;

    const double alphaT = 1.50e0; // ! declare as class private STATIC member
    int nMerged = 0;

    for (size_t i = 0; i < _sp.size(); i++)
    {
    lv1Loop:
        _pair = _par2gridNeighbor[i];
        if (_Dp[i] > (sqrt(2.0e0) * _sp[i]) && _isMerged[i] == false)
        {
            int _idxjMin = -1;
            double _rijMin = Pars::ly;
            for (size_t j = 0; j < _pair.size(); j++)
            {
                int _idxi = i;
                int _idxj = _pair[j];
                if (_idxi != _idxj && _isMerged[_idxj] == false)
                {
                    double _rij = std::sqrt(std::pow(_xp[_idxj] - _xp[_idxi], 2) + std::pow(_yp[_idxj] - _yp[_idxi], 2));

                    int _idxjTemp = _idxjMin;
                    double _rijTemp = _rijMin;
                    _rijMin = _rijTemp > _rij ? _rij : _rijTemp;
                    _idxjMin = _rijTemp > _rij ? _idxj : _idxjTemp;
                }
                else
                {
                    continue;
                }
            }

            if (_idxjMin != -1)
            {
                double _Li = _sp[i];
                double _Lj = _sp[_idxjMin];
                if (_rijMin < (alphaT * (_Li + _Lj) * 0.5) && _isMerged[_idxjMin] == false)
                {
                    _isMerged[i] = true;
                    _isMerged[_idxjMin] = true;
                    double _xpNew = (std::pow(_Li, 2) * _xp[i] + std::pow(_Lj, 2) * _xp[_idxjMin]) / (std::pow(_Li, 2) + std::pow(_Lj, 2));
                    double _ypNew = (std::pow(_Li, 2) * _yp[i] + std::pow(_Lj, 2) * _yp[_idxjMin]) / (std::pow(_Li, 2) + std::pow(_Lj, 2));
                    double _spNew = std::sqrt(std::pow(_Li, 2) + std::pow(_Lj, 2));

                    _xpMerge.push_back(_xpNew);
                    _ypMerge.push_back(_ypNew);
                    _spMerge.push_back(_spNew);

                    nMerged++;

                    if (i < _sp.size() - 1)
                    {
                        i++;
                        goto lv1Loop;
                    }
                    else
                    {
                        goto endloop;
                    }
                    // break;
                }
            }
            else
            {
                _sp[i] = std::sqrt(std::pow(_sp[i], 2) + std::pow(_sp[i], 2));
                continue;
            }
        }
        else
        {
            continue;
        }
    }
endloop:

    std::cout << "merged particles " << nMerged << " of " << _sp.size() << std::endl;
    xp.clear();
    yp.clear();
    sp.clear();

    for (size_t i = 0; i < _sp.size(); i++)
    {
        if (_isMerged[i] == false)
        {
            xp.push_back(_xp[i]);
            yp.push_back(_yp[i]);
            sp.push_back(_sp[i]);
        }
    }

    for (size_t i = 0; i < _spMerge.size(); i++)
    {
        xp.push_back(_xpMerge[i]);
        yp.push_back(_ypMerge[i]);
        sp.push_back(_spMerge[i]);
    }

    particle.num = xp.size();
    particle.x = xp;
    particle.y = yp;
    particle.s = sp;
    particle.u.resize(particle.num, 0.0e0);
    particle.v.resize(particle.num, 0.0e0);
    particle.gz.resize(particle.num, 0.0e0);
    particle.neighbor.resize(particle.num, std::vector<int>());
    particle.isActive.resize(particle.num, false);
    //// cout << "\n\nnumber of particles: " << xp.size() << endl;
}

void mps_adaption::particle_moved(Particle &particle)
{
    static const double C_m = 5.0e-2; // ! user defined constant
    std::vector<double> _ci = this->get_packing_ratio(particle);
    std::vector<double> _dx(particle.num, 0.0e0);
    std::vector<double> _dy(particle.num, 0.0e0);
    std::vector<double> _Cshift(particle.num, 0.0e0);
    std::vector<double> _R2(particle.num, 0.0e0);

    std::vector<std::vector<int>> neighborlist = _base_remeshing.inter_search(particle, particle); // neighbor : particle --> grid

    for (size_t i = 0; i < particle.num; i++)
    {
        std::vector<int> _neighbor = neighborlist[i];
        for (size_t j = 0; j < _neighbor.size(); j++)
        {
            int _idxi = i;
            int _idxj = _neighbor[j];
            double _xij = particle.x[_idxj] - particle.x[_idxi];
            double _yij = particle.y[_idxj] - particle.y[_idxi];
            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2));
            double _Rij = (particle.s[_idxi] + particle.s[_idxj]) * Pars::r_scale * 0.5;
            double _wij = this->weighting_function(_rij, _Rij);
            double _Pi = _ci[_idxi];
            double _Pj = _ci[_idxj];

            if (_rij > 0 && _Pi > 0)
            {
                _dx[i] += _wij * C_m * (_rij - _Rij) * (_Pj / _Pi) * (_xij / _rij);
                _dx[i] += _wij * C_m * (_rij - _Rij) * (_Pj / _Pi) * (_yij / _rij);
            }
        }
        particle.x[i] += _dx[i];
        particle.y[i] += _dy[i];
    }
    // for (size_t i = 0; i < particle.num; i++)
    // {
    //     std::vector<int> _neighbor = neighborlist[i];
    //     for (size_t j = 0; j < _neighbor.size(); j++)
    //     {
    //         int _idxi = i;
    //         int _idxj = _neighbor[j];
    //         double _xij = particle.x[_idxj] - particle.x[_idxi];
    //         double _yij = particle.y[_idxj] - particle.y[_idxi];
    //         double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2));
    //         double _Rij = (particle.s[_idxi] + particle.s[_idxj]) * Pars::r_scale * 0.5;
    //         double _wij = this->weighting_function(_rij, _Rij);
    //         double _Vj = Pars::pi * std::pow(particle.s[_idxj], 2);

    //         if (_idxi != _idxj)
    //         {
    //             _Cshift[i] += _wij * _Vj;
    //         }
    //     }
    //     _R2[i] = std::pow(particle.s[i], 2);
    // }

    // LSMPSa _lsmpsa;
    // _lsmpsa.set_LSMPS(particle.x, particle.y, particle.s, _Cshift, neighborlist);
    // _dx = _lsmpsa.get_ddx();
    // _dy = _lsmpsa.get_ddy();

    // const static double CSHIFT = 2.0e1;
    // for (size_t i = 0; i < particle.num; i++)
    // {
    //     particle.x[i] -= CSHIFT * _R2[i] * _dx[i];
    //     particle.y[i] -= CSHIFT * _R2[i] * _dy[i];
    // }
}

double mps_adaption::weighting_function(const double &rij, const double &Rij)
{
    if (rij <= Rij)
    {
        return std::pow(1 - rij / Rij, 2);
    }
    else
    {
        return 0.0e0;
    }
}
#pragma endregion

#pragma region private_methods_surface - detection
std::vector<double> mps_adaption::get_packing_ratio(const Particle &particle)
{
    const int R_eff = 2.1;  // ! effective radius, SHOULD be move to class' STATIC member
    const int R_surf = 2.1; // ! effective radius, SHOULD be move to class' STATIC member
    std::vector<double> ci(particle.num);
    // generating neighorhood // ! neighborscale: SHOULD be 2.1
    std::vector<std::vector<int>> neighborlist = _base_remeshing.inter_search(particle, particle); // neighbor : particle --> grid
    // std::vector<std::vector<int>> neighborlist = d_neighbor.link_list(particle.num, particle.s, particle.x, particle.y, R_surf);
    // std::vector<std::vector<int>> neighborlist = particle.neighbor;

    for (size_t i = 0; i < particle.num; i++)
    {
        std::vector<int> _neighbor = neighborlist[i];
        double _Ni = 0;
        double _Vr = 0;

        for (size_t j = 0; j < _neighbor.size(); j++)
        {
            int idxi = i;
            int idxj = _neighbor[j];
            double _xij = particle.x[idxj] - particle.x[idxi];
            double _yij = particle.y[idxj] - particle.y[idxi];
            double _rij = std::sqrt(std::pow(_xij, 2) + std::pow(_yij, 2));
            double _R = ((particle.s[idxj] + particle.s[idxi]) * 0.5) * R_eff;
            double _Lj = particle.s[idxj];
            double _Rj = particle.s[idxj] * R_eff;

            double _Wij = 0.0e0;
            if (_R < _rij - 0.5 * _Lj)
            {
                _Wij = 0.0e0;
            }
            if (_R >= (_rij - 0.5 * _Lj) && _R <= (_rij + 0.5 * _Lj))
            {
                _Wij = (_R - _rij + _Lj * 0.5) * (std::pow(_Lj, 2)) / _Lj;
            }
            if (_R > _rij + 0.5 * _Lj)
            {
                _Wij = (std::pow(_Lj, 2));
            }

            if (idxi == idxj)
            {
                continue;
            }
            else
            {
                _Ni = _Ni + _Wij;
                // _Vr = _Vr + Pars::pi * std::pow(_R, 2);
            }
        }
        _Vr = Pars::pi * std::pow(particle.s[i] * R_eff, 2);
        ci[i] = _Ni / _Vr;
    }

    return ci;
}
#pragma endregion