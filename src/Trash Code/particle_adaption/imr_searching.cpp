#include "particle_adaption.hpp"

// bool particle_adaption::sortcol( const vector<double>& v1, const vector<double>& v2 )
// {
//   return v1[1] < v2[1];
// }
// =====================
// =====================
void particle_adaption::_directFind(double xpi, double ypi, double hsmli,
                                    const int ntotal, const std::vector<double> &hsml,
                                    const std::vector<double> &xp, const std::vector<double> &yp, std::vector<int> &pair)
{
    // == clearing pair
    pair.clear();
    // internal variables
    double driac, mhsml;

    for (int j = 0; j < ntotal; j++)
    {
        driac = std::pow((xpi - xp[j]), 2) + std::pow((ypi - yp[j]), 2);
        if (std::sqrt(driac) < hsmli && std::sqrt(driac) > 1.0e-6)
        {
            pair.push_back(j);
        }
    }
}
// =====================
// =====================
void particle_adaption::_directFindMR(int i, const int ntotal, const std::vector<double> &hsml,
                                      const std::vector<double> &xp, const std::vector<double> &yp, const std::vector<bool> &isFused,
                                      std::vector<int> &pair)
{
    // == clearing pair
    pair.clear();
    // internal variables
    double driac, mhsml;

    for (int j = 0; j < ntotal; j++)
    {
        if ((isFused[j] == false) && (j != i))
        {
            driac = std::pow((xp[i] - xp[j]), 2) + std::pow((yp[i] - yp[j]), 2);
            mhsml = hsml[i] < hsml[j] ? hsml[i] : hsml[j];
            if (std::sqrt(driac) < mhsml)
            {
                pair.push_back(j);
            }
        }
    }
}
// =====================
// =====================
void particle_adaption::_directFindMRinterp(double xpi, double ypi, double hsmli,
                                            const int ntotal, const std::vector<double> &hsml,
                                            const std::vector<double> &xp, const std::vector<double> &yp, std::vector<int> &pair)
{
    // == clearing pairs
    pair.clear();
    // internal variables
    double driac, mhsml;
    // std::vector<int> _pair;
    // std::vector<double> _dr;

    for (int j = 0; j < ntotal; j++)
    {
        driac = std::sqrt(std::pow((xpi - xp[j]), 2) + std::pow((ypi - yp[j]), 2));
        mhsml = hsmli < hsml[j] ? hsmli : hsml[j];
        if (driac < mhsml /*&& driac>1e-6*/)
        {
            // _pair.push_back(j);
            // _dr.push_back(driac);
            pair.push_back(j);
        }
    }

    // std::vector<std::vector<double>> _sortpair(_pair.size(), std::vector<double>(2));
    // for (int i = 0; i < _pair.size(); i++)
    // {
    //   _sortpair[i][0] = _pair[i];
    //   _sortpair[i][1] = _dr[i];
    // }

    // std::sort(_sortpair.begin(), _sortpair.end(), customLess);

    // int _size = _pair.size() < 15/*2*Pars::l_2_0_2*/ ? _pair.size() : 15/*2*Pars::l_2_0_2*/;
    // pair.resize(_size);
    // for (int i = 0; i < _size; i++)
    // {
    //   pair[i] = _sortpair[i][0];
    //   // printf("%f ", _sortpair[i][1]);
    // }
    // // printf("\n\n");
}
// =====================
// =====================

void particle_adaption::_directFindMR(int ntotal, const std::vector<double> &hsml,
                                      const std::vector<double> &xp, const std::vector<double> &yp,
                                      std::vector<int> &pair_i, std::vector<int> &pair_j)
{
    /*************************************************************************
    **    Subroutine to calculate the smoothing funciton for each particle and
    **    the interaction parameters used by the SPH algorithm. Interaction 
    **    pairs are determined by directly comparing the particle distance 
    **    with the corresponding smoothing length.
    **
    **      itimestep : Current time step                                 [in]
    **      ntotal    : Number of particles                               [in]
    **      hsml      : Smoothing Length                                  [in]
    **      x         : Coordinates of all particles                      [in]
    **      niac      : Number of interaction pairs                      [out]
    **      pair_i    : List of first partner of interaction pair        [out]
    **      pair_j    : List of second partner of interaction pair       [out]
    **      countiac  : Number of neighboring particles                  [out]
    *************************************************************************/
    // clearing variables
    pair_i.clear();
    pair_j.clear();
    // internal variables
    // int sumiac, maxiac, miniac, noiac, maxp, minp;
    double driac, r, mhsml;
    // int countiac[ntotal]{0};

    for (int i = 0; i < ntotal - 1; i++)
    {
        for (int j = i + 1; j < ntotal; j++)
        {
            driac = std::pow((xp[i] - xp[j]), 2) + std::pow((yp[i] - yp[j]), 2);
            mhsml = hsml[i] < hsml[j] ? hsml[i] : hsml[j];
            if (std::sqrt(driac) < mhsml * Pars::r_scale)
            {
                pair_i.push_back(i);
                pair_j.push_back(j);
                // countiac[i] += 1;
                // countiac[j] += 1;
            }
        }
    }

#pragma region statistic
// // statisticas for the interaction
// sumiac = 0;
// maxiac = 0;
// miniac = 1000;
// noiac = 0;
// maxp = 0;
// minp = 0;
// for (int i = 0; i < ntotal; i++)
// {
//     sumiac = sumiac + countiac[i];
//     if (countiac[i] > maxiac)
//     {
//         maxiac = countiac[i];
//         maxp = i;
//     }
//     else if (countiac[i] < miniac)
//     {
//         miniac = countiac[i];
//         minp = i;
//     }
//     else if (countiac[i] == 0)
//     {
//         noiac = noiac + 1;
//     }
// }

// printf("\n >> Statistics: interactions per particle:");
// printf("\n**** Particle: %d maximal interactions: %d", maxp, maxiac);
// printf("\n**** Particle: %d minimal interactions: %d", minp, miniac);
// printf("\n**** Average : %f", float(sumiac) / float(ntotal) );
// printf("\n**** Total pairs : %d", niac);
// printf("\n**** Particles with no interactions: %d \n", noiac);
// printf("\n");
#pragma endregion
}
// =====================
// =====================