#include "particle_adaption.hpp"

// ===================
// ===================
void particle_adaption::removal()
{
	// == internal variables
	double dr, minDpq;
	std::vector<int> pair;
	std::vector<double> Dpq(2, 0.0e0);

	// == initialization
	int _nfused = 0;
	this->_np = this->_xp.size();

	this->_isFused.clear();
	this->_isFused.resize(this->_np, false);
	// == performing romval
	printf("fusing particles ... ");
	for (int i = 0; i < _np; i++)
	{
		if (_isFused[i] == false)
		{
			// _NeighborSearchAlgorithm.DirectFindMR(i, _np, _sp, _xp, _yp, _isFused, pair);
			_directFindMR(i, _np, _sp, _xp, _yp, _isFused, pair);
			for (int j = 0; j < pair.size(); j++)
			{
				dr = std::sqrt(std::pow((_xp[pair[j]] - _xp[i]), 2) + std::pow((_yp[pair[j]] - _yp[i]), 2));
				Dpq[0] = _sp[i] / Pars::r_scale;
				Dpq[1] = _sp[pair[j]] / Pars::r_scale;
				minDpq = *min_element(Dpq.begin(), Dpq.end());
				if (dr < minDpq * 0.51)
				{
					_nfused++;
					_isFused[pair[j]] = true;
				}
			}
			pair.clear(); // reset pair
		}
	}
	printf("%d particles has been fused\n", _nfused);
	// == update particle <- need to be edited !!!!!!!!!!, for single iteration still ok
	// -- storing temporary data
	_xptemp = _xp;
	_yptemp = _yp;
	_sptemp = _sp;

	// -- reset to be empty vector
	_xp.clear();
	_yp.clear();
	_sp.clear();
	// -- assign value
	for (int i = 0; i < this->_np; i++)
	{
		if (_isFused[i] == false)
		{
			_xp.push_back(this->_xptemp[i]);
			_yp.push_back(this->_yptemp[i]);
			_sp.push_back(this->_sptemp[i]);
		}
	}
}
// ===================
// ===================
void particle_adaption::insertion()
{
	// == internal variables
	double tetha, Dpq, _xadd, _yadd;
	int n_neighbor = 0;
	std::vector<int> pair; // neighbor of particle
	double d_random[2];	// distance of random particles to be added

	this->_np = this->_xp.size();	// current number of particles
	this->_npnew = this->_xp.size(); // initializing npnew.
	_isFused.resize(_np, false);

	int ninsert = 0;

	// int l = Pars::l_2_0_2;
	int l = 6; // !

	printf("inserting particles ... ");
	for (int i = 0; i < _npnew /*_np*/; i++)
	{
		// _NeighborSearchAlgorithm.DirectFindMR(i, _npnew, _sp, _xp, _yp, _isFused, pair);
		_directFindMR(i, _npnew, _sp, _xp, _yp, _isFused, pair);
		n_neighbor = pair.size();

		if (n_neighbor < l) // perform do-while so that at the end of this iteration each particle has at least l_2_0_2 neighbour
		{
			// int numadd = l - n_neighbor;
			int numadd = 2;
			for (int nadd = 0; nadd < numadd; nadd++)
			{
				Dpq = _sp[i] / Pars::r_scale;
				do
				{
					tetha = (double)(std::rand() % 360) / 180; // !!!!!
					d_random[0] = 1.0 * Dpq * std::cos(tetha * M_PI);
					d_random[1] = 1.0 * Dpq * std::sin(tetha * M_PI);
					_xadd = _xp[i] + d_random[0];
					_yadd = _yp[i] + d_random[1];
				} while (
					_xadd <= _xcorner[0] || _xadd >= _xcorner[1] ||
					_yadd <= _ycorner[0] || _yadd >= _ycorner[3]);

				_xp.push_back(_xadd);
				_yp.push_back(_yadd);
				_sp.push_back(_sp[i]);

				ninsert += 1;
				_isFused.push_back(false);
			}
		}
		_npnew = this->_xp.size(); // update number of particles
		pair.clear();			   // !!!!!
	}
	printf("%d particles has been inserted\n", ninsert);
}
// ===================
// ===================
