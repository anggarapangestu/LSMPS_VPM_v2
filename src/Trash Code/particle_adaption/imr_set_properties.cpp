#include "particle_adaption.hpp"

// ===================
// ===================
particle_adaption::particle_adaption()
{
	//
}
// ===================
// ===================
void particle_adaption::set_properties(const int np, const std::vector<double> &xp, const std::vector<double> &yp,
									   const std::vector<double> &sp, const std::vector<double> &fp)
{
	this->_xp.clear();
	this->_yp.clear();
	this->_sp.clear();
	this->_fp.clear();
	
	this->_xp = xp;
	this->_yp = yp;
	this->_sp = sp;
	this->_fp = fp;

	// -- domain corner boundary
	this->_xcorner.resize(4);
	this->_ycorner.resize(4);

	this->_xcorner[0] = *std::min_element(_xp.begin(), _xp.end());
	this->_xcorner[1] = *std::max_element(_xp.begin(), _xp.end());
	this->_xcorner[2] = _xcorner[1];
	this->_xcorner[3] = _xcorner[0];

	this->_ycorner[0] = *std::min_element(_yp.begin(), _yp.end());
	this->_ycorner[3] = *std::max_element(_yp.begin(), _yp.end());
	this->_ycorner[1] = _ycorner[0];
	this->_ycorner[2] = _ycorner[3];
}
// ===================
// ===================
particle_adaption::~particle_adaption()
{
	//
}
