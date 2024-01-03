#ifndef INCLUDED_SAVE_DATA_BASE
#define INCLUDED_SAVE_DATA_BASE
#include <iomanip> // std::setw, std::setfill 
#include <sstream> // std::stringstream
#include <fstream> // std::c_str()

using namespace std;

class base_save_data
{
	#define w20 std::setw(20) // spare width while saving data
private:
	std::vector<double> vdif;
	std::vector<double> udif;
public:
	void force_pen(int it, Particle p, double x_loc[2]);
};

#endif
