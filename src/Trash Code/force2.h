/*
 *  Created by Van Luc Nguyen
 *  Copyright 2020. All rights reserved.
 */


#ifndef FORCE2_H                 // header guards
#define FORCE2_H
#include<vector>
#include "global.h"

using namespace std;
                    
extern vector<double> u_jmin1;                                
extern vector<double> u_jmax1;	
extern vector<double> v_imin1;                                
extern vector<double> v_imax1;	

extern vector<double> u_jmin2;                                
extern vector<double> u_jmax2;	
extern vector<double> v_imin2;                                
extern vector<double> v_imax2;

extern vector<double> omegasumX2;
extern vector<double> omegasumY2;

class force2
{

	public: 
		force2();
		void dragForce2();
		

	protected:

	private:


};

#endif



