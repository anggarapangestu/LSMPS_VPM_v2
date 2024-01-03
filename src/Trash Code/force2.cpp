/*
 *  Created by Van Luc Nguyen
 *  Copyright 2020. All rights reserved.
 */

#include <omp.h>
#include "global.h"
#include <vector>
#include <math.h> 
#include <cmath> 
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <string>
#include "grid.h"
#include "force2.h"

using namespace std;


vector<double> u_jmin1(nx,0.0);                                  
vector<double> u_jmax1(nx,0.0);  	
vector<double> v_imin1(ny,0.0);                                 
vector<double> v_imax1(ny,0.0);  	

vector<double> u_jmin2(nx,0.0);                                  
vector<double> u_jmax2(nx,0.0);  	
vector<double> v_imin2(ny,0.0);                                 
vector<double> v_imax2(ny,0.0); 

vector<double> omegasumX2(Ninter,0.0);
vector<double> omegasumY2(Ninter,0.0);

force2::force2()
{

}

void force2::dragForce2()
{
	int imin, imax, jmin, jmax;
    double x_distance = 2.0;
    double y_distance = 2.0;
    double x_center1 = 0.0;
    double y_center1 = 0.0;

	for (int i=0;i<nx;i++){ 
		if (X[i]>=x_center1 -x_distance-dx && X[i]<=x_center1 -x_distance)            
            imin=i;
		else if (X[i]>=x_center1 +x_distance && X[i]<=x_center1 +x_distance +dx)      
            imax=i;
	}

	for (int j=0;j<ny;j++){
		if (Y[j] >= y_center1 -y_distance-dy  && Y[j] <= y_center1 -y_distance)             
            jmin=j;
		else if (Y[j] >= (y_center1 +y_distance) && Y[j] <= (y_center1 +y_distance+ dy) )   
            jmax=j;
	}

	double A1A2=0.0, A2A3=0.0, A3A4=0.0, A1A4=0.0;
	double B1B2=0.0, B2B3=0.0, B3B4=0.0, B1B4=0.0;

	double A1A2y=0.0, A2A3y=0.0, A3A4y=0.0, A1A4y=0.0;
	double B1B2y=0.0, B2B3y=0.0, B3B4y=0.0, B1B4y=0.0;

	double gammaX2=0;
	double gammaY2=0;
	
     for (int i=imin;i<=imax;i++)
     {
		for (int j=jmin;j<=jmax;j++){
			gammaX2=gammaX2+(ug[i][j]*dx*dy);
			gammaY2=gammaY2+(vg[i][j]*dx*dy);
		}
	}
	
	omegasumX2[loop+1]=gammaX2;
	omegasumY2[loop+1]=gammaY2;

    #pragma omp parallel for reduction (+:A1A2,A3A4,A1A2y,A3A4y)
	for (int i=imin;i<=imax;i++)
	{

		//----------------------------------------------for Fx------------------------------------------------------//
		double a1=ug[i][jmin]*vg[i][jmin]+vg[i][jmin]*vortg[i][jmin]*Y[jmin];
		double a2=-Y[jmin]*(ug[i][jmin]-u_jmin1[i])/(2.0*dt);
		double a3= 1.0/re*(2.0*(ug[i-1][jmin]-2.0*ug[i][jmin]+ug[i+1][jmin])/dx/dx+\
		      		          (ug[i][jmin-1]-2.0*ug[i][jmin]+ug[i][jmin+1])/dy/dy+\
			  				  (vg[i+1][jmin+1]-vg[i+1][jmin-1]-vg[i-1][jmin+1]+vg[i-1][jmin-1])/(4.0*dx*dy))*Y[jmin];
		double a4=-1.0/re*((ug[i][jmin+1]-ug[i][jmin-1])/(2.0*dy)+(vg[i+1][jmin]-vg[i-1][jmin])/(2.0*dx));
		A1A2=A1A2+(a1+a2+a3+a4)*dx;

		double b1=-ug[i][jmax]*vg[i][jmax]-vg[i][jmax]*vortg[i][jmax]*Y[jmax];
		double b2=Y[jmax]*(ug[i][jmax]-u_jmax1[i])/(2.0*dt);
		double b3=-1.0/re*(2.0*(ug[i-1][jmax]-2.0*ug[i][jmax]+ug[i+1][jmax])/(dx*dx)+\
		       			 	  (ug[i][jmax-1]-2.0*ug[i][jmax]+ug[i][jmax+1])/(dy*dy)+\
			   				  (vg[i+1][jmax+1]-vg[i+1][jmax-1]-vg[i-1][jmax+1]+vg[i-1][jmax-1])/(4.0*dx*dy))*Y[jmax];
		double b4=1.0/re*((ug[i][jmax+1]-ug[i][jmax-1])/(2.0*dy)+(vg[i+1][jmax]-vg[i-1][jmax])/(2.0*dx));  
		A3A4= A3A4+(b1+b2+b3+b4)*dx;

		//----------------------------------------------for Fy------------------------------------------------------//

		double a1y=0.5*(vg[i][jmin]*vg[i][jmin]-ug[i][jmin]*ug[i][jmin])-vg[i][jmin]*vortg[i][jmin]*X[i];
		double a2y=X[i]*(ug[i][jmin]-u_jmin1[i])/(2.0*dt);
		double a3y=-1.0/re*(2.0*(ug[i-1][jmin]-2.0*ug[i][jmin]+ug[i+1][jmin])/dx/dx+\
		      		          (ug[i][jmin-1]-2.0*ug[i][jmin]+ug[i][jmin+1])/dy/dy+\
			  				  (vg[i+1][jmin+1]-vg[i+1][jmin-1]-vg[i-1][jmin+1]+vg[i-1][jmin-1])/(4.0*dx*dy))*X[i];
		double a4y=-2.0/re*(vg[i][jmin+1]-vg[i][jmin-1])/(2.0*dy);
		A1A2y=A1A2y+(a1y+a2y+a3y+a4y)*dx;

		double b1y=-0.5*(vg[i][jmax]*vg[i][jmax]-ug[i][jmax]*ug[i][jmax])+vg[i][jmax]*vortg[i][jmax]*X[i];
		double b2y=-X[i]*(ug[i][jmax]-u_jmax1[i])/(2.0*dt);
		double b3y=1.0/re*(2.0*(ug[i-1][jmax]-2.0*ug[i][jmax]+ug[i+1][jmax])/dx/dx+\
		      		          (ug[i][jmax-1]-2.0*ug[i][jmax]+ug[i][jmax+1])/dy/dy+\
			  				  (vg[i+1][jmax+1]-vg[i+1][jmax-1]-vg[i-1][jmax+1]+vg[i-1][jmax-1])/(4.0*dx*dy))*X[i];
		double b4y=2.0/re*(vg[i][jmax+1]-vg[i][jmax-1])/(2.0*dy);
		A3A4y= A3A4y+(b1y+b2y+b3y+b4y)*dx;
	}

    #pragma omp parallel for reduction (+:A2A3,A1A4,A2A3y,A1A4y)
	for (int j=jmin;j<=jmax;j++)
	{
		//--------------------------------------------for Fx---------------------------------------------------------//
		double c1=0.5*(vg[imax][j]*vg[imax][j]-ug[imax][j]*ug[imax][j]);
		double c2=-ug[imax][j]*vortg[imax][j]*Y[j]-Y[j]*(vg[imax][j]-v_imax1[j])/(2.0*dt);
		double c3=1.0/re*(2.0*(vg[imax][j-1]-2.0*vg[imax][j]+vg[imax][j+1])/(dy*dy)+\
						     (vg[imax-1][j]-2.0*vg[imax][j]+vg[imax+1][j])/(dx*dx)+\
					         (ug[imax+1][j+1]-ug[imax+1][j-1]-ug[imax-1][j+1]+ug[imax-1][j-1])/(4.0*dx*dy) )*Y[j];
		double c4=2.0/re*(ug[imax+1][j]-ug[imax-1][j])/(2.0*dx);
		A2A3=A2A3+ (c1+c2+c3+c4)*dy;

		double d1=-0.5*(vg[imin][j]*vg[imin][j]-ug[imin][j]*ug[imin][j]);
		double d2=ug[imin][j]*vortg[imin][j]*Y[j]+Y[j]*(vg[imin][j]-v_imin1[j])/(2.0*dt);
		double d3=-1.0/re*(2.0*(vg[imin][j-1]-2.0*vg[imin][j]+vg[imin][j+1])/(dy*dy)+\
						      (vg[imin-1][j]-2.0*vg[imin][j]+vg[imin+1][j])/(dx*dx)+\
					          (ug[imin+1][j+1]-ug[imin+1][j-1]-ug[imin-1][j+1]+ug[imin-1][j-1])/(4.0*dx*dy))*Y[j];
		double d4=-2.0/re*(ug[imin+1][j]-ug[imin-1][j])/(2.0*dx);
		A1A4=A1A4+ (d1+d2+d3+d4)*dy;

		//-------------------------------------------for Fy----------------------------------------------------------//

		double c1y=-ug[imax][j]*vg[imax][j]+ug[imax][j]*vortg[imax][j]*X[imax];
		double c2y=X[imax]*(vg[imax][j]-v_imax1[j])/(2.0*dt);
		double c3y=-1.0/re*(2.0*(vg[imax][j-1]-2.0*vg[imax][j]+vg[imax][j+1])/(dy*dy)+\
						     (vg[imax-1][j]-2.0*vg[imax][j]+vg[imax+1][j])/(dx*dx)+\
					         (ug[imax+1][j+1]-ug[imax+1][j-1]-ug[imax-1][j+1]+ug[imax-1][j-1])/(4.0*dx*dy) )*X[imax];
		double c4y=1.0/re*((vg[imax+1][j]-vg[imax-1][j])/(2.0*dx)\
						 +(ug[imax][j+1]-ug[imax][j-1])/(2.0*dy));
		A2A3y=A2A3y+ (c1y+c2y+c3y+c4y)*dy;

		double d1y=ug[imin][j]*vg[imin][j]-ug[imin][j]*vortg[imin][j]*X[imin];
		double d2y=-X[imin]*(vg[imin][j]-v_imin1[j])/(2.0*dt);
		double d3y=1.0/re*(2.0*(vg[imin][j-1]-2.0*vg[imin][j]+vg[imin][j+1])/(dy*dy)+\
						     (vg[imin-1][j]-2.0*vg[imin][j]+vg[imin+1][j])/(dx*dx)+\
					         (ug[imin+1][j+1]-ug[imin+1][j-1]-ug[imin-1][j+1]+ug[imin-1][j-1])/(4.0*dx*dy) )*X[imin];
		double d4y=-1.0/re*((vg[imin+1][j]-vg[imin-1][j])/(2.0*dx)\
						 +(ug[imin][j+1]-ug[imin][j-1])/(2.0*dy));
		A1A4y=A1A4y+ (d1y+d2y+d3y+d4y)*dy;
	}


	for (int i=0;i<nx;i++)
	{
		u_jmin1[i]=u_jmin2[i];
		u_jmax1[i]=u_jmax2[i];
		u_jmin2[i]=ug[i][jmin];
		u_jmax2[i]=ug[i][jmax];
	}


	for (int j=0;j<ny;j++)
	{
		v_imin1[j]=v_imin2[j];
		v_imax1[j]=v_imax2[j];
		v_imin2[j]=vg[imin][j];
		v_imax2[j]=vg[imax][j];
	}


	double FX=A1A2+A2A3+A3A4+A1A4-(omegasumX2[loop+1]-omegasumX2[loop-1])/(2.0*dt);
	double FY=A1A2y+A2A3y+A3A4y+A1A4y-(omegasumY2[loop+1]-omegasumY2[loop-1])/(2.0*dt);
  	
	if (loop==0)
	{
		Cd2[loop] = 0.0;
		Cl2[loop] = 0.0;
	}
	else  
	{
		Cd2[loop] = 2.0*FX/((U_inf*U_inf)*D);
		Cl2[loop] = -2.0*FY/((U_inf*U_inf)*D);
	}

}










