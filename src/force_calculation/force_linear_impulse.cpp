#include "force_calc.hpp"

// TODO: Using Force Linear Impulse
// omega: vorticity, gamma vortex strength
// F/rho = - dI/dt
// ---> Fd= -rho * (Isumx^(it+1)-Isumx^(it-1))/(2*delta_t)
// but, I = 1/(N-1) *  integral of (x X omega dV)   (Rasmussen 2011) with dV = h^3,  but, omega*h^3=Gamma  (h=core size)
// --> Isumx = 1/2 * SUM(y*gamma_z-z*gamma_y), Isumx,Isumx in the same manner.
// or as introduce by KOUMOUTSAKOS 1995: d/dt  * SUM (i=1-->N) of (GAMMAi xi  X  e^_z)
// Cd=Fd/(1/2 *rho*U^2 *A), A: reference area, Sphere= pi*R^2 =pi*D^2/4  , FlatPlate with normal velocity: Length*Chord
// ----> Cd=-(Isumx^(it+1)-Isumx^(it-1))/(dt*U^2 * A)
// Input are: Particles location, number of Particles, Particles vortex Strength
// Output: Force and Coefficient of Forces
// =======================================================

void force_calculation::set_variables(int numberofparticles_x, int numberofparticles_y, int totalnumberofparticles){
    this->u_jmin1.resize(numberofparticles_x,Pars::u_inf);                                  
    this->u_jmax1.resize(numberofparticles_x,Pars::u_inf);  	
    this->v_imin1.resize(numberofparticles_y,Pars::v_inf);                                 
    this->v_imax1.resize(numberofparticles_y,Pars::v_inf);  	

    this->u_jmin2.resize(numberofparticles_x,Pars::u_inf);                                  
    this->u_jmax2.resize(numberofparticles_x,Pars::u_inf);  	
    this->v_imin2.resize(numberofparticles_y,Pars::v_inf);                                 
    this->v_imax2.resize(numberofparticles_y,Pars::v_inf);

    //this->_omegasumx.resize(totalnumberofparticles, 0.0);
    //this->_omegasumy.resize(totalnumberofparticles, 0.0);
}

void force_calculation::insert_data(const Particle &p, int i){
    if (DIM == 2){
        this->_p.x.push_back(p.x[i]);
        this->_p.y.push_back(p.y[i]);
        this->_p.s.push_back(p.s[i]);
        this->_p.u.push_back(p.u[i]);
        this->_p.v.push_back(p.v[i]);
        this->_p.gz.push_back(p.gz[i]);
        this->_p.neighbor.push_back(p.neighbor[i]); ///
        this->point_lists.push_back(i);
        this->_p.chi.push_back(p.chi[i]);
    }else{
        this->_p.x.push_back(p.x[i]);
        this->_p.y.push_back(p.y[i]);
        this->_p.z.push_back(p.z[i]); 
        this->_p.s.push_back(p.s[i]);
        this->_p.u.push_back(p.u[i]);
        this->_p.v.push_back(p.v[i]);
        this->_p.w.push_back(p.w[i]); 
        this->_p.gz.push_back(p.gz[i]);
        this->point_lists.push_back(i);
    }
}

void force_calculation::save_dat(Particle &p, int np, std::vector<int> &numbar, std::vector<int> &batas_bawah, std::vector<int> &batas_atas, 
				  std::vector<int> &batas_kiri, std::vector<int> &batas_kanan, std::string s, std::vector<double> &satu, std::vector<double> &dua,
				  std::vector<double> &tiga, std::vector<double> &empat, std::vector<double> &lima ){

    std::string name1;
    std::ofstream ofs;
    name1.append("Gridlike_");
    name1.append(s);
    name1.append(".csv");
	
    ofs.open(name1.c_str());
	ofs << "" << "x" << "," << "y" << "," << "u" << "," << "v" << "," << "size" << "," << "Strength" << "," << "neighbours" 
        << "," << "ddx" << "," << "ddy" << "," << "d2dx2" << "," << "d2dy2" << "," << "d2dxdy\n";

    for (int i = 0; i <np; i++){
        		ofs << "" << p.x[i]
					<< "," << p.y[i]
                    << "," << p.u[i]
                    << "," << p.v[i]
                    << "," << p.s[i]
                    << "," << p.gz[i]
                    << "," << p.neighbor[i].size()
                    << "," << satu[numbar[i]]
                    << "," << dua[numbar[i]]
                    << "," << tiga[numbar[i]]
                    << "," << empat[numbar[i]]
                    << "," << lima[numbar[i]]
					<< "\n";
    }

    ofs.close();
}


//=========================================================================================================================
void force_calculation::Force2(int iT, double A1, double A2, double A3, double A4, int nx, int ny, const Particle &p)
{
    //box size
	int imin, imax, jmin, jmax;
    double x_distance = 1.5; //coba 2 abis ini
    double y_distance = 2;
    double x_center1 = 0.0;
    double y_center1 = 0.0;
    double test1, test2;

    if (iT == 0){
        // test1 = omp_get_wtime();
		
        this->_omegasumx.clear();
        this->_omegasumy.clear();
        this->_omegasumx.push_back(0.0);
        this->_omegasumy.push_back(0.0);
        this->particles_inside.clear();

        int k = 0;
        for (int i = 0 ; i < p.num; i++){
            if (p.x[i] >= x_center1 - x_distance -  Pars::sigma && p.x[i] <= x_center1 + x_distance + Pars::sigma 
                && p.y[i] >= y_center1 - y_distance -  Pars::sigma && p.y[i] <= y_center1 + y_distance + Pars::sigma)
            {
                //konvensi
                //a1a2 = partikel bagian boundary bawah
                //a2a3 = partikel bagian boundary kanan
                //a3a4 = partikel bagian boundary atas
                //a4a1 = partikel bagian boundary kiri

                insert_data(p, i);
                this->particles_inside.push_back(k);

                if (p.y[i] >= y_center1 - y_distance - Pars::sigma && p.y[i] <= y_center1 - y_distance){
                    //this->a1a2.push_back(i);
                    this->a1a2.push_back(k);
                    
                    //kalo di bagian kanan atau kiri banget masukin juga ke a3a4/a2a3
                    if (p.x[i] <= x_center1 + x_distance + Pars::sigma && p.x[i] >= x_center1 + x_distance ){
                        this->a2a3.push_back(k); 
                    } 

                    if (p.x[i] >= x_center1 - x_distance - Pars::sigma && p.x[i] <= x_center1 - x_distance){
                        this->a4a1.push_back(k); 
                    }
                } else 

                //a3a4
                if (p.y[i] <= y_center1 + y_distance + Pars::sigma  && p.y[i] >= y_center1 + y_distance){
                    //this->a3a4.push_back(i);
                    this->a3a4.push_back(k);
                    
                    //kalo di bagian kanan atau kiri banget masukin juga ke a3a4/a2a3
                    if (p.x[i] <= x_center1 + x_distance + Pars::sigma && p.x[i] >= x_center1 + x_distance){
                        this->a2a3.push_back(k);
                    } 

                    if (p.x[i] >= x_center1 - x_distance - Pars::sigma && p.x[i] <= x_center1 - x_distance){
                        this->a4a1.push_back(k);
                    }
                } else

                //a2a3
                if (p.x[i] <= x_center1 + x_distance + Pars::sigma && p.x[i] >= x_center1 + x_distance){
                    this->a2a3.push_back(k);
                } else

                //a4a1
                if (p.x[i] >= x_center1 - x_distance - Pars::sigma && p.x[i] <= x_center1 - x_distance){
                    this->a4a1.push_back(k);
                }

                k++;
            }
        }

        //set the variables for calculation
        set_variables(this->a1a2.size(), this->a2a3.size(), this->particles_inside.size());
        // test2 = omp_get_wtime();
        printf("(Noca) step 1.0: %f s\n", test2-test1);
    }   
    
    if (iT != 0){
        // test1 = omp_get_wtime();
        //update velocity and particle strength
        #pragma omp parallel for
        for (size_t i = 0; i < this->point_lists.size();i++){
            this->_p.u[i] = p.u[point_lists[i]];
            this->_p.v[i] = p.v[point_lists[i]];
            this->_p.gz[i] = p.gz[point_lists[i]];
        }
        // test2 = omp_get_wtime();
        printf("(Noca) step 1.0: %f s\n", test2-test1);
    }
   
	double A1A2=0.0, A2A3=0.0, A3A4=0.0, A1A4=0.0;
	double B1B2=0.0, B2B3=0.0, B3B4=0.0, B1B4=0.0;

	double A1A2y=0.0, A2A3y=0.0, A3A4y=0.0, A1A4y=0.0;
	double B1B2y=0.0, B2B3y=0.0, B3B4y=0.0, B1B4y=0.0;

	double gammaX2=0.0;
	double gammaY2=0.0;

    // test1 = omp_get_wtime();
    #pragma omp parallel for reduction (+:gammaX2, gammaY2)
    for (size_t i = 0; i < this->particles_inside.size(); i++){
        int ii = i;
        gammaX2 = gammaX2 +  this->_p.u[ii] * std::pow(this->_p.s[ii],2); //dxdy using area define by multiresolution LSMPSVPM
        gammaY2 = gammaY2 +  this->_p.v[ii] * std::pow(this->_p.s[ii],2); 
    }
    // test2 = omp_get_wtime();
    printf("(Noca) step 2: %f s\n", test2-test1);
    
    this->_omegasumx.push_back(gammaX2);
	this->_omegasumy.push_back(gammaY2);

    //LSMPSA full domain
    LSMPSa lsmpsa_uu;
    LSMPSa lsmpsa_vv;
    // test1 = omp_get_wtime();
    //lsmpsa_uu.set_LSMPS(p.x, p.y, p.s, p.u, p.neighbor);
    //lsmpsa_vv.set_LSMPS(p.x, p.y, p.s, p.v, p.neighbor);
    lsmpsa_uu.set_LSMPS(this->_p.x, this->_p.y, this->_p.s, this->_p.u, p.x, p.y, p.s, p.u, this->_p.neighbor);
    lsmpsa_vv.set_LSMPS(this->_p.x, this->_p.y, this->_p.s, this->_p.v, p.x, p.y, p.s, p.v, this->_p.neighbor);
    // test2 = omp_get_wtime();
    printf("(Noca) step 3: %f s\n", test2-test1);
    
    std::vector<double> d2udx2 = lsmpsa_uu.get_d2d2x();
	std::vector<double> d2udy2 = lsmpsa_uu.get_d2d2y();
    std::vector<double> d2udxdy = lsmpsa_uu.get_d2dxdy();
    std::vector<double> dudy = lsmpsa_uu.get_ddy();
    std::vector<double> dudx = lsmpsa_uu.get_ddx();

    std::vector<double> d2vdx2 = lsmpsa_vv.get_d2d2x();
	std::vector<double> d2vdy2 = lsmpsa_vv.get_d2d2y();
    std::vector<double> d2vdxdy = lsmpsa_vv.get_d2dxdy();
    std::vector<double> dvdy = lsmpsa_vv.get_ddy();
    std::vector<double> dvdx = lsmpsa_vv.get_ddx();
    
    std::ofstream ofs;
   
    //#pragma omp parallel for reduction (+:A1A2,A3A4,A1A2y,A3A4y)
    /*for (size_t i = 0; i < this->a1a2.size(); i++)
    {
		//----------------------------------------------for Fx------------------------------------------------------//
        int i_b = this->a1a2[i]; //boundary bawah
        int i_a = this->a3a4[i]; //boundary atas
        
        double a1 = (this->_p.u[i_b] * this->_p.v[i_b]) + (this->_p.v[i_b] * (this->_p.gz[i_b] / std::pow(this->_p.s[i_b],2)) * this->_p.y[i_b]);
        double a2 = -(this->_p.y[i_b]) * (this->_p.u[i_b] - this->u_jmin1[i]) / (2 * Pars::dt);
        double a3 = 1.0 / Pars::Re * (2.0 * d2udx2[this->point_lists[i_b]] + d2udy2[this->point_lists[i_b]] + d2vdxdy[this->point_lists[i_b]]) * this->_p.y[i_b];
        double a4 = -1.0/ Pars::Re * (dudy[this->point_lists[i_b]] + dvdx[this->point_lists[i_b]]);
        A1A2 = A1A2+(a1+a2+a3+a4)*(this->_p.s[i_b]);

        double b1 = -this->_p.u[i_a] * this->_p.v[i_a] - this->_p.v[i_a] * (this->_p.gz[i_a] / std::pow(this->_p.s[i_a],2)) * this->_p.y[i_a];
        double b2 = this->_p.y[i_a] * ( this->_p.u[i_a] - u_jmax1[i]) / (2 * Pars::dt);
        double b3 = - 1.0 / Pars::Re * (2.0 * d2udx2[this->point_lists[i_a]] + d2udy2[this->point_lists[i_a]] + d2vdxdy[this->point_lists[i_a]]) * this->_p.y[i_a];
        double b4 = 1.0/ Pars::Re * (dudy[this->point_lists[i_a]] + dvdx[this->point_lists[i_a]]);
        A3A4= A3A4+(b1+b2+b3+b4) * (this->_p.s[i_a]);
		//----------------------------------------------for Fy------------------------------------------------------//

        double a1y = 0.5 * (this->_p.v[i_b] * this->_p.v[i_b] - this->_p.u[i_b] * this->_p.u[i_b]) - this->_p.v[i_b] * (this->_p.gz[i_b] / std::pow(this->_p.s[i_b],2)) * this->_p.x[i_b];
		double a2y = this->_p.x[i_b] * (this->_p.u[i_b] - u_jmin1[i])/(2.0 * Pars::dt);
		double a3y = -1.0/Pars::Re * (2.0 * d2udx2[this->point_lists[i_b]] + d2udy2[this->point_lists[i_b]] + d2vdxdy[this->point_lists[i_b]]) * this->_p.x[i_b];
		double a4y = -2.0/Pars::Re* (dvdy[this->point_lists[i_b]]);
		A1A2y = A1A2y+(a1y+a2y+a3y+a4y) * (this->_p.s[i_b]);

        double b1y = -0.5 * (this->_p.v[i_a] * this->_p.v[i_a] - this->_p.u[i_a] * this->_p.u[i_a]) + this->_p.v[i_a] * (this->_p.gz[i_a] / std::pow(this->_p.s[i_a],2)) * this->_p.x[i_a];
		double b2y = -this->_p.x[i_a]*(this->_p.u[i_a] - u_jmax1[i])/(2.0 * Pars::dt);
		double b3y = 1.0/Pars::Re * (2.0 * d2udx2[this->point_lists[i_a]] + d2udy2[this->point_lists[i_a]] + d2vdxdy[this->point_lists[i_a]])*this->_p.x[i_a];
		double b4y = 2.0/Pars::Re * (dvdy[this->point_lists[i_a]]);
		A3A4y = A3A4y+(b1y+b2y+b3y+b4y) * (this->_p.s[i_a] );
    }*/

    for (int i = 0; i < this->a1a2.size(); i++)
    {
		//----------------------------------------------for Fx------------------------------------------------------//
        int i_b = this->a1a2[i]; //boundary bawah
        int i_a = this->a3a4[i]; //boundary atas
        
        double a1 = (this->_p.u[i_b] * this->_p.v[i_b]) + (this->_p.v[i_b] * (this->_p.gz[i_b] / std::pow(this->_p.s[i_b],2)) * this->_p.y[i_b]);
        double a2 = -(this->_p.y[i_b]) * (this->_p.u[i_b] - this->u_jmin1[i]) / (2 * Pars::dt);
        double a3 = 1.0 / Pars::RE * (2.0 * d2udx2[i_b] + d2udy2[i_b] + d2vdxdy[i_b]) * this->_p.y[i_b];
        
        double a4 = -1.0/ Pars::RE * (dudy[i_b] + dvdx[i_b]);
        
        A1A2 = A1A2+(a1+a2+a3+a4)*(this->_p.s[i_b]);
        
        double b1 = -this->_p.u[i_a] * this->_p.v[i_a] - this->_p.v[i_a] * (this->_p.gz[i_a] / std::pow(this->_p.s[i_a],2)) * this->_p.y[i_a];
        double b2 = this->_p.y[i_a] * ( this->_p.u[i_a] - u_jmax1[i]) / (2 * Pars::dt);
        double b3 = - 1.0 / Pars::RE * (2.0 * d2udx2[i_a] + d2udy2[i_a] + d2vdxdy[i_a]) * this->_p.y[i_a];
        double b4 = 1.0/ Pars::RE * (dudy[i_a] + dvdx[i_a]);
        A3A4= A3A4+(b1+b2+b3+b4) * (this->_p.s[i_a]);
		//----------------------------------------------for Fy------------------------------------------------------//
        
        double a1y = 0.5 * (this->_p.v[i_b] * this->_p.v[i_b] - this->_p.u[i_b] * this->_p.u[i_b]) - this->_p.v[i_b] * (this->_p.gz[i_b] / std::pow(this->_p.s[i_b],2)) * this->_p.x[i_b];
		double a2y = this->_p.x[i_b] * (this->_p.u[i_b] - u_jmin1[i])/(2.0 * Pars::dt);
		double a3y = -1.0/Pars::RE * (2.0 * d2udx2[i_b] + d2udy2[i_b] + d2vdxdy[i_b]) * this->_p.x[i_b];
		double a4y = -2.0/Pars::RE* (dvdy[i_b]);
		A1A2y = A1A2y+(a1y+a2y+a3y+a4y) * (this->_p.s[i_b]);
        
        double b1y = -0.5 * (this->_p.v[i_a] * this->_p.v[i_a] - this->_p.u[i_a] * this->_p.u[i_a]) + this->_p.v[i_a] * (this->_p.gz[i_a] / std::pow(this->_p.s[i_a],2)) * this->_p.x[i_a];
		double b2y = -this->_p.x[i_a]*(this->_p.u[i_a] - u_jmax1[i])/(2.0 * Pars::dt);
		double b3y = 1.0/Pars::RE * (2.0 * d2udx2[i_a] + d2udy2[i_a] + d2vdxdy[i_a]) *this->_p.x[i_a];
		double b4y = 2.0/Pars::RE * (dvdy[i_a]);
		A3A4y = A3A4y+(b1y+b2y+b3y+b4y) * (this->_p.s[i_a] );
        
    }


    
    //printf("%d\n", this->a1a2.size());
    //#pragma omp parallel for reduction (+:A2A3,A1A4,A2A3y,A1A4y)
    /*for (int i = 0; i < this->a2a3.size(); i++)
    {     
        int i_kn = this->a2a3[i]; //boundary kanan
        int i_kr = this->a4a1[i]; //boundary kiri
        
		//--------------------------------------------for Fx---------------------------------------------------------//
		double c1=0.5 * ( this->_p.v[i_kn] * this->_p.v[i_kn] - this->_p.u[i_kn] * this->_p.u[i_kn] );
		double c2=-this->_p.u[i_kn] * (this->_p.gz[i_kn] / std::pow(this->_p.s[i_kn],2))* this->_p.y[i_kn] - this->_p.y[i_kn] * (this->_p.v[i_kn] - this->v_imax1[i])/(2 * Pars::dt);
		double c3=1.0/Pars::Re * ( 2.0 * d2vdy2[this->point_lists[i_kn]] + d2vdx2[this->point_lists[i_kn]] + d2udxdy[this->point_lists[i_kn]] ) * this->_p.y[i_kn];
		double c4=2.0/Pars::Re * ( dudx[this->point_lists[i_kn]] );
		A2A3 = A2A3+ (c1+c2+c3+c4) * (this->_p.s[i_kn]);

		double d1=-0.5 * ( this->_p.v[i_kr] * this->_p.v[i_kr] - this->_p.u[i_kr] * this->_p.u[i_kr] );
		double d2=this->_p.u[i_kr] * (this->_p.gz[i_kr] / std::pow(this->_p.s[i_kr],2)) * this->_p.y[i_kr] + this->_p.y[i_kr] * (this->_p.v[i_kr]- this->v_imin1[i])/(2 * Pars::dt);
		double d3=-1.0/Pars::Re * ( 2.0 * d2vdy2[this->point_lists[i_kr]] + d2vdx2[this->point_lists[i_kr]] + d2udxdy[this->point_lists[i_kr]] ) * this->_p.y[i_kr];
		double d4=-2.0/Pars::Re * (dudx[this->point_lists[i_kr]] );
		A1A4 = A1A4+ (d1+d2+d3+d4)* (this->_p.s[i_kr]);

		//-------------------------------------------for Fy----------------------------------------------------------//

		double c1y=-this->_p.u[i_kn]*this->_p.v[i_kn] + this->_p.u[i_kn] * (this->_p.gz[i_kn] / std::pow(this->_p.s[i_kn],2)) * this->_p.x[i_kn];
		double c2y= this->_p.x[i_kn]*(this->_p.v[i_kn]-v_imax1[i])/(2.0*Pars::dt);
		double c3y=-1.0/Pars::Re*( 2 * d2vdy2[this->point_lists[i_kn]] + d2vdx2[this->point_lists[i_kn]] + d2udxdy[this->point_lists[i_kn]])*p.x[i_kn];
		double c4y= 1.0/Pars::Re*(dvdx[this->point_lists[i_kn]] + dudy[this->point_lists[i_kn]]);
		A2A3y = A2A3y+ (c1y+c2y+c3y+c4y)* (this->_p.s[i_kn]);

        double d1y=this->_p.u[i_kr]*this->_p.v[i_kr] - this->_p.u[i_kr] * (this->_p.gz[i_kr] / std::pow(this->_p.s[i_kr],2)) * this->_p.x[i_kr];
		double d2y= -this->_p.x[i_kr]*(this->_p.v[i_kr]-v_imin1[i])/(2.0*Pars::dt);
		double d3y= 1.0/Pars::Re*( 2 * d2vdy2[this->point_lists[i_kr]] + d2vdx2[this->point_lists[i_kr]] + d2udxdy[this->point_lists[i_kr]]) * p.x[i_kr];
		double d4y= -1.0/Pars::Re*(dvdx[this->point_lists[i_kr]] + dudy[this->point_lists[i_kr]]);
		A1A4y = A1A4y+ (d1y+d2y+d3y+d4y) * (this->_p.s[i_kr]);
	}*/
    
    for (size_t i = 0; i < this->a2a3.size(); i++)
    {     
        int i_kn = this->a2a3[i]; //boundary kanan
        int i_kr = this->a4a1[i]; //boundary kiri
        
		//--------------------------------------------for Fx---------------------------------------------------------//
		double c1=0.5 * ( this->_p.v[i_kn] * this->_p.v[i_kn] - this->_p.u[i_kn] * this->_p.u[i_kn] );
		double c2=-this->_p.u[i_kn] * (this->_p.gz[i_kn] / std::pow(this->_p.s[i_kn],2))* this->_p.y[i_kn] - this->_p.y[i_kn] * (this->_p.v[i_kn] - this->v_imax1[i])/(2 * Pars::dt);
		double c3=1.0/Pars::RE * ( 2.0 * d2vdy2[i_kn] + d2vdx2[i_kn] + d2udxdy[i_kn] ) * this->_p.y[i_kn];
		double c4=2.0/Pars::RE * ( dudx[i_kn] );
		A2A3 = A2A3+ (c1+c2+c3+c4) * (this->_p.s[i_kn]);

		double d1=-0.5 * ( this->_p.v[i_kr] * this->_p.v[i_kr] - this->_p.u[i_kr] * this->_p.u[i_kr] );
		double d2=this->_p.u[i_kr] * (this->_p.gz[i_kr] / std::pow(this->_p.s[i_kr],2)) * this->_p.y[i_kr] + this->_p.y[i_kr] * (this->_p.v[i_kr]- this->v_imin1[i])/(2 * Pars::dt);
		double d3=-1.0/Pars::RE * ( 2.0 * d2vdy2[i_kr] + d2vdx2[i_kr] + d2udxdy[i_kr] ) * this->_p.y[i_kr];
		double d4=-2.0/Pars::RE * (dudx[i_kr] );
		A1A4 = A1A4+ (d1+d2+d3+d4)* (this->_p.s[i_kr]);

		//-------------------------------------------for Fy----------------------------------------------------------//

		double c1y=-this->_p.u[i_kn]*this->_p.v[i_kn] + this->_p.u[i_kn] * (this->_p.gz[i_kn] / std::pow(this->_p.s[i_kn],2)) * this->_p.x[i_kn];
		double c2y= this->_p.x[i_kn]*(this->_p.v[i_kn]-v_imax1[i])/(2.0*Pars::dt);
		double c3y=-1.0/Pars::RE*( 2 * d2vdy2[i_kn] + d2vdx2[i_kn] + d2udxdy[i_kn])*p.x[i_kn];
		double c4y= 1.0/Pars::RE*(dvdx[i_kn] + dudy[i_kn]);
		A2A3y = A2A3y+ (c1y+c2y+c3y+c4y)* (this->_p.s[i_kn]);

        double d1y=this->_p.u[i_kr]*this->_p.v[i_kr] - this->_p.u[i_kr] * (this->_p.gz[i_kr] / std::pow(this->_p.s[i_kr],2)) * this->_p.x[i_kr];
		double d2y= -this->_p.x[i_kr]*(this->_p.v[i_kr]-v_imin1[i])/(2.0*Pars::dt);
		double d3y= 1.0/Pars::RE*( 2 * d2vdy2[i_kr] + d2vdx2[i_kr] + d2udxdy[i_kr]) * p.x[i_kr];
		double d4y= -1.0/Pars::RE*(dvdx[i_kr] + dudy[i_kr]);
		A1A4y = A1A4y+ (d1y+d2y+d3y+d4y) * (this->_p.s[i_kr]);
	}


    //printf("%f, %f, %f ,%f \n", A1A2, A2A3, A3A4, A1A4);
	
    for (size_t i = 0; i < this->a1a2.size(); i++)
    {
        int i_b = this->a1a2[i]; 
        int i_a = this->a3a4[i];
		this->u_jmin1[i]= this->u_jmin2[i];
		this->u_jmax1[i]= this->u_jmax2[i];
		this->u_jmin2[i]= this->_p.u[i_b]; 
		this->u_jmax2[i]= this->_p.u[i_a]; 
	}

	//for (int j=0;j<ny;j++)
	for (size_t i = 0; i < this->a2a3.size(); i++)
    {
        int i_kn = this->a2a3[i]; //boundary kanan
        int i_kr = this->a4a1[i]; //boundary kiri
		this->v_imin1[i] = this->v_imin2[i];
		this->v_imax1[i] = this->v_imax2[i];
		this->v_imin2[i] = this->_p.v[i_kr];
		this->v_imax2[i] = this->_p.v[i_kn];
	}

    double FX;
  	double FY;


    FX = A1A2+A2A3+A3A4+A1A4-(this->_omegasumx[iT+1]-this->_omegasumx[iT-1])/(2 * Pars::dt);
    FY = A1A2y+A2A3y+A3A4y+A1A4y-(this->_omegasumy[iT+1]-this->_omegasumy[iT-1])/(2 * Pars::dt);  

    double Cd2, Cl2;
    if (iT==0 && iT==1  )
	{
		Cd2 = 0.0;
		Cl2 = 0.0;
	}
	else  
	{
		Cd2 = 2.0*FX/((Pars::U_inf*Pars::U_inf)*Pars::Df);
		Cl2 = 2.0*FY/((Pars::U_inf*Pars::U_inf)*Pars::Df);
	}

    //for EOM
    Pars::gaya = FY;


    //saving
    
    if(iT == 0){
        printf("saving FORCES2....\n");
        std::string nama;
        nama.append("FORCES.csv");

        ofs.open("FORCES.csv");
        ofs <<"Iter" << "," << "Cd" << "," << "Cl" << "," << "A1A2" 
            << "," << "A2A3"<< "," << "A3A4"<< "," << "A4A1" << "," <<"omega1" 
            << "," << "omega2" << "," << "A1A2y" 
            << "," << "A2A3y"<< "," << "A3A4y"<< "," << "A4A1y" << "," <<"omega1y" 
            << "," << "omega2y\n";
        ofs << iT << "," << Cd2 << "," << Cl2 << "," << A1A2 
            << "," << A2A3 << "," << A3A4 << "," << A1A4 << "," << 0.0
            << "," << 0.0 << "," << A1A2y 
            << "," << A2A3y << "," << A3A4y << "," << A1A4y << "," << 0.0
            << "," << 0.0 <<"\n";
        ofs.close();
    }else{
        printf("saving FORCES2....\n");
        ofs.open("FORCES.csv", std::ofstream::out | std::ofstream::app);
	    ofs << iT * Pars::dt << "," << Cd2 << "," << Cl2 << "," << A1A2 
            << "," << A2A3 << "," << A3A4 << "," << A1A4 << "," << this->_omegasumx[iT+1]
            << "," << this->_omegasumx[iT-1] << ","<< A1A2y 
            << "," << A2A3y << "," << A3A4y << "," << A1A4y << "," << this->_omegasumy[iT+1]
            << "," << this->_omegasumy[iT-1]<< "\n";
	    ofs.close();
    }
}

void force_calculation::sep(Particle p, std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, std::vector<double> e, std::string s){
    std::ofstream ofs;
    std::string nama;
    nama.append("particles_");
    nama.append(s);
    nama.append(".csv");
    //ofs.open("Particles_calculation.csv");
    ofs.open(nama.c_str());
    //sep(this->_p, _d2udx2, _d2udy2, _d2udxdy, _dudy, _dudx);
    ofs <<"xp" << "," << "yp" << "," << "d2dx2" << "," << "d2dy2"<< "," << "d2dxdy"<< "," << "ddy"<< "," << "ddx\n";
    
    for (size_t i = 0; i<p.x.size();i++){
         ofs << p.x[i] << "," << p.y[i] << "," << a[i] << "," << b[i] << "," << c[i] << "," << d[i] << "," << e[i] <<"\n";
    }

    ofs.close();
}
