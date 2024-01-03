#include "geometry.hpp"

// 
/* ===== Code Description ===== //
This code is performing "specified" body motion by updating each body particle
The motion follows the velocity from u_var method in this geomerty class
* The method u_var runs in the main code
*/
void geometry::moving_body(int it, Body &b)
{   
    // Internal variables 
    double _xcenter, _ycenter;
    double *_maxBody = new double[2];
    double *_minBody = new double[2];
    double new_distance;
    std::ofstream ofs;
    
    if (Pars::vib == 1)     // Vibration type 1
    {
        // Rigid body assumption
        if (Pars::m_star != 0)
        {
            // Chcek the current body center point
            _minBody[0] = *std::min_element(b.x.begin(), b.x.end());
            _minBody[1] = *std::min_element(b.y.begin(), b.y.end());
            _maxBody[0] = *std::max_element(b.x.begin(), b.x.end());
            _maxBody[1] = *std::max_element(b.y.begin(), b.y.end());
            _xcenter = (_maxBody[0] + _minBody[0]) / 2;
            _ycenter = (_maxBody[1] + _minBody[1]) / 2;

            // 4th order runge kutta method 
            double y, ydot;
            double dx1, dx2, dx3, dx4, dv1, dv2, dv3, dv4, a1, a2, a3, a4;

            y =_ycenter; // posisi body pada i saat ini; Kasus silinder = 0; kasus airfoil ?
            ydot = b.vT[1]; //kecepatan i ; Bisa begini karena asumsi rigid, 
                        //jadi tiap partikel body kecepatannya sama semua
                        //Bisa juga kecepatannya dibilang nol di awal. atau sama dengan kecepatan pada time step n                    
       
            if (it != 0 && it*Pars::dt >= 5.00){  //5detik karena katanya pas detik segitu udah gak ada instability
                                            // tes 15 detik
            // untuk k1
                //dx1 = y; 
                //dv1 = ydot;
                //a1 = (Pars::gaya - Pars::DamperConst * dv1 - Pars::SpringConst * dx1)/(Pars::mass - Pars::m_d);
            
                dx1 = Pars::dt * ydot; 
                dv1 = Pars::dt * ((Pars::gaya - Pars::DamperConst * ydot - Pars::SpringConst * y)/(Pars::mass - Pars::m_d));
            //untuk k2
                //dx2 = y + (dv1 * Pars::dt / 2.0);
                //dv2 = ydot + (a1 * Pars::dt / 2.0);
                //a2 = (Pars::gaya - Pars::DamperConst * dv2 - Pars::SpringConst * dx2) / (Pars::mass - Pars::m_d);
                dx2 = Pars::dt * (ydot+0.5*dv1); 
                dv2 = Pars::dt * ((Pars::gaya - Pars::DamperConst * (ydot+0.5*dv1) - Pars::SpringConst * (y+0.5*dx1))/(Pars::mass - Pars::m_d));
            //untuk k3
                //dx3 = y + (dv2 * Pars::dt / 2.0);
                //dv3 = ydot + (a2 * Pars::dt / 2.0);
                //a3 = (Pars::gaya - Pars::DamperConst * dv3 - Pars::SpringConst * dx3) / (Pars::mass - Pars::m_d);
                dx3 = Pars::dt * (ydot+0.5*dv2); 
                dv3 = Pars::dt * ((Pars::gaya - Pars::DamperConst * (ydot+0.5*dv2) - Pars::SpringConst * (y+0.5*dx2))/(Pars::mass - Pars::m_d));
            //untuk k4
                //dx4 = y + (dv3 * Pars::dt);
                //dv4 = ydot + (a3 * Pars::dt);
                //a4 = (Pars::gaya - Pars::DamperConst * dv4 - Pars::SpringConst * dx4) / (Pars::mass - Pars::m_d);
                dx4 = Pars::dt * (ydot+dv2); 
                dv4 = Pars::dt * ((Pars::gaya - Pars::DamperConst * (ydot+dv2) - Pars::SpringConst * (y+dv2))/(Pars::mass - Pars::m_d));
            //Calculating new value 
                //y = y + Pars::dt * ((dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4) / 6.0); //lokasi baru nih
                //new_distance = Pars::dt * ((dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4) / 6.0); 
                //ydot += Pars::dt * ((a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0); // kecepatan baru (?)
                y = y + ((dx1+2*dx2+2*dx3+dx4)/6);
                ydot = ydot + ((dv1+2*dv2+2*dv3+dv4)/6);
                new_distance = ((dx1+2*dx2+2*dx3+dx4)/6);
            }
    
            printf("Kecepatan saat ini: %f \n", ydot);

            printf("Body is moving ... \n");
            // -- Current position of body-points
            for (size_t i = 0; i < b.x.size(); i++)
            {
                b.x[i] += b.uT[it] * Pars::dt;
                //b.y[i] += b.vT[it] * Pars::dt;
                b.vT[i] = ydot;
                //b.y[i] += b.vT[i] * Pars::dt;
                b.y[i] += new_distance;
            }
    
            //Update center body
            _minBody[0] = *std::min_element(b.x.begin(), b.x.end());
            _minBody[1] = *std::min_element(b.y.begin(), b.y.end());
            _maxBody[0] = *std::max_element(b.x.begin(), b.x.end());
            _maxBody[1] = *std::max_element(b.y.begin(), b.y.end());
            _xcenter = (_maxBody[0] + _minBody[0]) / 2;
            _ycenter = (_maxBody[1] + _minBody[1]) / 2;
    
            printf("\nBody center is: [%f, %f]\n", _xcenter, _ycenter);
            // saving data
        
	        if (it==0){
            ofs.open("Vibration data.csv");
            ofs << "" << "it" << "," << "y_center" << "," << "y/d" << "\n";
            ofs << "" << it*Pars::dt << "," << y << "," << double(y)/Pars::Df << "\n";
            }else{
            double kk ;
            kk = it*Pars::dt;
        
            ofs.open ("Vibration data.csv", std::ios_base::app);    
            ofs << "" << kk << "," << y << "," << double(y)/Pars::Df << "\n";
            }
        }   
        else if (Pars::m_star == 0)
        {
            //Asumsi rigid.
            //cek center dari body pada saat ini. 
            _minBody[0] = *std::min_element(b.x.begin(), b.x.end());
            _minBody[1] = *std::min_element(b.y.begin(), b.y.end());
            _maxBody[0] = *std::max_element(b.x.begin(), b.x.end());
            _maxBody[1] = *std::max_element(b.y.begin(), b.y.end());
            _xcenter = (_maxBody[0] + _minBody[0]) / 2;
            _ycenter = (_maxBody[1] + _minBody[1]) / 2;

        //===============================================================================================================   
        //4th order runge kutta method 
            double y, ydot;
            double dx1, dx2, dx3, dx4, a1, a2, a3, a4;

            y =_ycenter; // posisi body pada i saat ini; Kasus silinder = 0; kasus airfoil ?
            ydot = b.vT[1]; //kecepatan i ; Bisa begini karena asumsi rigid, 
                        //jadi tiap partikel body kecepatannya sama semua
                        //Bisa juga kecepatannya dibilang nol di awal. atau sama dengan kecepatan pada time step n                    
       
            if (it != 0 && it*Pars::dt >= 5.00){  //5detik karena katanya pas detik segitu udah gak ada instability
                                            // tes 15 detik
                // untuk k1
                dx1 = Pars::dt * (Pars::gaya -Pars::SpringConst * y) / (Pars::DamperConst);
            
                //untuk k2
                dx2 = Pars::dt * (Pars::gaya -Pars::SpringConst * ( y + (dx1 / 2.0))) / (Pars::DamperConst);
                      
                //untuk k3
                dx3 = Pars::dt * (Pars::gaya -Pars::SpringConst * ( y + (dx2 / 2.0))) / (Pars::DamperConst);
                            
                //untuk k4
                dx4 = Pars::dt * (Pars::gaya -Pars::SpringConst * ( y + dx3))/ (Pars::DamperConst);
                
                //Calculating new value 
                new_distance = y + ((dx1 + 2.0 * dx2 + 2.0 * dx3 + dx4) / 6.0); //lokasi baru nih
                ydot  = (new_distance - y) / Pars::dt; // kecepatan baru (?)
            }
    
            printf("Kecepatan saat ini: %f \n", ydot);

            printf("Body is moving ... \n");
            // -- Current position of body-points
            for (size_t i = 0; i < b.x.size(); i++)
            {
                b.x[i] += b.uT[it] * Pars::dt; //nol terus
                //b.y[i] += b.vT[it] * Pars::dt;
                b.vT[i] = ydot;
                //b.y[i] += b.vT[i] * Pars::dt;
                b.y[i] += new_distance;
            }
    
            //Update center body
            _minBody[0] = *std::min_element(b.x.begin(), b.x.end());
            _minBody[1] = *std::min_element(b.y.begin(), b.y.end());
            _maxBody[0] = *std::max_element(b.x.begin(), b.x.end());
            _maxBody[1] = *std::max_element(b.y.begin(), b.y.end());
            _xcenter = (_maxBody[0] + _minBody[0]) / 2;
            _ycenter = (_maxBody[1] + _minBody[1]) / 2;
    
            printf("\nBody center is: [%f, %f]\n", _xcenter, _ycenter);
            // saving data
        
        
	        if (it==0){
            ofs.open("Vibration data.csv");
            ofs << "" << "it" << "," << "y_center" << "," << "y/d" << "\n";
            ofs << "" << it*Pars::dt << "," << y << "," << double(y)/Pars::Df << "\n";
            }else{
        
            double kk ;
            kk = it*Pars::dt;
        
            ofs.open ("Vibration data.csv", std::ios_base::app);    
            ofs << "" << kk << "," << y << "," << double(y)/Pars::Df << "\n";
            }  
        }
    }
    else if(Pars::vib == 2) // Vibration type 2
    {
        // Variabel untuk Runge-Kutta metode;
        double y, ydot;
        double dx1, dx2, dx3, dx4, dv1, dv2, dv3, dv4, a1, a2, a3, a4;

        y = Pars::tetha - Pars::tetha_nol;   //Angle saat ini (dalam deg)
        ydot = Pars::tetha_dot; //kecepatan angular saat ini; 
                        //jadi tiap partikel body kecepatannya sama semua
                        //Bisa juga kecepatannya dibilang nol di awal. atau sama dengan kecepatan pada time step n                    
       
        if (it != 0 && it*Pars::dt >= 5.00){  //5detik karena katanya pas detik segitu udah gak ada instability
                                            // tes 15 detik
            // untuk k1
            dx1 = y; 
            dv1 = ydot;
            a1 = (Pars::momen - Pars::DamperConst * dv1 - Pars::SpringConst * dx1)/(Pars::inertia);
            
            //untuk k2
            dx2 = y + (dv1 * Pars::dt / 2.0);
            dv2 = ydot + (a1 * Pars::dt / 2.0);
            a2 = (Pars::momen - Pars::DamperConst * dv2 - Pars::SpringConst * dx2) / (Pars::inertia);
            
            //untuk k3
            dx3 = y + (dv2 * Pars::dt / 2.0);
            dv3 = ydot + (a2 * Pars::dt / 2.0);
            a3 = (Pars::momen - Pars::DamperConst * dv3 - Pars::SpringConst * dx3) / (Pars::inertia);
            
            //untuk k4
            dx4 = y + (dv3 * Pars::dt);
            dv4 = ydot + (a3 * Pars::dt);
            a4 = (Pars::momen - Pars::DamperConst * dv4 - Pars::SpringConst * dx4) / (Pars::inertia);

            //Calculating new value 
            y = y + Pars::dt * ((dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4) / 6.0); //angle baru
            new_distance = Pars::dt * ((dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4) / 6.0); 
            ydot += Pars::dt * ((a1 + 2.0 * a2 + 2.0 * a3 + a4) / 6.0); // kecepatan angular baru (?)
            
            //changing new angle and angular velocity
            double change;
            change = (y + Pars::tetha_nol) -Pars::tetha; //deg
            Pars::tetha_dot = ydot; //deg/s
            Pars::tetha = y+ Pars::tetha_nol;
            printf("Body is rotating: %f \n", change);
        
        double Rotz [2][2];	

        Rotz[0][0] =  std::cos(change*M_PI/180.0e0);
	    Rotz[0][1] =  std::sin(change*M_PI/180.0e0);
	    Rotz[1][0] = -std::sin(change*M_PI/180.0e0);
	    Rotz[1][1] =  std::cos(change*M_PI/180.0e0);

       //printf("Kecepatan saat ini: %f \n", ydot);
        // -- Current position of body-points
        double xx,yy;
        for (size_t i = 0; i < b.x.size(); i++)
        {
            xx = Rotz[0][0]*b.x[i] + Rotz[0][1]*b.y[i];
		    yy = Rotz[1][0]*b.x[i] + Rotz[1][1]*b.y[i];
            b.x[i] = xx;
		    b.y[i] = yy;
        }
        }

        printf("\nBody center is: [%f, %f]\n", _xcenter, _ycenter);
        // saving data 
        
	    if (it==0){
        
            ofs.open("Vibration data.csv");
            ofs << "" << "it" << "," << "theta" << "," << "thetadot" << "\n";
            ofs << "" << it*Pars::dt << "," << Pars::tetha << "," << Pars::tetha_dot << "\n";
        
        }else{
    
            double kk ;
            kk = it*Pars::dt;
            ofs.open ("Vibration data.csv", std::ios_base::app);    
            ofs << "" << kk << "," << Pars::tetha << "," << Pars::tetha_dot << "\n";
        
        }
    }
    else                    // No vibration
    {
        // Rigid body no deformation
        printf("Body is moving ... \n");
        
        // Update the body-points position
        for (size_t i = 0; i < b.x.size(); i++)
        {
            b.x[i] += b.uT[it] * Pars::dt;
            b.y[i] += b.vT[it] * Pars::dt;
        }
    
        // Update the body center point (Not really neccesarry) [ND]
        _minBody[0] = *std::min_element(b.x.begin(), b.x.end());
        _minBody[1] = *std::min_element(b.y.begin(), b.y.end());
        _maxBody[0] = *std::max_element(b.x.begin(), b.x.end());
        _maxBody[1] = *std::max_element(b.y.begin(), b.y.end());
        _xcenter = (_maxBody[0] + _minBody[0]) / 2;
        _ycenter = (_maxBody[1] + _minBody[1]) / 2;
    
        printf("\nBody center is: [%f, %f]\n", _xcenter, _ycenter);
    }
	
    // Deallocating memory
    delete[] _maxBody;
    delete[] _minBody;

    // Write the final position of each body particle
    ofs.open("output/2Dbody.dat");
    for (int i = 0; i < Pars::n_a + 1; i++)
    {
        ofs << b.x[i] << " " << b.y[i] << "\n";
    }
    ofs.close();
}
