#include "penalization.hpp"
#include "../geometry/geometry.hpp"

/**
 *  @brief Calculate the penalization mark (chi).
 *  
 *  @param	_particle  The particle data for chi to be calculated.
 *  @param  _bodyList  The body data list for calculation tools.
 */
void penalization::calculate_chi(Particle &_p, const std::vector<Body> &bList) const
{
    // Calculate the nearest particle distance to body surface
    
    // Internal variable and tool
    geometry geom_tool;
    std::vector<double> nearestDist(_p.num, 0.0e0);

    // Calculate the nearest distance of particle to the body
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++){
        // Create the position array of particle
        std::vector<double> _pos(DIM);
        _pos[0] = _p.x[i];
        _pos[1] = _p.y[i];
        if (DIM > 2)
        _pos[2] = _p.z[i];
        
        // Aliasing the body and calculate the minimum distance (relative to body surface normal)
        const Body &b = bList[_p.bodyPart[i]];
        nearestDist[i] = geom_tool.distance_calc(_pos, b, false);
    }
    
    // Adjust the mollification parameter
    printf("<+> Calculate chi parameter\n");
    _p.chi.resize(_p.num, 0.0e0);
    #pragma omp parallel for
    for (int i = 0; i < _p.num; i++){
        // Chi calculation and evaluation
        _p.chi[i] = this->get_chi(nearestDist[i]);
    }
    
    // // Special mask function for semi-implicit scheme. [Rasmussen 2011]
    // if (Pars::opt_kaipen == 1){
    //     // Recalculating chi
    //     #pragma omp parallel for
    //     for (int i = 0; i < _p.num; i++){
    //         double kai_tilde;
            
    //         // Recalculation kai for semi-implicit brinkmann penalization (Combine the implicit and explicit)
    //         if (_p.chi[i] < 1.0e0){
    //             kai_tilde = _p.chi[i] / (Pars::dt * lambda * (1.0e0 - _p.chi[i]));
    //         }
    //         else{
    //             kai_tilde = _p.chi[i];
    //         }

    //         // Minimum between 1.0 and kai_tilde
    //         _p.chi[i] = kai_tilde < 1.0 ? kai_tilde : 1.0;
    //     }
    // }
}

/**
 *  @brief  Calculate the penalization mark value.
 *         
 *  @param  _normalDist  Normal distance of particle toward near body surface.
 * 
 *  @return The value of penalization mark (χ).
*/
double penalization::get_chi(const double _normalDist) const
{   
    // Adjust the mollification parameter
	double chi = 0.0;
    
	if (_normalDist > Pars::hmollif){                 // Outside the body and mollification region
		chi = 0.0e0;
	}
	else if (std::abs(_normalDist) <= Pars::hmollif){ // Inside the mollification region (transition region)
		chi = -0.5e0 * (-1.0e0 + _normalDist / Pars::hmollif + 
              (1.0e0 / M_PI) * sin(M_PI * _normalDist / Pars::hmollif));
	}
	else if (_normalDist < -Pars::hmollif){           // Inside the body region only
		chi = 1.0e0;
	}


    // Special mask function for SEMI_IMPLICIT scheme. [Rasmussen 2011]
    if (Pars::opt_kaipen == 1){
        // Recalculation chi.
        double kai_tilde = 1.0;
                
        // Recalculation kai for SEMI_IMPLICIT brinkmann penalization (Combine the implicit and explicit)
        if (chi < 1.0e0){   // Preventing division by 0
            kai_tilde = chi / (Pars::dt * lambda * (1.0e0 - chi));
        }

        // Minimum between 1.0 and kai_tilde
        chi = kai_tilde < 1.0 ? kai_tilde : 1.0;
    }

	return chi;
}

/**
 *  @brief Define the value of penalization constant (lambda).
*/
void penalization::lambda_def(){
    // Update the value of lambda; References [Rasmussen 2011]
    
    // IMPLICIT Scheme
    if (Pars::opt_pen == 1){
        lambda = Pars::lambda;
    }
    // SEMI-IMPLICIT Scheme
    else if (Pars::opt_pen == 2){
        Pars::opt_kaipen = 1;
        lambda = Pars::lambda;
    }
    // EXPLICIT Scheme
    else if (Pars::opt_pen == 3){
        lambda = 1.0e0 / Pars::dt;
    }

    // Update the value of global lambda
    Pars::lambda = lambda;
}