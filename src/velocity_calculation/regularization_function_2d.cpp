#include "velocity_biot_savart.hpp"


void velocity_biot_savart::gaussian_function(double rij, double sij, double &gaussij)
{
    // Gaussian for Diffustion form 2D
    // ------------ 2nd order -----------
    gaussij = std::exp(-rij * rij / (2.0e0 * pow(sij, 2))) / (pow(sij, 2) * (2.0e0 * M_PI));
}

void velocity_biot_savart::regul_func_2d(double rij, double sij, double &q)
{
    double qt, qte, rho, rho2, rij2, sij2, sij3;

    rij2 = pow(rij, 2);
    sij2 = pow(sij, 2);
    sij3 = std::sqrt(pow(sij2, 3));
    rho2 = rij2 / sij2;
    rho = std::sqrt(rho2);

    // icutoff =  0.singular ; 1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian
    if (Pars::icutoff == 0){
        q = 1.0e0 / (2.0e0 * M_PI);
    }
    else if (Pars::icutoff == 1){
        q = ((rho2 * (rho2 + 2.0e0)) / pow((rho2 + 1.0e0), 2)) / (2.0e0 * M_PI);
    }
    else if (Pars::icutoff == 2){
        qt = -rho2 / (2.0e0);
        q = (1.0e0 - std::exp(qt)) / (2.0e0 * M_PI);
    }
    else if (Pars::icutoff >= 3){
        qt = -rho2 / (2.0e0);
        qte = (1 - qt) * std::exp(qt);
        q = (1.0e0 - qte) / (2.0e0 * M_PI);
    }
}
