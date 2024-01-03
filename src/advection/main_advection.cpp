#include "advection.hpp"

#define ADVECTION_TYPE 1    // The type of advection calculation

// The basic advection
void advection::main_advection(Particle &p){
    // Advection prompt
    printf("\nCalculating advection ...\n");

    // Advection computational time manager
    double _time = omp_get_wtime();

    if (ADVECTION_TYPE == 1){
        printf("%s<+> Type 1: Euler type %s\n", FONT_CYAN, FONT_RESET);
        this->advection_euler(p);
    }
    else if (ADVECTION_TYPE == 2){
        printf("%s<+> Type 2: Runge Kutta 2nd order type %s\n", FONT_CYAN, FONT_RESET);
        this->advection_rk2(p);
    }
    
    // Display computational time
    _time = omp_get_wtime() - _time;
    printf("<-> Advection computational time:      [%f s]\n", _time);

    return;
}

// The euler advection
void advection::advection_euler(Particle &p)
{   
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        p.x[i] += p.u[i] * Pars::dt;
        p.y[i] += p.v[i] * Pars::dt;
    }}
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        p.x[i] += p.u[i] * Pars::dt;
        p.y[i] += p.v[i] * Pars::dt;
        p.z[i] += p.w[i] * Pars::dt;
    }}
    
    return;
}

// Advection runge kutta 2nd order method
void advection::advection_rk2(Particle &p)
{   
    // Advection computation 2D
    if (DIM == 2){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        double u1 = (p.u[i] / 2);
        double u2 = u1 / 2;
        double u3 = u2 ;
        double v1 = (p.v[i] / 2);
        double v2 = v1 / 2;
        double v3 = v2;

        p.x[i] += (p.u[i] + 2*u1 + 2*u2 + u3)*Pars::dt/6;
        p.y[i] += (p.v[i] + 2*v1 + 2*v2 + v3)*Pars::dt/6;
    }}
    
    // Advection computation 3D
    else if (DIM == 3){
    #pragma omp parallel for
    for (int i = 0; i < p.num; i++)
    {
        double u1 = (p.u[i] / 2);
        double u2 = u1 / 2;
        double u3 = u2 ;
        double v1 = (p.v[i] / 2);
        double v2 = v1 / 2;
        double v3 = v2;
        double w1 = (p.w[i] / 2);
        double w2 = w1 / 2;
        double w3 = w2;

        p.x[i] += (p.u[i] + 2*u1 + 2*u2 + u3)*Pars::dt/6;
        p.y[i] += (p.v[i] + 2*v1 + 2*v2 + v3)*Pars::dt/6;
        p.z[i] += (p.w[i] + 2*w1 + 2*w2 + w3)*Pars::dt/6;
    }}
    return;
}