#include "remeshing.hpp"

void remeshing::create_extra_domain(Particle &p, double x1, double y1, double x2, double y2, double x_1, double y_1, double x_2, double y_2)
{
    double lenx = abs(x1 - x2);
    double leny = abs(y1 - y2);
    double s = Pars::D_0 * 2; //Changeable uhuy
    int nx = lenx / s;
    int ny = leny / s;
    
    for (int i = 0; i <= nx; i++){
        for (int j = 0; j <= ny; j++){
            double xx = x1 + i*s;
            double yy = y1 + j*s;
            if ((xx < x_1 || xx > x_2) || 
                ( yy < y_1 || yy > y_2)){
                p.x.push_back(xx);
                p.y.push_back(yy);
                p.s.push_back(s);
            }
        }
    }

    p.num = p.x.size();
    p.u.resize(p.num, 0.0e0);
    p.v.resize(p.num, 0.0e0);
    p.gz.resize(p.num, 0.0e0);
    p.neighbor.resize(p.num, std::vector<int>());
    p.isActive.resize(p.num, false);
    p.chi.resize(p.num, 0.0);
}