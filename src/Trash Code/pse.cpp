#include "diffusion.hpp"

std::vector<double> diffusion::pse(int np, const std::vector<double> &xp, const std::vector<double> &yp,
                                   const std::vector<double> &sp, const std::vector<double> &gpz,
                                   const std::vector<int> &pair_i, const std::vector<int> &pair_j)
{
    // internal variables
    double *ar = new double[np]; // area of particles
    std::vector<double> dgdt_df_z(np);
    double dxij, dyij, rij, rij2, sij, sij2, sij3, spij2, tempzi, tempzj, xnij, gaussij;

    for (size_t i = 0; i < np; i++)
    {
        ar[i] = pow(sp[i], 2);
    }

    for (size_t k = 0; k < pair_i.size(); k++)
    {
        int i = pair_i[k];
        int j = pair_j[k];

        dxij = xp[i] - xp[j];
        dyij = yp[i] - yp[j];

        rij2 = pow(dxij, 2) + pow(dyij, 2);
        spij2 = (pow(sp[i], 2) + pow(sp[j], 2)) / 4.0e0;

        rij = std::sqrt(rij2);
        sij = std::sqrt(spij2);
        sij3 = sij * spij2;

        d_base_poisson.gaussian_function(rij, sij, gaussij);
        xnij = gaussij;

        tempzi = (2.0e0 * Pars::NU / sij3) * ((sp[i] + sp[j]) * 0.5) * (ar[i] * gpz[j] - ar[j] * gpz[i]) * xnij;
        dgdt_df_z[i] = dgdt_df_z[i] + tempzi;
        tempzj = (2.0e0 * Pars::NU / sij3) * ((sp[i] + sp[j]) * 0.5) * (ar[j] * gpz[i] - ar[i] * gpz[j]) * xnij;
        dgdt_df_z[j] = dgdt_df_z[j] + tempzj;
    }
    // -- delete internal variables
    delete[] ar;
    // -- return value
    return dgdt_df_z;
}
