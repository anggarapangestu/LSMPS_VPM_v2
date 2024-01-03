#ifndef INCLUDED_FORTRAN_UTILS
#define INCLUDED_FORTRAN_UTILS
#include <iostream>
using namespace std;

namespace FortranUtils
{
extern "C"
{
    void biotsavart_fmm_2d_(int *n0, int *n1, int *npi, double *xi, double *yi, double *si, double *ui, double *vi,
                            int *n2, int *n3, int *npj, double *xj, double *yj, double *sj, double *gj,
                            int *icutoff, int *n_s, int *n_inter, int *ndp);
}

} // namespace FortranUtils

#endif