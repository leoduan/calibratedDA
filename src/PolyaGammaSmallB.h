// -*- mode: c++; -*-

////////////////////////////////////////////////////////////////////////////////
// Scott Linderman
//
// Polya-gamma sampling in the small shape parameter regime.
////////////////////////////////////////////////////////////////////////////////

#ifndef __POLYAGAMMASMALLB__
#define __POLYAGAMMASMALLB__

#include "RNG.h"
#include <cmath>

#define _MAXITER 1000

class PolyaGammaSmallB
{
    public:

    // Constructors.
    PolyaGammaSmallB();

    // Draw.
    double draw(double b, double z, RNG& r);

    private:

    double draw_invgauss_rej(double b, double z, RNG& r);
    double draw_invgamma_rej(double b, RNG& r);

    // Helper.
    inline double one_minus_psi(double x, double b);

};

#endif
