#include "G4KM2A_RecE.h"
#include "TMath.h"
#include <cmath>


double getNKGdensity(double age, double size, double r)
{
    double cs, r0;
    double rm  = 130;
    cs = TMath::Gamma(4.5 -age)/( TMath::Gamma(age - 0.5) * TMath::Gamma(5 - 2 * age));
    cs = cs / (2 * TMath::Pi() * rm * rm);
    r0 = cs * size * pow(r/rm, age - 2.5) * pow(1 + r/rm, age -4.5);
    return r0;
}



/*
    p0, p1, p2 : Parameter Depends on Zenith Angle
    Age Must in 0.6 ~ 2.4
    Zenith Angle (deg) : 0 ~ 50

*/
double recer50new3(double age, double size, double theta)
{
    double p2[12]={0.02503, 0.02361, 0.02489, 0.02357, 0.02149, 0.02303, 0.01852, 0.01409, 0.01289, 0.02556, 0.03406, 0.03635,};
    double p1[12]={0.96475, 0.94221, 0.91961, 0.90142, 0.88660, 0.86839, 0.86159, 0.85272, 0.84218, 0.82113, 0.79659, 0.77589};
    double p0[12]={1.74558, 1.75789, 1.77675, 1.80182, 1.83308, 1.86808, 1.91001, 1.95468, 1.99553, 2.04191, 2.09133, 2.14474};
    double r50, re;
    int k ;
    if( age < 0.6 || age > 2.4)
    {
        return -1;
    }
    if( theta < 0 || (theta * 57.3 > 51))
    {
        return -1;
    }

    k = int ((1./(cos(theta) - 1)) / 0.05);
    if( k > 11)
    {
        k = 11;
    }
    if( k < 0)
    {
        k = 0;
    }
    r50 = log10(getNKGdensity(age, size, 50.));
    return pow(r50, 2) * p2[k] + r50 * p1[k] + p0[k];

}


double angle_between (double azimuth1, double altitude1, double azimuth2, double altitude2)
{
   double ax1 = cos(azimuth1)*cos(altitude1);
   double ay1 = sin(-azimuth1)*cos(altitude1);
   double az1 = sin(altitude1);
   double ax2 = cos(azimuth2)*cos(altitude2);
   double ay2 = sin(-azimuth2)*cos(altitude2);
   double az2 = sin(altitude2);
   double cos_ang = ax1*ax2 + ay1*ay2 + az1*az2;
   /* Check for rounding errors pushing us outside the valid range. */
   if ( cos_ang <= -1. )
      return M_PI;
   else if ( cos_ang >= 1. )
      return 0.;
   else
      return acos(cos_ang);
}