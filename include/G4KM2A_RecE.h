#ifndef _KM2A_RECE_H
#define _KM2A_RECE_H


#include "TMath.h"


double getNKGdensity(double age, double size, double r);
double getzenithcorrect(double e);
double recer50new3(double age, double size, double theta);
double angle_between(double azimuth1, double altitude1, double azimuth2, double altitude2);












#endif