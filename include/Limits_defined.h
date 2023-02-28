#ifndef _Limits_
#define _Limits_


#ifndef LACT_MAXTEL
#define LACT_MAXTEL 39
#endif

#ifndef LACT_MAXPIXELS
#define LACT_MAXPIXELS 4000
#endif

#ifndef LACT_SUMWINDOW
#define LACT_SUMWINDOW 128
#endif

#ifndef LACT_MAX_TIMELEVELS
#define LACT_MAX_TIMELEVELS 10
#endif
#ifndef MAX_NEIGHBOR
#define MAX_NEIGHBOR 8
#endif
#define HI_GAIN 0
#define LOW_GAIN 1
/** The factor needed to transform from mean p.e. units to units of the single-p.e. peak:
    Depends on the collection efficiency, the asymmetry of the single p.e. amplitude
    distribution and the electronic noise added to the signals. */
#define CALIB_SCALE 0.92 









#endif
