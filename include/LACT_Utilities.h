#ifndef _LACT_UTILITIES_H
#define _LACT_UTILITIES_H






#include <cmath>
#include "TCanvas.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TEllipse.h"
#include "TLine.h"

namespace Utilities 
{
        inline double line_point_distance (double xp1, double yp1, double zp1, 
                    double cx, double cy, double cz,
                double x, double y, double z)
    {
        double a, a1, a2, a3, b;
   
        a1 = (y-yp1)*cz - (z-zp1)*cy;
        a2 = (z-zp1)*cx - (x-xp1)*cz;
        a3 = (x-xp1)*cy - (y-yp1)*cx;
        a  = a1*a1 + a2*a2 + a3*a3;
        b = cx*cx + cy*cy + cz*cz;
        if ( a<0. || b<= 0. )
            return -1;
        return sqrt(a/b);
    }

        inline void angles_to_offset (double obj_azimuth, double obj_altitude,
   double azimuth, double altitude, double focal_length, 
   double *xoff, double *yoff)
    {
        double daz = obj_azimuth - azimuth;
        double coa = cos(obj_altitude);

        double xp0 = -cos(daz) * coa;
        double yp0 = sin(daz) * coa;
        double zp0 = sin(obj_altitude);

        double cx = sin(altitude);
        double sx = cos(altitude);

        double xp1 = cx*xp0 + sx*zp0;
        double yp1 = yp0;
        double zp1 = -sx*xp0 + cx*zp0;
        
        if ( xp1 == 0 && yp1 == 0 ) /* On-axis ? */
        {
            *xoff = *yoff = 0.;
            return;
        }
        #if 0 /* The more intuitive way but taking more CPU time */
        double q = acos(zp1); // off-axis angle
        double d = tan(q) * focal_length; // distance to camera center
        double alpha = atan2(yp1,xp1); // orientation

        *xoff = d * cos(alpha);
        *yoff = d * sin(alpha);
        #else /* Actually we don't need any further function calls */
        double s = focal_length / zp1;

        *xoff = s * xp1;
        *yoff = s * yp1;
        #endif
    } 

        inline void offset_to_angles (double xoff, double yoff, 
    double azimuth, double altitude, double focal_length, 
    double *obj_azimuth, double *obj_altitude)
    {
        if ( xoff == 0. && yoff == 0. ) /* Avoid division by zero */
        {
            *obj_azimuth = azimuth;
            *obj_altitude = altitude;
            return;
        }
        else
        {
            double d = sqrt(xoff*xoff+yoff*yoff);
            double q = atan(d/focal_length);

            double sq = sin(q);
            double xp1 = xoff * (sq/d);
            double yp1 = yoff * (sq/d);
            double zp1 = cos(q);

            double cx = sin(altitude);
            double sx = cos(altitude);

            double xp0 = cx*xp1 - sx*zp1;
            double yp0 = yp1;
            double zp0 = sx*xp1 + cx*zp1;

            *obj_altitude = asin(zp0);
            *obj_azimuth  = atan2(yp0,-xp0) + azimuth;
            if ( *obj_azimuth < 0. )
                *obj_azimuth += 2.*M_PI;
            else if ( *obj_azimuth >= (2.*M_PI ) )
                *obj_azimuth -= 2.*M_PI;
        }
    }

        inline void get_shower_trans_matrix (double azimuth, double altitude,
   double trans[][3])
    {
        double cos_z = sin(altitude);
        double sin_z = cos(altitude);
        double cos_az = cos(azimuth);
        double sin_az = sin(azimuth);
    
        trans[0][0] = cos_z*cos_az;
        trans[1][0] = sin_az;
        trans[2][0] = sin_z*cos_az;
        
        trans[0][1] = -cos_z*sin_az;
        trans[1][1] = cos_az;
        trans[2][1] = -sin_z*sin_az;
        
        trans[0][2] = -sin_z;
        trans[1][2] = 0.;
        trans[2][2] = cos_z;
    }

        inline void cam_to_ref (double ximg, double yimg, double phi,
   double ref_azimuth, double ref_altitude, double cam_rot,
   double azimuth, double altitude, double focal_length,
   double *axref, double *ayref, double *phiref)
    {
        double s, c;
        double ximg_rot, yimg_rot;
        double azm_img, alt_img, dphi1, dphi2;
        
        if ( fabs(cam_rot) > 1e-14 )
        {
            c = cos(cam_rot);
            s = sin(cam_rot);
            ximg_rot = ximg*c + yimg*s;
            yimg_rot = yimg*c - ximg*s;
        }
        else
        {
            ximg_rot = ximg;
            yimg_rot = yimg;
        }

        offset_to_angles(ximg_rot, yimg_rot, azimuth, altitude, focal_length,
            &azm_img, &alt_img);
        dphi1 = -atan(tan(azm_img-azimuth)*sin(alt_img));
        angles_to_offset(azm_img, alt_img, ref_azimuth, ref_altitude, 1.0,
            axref, ayref);
        dphi2 = -atan(tan(azm_img-ref_azimuth)*sin(alt_img));

        *phiref = phi + cam_rot + (dphi2-dphi1); 
    }

        inline int intersect_lines (double xp1, double yp1, double phi1,
   double xp2, double yp2, double phi2, double *xs, double *ys, double *sang)
    {
        double A1, B1, C1;
        double A2, B2, C2;
        double detAB, detBC, detCA;
        double s1, c1, s2, c2;
        
        /* Hesse normal form for line 1 */
        s1 = sin(phi1);
        c1 = cos(phi1);
        A1 = s1;
        B1 = -c1;
        C1 = yp1*c1 - xp1*s1;

        /* Hesse normal form for line 2 */
        s2 = sin(phi2);
        c2 = cos(phi2);
        A2 = s2;
        B2 = -c2;
        C2 = yp2*c2 - xp2*s2;

        detAB = (A1*B2-A2*B1);
        detBC = (B1*C2-B2*C1);
        detCA = (C1*A2-C2*A1);
        
        if ( fabs(detAB) < 1e-14 ) /* parallel */
        {
            if ( sang )
                *sang = 0.;
            if ( fabs(detBC) < 1e-14 && fabs(detCA) < 1e-14 ) /* same lines */
            {
                /* We could take any point on the line but use the middle here. */
                *xs = 0.5*(xp1+xp2);
                *ys = 0.5*(yp1+yp2);
                return 2;
            }
            *xs = *ys = 0.;
            return 0;
        }
        
        *xs = detBC / detAB;
        *ys = detCA / detAB;
        
        if ( sang != NULL )
        {
            double dx1 = (*xs-xp1);
            double dx2 = (*xs-xp2);
            double dy1 = (*ys-yp1);
            double dy2 = (*ys-yp2);
            double dr1 = sqrt(dx1*dx1+dy1*dy1);
            double dr2 = sqrt(dx2*dx2+dy2*dy2);
            double cos_ang;
            if ( dr1*dr2 == 0. )
                *sang = 0.;
            else
            {
                cos_ang = (dx1*dx2+dy1*dy2) / (dr1*dr2);
                if ( cos_ang >= 1. )
                    *sang = 0.;
                else if ( cos_ang <= -1. )
                    *sang = M_PI;
                else
                    *sang = acos(cos_ang);
            }
        }
        
        return 1;
    }
    inline double  angle_between(double azimuth1, double altitude1, double azimuth2, double altitude2)
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

}













#endif