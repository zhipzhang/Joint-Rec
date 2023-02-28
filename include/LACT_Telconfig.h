#ifndef _LACT_TEL_CONFIG_H
#define _LACT_TEL_CONFIG_H

#include "Limits_defined.h"
#include <vector>
#include <map>
#include <iostream>
class LACT_TelConfig
{
    private:
        double tel_position[3];
        double focal_length;
        double effective_focal_length;

        int npix;                                   // Be Careful ! I assum all pixels are same
        int pixel_shape;                            // Pixel Shape
        double pixel_size;                          // Pixel Size (m)

        double x_pix[LACT_MAXPIXELS]{0};
        double y_pix[LACT_MAXPIXELS]{0};
        std::map<int, std::vector<int> > pixel_neighbors;
    public:
        LACT_TelConfig();
        ~LACT_TelConfig();
        int SetTelPos(int i, float* pos)
        {
            if( i != 3)
            {
                return -1;
            }
            for( int j = 0; j < i; j++)
            {
                tel_position[j] = pos[j];
            }
            return 1;
        } 
        int SetNpix(int n)
        {
            if( n > LACT_MAXPIXELS)
            {
                return -1;
            }
            npix = n;
            return 1;
        }
        void SetPix_Shape_Size(float size, int shape)
        {
            pixel_size = size / 1000;
            pixel_shape = shape;
        }
        int SetPixPos(float* x, float* y)
        {
            if( npix > 0 && npix < LACT_MAXPIXELS)
            {
                for( int i = 0; i < npix; i++)
                {
                    x_pix[i] = x[i] / 1000;
                    y_pix[i] = y[i] / 1000;
                }
            }
            else 
            {
                return -1;
            }
            return 1;
        }
        double GetFocalLength()
        {
            return focal_length;
        }
        void SetFocalLength(int f)
        {
            focal_length = f;
        }
        double GetPixSize()
        {
            return pixel_size;
        }
        void ComputeNeighbor();
        double GetPixX(int ipix)
        {
            return x_pix[ipix];
        }
        double GetPixY(int ipix)
        {
            return y_pix[ipix];
        }
        double GetTelpos(int i)
        {
            return tel_position[i];
        }
        int GetNpix()
        {
            return npix;
        }
        std::map<int, std::vector<int> >& GetNeighbor()
        {
            return pixel_neighbors;
        }
};




























#endif