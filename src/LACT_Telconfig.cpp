#include "LACT_Telconfig.h"
#include "Limits_defined.h"
#include <algorithm>
#include <cmath>
void LACT_TelConfig::ComputeNeighbor()
{
    double x_i, y_i;
    double x_j, y_j;
    for(int i = 0; i < npix; i++)
    {
        x_i = x_pix[i];
        y_i = y_pix[i];
        for( int j = 0; j < npix; j++)
        {
            if( i != j)
            {
                x_j = x_pix[j];
                y_j = y_pix[j];
            }
            if( sqrt( pow(x_i - x_j, 2) + pow(y_i - y_j, 2)) < 1.7 * pixel_size)
            {
                if( pixel_neighbors[i].size() < MAX_NEIGHBOR)
                    pixel_neighbors[i].push_back(j);
                else
                    break;
            }
        }
        
        
    }
}

LACT_TelConfig::LACT_TelConfig()
{
    npix = 0;
    focal_length = effective_focal_length = 0;
    pixel_shape = -1;
    pixel_size = 0;

}