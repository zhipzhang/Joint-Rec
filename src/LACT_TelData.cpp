#include "LACT_TelData.h"
#include <algorithm>
#include <array>
#include <map>
#include <vector>
#include <queue>

LACT_TelData::LACT_TelData(int n)
{
    npix = n;
    pe = new float[npix];
    itel = -1;
    flag = false;
    overflow = 0;
    hottest = 0;
}

LACT_TelData::~LACT_TelData()
{
    delete []pe;
    image_pixel.clear();
}

bool LACT_TelData::ImageClean(double tail1, double tail2, std::map<int, std::vector<int> > &pixel_neighbor)
{
    std::queue<int> tmp;
    bool* res = new bool[npix];
    for( int i = 0; i < npix; i++)
    {
        res[i] = 0;
    }
    int max_pos = std::max_element(pe, pe + npix) - pe;
    if( pe[max_pos] < 10)
    {
        return false;
    }
    hottest = pe[max_pos];
    tmp.push(max_pos);
    while ( !tmp.empty() )
    {
        int tmp_p = tmp.front();
        res[tmp_p] = 1;
        image_pixel.push_back(tmp_p);
        tmp.pop();
        if( pe[tmp_p] > tail2)
        {
            for( auto k : pixel_neighbor[tmp_p])
            {
                if( res[k] )
                {
                    continue;
                }
                else 
                {
                    if(pe[k] >tail1)
                    {
                        res[k] = 1;
                        tmp.push(k);
                    }
                }
            }
        }
        else if (pe[tmp_p] > tail1) 
        {
            for( auto k : pixel_neighbor[tmp_p])
            {
                if( res[k])
                {
                    continue;
                }
                else 
                {
                    if( pe[k] > tail2)
                    {
                        res[k] = 1;
                        tmp.push(k);
                    }
                }
            }
        }
    }
    delete [] res;
    if( image_pixel.size() < 4)
    {
        flag = false;
    }
    else 
    {
        flag = true;    
    }
    return true;
}
