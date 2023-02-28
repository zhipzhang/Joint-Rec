#ifndef _LACT_TELDATA_H
#define _LACT_TELDATA_H

#include <vector>
#include <map>


class LACT_TelData
{
    private:
        double tel_diretion[2];               //  [0] : Azimuth  [1] : Altitude
        std::vector<int > image_pixel;
        int itel;
        int npix;
        float* pe;
        float hottest;
        int overflow;
        bool  flag;                          // flag whether used for image

    public:
        LACT_TelData(int n);
        ~LACT_TelData();
        void SetTelid(int i)
        {
            itel = i;
        }
        int GetTelid()
        {
            return itel;
        }
        float GetHottest()
        {
            return hottest;
        }
        int GetOverFlow()
        {
            return overflow;
        }
        float* GetPe()
        {
            return pe;
        }
        float GetPe(int ipix)
        {
            return pe[ipix];
        }
        bool GetFlag()
        {
            return flag;
        }
        void SetPointDirection(double az, double al)
        {
            tel_diretion[0] = az;
            tel_diretion[1] = al;
        }
        double GetTelAzimuth()
        {
            return tel_diretion[0];
        }
        double GetTelAltitude()
        {
            return tel_diretion[1];
        }
        int GetImagePixnum()
        {
            return image_pixel.size();
            
        }
        int GetImagePixId(int i)
        {
            if( i < npix)
                return image_pixel[i];
            else
                return -1;
        }
        void AddOverFlow()
        {
            overflow++;
        }
        //Use Tail-Cuts to do the Image Cleaning
        bool ImageClean(double tail1, double tail2, std::map<int, std::vector< int> > &pixel_nerighbor);
};


























#endif