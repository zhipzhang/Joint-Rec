#ifndef _LACT_EVENT_H
#define _LACT_EVENT_H

#include "LACT_TelData.h"
#include <vector>
#include "Limits_defined.h"
#include "TTree.h"
#include "LACT_Telconfig.h"


class LACTEvent
{   
    private:
        std::vector<LACT_TelData*> tel_data;
        std::vector<int> tel_pass_cuts;

        // Mc Data
        unsigned int    Runnumber;
        unsigned int    Eventnumber;
        unsigned short int    particle_id;
        float  MCenergy;
        float  MCxcore;
        float  MCycore;
        float  MCaz;
        float  weight;
        float  MCze;
        float  MCal;
        float  xmax;
        float  hmax;
        float  Point_az;                          // Not Considered the Pointing Error !
        float  Point_ze;                          // Same as Above
        float  Point_al;

        float tel_az[LACT_MAXTEL];
        float tel_al[LACT_MAXTEL];

        int flag;
        unsigned int ntel;
        int npix[LACT_MAXTEL];

        unsigned int Trig_list[LACT_MAXTEL];
        unsigned int ntrig;
        std::map<int, int> map_telid_index;
        unsigned short int  pe_list[LACT_MAXTEL][LACT_MAXPIXELS];

        
    public:
        LACTEvent();
        ~LACTEvent();
        void Init(TTree *);
        bool GetData();
        void SetTelpix(int i, int n)
        {
            npix[i] = n;
        }
        void Reset();
        float GetAzimuth()
        {
            return MCaz;
        }
        float GetNTrig()
        {
            return ntrig;
        }
        float GetAltitude()
        {
            return MCal;
        }
        int GetTelnum()
        {
            return tel_data.size();
        }
        void ReWeight(float i)
        {
            weight = weight * i;
        }
        void SetWeight(float w)
        {
            weight = w;
        }
        std::vector<LACT_TelData*>& GetTelData( )
        {
            return tel_data;
            
        }
        LACT_TelData* GetTelData(int i)
        {
            return tel_data[i];
        }
        float GetMCenergy()
        {
            return MCenergy;
        }
        float GetMCxcore()
        {
            return MCxcore;
        }
        float GetMCycore()
        {
            return MCycore;
        }
        float GetMCal()
        {
            return MCal;
        }
        float GetMCaz()
        {
            return MCaz;
        }
        int GetParticleid()
        {
            return particle_id;
        }
        float GetWeight()
        {
            return weight;
        }
        float GetPointAz()
        {
            return Point_az;
        }
        float GetPointAl()
        {
            return Point_al;
        }
        int IsTrigger()
        {
            return flag;
        }
        int GetTelIndex(int tel_id)
        {
            return map_telid_index[tel_id];
        }
        float GetXmax()
        {
            return xmax;
        }
        float GetHmax()
        {
            return hmax;
        }




};


















#endif