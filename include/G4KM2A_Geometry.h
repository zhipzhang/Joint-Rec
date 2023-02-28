// -----------------------------------------------------------------------
//                                G4KM2A
//                       --- G4KM2A_Geometry ---
//                           class header file
// -----------------------------------------------------------------------
//  Description:
//
//    Managers the LHAASO detector position
//
// -----------------------------------------------------------------------
// History:
//    Songzhan Chen    -- 2017.3 created
//           
//========================================================================
//
#ifndef __G4KM2AGEOMETRY_HH__
#define __G4KM2AGEOMETRY_HH__

//#include "globals.hh"

//#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#define CORE_X  0.0 //m    
#define CORE_Y  0.0 //m   
#define CORE_R  575.0 //m  
#define CORE_R2  637.0 //m 

class G4KM2A_Geometry
{
  public:
     static G4KM2A_Geometry * GetInstance(int tFlag)
      {
        
       if (m_myself == 0) m_myself = new G4KM2A_Geometry(tFlag);
       return m_myself;
      }
     ~G4KM2A_Geometry();

  private:
     G4KM2A_Geometry(int tFlag);
     static G4KM2A_Geometry * m_myself;
     static int Flag;
  public:
     void ReadFileED(const char filename[100]);
     void ReadFileMD(const char filename[100]);
     void ReadFileWCDA(const char filename[100]);
     void ReadFile();

     int GetNED(){return fNED;};
     int GetNMD(){return fNMD;};
     int GetNWCDA(){return fNWCDA;};
     //from  the id of array to the real id 
     int GetEDId(int n){return fEDId[n];};
     int GetMDId(int n){return fMDId[n];}; 
     int GetWCDAId(int n){return fWCDAId[n];}; 
     //from real id to the id of array
     int GetEDId2(int n);
     int GetMDId2(int n);
     int GetWCDAId2(int n);

     int GetEDxyz(int id, double &x,double &y, double &z, int flag);
     int GetMDxyz(int id, double &x,double &y, double &z, int flag);
     int GetWCDAxyz(int id, double &x,double &y, double &z, int flag);

     int Getxyz(int id, double &x,double &y, double &z, int flag,const char style[20]);

     double Get2EDr(int id1,int id2);  
     double Get2MDr(int id1,int id2);
     double Get2WCDAr(int id1,int id2);
     void GetRangeX(double &minx,double &maxx);
     void GetRangeY(double &miny,double &maxy);
     void GetRangeZ(double &minz,double &maxz);
     double GetRangeX();
     double GetRangeY();
     double GetRangeZ();
     double GetRotationED( ) { return fRotationED;};
     void SetArrayRange();
     int CheckCore(double x,double y,int nd);
     double GetCoreEdge(double dx,double dy);
   private:
     //string EDfilename,MDfilename;
     std::vector<int>    fEDId;
     std::vector<int>    fEDId2;
     std::vector<double> fEDx;
     std::vector<double> fEDy;
     std::vector<double> fEDz;
     std::vector<int>    fMDId;
     std::vector<int>    fMDId2;
     std::vector<double> fMDx;
     std::vector<double> fMDy;
     std::vector<double> fMDz;
     std::vector<int>    fWCDAId;
     std::vector<int>    fWCDAId2;
     std::vector<double> fWCDAx;
     std::vector<double> fWCDAy;
     std::vector<double> fWCDAz;


     int    fNED,fNMD,fNWCDA; 
     int    fEDIdmax,fMDIdmax,fWCDAIdmax; 
     double fMinx,fMaxx,fMiny,fMaxy,fMinz,fMaxz;  
     double fRotationED; //degree
     int fXX[140][2];
     int fYY[140][2];
};

#endif  /* __AG4GEOMETRYCONFIG_HH__ */
