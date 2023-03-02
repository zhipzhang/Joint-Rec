/*  main reconstruction code for LHAASOFM
 *  If you find any bug please send email to
 *  chensz@ihep.ac.cn
 *    2016-04-18 
 */
#ifndef __G4KM2A_RECONSTRUCTION_HH__
#define __G4KM2A_RECONSTRUCTION_HH__

#include "G4KM2A_Geometry.h"
#include "TClonesArray.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "math.h"
#include "LHEvent.h"
#include "KM2ARecEvent.h"

#define PI 3.14159265358979312
#include "TMinuit.h"
#include "TMath.h"

class G4KM2A_Reconstruction
{
  public:
     static G4KM2A_Reconstruction * GetInstance(int tFlag)
      {

       if (m_myself == 0) m_myself = new G4KM2A_Reconstruction(tFlag);
       return m_myself;
      }
     ~G4KM2A_Reconstruction();

  private:
     G4KM2A_Reconstruction(int tFlag);
     static G4KM2A_Reconstruction * m_myself;
     static TMinuit *minuit;                 //=new TMinuit(4) Mean 4 Parameters
     static double _np,_dirl,_dirm,_dirn;
     static char _style[20];
     static float nED[6000]; 
     static float nMD[2000];
     static float nWCDA[4000];
     std::vector<int> events;

     static G4KM2A_Geometry *geom;// = G4KM2A_Geometry::GetInstance(ARRAY); 
  public:
     int    Init(KM2ARecEvent *trec);
     void   setnparticle(TClonesArray &tHits,int np, double pe);
     int    trigger(TClonesArray &tHits,int np,int twind,const char *style,int th); 
     int    timefilter(TClonesArray &tHits,int np,int twind,const char *style, int flag);
     int    spacetimefilter(TClonesArray &tHits,int np,int twind,int rwind,const char *style, int flag);
     int    planarfit(TClonesArray &tHits,int np,KM2ARecEvent *trec,const char *style);
     int    conicalfit(TClonesArray &tHits,int np,KM2ARecEvent *trec,float alpha,const char *style);
     void   core_centre(TClonesArray &tHits,int np,KM2ARecEvent *trec,const char *style);
     int    core_centre2(TClonesArray &tHits,int np,int rp,KM2ARecEvent *trec,const char *style);
     double compact(TClonesArray &tHits,int np,int rp,KM2ARecEvent *trec,const char *style);
     int    noisefilterPlanar(TClonesArray &tHits,int np,int twind,int twinu,KM2ARecEvent *trec,const char *style);
     int    noisefilter(TClonesArray &tHits,int np,int twind,int twinu,int rwin,KM2ARecEvent *trec,const char *style);
     float  getnp(TClonesArray &tHits,int np,int rcut,KM2ARecEvent *trec,const char *style);
     float  getdensity(TClonesArray &tHits,int np,int rcut1,int rcut2,KM2ARecEvent *trec,const char *style);
     float  getdensityE(TClonesArray &tHits, int np, int rcut1, int rcut2, KM2ARecEvent *trec, const char *style);
     float  getdensityM(TClonesArray &tHits, int np, int rcut1, int rcut2, KM2ARecEvent *trec, const char *style);
     float  getmuM(TClonesArray &tHits,int np,int rcut,KM2ARecEvent *trec);
     float  getmuW(TClonesArray &tHits,int np,int cut,KM2ARecEvent *trec);
     void   Draw(LHEvent* tevent, int i, KM2ARecEvent* trec);

     //NKG function1 for likelihood1 taken into account the detectors with no hit
     static void NKGfunction1(int &npar,double *gin,double &f,double *par,int iflag);
     //NKG function2 for likelihood2 only taken into account the detectors with  hit 
     static void NKGfunction2(int &npar,double *gin,double &f,double *par,int iflag);
     //Using TMinuit do the Maximum Likelihood Methods
     int core_likelihood(TClonesArray &tHits,int np,KM2ARecEvent *trec,const char *style,int method);

     int eventrecline(LHEvent *tevent,KM2ARecEvent *trec);
     void SetDrawEvents(std::vector<int> e)
     {
         events = e;
     }
     std::vector<int> GetDrawEvents()
     {
         return events;
     }

};

#endif /* __G4KM2A_RECONSTRUCTION_HH__ */ 
