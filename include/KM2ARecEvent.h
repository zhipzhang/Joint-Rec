#ifndef __LHRECEVENT_HH__
#define __LHRECEVENT_HH__
#include "RtypesCore.h"
#include "TObject.h"
#include "TClonesArray.h"
#define _IS_MC 1

class KM2ARecEvent : public TObject {
 public:
  Int_t ev_n;
  Double_t mjd;
  Double_t dt;
#ifdef _IS_MC
  Int_t   id;
  Float_t E;
  Float_t theta;
  Float_t phi;
  Float_t corex;
  Float_t corey;
  //Float_t corez;
  Int_t pNpE;  //primary number of particles of ED
  Int_t pNpW;  //primary number of particles of WCDA
  Int_t pNuM;  //primary number of muons of MD
  Int_t pNuW;  //primary number of muons of WCDA
  Int_t pNuW2; //primary number of muons of WCDA above 1GeV
#endif
  Float_t rec_x;
  Float_t rec_y;
  //Float_t rec_z;
  Float_t rec_theta;  
  Float_t rec_phi;
  Float_t rec_a;
  Float_t rec_c0;
  Float_t rec_sigma; 

  Float_t rec_Etheta_p; //planar fit 
  Float_t rec_Ephi_p;
  Float_t rec_Ec0_p;
  Float_t rec_Esigma_p;
  Float_t rec_Wtheta_p;
  Float_t rec_Wphi_p;
  Float_t rec_Wc0_p;
  Float_t rec_Wsigma_p;
  Float_t rec_Etheta_c; //conical fit 
  Float_t rec_Ephi_c;
  Float_t rec_Ec0_c;
  Float_t rec_Esigma_c;
  Float_t rec_Ea;
  Float_t rec_Wtheta_c;
  Float_t rec_Wphi_c;
  Float_t rec_Wc0_c;
  Float_t rec_Wsigma_c;
  Float_t rec_Wa;

  Float_t rec_Esize;
  Float_t rec_Wsize;
  Float_t rec_Eage;
  Float_t rec_Wage;

  Float_t rec_Ex;
  Float_t rec_Ey;
  Float_t rec_Ez;
  Float_t rec_Wx;
  Float_t rec_Wy;
  Float_t rec_Wz;

  Int_t NhitE;  //Nhit of ED
  Int_t NhitM;  //Nhit of MD
  Int_t NhitW;  //Nhit of WCDA

  Int_t NtrigE;  //Ntrig of ED
  Int_t NtrigW;  //Ntrig of WCDA
  Int_t NfiltE;  //Nhit of ED after NoiseFilter Function with r : 0m ~ 200m, t: -50ns ~ 100ns
  Int_t NfiltE2; //Nhit After TWind and RWind
  Int_t NfiltM;  //Nhit of MD after filter out noise
  Int_t NfiltW;  //Nhit of WCDA after filter out noise

  Float_t NpE1;  //number of particles of ED  with r : 0m  ~ 100m, t : -30ns ~ 50ns
  Float_t NpE2;  //number of particles of ED  with r : 40m ~ 100m, t : -30ns ~ 50ns
  Float_t NpE3;  //number of particles of ED  with r : 40m ~ 200m, t : -30ns ~ 50ns
  Float_t NpE4;  //number of particles of ED  with r : 0m  ~ 400m, t : -30ns ~ 50ns
  Float_t NpE5;  //number of particles of ED  with r : 0m  ~ 600m, t : -30ns ~ 50ns
  Float_t NpE6;  //number of particles of ED  with r : 0m  ~ 800m, t : -30ns ~ 50ns
  Float_t DenE1;
  Float_t DenE2;
  Float_t DenE3;

  Float_t NpW;  //number of particles of WCDA
  Float_t NuM1;  //number of muons of MD with r : 15m ~ 200m, t : -30s ~ 50s
  Float_t NuM2;  //number of muons of MD with r : 40m ~ 200m, t : -30s ~ 50s
  Float_t NuM3;  //number of muons of MD with r : 15m ~ 400m, t : -30s ~ 50s 
  Float_t NuM4;  //number of muons of MD with r : 40m ~ 400m, t : -30s ~ 50s
  Float_t NuM5;  //number of muons of MD with r : 15m ~ 600m, t : -30s ~ 50s
  Float_t NuM6;  //number of muons of MD with r : 15m ~ 800m, t : -30s ~ 50s
  Float_t NuW1;  //number of muons of WCDA
  Float_t NuW2;  //number of muons of WCDA
  Float_t NuW3;  //number of muons of WCDA 

  Float_t Rec_Erho;  // Reconstruction Energy Use rho50
  Float_t Direction_Error;
  KM2ARecEvent();
  virtual ~KM2ARecEvent();
  ClassDef(KM2ARecEvent,1)
};

#endif
