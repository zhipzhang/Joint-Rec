//////////////////////////////////////////////////////////////////////////
// Event                                                                //
// Description of the event and hit parameters                        //
//////////////////////////////////////////////////////////////////////////
#ifndef __LHEVENT_HH__
#define __LHEVENT_HH__

#include "TObject.h"
#include "TClonesArray.h"
#define _IS_MC 1

class LHEvent : public TObject {

 private:
  Int_t ev_n; 
  Double_t mjd;
  Double_t dt;  //the difference of referee time for KM2A and WCDA (T_km2a-T_wcda)
  Int_t NhitE;  //Nhit of ED
  Int_t NhitM;  //Nhit of MD
  Int_t NhitW;  //Nhit of WCDA
  Int_t NwaveE;  //Nhit of EDWave
  Int_t NwaveM;  //Nhit of MDWave
  Int_t NtrigE;  //Ntrig of ED;
  Int_t NtrigW;  //Ntrig of WCDA;
  Int_t TrigEt;
#ifdef _IS_MC
  Double_t E;
  Int_t    Id;   //the primary type of particle
  Int_t    NpE;  //number of e/g/u in MD
  Int_t    NpW;  //number of e/g/u in WCDA
  Int_t    NuM;  //number of muon in MD
  Int_t    NuW;  //number of muon in WCDA
  Int_t    NuW2; //number of muon in WCDA E>1GeV
  Double_t Theta;
  Double_t Phi;
  Double_t Corex; 
  Double_t Corey;
#endif

  TClonesArray * HitsE;
  TClonesArray * HitsM;
  TClonesArray * HitsW;
  TClonesArray * WaveE;
  TClonesArray * WaveM;

 public:
  LHEvent();
  virtual ~LHEvent();

  // Set Methods
  void     SetEvN(Int_t s_ev_n)        { ev_n    = s_ev_n  ; } ; 
  void	   SetMjd(Double_t s_mjd)      { mjd     = s_mjd   ; } ;
  void     SetDt(Double_t s_dt)      { dt     = s_dt   ; } ;
  void	   SetNhitE(Int_t s_nhite)     { NhitE   = s_nhite ; } ;
  void     SetNtrigE(Int_t s_nhite)    { NtrigE  = s_nhite ; } ;
  void     SetTrigEt(Int_t s_nhite)    { TrigEt  = s_nhite ; } ;
  void     SetNhitM(Int_t s_nhitm)     { NhitM   = s_nhitm ; } ;
  void     SetNhitW(Int_t s_nhitw)     { NhitW   = s_nhitw ; } ;
  void     SetNtrigW(Int_t s_nhite)    { NtrigW  = s_nhite ; } ;
#ifdef _IS_MC
  void     SetE(Double_t s_E)         { E     = s_E     ; } ;
  void     SetId(Int_t s_id)          { Id    = s_id    ; } ;
  void     SetNpE(Int_t s_num)        { NpE   = s_num   ; } ;
  void     SetNpW(Int_t s_nuw)        { NpW   = s_nuw   ; } ;
  void     SetNuM(Int_t s_num)        { NuM   = s_num   ; } ;
  void     SetNuW(Int_t s_nuw)        { NuW   = s_nuw   ; } ;
  void     SetNuW2(Int_t s_nuw2)      { NuW2  = s_nuw2  ; } ;
  void     SetTheta(Double_t s_Theta) { Theta = s_Theta ; } ;
  void     SetPhi(Double_t s_Phi)     { Phi   = s_Phi   ; } ;
  void     SetCorex(Double_t s_Corex) { Corex = s_Corex ; } ;
  void     SetCorey(Double_t s_Corey) { Corey = s_Corey ; } ;
#endif 
 
  // Get Methods
  Int_t    GetEvN()     const { return ev_n  ; } ;
  Double_t GetMjd()     const { return mjd   ; } ; 
  Double_t GetDt()      const { return dt   ; } ;
  Int_t    GetNhitE()   const { return NhitE ; } ;
  Int_t    GetNhitM()   const { return NhitM ; } ;
  Int_t    GetNhitW()   const { return NhitW ; } ;
  Int_t    GetNwaveE()  const { return NwaveE ; } ;
  Int_t    GetNwaveM()  const { return NwaveM ; } ;
  Int_t    GetNtrigE()  const { return NtrigE ; } ;
  Int_t    GetNtrigW()  const { return NtrigW ; } ;
  Int_t    GetTrigEt()  const { return TrigEt ; } ;
#ifdef _IS_MC
  Double_t GetE()       const { return E     ; } ;
  Int_t    GetId()      const { return Id    ; } ;  
  Int_t    GetNpE()     const { return NpE   ; } ;
  Int_t    GetNpW()     const { return NpW   ; } ;
  Int_t    GetNuM()     const { return NuM   ; } ; 
  Int_t    GetNuW()     const { return NuW   ; } ;
  Int_t    GetNuW2()    const { return NuW2  ; } ; 
  Double_t GetTheta()   const { return Theta ; } ;
  Double_t GetPhi()     const { return Phi   ; } ;
  Double_t GetCorex()   const { return Corex ; } ;
  Double_t GetCorey()   const { return Corey ; } ;
#endif

  TClonesArray * GetHitsE() const  { return HitsE    ; }  ;
  TClonesArray * GetHitsM() const  { return HitsM    ; }  ;
  TClonesArray * GetHitsW() const  { return HitsW    ; }  ;
  TClonesArray * GetWaveE() const  { return WaveE    ; }  ;
  TClonesArray * GetWaveM() const  { return WaveM    ; }  ;

  // To add a hit ...
  void AddHitE(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np);
  void AddHitM(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np);
  void AddHitW(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np);
  void AddWaveE(Double_t dt, Double_t da);
  void AddWaveM(Double_t dt, Double_t da);

  void AddHitE(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e);
  void AddHitM(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e);
  void AddHitW(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e);

  void Initcsz()
  //void Clear()
  {
    HitsE->Clear() ;  NhitE = 0;
    HitsM->Clear() ;  NhitM = 0;
    HitsW->Clear() ;  NhitW = 0;
    WaveE->Clear() ;  NwaveE = 0;
    WaveM->Clear() ;  NwaveM = 0; 
  };

  ClassDef(LHEvent,1) 

 };
/*
class LHHitColl : public TObject {
     private:
       Int_t Nhit;
       TClonesArray * Hits;

     public:
       LHHitColl();
       virtual ~LHHitColl();

       Int_t          GetNhit() const  { return Nhit ; } ;
       TClonesArray * GetHits() const  { return Hits ; }  ;

     void AddHit(Int_t p_id, Double_t p_time);
     //void AddHit(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e);
     //void AddHit(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np);
     void Initcsz()
     {
         Hits->Clear() ;  Nhit = 0;
     }
     ClassDef(LHHitColl,1)
 };
*/
class LHHit : public TObject {
  
 private:
   Int_t    id; 
   Int_t    status; //5 is ok,2 is noise outof circle, 1 is noise after planar fit, 0 is noise hit after time filter, -1 is bad detector 
   Double_t time;
   Double_t pe; // number of pe or TDC count
   Int_t    np; //number of primary secondary particles
  
 public:
  
   // constructor and distructor
   LHHit() {};
  
   LHHit(Int_t p_id, Double_t p_time, Double_t p_pe, Int_t p_np){
       id = p_id; time = p_time; pe = p_pe;  np=p_np; status=5;
   };

   LHHit(Int_t p_id, Double_t p_time, Double_t p_pe, Int_t p_np, Int_t p_e){
       id = p_id; time = p_time; pe = p_pe;  np=p_np; status=p_e;
   };

   virtual ~LHHit() {} ;

   // Set Methods
   void     SetId(Int_t p_id)         { id     = p_id   ; } ;
   void     SetStatus(Int_t p_s)      { status = p_s    ; } ;
   void     SetTime(Double_t p_time)  { time   = p_time ; } ;
   void     SetPe(Double_t p_pe)      { pe     = p_pe   ; } ;
   void     SetNp(Int_t p_np)         { np     = p_np   ; } ;

   // Get Methods
   Int_t    GetId()   const  { return  id    ; }  ;
   Int_t    GetStatus()const { return  status; }  ;
   Double_t GetTime() const  { return  time  ; }  ;
   Double_t GetPe()   const  { return  pe    ; }  ;
   Int_t    GetNp()   const  { return  np    ; }  ;

   ClassDef(LHHit,1)
 
 };
/*
  //class to store single photon
class LHPhoton : public TObject {

 private:
   Int_t    id; // id
   Double_t time; //time  

 public:
   LHPhoton() {};

   LHPhoton(Int_t p_id, Double_t p_time){
       id = p_id; time = p_time; 
   };
   
   virtual ~LHPhoton() {} ;
   // Set Methods
   void     SetId(Int_t p_id)         { id     = p_id   ; } ;
   void     SetTime(Double_t p_time)  { time   = p_time ; } ;

   // Get Methods
   Int_t    GetId()   const  { return  id    ; }  ;
   Double_t GetTime() const  { return  time  ; }  ;
   ClassDef(LHPhoton,1)
};
*/
  //class to store ED and MD wave
class LHWave : public TObject {

 private:
   Double_t t; //time 
   Double_t a; //ampl 

 public:
   LHWave() {};

   LHWave(Double_t dt, Double_t da){
      t =dt; a = da;
   };

   virtual ~LHWave() {} ;
   // Set Methods
   void     SetT(Double_t dt)  { t   = dt   ; } ;
   void     SetA(Double_t da)  { a   = da ; } ;
   // Get Methods
   Int_t    GetT() const  { return  t    ; }  ;
   Double_t GetA() const  { return  a  ; }  ;
   ClassDef(LHWave,1)
};

#endif
