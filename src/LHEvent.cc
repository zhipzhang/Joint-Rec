// -----------------------------------------------------------------------
// //                                KM2A
// //                            --- LHEvent ---
// //                       class implementation file
// // -----------------------------------------------------------------------
// //  Description:
// //
// //    A event class to store LHAASO event include ED, MD and WCDA
// //
// // =======================================================================
// // History:
// //
// //    Songzhan Chen  -- 2016.1 Created for fast simulation
// //    Songzhan Chen  -- 2016.11  adopt it for LHAASO-KM2A simulation
// //                   -- 2017.3  add LHHitColl for hit collection         
// //========================================================================
// //
#include "LHEvent.h"

ClassImp(LHEvent)
ClassImp(LHHit)
//ClassImp(LHHitColl)
//ClassImp(LHPhoton)
ClassImp(LHWave)

LHEvent::LHEvent()
{
  NhitE  =  0  ;
  NhitM  =  0  ;
  NhitW  =  0  ;
  NwaveE =  0  ;
  NwaveM =  0  ;
  HitsE  =  new TClonesArray("LHHit", 100) ;
  HitsM  =  new TClonesArray("LHHit", 100) ;
  HitsW  =  new TClonesArray("LHHit", 100) ;
  WaveE  =  new TClonesArray("LHWave", 100) ;
  WaveM  =  new TClonesArray("LHWave", 100) ;
}

LHEvent::~LHEvent(){
  HitsE->Clear() ; 
  delete  HitsE  ;
  HitsM->Clear() ;
  delete  HitsM  ;
  HitsW->Clear() ;
  delete  HitsW  ;
  WaveE->Clear() ;
  delete  WaveE  ;
  WaveM->Clear() ;
  delete  WaveM  ;
}

void LHEvent::AddHitE(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np){
  new((*HitsE)[NhitE++]) LHHit(p_id,p_time,p_pe,p_np) ;
}
void LHEvent::AddHitM(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np){
  new((*HitsM)[NhitM++]) LHHit(p_id,p_time,p_pe,p_np) ;
}
void LHEvent::AddHitW(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np){
  new((*HitsW)[NhitW++]) LHHit(p_id,p_time,p_pe,p_np) ;
}
void LHEvent::AddWaveE(Double_t dt, Double_t da){
  new((*WaveE)[NwaveE++]) LHWave(dt,da) ;
}
void LHEvent::AddWaveM(Double_t dt, Double_t da){
  new((*WaveM)[NwaveM++]) LHWave(dt,da) ;
}
void LHEvent::AddHitE(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e){
  new((*HitsE)[NhitE++]) LHHit(p_id,p_time,p_pe,p_np,p_e) ;
}
void LHEvent::AddHitM(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np,Int_t p_e){
  new((*HitsM)[NhitM++]) LHHit(p_id,p_time,p_pe,p_np,p_e) ;
}
void LHEvent::AddHitW(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np,Int_t p_e){
  new((*HitsW)[NhitW++]) LHHit(p_id,p_time,p_pe,p_np,p_e) ;
}
/*
LHHitColl::LHHitColl()
{
  Nhit  =  0  ;
  //Hits  =  new TClonesArray("LHHit", 100) ;
  Hits  =  new TClonesArray("LHPhoton", 100) ;
}

LHHitColl::~LHHitColl(){
  Hits->Clear() ;
  delete  Hits  ;
}
void LHHitColl::AddHit(Int_t p_id, Double_t p_time){
  new((*Hits)[Nhit++]) LHPhoton(p_id,p_time) ;
}
*/
/*
void LHHitColl::AddHit(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np, Int_t p_e){
  new((*Hits)[Nhit++]) LHHit(p_id,p_time,p_pe,p_np,p_e) ;
}
void LHHitColl::AddHit(Int_t p_id, Double_t p_time, Double_t p_pe,Int_t p_np){
  new((*Hits)[Nhit++]) LHHit(p_id,p_time,p_pe,p_np) ;
}
*/
