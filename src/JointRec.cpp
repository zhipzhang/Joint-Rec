/*
        Main Program to reconstruction the LACT and KM2A
        First Version (2023.2.7): Both Rec Shower individually.
        Usage : 
*/


#include "G4KM2A_Geometry.h"
#include "G4KM2A_Reconstruction.h"
#include "LACTree.h"
#include "Detect_config.h"
#include "LACT_Reconstruction.h"
#include "LACTEvent.h"
#include "LACTRecEvent.h"
#include <cstring>
#include <string>
#include <vector>
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "straux.h"
#include "LACT_RunPara.h"
#include "KM2ARecEvent.h"
#include "LHEvent.h"
#include "TRandom2.h"

static double Tresolution = 0.2; // Time Resolution of KM2A 0.2ns
static   int    array_flag  = 6  ; // KM2A Full Array Flag

void SetKM2AStatus(LHEvent* km2aevent, G4KM2A_Geometry* geom);

int main(int argc, char** argv)
{
        LACT_RUNPARA* LACT_Runpara = new LACT_RUNPARA(); 

        LACTRecEvent* lactrec  = new LACTRecEvent();
        LACTEvent*   lactevent = new LACTEvent();
        KM2ARecEvent* km2arec  = new KM2ARecEvent();
        LHEvent*     km2aevent = new LHEvent();
        LACT_Reconstruction* Lact_Reconstruction = new LACT_Reconstruction();
        LACT_Runpara->ProcessCommandLine(argc, argv);
        TFile* out_root  = TFile::Open(LACT_Runpara->out_file.c_str(), "recreate");
        TTree* joint_rectree = new TTree("JointRectree", "LACT_KM2ARec");
        joint_rectree->Branch("LactRecEvent", &lactrec);
        joint_rectree->Branch("Km2aRecEvent", &km2arec);

        G4KM2A_Reconstruction* Km2a_Reconstruction = G4KM2A_Reconstruction::GetInstance(array_flag);
        G4KM2A_Geometry* geom = G4KM2A_Geometry::GetInstance(array_flag);

        Lact_Reconstruction->GetCommandConfig(LACT_Runpara);
        for( auto input : LACT_Runpara->input_file)
        {
                TFile* input_root = TFile::Open(input.c_str(),"read");
                if( input_root->IsZombie())
                {
                        std::cout << "File " << input << " Failed to Open !" << std::endl;
                        continue;
                }
                TTree* events = (TTree*) input_root->Get("match_events");
                TTree* config_tree = (TTree*) input_root->Get("config_tree");
                Lact_Reconstruction->GetTelConfig(config_tree);

                events->SetBranchAddress("km2a_events", &km2aevent);
                lactevent->Init(events);
                int n = events->GetEntries();

                if(Lact_Reconstruction->OnlyDraw())
                {
                        Lact_Reconstruction->Draw(events, lactevent);    
                        Km2a_Reconstruction->SetDrawEvents(Lact_Reconstruction->GetDrawEvents());
                        for (auto i : Km2a_Reconstruction->GetDrawEvents())
                        {
                                events->GetEntry(i);
                                SetKM2AStatus(km2aevent, geom);
                                Km2a_Reconstruction->Draw(km2aevent, i, km2arec);
                        }
                        input_root->Close();
                        break;

                }
                Lact_Reconstruction->SetEventPix(lactevent);
                for( int i = 0; i < events->GetEntries(); i++)
                {
                        events->GetEntry(i);
                        lactevent->GetData();
                        SetKM2AStatus(km2aevent, geom);
                        if( km2aevent->GetNtrigE() > 0 )
                        {
                                Km2a_Reconstruction->eventrecline(km2aevent, km2arec);
                        }

                        if(lactevent->IsTrigger())
                                Lact_Reconstruction->EventRec(lactevent, lactrec);
                        else
                        {
                                Lact_Reconstruction->Weight(lactevent);
                                lactrec->GetMCData(lactevent);
                                        
                        }
                        joint_rectree->Fill();
                        lactevent->Reset();
                        lactrec->Reset();
                }
                input_root->Close();


        }
        out_root->cd();
        joint_rectree->Write();
        out_root->Write();
        out_root->Close();


        
}


void SetKM2AStatus(LHEvent* km2aevent, G4KM2A_Geometry* geom)
{
        TRandom2* rd2 = new TRandom2();
        TClonesArray* HitsE, *HitsM;
        LHHit* Hit;
        HitsE = km2aevent->GetHitsE();
        HitsM = km2aevent->GetHitsM();
        int nhite = km2aevent->GetNhitE();
        int nhitm = km2aevent->GetNhitM();

        for( int i = 0; i < nhite; i++)
        {
                Hit = (LHHit*) (*HitsE)[i];
                int id = Hit->GetId();
                double x, y, z, dt;
                if( geom->GetEDxyz(id, x, y, z, 1) < 0)
                {
                        Hit->SetStatus(-1);
                }
                else 
                {
                        Hit->SetStatus(5);
                }
                dt  = Hit->GetTime() + rd2->Gaus(0, Tresolution);
                Hit->SetTime(dt);

        }
        for( int j = 0; j < nhitm ; j++)
        {
                Hit = (LHHit *)(*HitsM)[j];
                int id = Hit->GetId();
                double x,y,z, dt;
                if( geom->GetMDxyz(id, x, y, z, 1) < 0)
                {
                        Hit->SetStatus(-1); // -1 Mean Cannot Read Coordinates
                }
                else 
                {
                        Hit->SetStatus(5);
                }
                dt = Hit->GetTime() + rd2->Gaus(0, Tresolution);
                Hit->SetTime(dt);
        }

}