/*
        Main Program to reconstruction the LACT and KM2A
        First Version (2023.2.7): Both Rec Shower individually.
        Usage : 
*/


#include "LACTree.h"
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
#include "TRandom2.h"
#include "TROOT.h"

static double Tresolution = 0.2; // Time Resolution of KM2A 0.2ns
static   int    array_flag  = 6  ; // KM2A Full Array Flag


int main(int argc, char** argv)
{
        gROOT->ProcessLine("#include <vector>");
        LACT_RUNPARA* LACT_Runpara = new LACT_RUNPARA(); 

        LACTRecEvent* lactrec  = new LACTRecEvent();
        LACTEvent*   lactevent = new LACTEvent();
        LACT_Reconstruction* Lact_Reconstruction = new LACT_Reconstruction();
        LACT_Runpara->ProcessCommandLine(argc, argv);
        TFile* out_root  = TFile::Open(LACT_Runpara->out_file.c_str(), "recreate");
        out_root->cd();
        TTree* lact_rectree = new TTree("LactRectree", "LACT_KM2ARec");
        lact_rectree->Branch("LactRecEvent", &lactrec) ;


        Lact_Reconstruction->GetCommandConfig(LACT_Runpara);
        for( auto input : LACT_Runpara->input_file)
        {
                TFile* input_root = TFile::Open(input.c_str(),"read");
                if( input_root->IsZombie())
                {
                        std::cout << "File " << input << " Failed to Open !" << std::endl;
                        continue;
                }
                TTree* events = (TTree*) input_root->Get("event");
                TTree* config_tree = (TTree*) input_root->Get("config_tree");
                Lact_Reconstruction->GetTelConfig(config_tree);

                lactevent->Init(events);
                int n = events->GetEntries();

                if(Lact_Reconstruction->OnlyDraw())
                {
                        Lact_Reconstruction->Draw(events, lactevent);    
                        out_root->Close();
                        return 0;

                }
                Lact_Reconstruction->SetEventPix(lactevent);
                for( int i = 0; i < events->GetEntries(); i++)
                {
                        events->GetEntry(i);
                        std::cout << "Processing Event" <<i <<std::endl;
                        lactevent->GetData();
                        if(lactevent->IsTrigger())
                        {
                                Lact_Reconstruction->EventRec(lactevent, lactrec);
                        }
                        else 
                        {
                                Lact_Reconstruction->Weight(lactevent);
                                lactrec->GetMCData(lactevent);
                        }
                        lact_rectree->Fill();
                        lactevent->Reset();
                        lactrec->Reset();
                }
                input_root->Close();


        }
        out_root->cd();
        lact_rectree->Write();
        out_root->Write();
        out_root->Close();
        
}

