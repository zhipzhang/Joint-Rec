/*
    Because of the Large Data, Analysis Must be Done by Program not in interactive ROOR
*/
#include "LACT_RunPara.h"
#include "KM2ARecEvent.h"
#include "LACTRecEvent.h"
#include "TFile.h"
#include "TTree.h"



int main(int argc, char** argv)
{
    LACT_RUNPARA* lact_run = new LACT_RUNPARA();
    lact_run->ProcessCommandLine(argc, argv);
    TFile* out_root = TFile::Open(lact_run->out_file.c_str(), "recreate");
    TH2D* h1 = new TH2D("h1", "distribution of R versus Energy", 20, 1, 3, 100, -8, 2);
    TH2D* h2 = new TH2D("h2", "Angular Resolution", 10, 1, 3, 1000, 0, 1);

    KM2ARecEvent* km2arec = new KM2ARecEvent();
    LACTRecEvent* lactrec = new LACTRecEvent();

    for( auto input : lact_run->input_file)
    {
        TFile* input_root = TFile::Open(input.c_str(), "read");
        if( input_root->IsZombie() )
        {
            std::cout << "File " << input << "Failed to Open !" << std::endl;
            continue;
        }
        TTree* jointrec = (TTree*) input_root->Get("JointRectree");
        jointrec->SetBranchAddress("Km2aRecEvent", &km2arec);
        jointrec->SetBranchAddress("LactRecEvent", &lactrec);
        for( int i = 0 ; i < jointrec->GetEntries(); i++)
        {
            double flag = 1.0;
            jointrec->GetEntry(i);
            if(lactrec->MCenergy < 100)
            {
                flag = flag * lact_run->num_weight[2];
            }
            else 
            {
                flag = flag * lact_run->num_weight[3];
            }
            if( km2arec->rec_Eage > 0.6 && km2arec->rec_Eage < 2.4 && km2arec->Rec_Erho > 0 && (km2arec->NpE1-km2arec->NpE2) > km2arec->NpE2 && km2arec->NpE2 > 10)
            {
                h1->Fill(log10(km2arec->Rec_Erho), log10((km2arec->NuM3 + 0.0001) / km2arec->NpE4), lactrec->weight * flag);
                if( km2arec->Direction_Error > 0)
                {
                    h2->Fill(log10(km2arec->Rec_Erho), km2arec->Direction_Error, flag * lactrec->weight);
                }
            }
        }
        input_root->Close();
    }

    out_root->cd();
    h1->Write();
    h2->Write();
    out_root->Close();
}