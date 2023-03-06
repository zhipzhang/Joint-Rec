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
    TH2D* h3 = new TH2D("h3", "Angular Resolution ( Core Dist < 250m)", 10, 1, 3, 1000, 0, 1);
    TH2D* h4 = new TH2D("h4", "Angular Resolution ( Core Dist < 400m)", 10, 1, 3, 1000, 0, 1);
    TH2D* h5 = new TH2D("h5", "Angular Resolution ( Core Dist < 550m)", 10, 1, 3, 1000, 0, 1);
    TH2D* h6 = new TH2D("h6", "Angular Resolution ( Core Dist > 550m)", 10, 1, 3, 1000, 0, 1);
    TH2D* h7 = new TH2D("h7", "Core Resolution ( Core Dist < 250m)", 10, 1, 3, 100, 0, 100);
    TH2D* h8 = new TH2D("h8", "Core Resolution ( Core Dist < 400m)", 10, 1, 3, 100, 0, 100);
    TH2D* h9 = new TH2D("h9", "Core Resolution ( Core Dist < 550m)", 10, 1, 3, 100, 0, 100);
    TH2D* h10 = new TH2D("h10", "Core Resolution ( Core Dist > 550m)", 10, 1, 3, 100, 0, 100);
    TH2D* h11 = new TH2D("h11", "Core_Pos not Pass the Condition 1 ", 320, -800, 800,  320, -800, 800);
    TH2D* h12 = new TH2D("h12", "Core_Dist Versus Energy", 20, 1, 3, 100, 0, 900);
    TH2D* h13 = new TH2D("h13", "pass age cut", 10, 1, 3, 1000, 0, 1);
    TH2D* h14 = new TH2D("h14", "pass age cut", 10, 1,3, 1000, 0, 1);
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
            double core_dist = sqrt(pow(km2arec->corex,2) +  pow(km2arec->corey,2));
            double core_err = sqrt(pow(km2arec->rec_x - km2arec->corex, 2) + pow( km2arec->rec_y - km2arec->corey, 2));
            if(core_dist < 250)
            {
                h3->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                h7->Fill(km2arec->Rec_Erho, core_err, flag * lactrec->weight);
            }
            else if ( core_dist < 400 )
            {
                h4->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                h8->Fill(km2arec->Rec_Erho, core_err, flag * lactrec->weight);
                if( km2arec->rec_Eage > 0.6 && km2arec->rec_Eage < 2.4 && (km2arec->NpE1- km2arec->NpE2 ) > km2arec->NpE2)
                {
                    h13->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                }
            }
            else if ( core_dist < 550)
            {
                h5->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                h9->Fill( km2arec->Rec_Erho, core_err, flag * lactrec->weight);
                if( km2arec->rec_Eage > 0.6 && km2arec->rec_Eage < 2.4 && (km2arec->NpE1- km2arec->NpE2 ) > km2arec->NpE2)
                {
                    h14->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                }
            }
            else 
            {
                h6->Fill(km2arec->Rec_Erho, km2arec->Direction_Error, flag * lactrec->weight);
                h10->Fill( km2arec->Rec_Erho, core_err, flag * lactrec->weight);
            }
            if( km2arec->rec_Eage > 0.6 && km2arec->rec_Eage < 2.4 && km2arec->Rec_Erho > 0 && (km2arec->NpE1-km2arec->NpE2) > km2arec->NpE2 && km2arec->NpE2 > 10)
            {
                h1->Fill(log10(km2arec->Rec_Erho), log10((km2arec->NuM3 + 0.0001) / km2arec->NpE4), lactrec->weight * flag);
                if( km2arec->Direction_Error > 0)
                {
                    h2->Fill(log10(km2arec->Rec_Erho), km2arec->Direction_Error, flag * lactrec->weight);
                }
            }
            if( km2arec->rec_Eage > 0.6 && km2arec->rec_Eage < 2.4 && (km2arec->NpE1 - km2arec->NpE2) < km2arec->NpE2  && km2arec->NpE1 > 10)
            {
                h11->Fill(km2arec->corex, km2arec->corey);
                h12->Fill(log10(lactrec->MCenergy), core_dist, flag * lactrec->weight);
            }
        }
        input_root->Close();
    }

    out_root->cd();
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
    h7->Write();
    h8->Write();
    h9->Write();
    h10->Write();
    h11->Write();
    h12->Write();
    h13->Write();
    h14->Write();
    out_root->Close();
}