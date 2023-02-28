#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include "LACTRecEvent.h"
#include "LACT_RunPara.h"
#include "TMath.h"



int main(int argc, char** argv)
{
    std::string out_file = "dst.root";
    std::string in_file = ""; 
    LACT_RUNPARA* lactrunpara = new LACT_RUNPARA();
    lactrunpara->ProcessCommandLine(argc, argv);

    TH1F* h1[50][50];
    TH1F* h2[50][50];
    TH1F* h3[50][50];
    for( int i = 0; i < 50; i++)
    {

        for( int j = 0; j < 50; j++)
        {
            h1[i][j] = new TH1F(Form("h%d", i * 100 + j), Form("h%d", i * 100 +j), 5000, 0, 0.05);
            h2[i][j] = new TH1F(Form("h2%d", i * 100 + j), Form("h%d", i * 100 +j), 5000, 0, 0.05);
            h3[i][j] = new TH1F(Form("h3%d", i * 100 + j), Form("h%d", i * 100 +j), 5000, 0, 0.05);

        }
    }
    TFile* out_root = TFile::Open(lactrunpara->out_file.c_str(), "recreate");
    out_root->cd();
    TH2D* mean_l = new TH2D("ml", "mean_l", 50, 2, 7, 50, 0, 750);
    TH2D* mean_AE = new TH2D("mae", "mean_ae", 50, 2, 7, 50, 0, 750);
    TH2D* sigma_l = new TH2D("sl", "sigma_l", 50, 2, 7, 50, 0, 750);
    TH2D* sigma_AE = new TH2D("sae", "sigma_ae", 50, 2, 7, 50, 0, 750);
    TH2D* mean_w  = new TH2D("mw", "mean_w",50, 2, 7, 50, 0, 750);
    TH2D* sigma_w  = new TH2D("sw", "sgima_w",50, 2, 7, 50, 0, 750);
    for( auto input :lactrunpara->input_file)
    {
        TFile* input_root = TFile::Open(input.c_str(), "read");
        LACTRecEvent* lactrec = new LACTRecEvent();
        if(input_root->IsZombie())
        {
            continue;
        }
        TTree* lactrectree = (TTree*) input_root->Get("LactRectree");
        lactrectree->SetBranchAddress("LactRecEvent", &lactrec);
        for( int i = 0; i < lactrectree->GetEntries(); i++)
        {
            lactrectree->GetEntry(i);
            if(lactrec->direction_error > 0)
            {
                for( int i = 0; i < lactrec->GetNtel(); i++)
                {
                    if( lactrec->good_image[i] >= 2)
                    {
                        int ir = mean_l->GetYaxis()->FindBin(lactrec->rec_rp[i]) - 1;
                        int is = mean_l->GetXaxis()->FindBin(log10(lactrec->size[i])) - 1;
                        if ( ir >= 50 || is >=50 || ir < 0 || is <0)
                        {
                            continue;
                        }
                        else 
                        {
                            h1[is][ir]->Fill(lactrec->length[i]);
                            h2[is][ir]->Fill(lactrec->width[i]);
                            h3[is][ir]->Fill(lactrec->size[i]/lactrec->MCenergy);
                        }
                    }
                }
            }
            std::cout << "End Entry "<< i <<std::endl;


        } 


        input_root->Close();

    }
    for( int i = 0 ; i < 50; i++)
    {
        for( int j = 0; j < 50; j++)
        {
            if(h1[i][j]->GetEntries() > 20)
            {
                double ia[] ={0,0,0};
                double ic[] = {0,0,0};
                double id[] = {0,0,0};
                double ib[] ={0.16, 0.5, 0.84};
                h1[i][j]->GetQuantiles(3, ia, ib);
                h2[i][j]->GetQuantiles(3, ic, ib);
                h3[i][j]->GetQuantiles(3, id, ib);
                mean_l->SetBinContent(i+1, j+1, ia[1]);
                mean_w->SetBinContent(i+1, j+1, ic[1]);
                mean_AE->SetBinContent(i+1, j+1, id[1]);
                sigma_l->SetBinContent(i+1, j+1, 0.5*(ia[2] - ia[0]));
                sigma_w->SetBinContent(i+1, j+1, 0.5*(ic[2] - ic[0]));
                sigma_AE->SetBinContent(i+1, j+1, 0.5 *(id[2] - id[0]));
            }
        }
    }
    out_root->cd();
   for( int i = 0; i < 50; i++)
   {
        for( int j = 0; j < 50; j++)
        {
            h1[i][j]->Write();
            h2[i][j]->Write();
            h3[i][j]->Write();
        }
   }
   mean_l->Write();
   sigma_l->Write();
   sigma_w->Write();
   mean_w->Write();
   sigma_AE->Write();
   mean_AE->Write();
   out_root->Write();
   return 0;
}