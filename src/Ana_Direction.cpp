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

    TH1F* h1[10][20];
    TH1F* h2[10][20];
    TH1F* h3[20];
    TH1F* h4[20];
    for( int i = 0; i < 10; i++)
    {

        for( int j = 0; j < 20; j++)
        {
            h1[i][j] = new TH1F(Form("h%d", i * 20 + j), Form("h%d", i * 20 +j), 2000, 0, 2);
            h2[i][j] = new TH1F(Form("hcore%d", i * 20 + j), Form("h%d", i * 20 +j), 300, 0, 300);

        }
    }
    for( int i = 0; i < 20; i++)
    {
        
        h3[i] = new TH1F(Form("h300%d", i), Form("h3%d", i), 2000, 0, 2);
        h4[i] = new TH1F(Form("h400%d", i), Form("h4%d", i), 300, 0, 300);
    }
    TFile* out_root = TFile::Open(lactrunpara->out_file.c_str(), "recreate");
    out_root->cd();
    TH2D* angular = new TH2D("angular", "angular_resolution", 20, -1, 3, 10, 0, 800);
    TH1D* nevent  = new TH1D("event", "event", 20, -1, 3);
    TH1D* tmp_dist = new TH1D("tmp", "tmp_dist", 10, 0, 800);
    TH1D* tmp_size = new TH1D("tmpsize", "tmp_size", 20, -1, 3);
    std::vector<TH1D*> nev;
    std::vector<TH1D*> nev2;
    TH1D* nev3;
    TH1D* nev4;

    nev3 = new TH1D("nev3", "nev3", 20, -1, 3);
    nev4 = new TH1D("nev4", "nev4", 20, -1, 3);
    std::vector<double> x[11];
    std::vector<double> x2[11];
    std::vector<double> y[11];
    std::vector<double> energy_bin;
    std::vector<double> bili[10];
    
    for( int i = 1; i <= 10; i++)
    {
        nev.push_back(new TH1D(Form("nev %d", i), "nev", 20 , -1, 3));
        nev2.push_back(new TH1D(Form("nev2 %d", i), "nev", 20 , -1, 3));
    }
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
            std::cout << "energy is "<< lactrec->MCenergy << std::endl;
            double core_dist = sqrt(pow(lactrec->GetMCCoreX(), 2) + pow(lactrec->GetMCCoreY(), 2)); 
            double core_error = sqrt(pow(lactrec->GetMCCoreX() - lactrec->GetRecCoreX(), 2) + pow(lactrec->GetMCCoreY() - lactrec->GetMCCoreY(), 2));
            int idist = tmp_dist->FindFixBin(core_dist) - 1;
            int ienergy = tmp_size->FindFixBin(log10(lactrec->MCenergy)) - 1;
            if ( idist >= 10)
            {
                idist = 9;
            }
            nev[idist]->Fill(log10(lactrec->MCenergy));
            nev3->Fill(log10(lactrec->MCenergy));
            if( lactrec->direction_error > 0)
            {
                nev4->Fill(log10(lactrec->MCenergy));
                nev2[idist]->Fill(log10(lactrec->MCenergy));
                h1[idist][ienergy]->Fill(lactrec->direction_error);
                h2[idist][ienergy]->Fill(core_error);
                h3[ienergy]->Fill( lactrec->direction_error);
                h4[ienergy]->Fill(core_error);
            }


        } 
        input_root->Close();

    }

    out_root->cd();
   for( int i = 0; i < 10; i++)
   {
        for( int j = 0; j < 20; j++)
        {
            h1[i][j]->Write();
            h2[i][j]->Write();
        }
   }
   for( int j = 0; j  < 20; j++)
   {
        h3[j]->Write();
        h4[j]->Write();
   }
   for( int i = 0; i < 10; i++)
   {
        nev[i]->Write();
        nev2[i]->Write();
   }
   nev3->Write();
   nev4->Write();
   out_root->Write();
   return 0;
}