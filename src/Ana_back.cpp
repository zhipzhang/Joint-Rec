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
            h2[i][j] = new TH1F(Form("hcore%d", i * 20 + j), Form("h%d", i * 20 +j), 2000, 0, 2);

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
    TGraph* g[11];
    TGraph* core[11];
    TGraph* rat[10];
    TGraph* effective_area;

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
            nev[idist]->Fill(lactrec->MCenergy);
            nev3->Fill(log10(lactrec->MCenergy));
            if( lactrec->direction_error > 0)
            {
                nev4->Fill(log10(lactrec->MCenergy));
                nev2[idist]->Fill(lactrec->MCenergy);
                h1[idist][ienergy]->Fill(lactrec->direction_error);
                h2[idist][ienergy]->Fill(core_error);
                h3[ienergy]->Fill( lactrec->direction_error);
                h4[ienergy]->Fill(core_error);
            }


        } 
        input_root->Close();

    }

    for( int i = 0 ; i < 10; i++)
    {
        for( int j = 0; j < 20; j++)
        {
            if( h1[i][j]->GetEntries() < 20 || h2[i][j]->GetEntries() < 20)
            {
                continue;
            }
            else 
            {
                double r68[] ={};
                double r682[] = {};
                double ia[]  = {0.68};
                h1[i][j]->GetQuantiles(1, r68, ia);
                h2[i][j]->GetQuantiles(1, r682, ia);
                x[i].push_back(*r68);
                x2[i].push_back(*r682);
                y[i].push_back(tmp_size->GetBinCenter(j + 1));

            }
        }
        g[i] = new TGraph(x[i].size(), &y[i][0], &x[i][0]);
        core[i] = new TGraph(x2[i].size(), &y[i][0], &x2[i][0]);
        core[i]->SetName(Form("core%d", i));
        g[i]->SetName(Form("g%d", i));
        g[i]->SetTitle(Form("%fm - %fm; Energy log10(E/TeV); Angular Resolution", tmp_dist->GetXaxis()->GetBinUpEdge(i), tmp_dist->GetXaxis()->GetBinUpEdge(i + 1)));
        core[i]->SetTitle(Form("%fm - %fm; Energy log10(E/TeV); Core Resolution", tmp_dist->GetXaxis()->GetBinUpEdge(i), tmp_dist->GetXaxis()->GetBinUpEdge(i + 1)));
    }


    for( int i = 0; i < 10; i++)
    {
        for( int k = 0; k < 20; k++)
        {
            double ratio = nev[i]->GetBinContent(k + 1)/ nev2[i]->GetBinContent(k + 1);
            bili[i].push_back(ratio);
        }
        
    }
    for( int i = 0; i < 10; i++)
    {
        rat[i] = new TGraph(bili[i].size(), &energy_bin[0], &bili[i][0]);
        rat[i]->SetName(Form("rat%d", i));
        rat[i]->SetTitle(Form("%fm - %fm; Energy log10(E/TeV); Ratio", tmp_dist->GetXaxis()->GetBinUpEdge(i), tmp_dist->GetXaxis()->GetBinUpEdge(i + 1)));

    }

    
    for( int i = 0; i < 20; i++)
    {
        double r68[]={};
        double ia[] = {0.68};
        double r683[] = {};
         if( h3[i]->GetEntries() < 20 || h4[i]->GetEntries() < 20)
          {
             continue;
         }
         else{
        h3[i]->GetQuantiles(1, r68, ia);
        double ia2[] = {0.68};
        h4[i]->GetQuantiles(1, r683, ia2);
        x[10].push_back(*r68);
        x2[10].push_back(*r683);
        energy_bin.push_back(tmp_size->GetXaxis()->GetBinCenter(i + 1));
         }
    }

   g[10] = new TGraph(x[10].size(),  &energy_bin[0], &x[10][0]);
   g[10]->SetName("g10");
   g[10]->SetTitle("Angular Resolution");
   core[10] = new TGraph(x2[10].size(), &energy_bin[0], &x2[10][0]);
   core[10]->SetName("core10");
   core[10]->SetTitle("Core Resolution");

    std::vector<double> ratio;
    std::vector<double> energy_bin2;
   for( int i = 0; i < 20; i++)
   {
        if(nev3->GetBinContent(i +1 ) > 0)
        {
            ratio.push_back(nev4->GetBinContent(i+1)/ nev3->GetBinContent(i+1) * TMath::Pi() * 1200 * 1200);
        
            energy_bin2.push_back(nev3->GetBinCenter(i+1 ));
        }
   }

   effective_area = new TGraph(energy_bin.size(), &energy_bin2[0], &ratio[0]);
   effective_area->SetName("effectivearea");
   effective_area->SetTitle("EffectiveArea");

    out_root->cd();
   for( int i = 0; i < 10; i++)
   {
        g[i]->Write();
        core[i]->Write();
        rat[i]->Write();
   }
   effective_area->Write();
   g[10]->Write();
   core[10]->Write();
   out_root->Write();
   return 0;
}