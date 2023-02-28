#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"



int main(int argc, char** argv)
{
    TFile* outroot = TFile::Open("dst.root", "recreate");
    TFile* inputroot  = TFile::Open("ball.root", "read");
    TH1F* h1[10][20];
    TH1F* h2[10][20];
    TH1F* h3[20];
    TH1F* h4[20];

    TH1D* tmp_dist = new TH1D("tmp", "tmp_dist", 10, 0, 800);
    TH1D* tmp_size = new TH1D("tmpsize", "tmp_size", 20, -1, 3);
    for( int i = 0; i < 10; i++)
    {
        for( int j = 0; j < 20 ; j++)
        {
            h1[i][j] = (TH1F*)inputroot->Get(Form("h%d", i * 20 +j));
            h2[i][j] = (TH1F*)inputroot->Get(Form("hcore%d", i * 20 +j))    ;
        }
    }
    for( int i = 0; i < 20 ; i++)
    {
        h3[i] = (TH1F*)inputroot->Get(Form("h300%d", i));
        h4[i] = (TH1F*)inputroot->Get(Form("h400%d", i));
    }
    TH1D* nev[10];
    TH1D* nev2[10];
    TH1D* nev3;
    TH1D* nev4;

    nev3 = (TH1D*)inputroot->Get("nev3");
    nev4 = (TH1D*)inputroot->Get("nev4");
    for( int i = 0; i < 10; i ++)
    {
        nev[i] = (TH1D*)inputroot->Get(Form("nev %d", i+1));
        nev2[i] = (TH1D*)inputroot->Get(Form("nev2 %d", i+1));
    }
    std::vector<double> x[11];
    std::vector<double> x2[11];
    std::vector<double> y[11];
    std::vector<double> energy_bin;
    std::vector<double> bili[10];
     TGraph* g[11];
    TGraph* core[11];
    TGraph* rat[10];
    TGraph* effective_area;
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
    std::vector<double>e_bins[10];
     for( int i = 0; i < 10; i++)
    {
        for( int k = 0; k < 20; k++)
        {
                if(nev[i]->GetBinContent(k + 1) >0)
            {
                double ratio = nev2[i]->GetBinContent(k + 1)/ nev[i]->GetBinContent(k + 1);
                bili[i].push_back(ratio);
                e_bins[i].push_back(tmp_size->GetBinCenter(k +1));
            }
        }
        
    }
    for( int i = 0; i < 10; i++)
    {
        rat[i] = new TGraph(bili[i].size(), &e_bins[i][0], &bili[i][0]);
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

    outroot->cd();
   for( int i = 0; i < 10; i++)
   {
        g[i]->Write();
        core[i]->Write();
        rat[i]->Write();
   }
   effective_area->Write();
   g[10]->Write();
   core[10]->Write();
   outroot->Write();
   return 0;
}