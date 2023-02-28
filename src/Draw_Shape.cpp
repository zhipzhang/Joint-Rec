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

    TH1D* MRSW = new TH1D("MRSW", "MRSW", 400, -10, 10);
    TH1D* MRSL = new TH1D("MRSL", "MRSL", 400, -10, 10);
    TFile* out_root = TFile::Open(lactrunpara->out_file.c_str(), "recreate");
    out_root->cd();
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
            if( lactrec->mrsl > -30 && lactrec->mrsw > -30)
            {
                MRSW->Fill(lactrec->mrsw, lactrec->weight);
                MRSL->Fill(lactrec->mrsl, lactrec->weight);
            }


        } 
        input_root->Close();

    }
    double sum_w = 0;
    double sum_l = 0;
    for( int i = 1; i <= MRSW->GetNbinsX(); i++)
    {
        sum_w += MRSW->GetBinContent(i);
        sum_l += MRSL->GetBinContent(i);
    }
    std::vector<double> width;
    std::vector<double> length;
    std::vector<double> wx;
    std::vector<double> wy;
    for( int i = 1; i <=MRSW->GetNbinsX(); i++)
    {
        width.push_back(MRSW->GetBinContent(i)/sum_w);
        wx.push_back(MRSW->GetBinCenter(i));
        length.push_back(MRSL->GetBinContent(i)/sum_l);
        wy.push_back(MRSL->GetBinCenter(i));
    }
    TGraph* g1 = new TGraph(wx.size(), &wx[0], &width[0]);
    g1->SetTitle("MRSW");
    TGraph* g2 = new TGraph(wy.size(), &wy[0], &length[0]);
    g2->SetTitle("MRSL");
    out_root->cd();
    MRSL->Write();
    MRSW->Write();
    g1->Write();
    g2->Write();
    out_root->Write();
   return 0;
}