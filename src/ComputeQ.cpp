#include "LACTRecEvent.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include <string>
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMultiGraph.h"



int main( int argc, char ** argv)
{
    std::string out_file = "q.root";

    TFile* outroot = TFile::Open(out_file.c_str(), "recreate");
    TFile* gamma_root = TFile::Open(argv[1], "read");
    TFile* proton_root = TFile::Open(argv[2], "read");
    LACTRecEvent* grec = new LACTRecEvent();
    LACTRecEvent* prec = new LACTRecEvent();
    TTree* gamma_tree = (TTree*) gamma_root->Get("LactRectree");
    TTree* proton_tree = (TTree*)proton_root->Get("LactRectree");

    gamma_tree->SetBranchAddress("LactRecEvent", &grec);
    proton_tree->SetBranchAddress("LactRecEvent", &prec);
    TH1D* gamma_all = new TH1D("gamma_all", "gamma_all", 20, -1, 3);
    TH1D* gamma_event = new TH1D("gamma_e", "gamma_e", 20, -1, 3);
    TH1D* gamma_pass  = new TH1D("gamma_pass", "gamma_pass", 20, -1, 3);
    TH1D* proton_all = new TH1D("proton_all", "proton_all", 20, -1, 3);
    TH1D* proton_event = new TH1D("proton_e", "proton_e", 20, -1, 3);
    TH1D* proton_pass = new TH1D("proton_pass", "proton_pass", 20, -1, 3);

    for( int i = 0 ; i < gamma_tree->GetEntries(); i++ )
    {
        gamma_tree->GetEntry(i);
        gamma_all->Fill(log10(grec->MCenergy), grec->weight);
        if( grec->direction_error > 0)
        {
            gamma_event->Fill(log10(grec->MCenergy), grec->weight);
            if( grec->mrsl > -1.3 && grec->mrsl < 1.3 && grec->mrsw > -1.9 && grec->mrsw <0.6)
            {
                gamma_pass->Fill(log10(grec->MCenergy), grec->weight);
            }
        }
    }
    float x[20];
    float y[20];
    float y2[20];
    float y3[20];
    float z[20];

    for( int i = 1; i <= 20 ; i++)
    {
        x[i - 1] = gamma_event->GetBinCenter(i);
        y[i - 1] = gamma_pass->GetBinContent(i) / gamma_event->GetBinContent(i);
        y2[i - 1] = gamma_event->GetBinContent(i) / gamma_all->GetBinContent(i) * TMath::Pi() * 700 *  700;
        y3[i - 1] = gamma_pass->GetBinContent(i) / gamma_all->GetBinContent(i) * TMath::Pi() * 700 * 700;
        z[i - 1] = 0.1;
    }
    for( int i = 0 ; i < proton_tree->GetEntries(); i++ )
    {
        proton_tree->GetEntry(i);
        proton_all->Fill(log10(prec->MCenergy), prec->weight);
        proton_event->Fill(log10(prec->MCenergy), prec->weight);
        if( prec->direction_error > 0)
        {
            if( prec->mrsl > -1.3 && prec->mrsl < 1.3 && prec->mrsw > -1.9 && prec->mrsw <0.6)
            {
                proton_pass->Fill(log10(prec->MCenergy), prec->weight);
            }
        }
    }
    float px[20];
    float py[20];
    float py2[20];
    float py3[20];
    float pz[20];

    for( int i = 1; i <= 20 ; i++)
    {
        px[i - 1] = proton_event->GetBinCenter(i);
        py[i - 1] = proton_pass->GetBinContent(i) / proton_event->GetBinContent(i);
        py2[i - 1] = proton_event->GetBinContent(i) / proton_all->GetBinContent(i);
        py3[i - 1] = proton_pass->GetBinContent(i) / proton_all->GetBinContent(i);
        pz[i - 1] = 0.1;
    }
    outroot->cd();
    TCanvas* q = new TCanvas("qfactor", "Fraction");
    q->cd();
    TMultiGraph* mg2 = new TMultiGraph("mg", "Fraction");
    for( int i = 0; i< 20; i++)
    {
        std::cout << y[i] << "Fraction" << std::endl;
        std::cout << py[i] << "Fraction" << std::endl;
    }
    auto g1 = new TGraphErrors(20, x, y, z);
    auto g2 = new TGraphErrors(20, px, py, pz);
    mg2->Add(g1);
    mg2->Add(g2);
    

    mg2->Draw();
    outroot->cd();
    q->Write();
/*
    TCanvas* c2 = new TCanvas("c2", "EffectiveArea");
    c2->cd();
    TMultiGraph* mg = new TMultiGraph("mg", "Effective Area") ;
    auto  g3 = new TGraph(20, x, y2);
    auto  g4 = new TGraph(20, x, y3 );
    mg->Add(g3);
    mg->Add(g4);
    mg->Draw();

    c2->Write();
*/

    outroot->Write();
}