#ifndef _TABLE_CALCULATE_H
#define _TABLE_CALCULATE_H

#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
//#include "LACT_MedianCalculator.h"
#include "TFile.h"


using namespace std;

class LACTTableCalculator
{
    public:

        LACTTableCalculator();


        vector<vector< TH1F*> > Oh;
        //vector<vector< MedianCalculator*> > OMedian;
        TH2F* hMedian;
        TProfile2D* hmean;
        TFile* lookupfile;

        bool fEnergy = false;
        unsigned int fBinning1DxbinsN;
        int NumSize;
        int NumDist;
        float* fBinning1Dxbins;
        float MinShowerPerBin;
        bool Write1DHistogram;
        string parameter;



        LACTTableCalculator(bool iEnergy, string para);
        ~LACTTableCalculator();
        void SetCalEnergy(bool flag)
        {
            fEnergy = flag;
        }
        void FillLookupTable(int Ntel, vector<int>& goodimage, vector<double>& rp, vector<double>& size, vector<double>& fill_para, double weight);
        bool create1DHistogram( int i , int j);
        //bool createMedianApprox( int i, int j);
        double interpolate(TH2F* h, double x, double y, bool iError);
        void terminate(TFile* );
        

};




















#endif