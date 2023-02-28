#include "LACT_TableCalculator.h"
#include "LACT_MedianCalculator.h"
#include "LACTStatistics.h"
#include "TH3.h"
#include "TAxis.h"
#include <cstdio>

LACTTableCalculator::LACTTableCalculator(bool iEnergy, string para)
{
    parameter = para;
    NumSize = 55;
    NumDist = 50;
    char hname[200];
    char htitle[200];
    sprintf(hname, "%s_median", para.c_str());
    sprintf(htitle, "%s vs log10(size) vs distance", para.c_str());
    hMedian = new TH2F(hname, htitle, 55, 1.5, 7, 50, 0, 750);
    hMedian->SetXTitle(" log_{10} size");
    hMedian->SetYTitle("distance[m]");
    if( !iEnergy )
    {
        sprintf(htitle, "%s (median) [deg]", para.c_str());
    }
    else 
    {
        sprintf(htitle, "%s (median) [pe/TeV]", para.c_str());
    }
    hMedian->SetZTitle(htitle);

    sprintf(hname, "%s_mean.", para.c_str());
    sprintf(htitle, "%s vs log10(size) vs distance", para.c_str());
    hmean = new TProfile2D(hname, htitle, 55, 1.5, 7, 50, 0, 750);
    hmean->SetXTitle("log_{10} size");
    hmean->SetYTitle("distance [m]");
    if( !iEnergy )
    {
        sprintf(htitle, "%s (mean)[deg]", para.c_str());
    }
    else 
    {
        sprintf(htitle, "%s (mean) [TeV]", para.c_str());
    }
    hmean->SetZTitle(htitle);
    SetCalEnergy(iEnergy);
    Oh.resize(NumSize);
    for( int i = 0; i < NumSize; i++)
    {
        Oh[i].resize(NumDist);
    }
}

bool LACTTableCalculator::create1DHistogram(int i, int j)
{
    if( i >= 0  && j >= 0 && i < (int) Oh.size() && j < (int)Oh[i].size() && !Oh[i][j])
    {
        char hisname[200];
        char histitle[200];
        int id = i * 1000 +j;

        sprintf(hisname, "h%s%d",parameter.c_str(), id);
        double is1 = hMedian->GetXaxis()->GetBinLowEdge(i + 1);
        double is2 = hMedian->GetXaxis()->GetBinLowEdge(i + 1) + hMedian->GetXaxis()->GetBinWidth(i + 1);
        double id1 = hMedian->GetYaxis()->GetBinLowEdge(j + 1);
        double id2 = hMedian->GetYaxis()->GetBinLowEdge(j + 1) + hMedian->GetYaxis()->GetBinWidth(j + 1);

        sprintf( histitle, "%.2f < log10 size < %.2f, %.1f < r < %.1f ", is1, is2, id1, id2);
        
        Oh[i][j] = new TH1F(hisname, histitle, 1000, 0, 0.03);
       // Oh[i][j]->GetXaxis()->SetCanExtend(true);

    }
    else 
    {
        return false;
    }
    return true;

}
/*
bool LACTTableCalculator::createMedianApprox(int i, int j)
{
    if( i >= 0  && j >= 0 && i < (int) OMedian.size() && j < (int)OMedian[i].size() && !OMedian[i][j])
    {
        OMedian[i][j] = new MedianCalculator();
    } 
    else 
    {
        return false;
    }
    return true;

}
*/


/*
    Use For Reading !
*/
double LACTTableCalculator::interpolate( TH2F* h, double x, double y, bool iError)
{
    if( !h )
    {
        return -999;
    }

    int i_x = h->GetXaxis()->FindFixBin( x );
    int i_y = h->GetYaxis()->FindBin(y);

    if( i_x == 0 || i_y == 0 || i_x == h->GetNbinsX() || i_y == h->GetNbinsY())
    {
        if(iError)
        {
            return h->GetBinError(i_x, i_y);
        }
        else 
        {
            return h->GetBinContent(i_x, i_y);
        }
    }
    if( x < h->GetXaxis()->GetBinCenter(i_x))
    {
        i_x --;
    }
    if( y < h->GetYaxis()->GetBinCenter(i_y))
    {
        i_y --;
    }
    double e1 = 0.;
    double e2 = 0.;
    double v = 0.;
    
    // first interpolate on distance axis, then on size axis
    if( !iError )
    {
        e1 = Statistics::interpolate( h->GetBinContent( i_x, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinContent( i_x, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        e2 = Statistics::interpolate( h->GetBinContent( i_x + 1, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinContent( i_x + 1, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        v = Statistics::interpolate( e1, h->GetXaxis()->GetBinCenter( i_x ),
                                      e2, h->GetXaxis()->GetBinCenter( i_x + 1 ),
                                      x, false, 0.5, 1.e-5 );
    }
    else
    {
        e1 = Statistics::interpolate( h->GetBinError( i_x, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinError( i_x, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
        e2 = Statistics::interpolate( h->GetBinError( i_x + 1, i_y ), h->GetYaxis()->GetBinCenter( i_y ),
                                       h->GetBinError( i_x + 1, i_y + 1 ), h->GetYaxis()->GetBinCenter( i_y + 1 ),
                                       y, false, 0.5, 1.e-5 );
                                       
        v = Statistics::interpolate( e1, h->GetXaxis()->GetBinCenter( i_x ), e2, h->GetXaxis()->GetBinCenter( i_x + 1 ), x, false, 0.5, 1.e-5 );
    }
    // final check on consistency of results
    // (don't expect to reconstruct anything below 1 GeV)
    if( e1 > 1.e-3 && e2 < 1.e-3 )
    {
        return e1;
    }
    if( e1 < 1.e-3 && e2 > 1.e-3 )
    {
        return e2;
    }
    
    return v;
}

/*
    Write the Histogram to Disk

*/ 
void LACTTableCalculator::terminate(TFile* out_file)
{
    if(!out_file->cd())
    {
        cout << "Error :Unable to cd out_file" << endl;
        return;
    }
    float med = 0.;
    float sigma = 0.;
    int   nevents = 0;
    double i_a[] = {0.16, 0.5, 0.84};
    double i_b[] = {0.0, 0.0, 0.0};
    for( int i = 0; i < NumSize; i++)
    {
        for( int j = 0; j < NumDist; j++)
        {
            if( Oh[i][j] && Oh[i][j]->GetEntries() > MinShowerPerBin)
            {
                Oh[i][j]->GetQuantiles(3, i_b, i_a);
                med   = i_b[1];
                sigma = i_b[2] - i_b[0];
                nevents = Oh[i][j]->GetEntries();
            }
            else 
            {
                med = 0.;
                sigma = 0.;
                nevents = 0.;
            }
            hMedian->SetBinContent(i + 1, j + 1, med);
            hMedian->SetBinError(i + 1, j + 1, sigma);
            if( Oh[i][j] && Oh[i][j]->GetEntries() > 0)
            {
                out_file->cd();
                Oh[i][j]->Write();
                delete Oh[i][j];
            }
        }
    }
        out_file->Write();
}

void LACTTableCalculator::FillLookupTable(int Ntel, vector<int> &goodimage, vector<double> &rp, vector<double> &size, vector<double> &fill_para, double weight)
{
    for( int i = 0; i < Ntel; i++)
    {
        if(goodimage[i] >= 2)
        {
            int ir = hMedian->GetYaxis()->FindFixBin(rp[i]);
            if( ir >= NumDist)
            {
                continue;
            }
            double ilogsize = log10(size[i]);
            int is = hMedian->GetXaxis()->FindFixBin(ilogsize);
            if( is >= NumSize)
            {
                continue;
            }
            if( ir >=0 && is >= 0 )
            {
                if(!Oh[is][ir] && !create1DHistogram(is, ir))
                {   
                    continue;
                }
                Oh[is][ir]->Fill(fill_para[i]);
            }

           /*
            if( !OMedian[is][ir])
            {
                if( !createMedianApprox(is, ir))
                {
                    continue;
                }

            }
            OMedian[is][ir]->fill(fill_para[i]);
            */

            hmean->Fill(ilogsize, rp[i], fill_para[i], weight);

        }
    }
}