
#ifndef _Median_H
#define _Median_H

#include "TMath.h"
#include <vector>
#include <iostream>



class MedianCalculator
{
    private:
        int nDim_exact;
        std::vector<double> x;
        int n_conter;

        float eta;
        float mean_x;
        float mean_xx;
        float quantiles[3];       // 0.16, 0.5, 0.84
        float prob[3];
    public:
        MedianCalculator();
        ~MedianCalculator(){};

        void fill(double x);
        double getMean();
        double getMedian(float & medianwidth, int& N);
        double getMedianWidth();
        int getN()
        {
            return n_conter;
        }      
        double getRMS();
        void reset();
        void setEta(double iEta = 0.01)
        {
            eta  = iEta;
        }
        void setNExact( int n = 1000)
        {
            nDim_exact = n;
        }
};


























#endif