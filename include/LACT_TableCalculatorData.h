#ifndef _TABLE_LOOKUP_H
#define _TABLE_LOOKUP_H

#include "TFile.h"
#include "TMath.h"

#include "LACTStatistics.h"
#include "LACT_TableCalculator.h"

#include <string>
#include <vector>

using namespace std;
class LACT_TableCalculatorData
{
    public:
        bool fEnergy;
        string fFillVariable;
        LACTTableCalculator* fTable;

        LACT_TableCalculatorData();
        ~LACT_TableCalculatorData(){};

        void InitTableCalculator();
        void terminate(TFile* );

};
























#endif