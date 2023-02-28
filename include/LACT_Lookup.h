#include <string>
#include "TTree.h"
#include "TFile.h"
#include "LACTRecEvent.h"
#include "LACT_TableCalculator.h"
#include "LACT_TableCalculatorData.h"
#include "LACT_MedianCalculator.h"


class LACTLookup
{
    public:
    enum Value{MRSW, MRSL, EREC};
    LACTRecEvent* lactrec;

    LACTTableCalculator* TableCalculator;
    vector<LACT_TableCalculatorData*> TableCalculatorData;

    TFile* lookupfile;
    TTree* LactRecTree;
    LACTLookup(std::string name);
    ~LACTLookup();
    void GetData(TTree* );
    void InitLookupTableData();
    void Loop();
    void terminate();

};