#include "LACT_TableCalculatorData.h"
#include "LACT_TableCalculator.h"


LACT_TableCalculatorData::LACT_TableCalculatorData()
{
    fFillVariable = "";
    fEnergy = false;
}

//Be Careful! Now I only Considered the zenith bin! Actually For Now There only exit 1 bin 20deg
void LACT_TableCalculatorData::terminate(TFile* out_file)
{
    fTable->terminate(out_file);
}

void LACT_TableCalculatorData::InitTableCalculator()
{
    fTable = new LACTTableCalculator(fEnergy, fFillVariable);
}