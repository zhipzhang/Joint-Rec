#include "LACT_Lookup.h"
#include "LACT_TableCalculatorData.h"
#include "TDirectory.h"
#include <vector>
LACTLookup::LACTLookup(std::string name)
{
    lookupfile = TFile::Open(name.c_str(), "recreate");
}
void LACTLookup::GetData(TTree* RecTree)
{
    LactRecTree = RecTree;
    lactrec = new LACTRecEvent();
    RecTree->SetBranchAddress("LactRecEvent", &lactrec);

}

void LACTLookup::InitLookupTableData()
{

    TableCalculatorData.push_back(new LACT_TableCalculatorData());
    TableCalculatorData[MRSW]->fEnergy = false;
    TableCalculatorData[MRSW]->fFillVariable = "width";
    TableCalculatorData[MRSW]->InitTableCalculator();

    TableCalculatorData.push_back(new LACT_TableCalculatorData());
    TableCalculatorData[MRSL]->fEnergy = false;
    TableCalculatorData[MRSL]->fFillVariable = "length";
    TableCalculatorData[MRSL]->InitTableCalculator();

    TableCalculatorData.push_back(new LACT_TableCalculatorData());
    TableCalculatorData[EREC]->fEnergy = true;
    TableCalculatorData[EREC]->fFillVariable = "size/Energy";
    TableCalculatorData[EREC]->InitTableCalculator();
    int zenithbin = 0;

}


void LACTLookup::terminate()
{
    for(auto Data: TableCalculatorData)
        Data->terminate(lookupfile);
}

void LACTLookup::Loop()
{
    if( LactRecTree)
    {
        std::cout << "Enrties is " << LactRecTree->GetEntries() << std::endl; 
        for( int i = 0; i < LactRecTree->GetEntries(); i++)
        {
            LactRecTree->GetEntry(i);
            std::vector<double> size_E;
            if( lactrec->direction_error < 0)
            {
                continue;
            }
            for( int j = 0; j < lactrec->GetNtel(); j++)
            {
                size_E.push_back(lactrec->GetTelSize(j)/ lactrec->MCenergy);
            }
            if( lactrec->GetNtel() > 0)
            {
                    TableCalculatorData[MRSW]->fTable->FillLookupTable(lactrec->GetNtel(), lactrec->good_image, lactrec->rec_rp, lactrec->size, lactrec->length, lactrec->weight);
                    TableCalculatorData[MRSL]->fTable->FillLookupTable(lactrec->GetNtel(), lactrec->good_image, lactrec->rec_rp, lactrec->size, lactrec->width, lactrec->weight);
                   // TableCalculatorData[EREC]->fTable->FillLookupTable(lactrec->GetNtel(), lactrec->good_image, lactrec->rec_rp, lactrec->size, size_E, lactrec->weight);
            }
            size_E.clear();
            std::cout << "End Entry " << i <<std::endl;


        }
    }
}