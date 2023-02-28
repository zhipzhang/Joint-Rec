/*

            Main Program Generate Lookup Tables

*/
#include "LACT_Lookup.h"
#include "LACT_RunPara.h"
#include "TFile.h"
#include "TTree.h"
#include "LACTRecEvent.h"
#include <cstdlib>
#include <iostream>
int main(int argc, char** argv)
{
    LACT_RUNPARA* lact_runpara = new LACT_RUNPARA();
    lact_runpara->ProcessCommandLine(argc, argv, "true");
    LACTLookup* lact_lookup = new LACTLookup(lact_runpara->GetLookupName());

    for( auto input : lact_runpara->input_file)
    {
        std::cout <<"Open File "<<input <<std::endl;
        TFile* input_root = TFile::Open(input.c_str(), "read");
        if(input_root->IsZombie() || !input_root)
        {
            exit(EXIT_FAILURE);
        }
        TTree* input_tree = (TTree*) input_root->Get("LactRectree");
        if( !input_tree)
        {
            exit(EXIT_FAILURE);
        }
        lact_lookup->GetData(input_tree);
        lact_lookup->InitLookupTableData();
        lact_lookup->Loop();
        input_root->Close();
    }
    lact_lookup->terminate();
}

