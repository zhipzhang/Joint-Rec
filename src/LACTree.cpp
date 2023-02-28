/*
    define the method in LACTree.h
*/

#include "LACTree.h"
#include "TTree.h"
LACTree::LACTree()
{
    Mctree = 0;
    EventTree = 0;
    LTrig = 0;
    Ntrig = 0;
    ntel = 0;
    ntel_data = 0;
    runNumber = 0;
    eventNumber = 0;
    primary = 0;
    energy = 0;
    xcore = 0;
    ycore = 0;
    flag = 0;
    fadc_read_write = 0;
    fillPeLeaf = 0;
    resetDataVectors();
}
bool LACTree::initMctree()
{
    return true;
}
bool LACTree::initEventTree()
{
    EventTree = new TTree("event","event data tree");
    EventTree->SetMaxTreeSize(1000 * Long64_t(200000000));
    EventTree->Branch("ntel_data", &ntel_data, "ntel_data/i") ;
    EventTree->Branch("tel_data", tel_data,"tel_data[ntel_data]/i");
    EventTree->Branch("ntel ", &ntel,"ntel/i");
    EventTree->Branch("Paz", Point_Az, "Paz[ntel]/F");
    EventTree->Branch("Pal", Point_Al, "Pal[ntel]/F");
    EventTree->Branch("tel_az", &Tel_az, "tel_az/F");
    EventTree->Branch("tel_al", &Tel_al, "tel_al/F");
    EventTree->Branch("flag", &flag, "flag/I") ;
    char tname[100];
    EventTree->Branch("ltrig", &LTrig,"ltrig/i");
    EventTree->Branch("ntrig", &Ntrig,"ntrig/i");
    EventTree->Branch("ltrig_list", LTrig_list,"ltrig_list[ntrig]/i");
    EventTree->Branch( "runNumber", &runNumber, "runNumber/i" );
    EventTree->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    EventTree->Branch( "MCprim", &primary, "MCprimary/s" );
    EventTree->Branch( "MCe0", &energy, "MCenergy/F" );
    EventTree->Branch( "MCxcore", &xcore, "MCxcore/F" );
    EventTree->Branch( "MCycore", &ycore, "MCycore/F" );
    EventTree->Branch( "MCze", &ze, "MCze/F" );
    EventTree->Branch( "MCaz", &az, "MCaz/F" );
    EventTree->Branch("xmax", &xmax, "xmax/F");
    EventTree->Branch("hmax", &hmax, "hmax/F");
    EventTree->Branch("emax", &emax, "emax/F");
    EventTree->Branch("cmax", &cmax, "cmax/F");
    EventTree->Branch( "weight", &weight,"weight/F");
    if(fillPeLeaf)
    {
        sprintf(tname,"Pe[ntel_data][%d]/s",LACT_MAXPIXELS);
        EventTree->Branch("Pe",pe_list,tname);
    }
    EventTree->Branch( "Write_fadc", &fadc_read_write, "Write_fadc/B");
    if(fadc_read_write)
    {
        sprintf(tname,"fadc_HG[ntel_data][%d]/S",LACT_MAXPIXELS);
        EventTree->Branch("fadc_HG", fadc_HG, tname);
        sprintf(tname,"fadc_sum[ntel_data][%d]/F",LACT_MAXPIXELS);
        EventTree->Branch("fadc_sum", fadc_sums, tname);
        sprintf(tname,"fadc_num_samples[ntel_data]/S");
        EventTree->Branch("fadc_num_samples", fadc_num_samples, tname);
        sprintf(tname,"fadc_trace[ntel_data][%d][%d]/S",LACT_SUMWINDOW,LACT_MAXPIXELS);
        EventTree->Branch("fadc_traces", fadc_trace, tname);
        sprintf(tname,"fadc_pedestal[ntel_data][%d]/F",LACT_MAXPIXELS);
        EventTree->Branch("fadc_pedestal", fadc_pedestal, tname);
    }

    
    resetDataVectors();
    return true;
}

bool LACTree::initMatchTree(TTree* match)
{
    match->SetMaxTreeSize(1000 * Long64_t(200000000));
    match->Branch("ntel_data", &ntel_data, "ntel_data/i") ;
    match->Branch("tel_data", tel_data,"tel_data[ntel_data]/i");
    match->Branch("ntel ", &ntel,"ntel/i");
    match->Branch("Paz", Point_Az, "Paz[ntel]/F");
    match->Branch("Pal", Point_Al, "Pal[ntel]/F");
    match->Branch("tel_az", &Tel_az, "tel_az/F");
    match->Branch("tel_al", &Tel_al, "tel_al/F");
    match->Branch("flag", &flag, "flag/I") ;
    char tname[100];
    match->Branch("ltrig", &LTrig,"ltrig/i");
    match->Branch("ntrig", &Ntrig,"ntrig/i");
    match->Branch("ltrig_list", LTrig_list,"ltrig_list[ntrig]/i");
    match->Branch( "runNumber", &runNumber, "runNumber/i" );
    match->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    match->Branch( "MCprim", &primary, "MCprimary/s" );
    match->Branch( "MCe0", &energy, "MCenergy/F" );
    match->Branch( "MCxcore", &xcore, "MCxcore/F" );
    match->Branch( "MCycore", &ycore, "MCycore/F" );
    match->Branch( "MCze", &ze, "MCze/F" );
    match->Branch( "MCaz", &az, "MCaz/F" );
    match->Branch("xmax", &xmax, "xmax/F");
    match->Branch("hmax", &hmax, "hmax/F");
    match->Branch("emax", &emax, "emax/F");
    match->Branch("cmax", &cmax, "cmax/F");
    match->Branch( "weight", &weight,"weight/F");
    if( fillPeLeaf )
    {
        sprintf(tname,"Pe[ntel_data][%d]/s",LACT_MAXPIXELS);
        match->Branch("Pe",pe_list,tname);
    }

}



void LACTree::resetDataVectors(unsigned int iNMAX_TEL, unsigned int iNMAX_PIXELS)
{   
    if(iNMAX_TEL>=LACT_MAXTEL)
    {
        iNMAX_TEL=LACT_MAXTEL;
    }
    if(iNMAX_PIXELS>=LACT_MAXPIXELS)
    {
        iNMAX_PIXELS=LACT_MAXPIXELS;
    }

    flag = 0;
    LTrig = 0;
    Ntrig = 0;
    for(unsigned int i=0; i<iNMAX_TEL; i++)
    {
        LTrig_list[i] = 0;
        LTime[i] = 0.;
        Point_Al[i] = 0.;
        Point_Az[i] = 0.;
    }
    memset(pe_list, 0, LACT_MAXTEL*LACT_MAXPIXELS*sizeof(pe_list[0][0]));
    memset(fadc_pedestal, 0., LACT_MAXPIXELS*LACT_MAXTEL*sizeof(fadc_pedestal[0][0]));
    memset(fadc_sums, 0., LACT_MAXTEL * LACT_MAXPIXELS * sizeof(fadc_sums[0][0]));
    memset(fadc_num_samples, 0, LACT_MAXTEL* sizeof(fadc_num_samples[0]));
    memset(fadc_trace, 0, LACT_MAXTEL * LACT_SUMWINDOW *LACT_MAXPIXELS * sizeof(fadc_trace[0][0][0]));
}

bool LACTree::initEventTree(TTree *t)
{
    EventTree = t;
    if(EventTree && EventTree->GetEntries() == 0)
    {
        std::cout << "Event Tree Have no Entries" << std::endl;
    }

    if(EventTree->GetBranch("MCe0"))
    {
        fMC = true;
    }

    EventTree->SetBranchAddress("runNumber", &runNumber);
    EventTree->SetBranchAddress("eventNumber", &eventNumber);
    EventTree->SetBranchAddress("ntel", &ntel);
    
    EventTree->SetBranchAddress("Paz", Point_Az);
    EventTree->SetBranchAddress("Pal", Point_Al);
    EventTree->SetBranchAddress("ltrig", &LTrig);
    EventTree->SetBranchAddress("ntrig", &Ntrig);
    EventTree->SetBranchAddress("ltrig_list", &LTrig_list);
    EventTree->SetBranchAddress("ntel_data", &ntel_data);
    EventTree->SetBranchAddress("tel_data", tel_data);
    EventTree->SetBranchAddress("Write_fadc", &fadc_read_write);
    EventTree->SetBranchAddress("weight", &weight);
    EventTree->SetBranchAddress("tel_az", &Tel_az);
    EventTree->SetBranchAddress("tel_al", &Tel_al);
    EventTree->SetBranchAddress("flag", &flag);

    if(EventTree->GetBranch("Pe"))
    {
        EventTree->SetBranchAddress("Pe", pe_list);
        setFillPEleaf(true);
    }
    
    if( fMC )
    {
        EventTree->SetBranchAddress( "MCprim", &primary );
        EventTree->SetBranchAddress( "MCe0", &energy );
        EventTree->SetBranchAddress( "MCxcore", &xcore );
        EventTree->SetBranchAddress( "MCycore", &ycore );
        EventTree->SetBranchAddress( "MCze", &ze );
        EventTree->SetBranchAddress( "MCaz", &az );
        EventTree->SetBranchAddress("hmax", &hmax);
        EventTree->SetBranchAddress("xmax", &xmax);
    }



}
