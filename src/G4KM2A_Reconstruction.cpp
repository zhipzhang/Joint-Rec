#include "G4KM2A_Reconstruction.h"
#include "G4KM2A_Geometry.h"
#include "G4KM2A_RecE.h"
#include "KM2ARecEvent.h"
#include "LHEvent.h"
#include "RtypesCore.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TDecompSVD.h"
#include "TGraph.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TMatrixDfwd.h"
#include "TMinuit.h"
#include "TVectorDfwd.h"
#include <cmath>
#include <cstring>
#include "TAxis.h"
#include "TLatex.h"
#include "TPaveText.h"

#define PETH 0.2
G4KM2A_Geometry * G4KM2A_Reconstruction::geom = 0;
TMinuit* G4KM2A_Reconstruction::minuit = new TMinuit(4);
G4KM2A_Reconstruction* G4KM2A_Reconstruction::m_myself = 0;

float  G4KM2A_Reconstruction::nED[] = {0};
float  G4KM2A_Reconstruction::nMD[] = {0};
float  G4KM2A_Reconstruction::nWCDA[] = {0};
double G4KM2A_Reconstruction::_np   = 0;
double G4KM2A_Reconstruction::_dirl = 0;
double G4KM2A_Reconstruction::_dirm = 0;
double G4KM2A_Reconstruction::_dirn = 0;
char G4KM2A_Reconstruction::_style[20]="ED";

G4KM2A_Reconstruction::G4KM2A_Reconstruction(int tFlag)
{
    int ned, nwc, i, nmd;
    geom = G4KM2A_Geometry::GetInstance(tFlag);
    ned = geom->GetNED();
    nmd = geom->GetNMD();
 
    for( i = 0; i < ned; i++)  
        nED[i] = 0;
    for( i = 0; i < nmd; i++)
        nMD[i] = 0;
    
}




int G4KM2A_Reconstruction::Init(KM2ARecEvent *trec){
    trec->rec_theta=-10;trec->rec_phi=-10;
    trec->rec_Etheta_p=-10;  trec->rec_Ephi_p=-10;
    trec->rec_Wtheta_p=-10;  trec->rec_Wphi_p=-10;
    trec->rec_Etheta_c=-10;  trec->rec_Ephi_c=-10;
    trec->rec_Wtheta_c=-10;  trec->rec_Wphi_c=-10;

    trec->rec_x=-10000  ; trec->rec_y=-10000;
    trec->rec_Ex=-10000 ; trec->rec_Ey=-10000;  trec->rec_Ez=-10000;
    trec->rec_Wx=-10000 ; trec->rec_Wy=-10000;  trec->rec_Wz=-10000;

    trec->rec_Esize=-10; trec->rec_Wsize=-10;
    trec->rec_Eage=-10;  trec->rec_Wage=-10;
    trec->rec_sigma=-1;
    trec->rec_a=0; 

    trec->NtrigE=-1 ;  trec->NtrigW=-1;
    trec->NfiltE=-1 ;  trec->NfiltM=-1;  trec->NfiltW=-1;
    trec->NfiltE2 = -1;
    trec->NpW=-1;
    trec->NpE1 = trec->NpE2 = trec->NpE3 = trec->NpE4 = trec->NpE5 = trec->NpE6 = trec->NpE6 = trec->NpE6 = trec->NpE6 = trec->NpE6 = trec->NpE6 = -1;
    trec->NuM1 =  trec->NuM2 = trec->NuM3 = trec->NuM4 = trec->NuM5 = trec->NuM6 = -1;
    trec->NuW1=-1; trec->NuW3=-1;trec->NuW3=-1;
    trec->NfiltE =-1; trec->NfiltM =-1; trec->NfiltW =-1;
    trec->Direction_Error = -1;
   return 0;
}
/*
    Convert Pe to particle number
    ED 40pe : 1 electrons
    MD 75pe : 1 Muons
*/
void G4KM2A_Reconstruction::setnparticle(TClonesArray &tHits, int np, double pe)
{
    LHHit *tHit;
    double ne;
    for( int i = 0; i < np ; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if( tHit->GetStatus() > -1 )
        {
            ne = tHit->GetPe()/pe;
            if ( ne < PETH)
            {
                ne = 0;
            }
            tHit->SetPe(ne);
        }
    }
}



/*
    Use Twind and rwind to filter noise
    Not Passed Cut will Set Status 0
*/
int G4KM2A_Reconstruction::spacetimefilter(TClonesArray &tHits, int np, int twind, int rwind, const char *style, int flag)
{
    LHHit *tHit;
    int t;
    int Nall = 0;
    int Nbin = 10000;
    int trigger_num[10000]{0};
    for(int  i = 0; i < np; i++)
    {
        tHit = (LHHit*)((tHits)[i]);
        if(tHit->GetStatus() > -1)
        {
            if( tHit->GetPe() < 1.e-3)
            {
                tHit->SetStatus(0); //  Status 0: Mean no particle 
            }
            else 
            {
                tHit->SetStatus(5); //  Status 5: Mean good
                t = (int) tHit->GetTime();
                if(t > 0 && t < Nbin)
                {
                    trigger_num[t]++;
                    Nall++;
                }

            }
        }
    }

    if( Nall < 4 )
    {
        return 0;                      // No Trigger
    }

    double nums = 0;
    int    num_ed  = 0;
    double max_nums = 0;
    int max_id = 0;
    int max_t = 0; 
    double x1,y1,z1;
    double x2,y2,z2;
    double r;
    int    id;
    for( int tflag = 0; tflag <= Nbin-twind; tflag++)
    {
        if( trigger_num[tflag] == 0 && tflag != (Nbin - twind))
        {
            continue;
        }
        nums = 0;
        for( int j =0 ;j < np; j++)
        {
            tHit = (LHHit *)((tHits)[j]);
            if( tHit->GetStatus() > 0 && tHit->GetTime() > tflag && tHit->GetTime() < tflag + twind)
            {
                if(geom->Getxyz(tHit->GetId(), x1, y1, z1, 1, style) < 0)
                {
                    continue;
                }
            }
            id = tHit->GetId();
            nums = 0;
            num_ed  = 0;
            for ( int k = 0; k < np; k++)
            {
                tHit = (LHHit *)((tHits)[k]);
                if( tHit->GetStatus() > 0 && tHit->GetTime() > tflag && tHit->GetTime() < tflag + twind)
                {
                    if(geom->Getxyz(tHit->GetId(), x2, y2, z2, 1, style) < 0)
                    {
                        continue;
                    }
                    r = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
                    if ( r < rwind )
                    {
                        nums += tHit->GetPe();
                        num_ed++;

                    }
                }
            }
            if( nums > max_nums && num_ed > 5)
            {
                max_nums = nums;
                max_id = id;
                max_t = tflag;
            }
            if( num_ed > max_nums)
            {
                max_nums = num_ed;
                max_id   = id;
                max_t    = tflag;
            }

        }

    }
    double x,y,z;
    double tx,ty,tz;
    if( max_id > 0 && geom->Getxyz(max_id, x, y, z, 1, style) > 0)
    {
        for( int i = 0; i < np; i++)
        {
            tHit = (LHHit*) ((tHits)[i]);
            if( tHit->GetStatus() > 0)
            {
                if( tHit->GetTime() < max_t || tHit->GetTime() >= max_t + twind)
                {
                    tHit->SetStatus(0);
                }
                else 
                {
                    if(geom->Getxyz(tHit->GetId(), tx, ty, tz, 1, style) > 0)
                    {
                        r = sqrt(pow(tx -x ,2) + pow( ty - y, 2));
                        if( r > rwind)
                        {
                            tHit->SetStatus(0);
                        }

                    }
                }
            }
        }

    }
    return (int) max_nums;

}

/*
    np     : The number of ED Detector which have data
    twind  : time window for trigger
    style  : May Can be changed to Trigger WCDA
    th     : Threshold in the twind
*/
int G4KM2A_Reconstruction::trigger(TClonesArray &tHits, int np, int twind, const char *style, int th)
{
    LHHit* tHit;
    int Nbin = 10000;
    int trigger_num[10000]{0};
    int  t;

    for( int i = 0; i < np; i++)
    {
        tHit = (LHHit*)((tHits[i]));
        if( tHit->GetStatus() > -1) // Status -1 means detector not work well
        {
            t = int(tHit->GetTime());
            if( t > 0 && t < 10000)
            {
                trigger_num[t]++;
            }
        }
    }

    // Total Time is 10000, test whether it can trigger in any twind 
    int nums = 0; // record the number in twind
    int max_nums = 0;
    for( int tflag = 0; tflag <= Nbin - twind; tflag++)
    {
        if( trigger_num[tflag] == 0 && tflag != (Nbin - twind))
        {
            continue;
        }
        nums = 0;
        for( int  j = tflag; j < (tflag + twind); j++)
        {
            nums += trigger_num[j];
        }
        if( nums > max_nums)
        {
            max_nums = nums;
        }

    }

    return max_nums;
}

int G4KM2A_Reconstruction::planarfit(TClonesArray &tHits, int np, KM2ARecEvent *trec, const char *style)
{
    if( np < 4)
    {
        return -1;
    }
    TMatrixD a(3,3);
    TVectorD b(3), r(3);
    LHHit *tHit;
    int nus;
    double x, y, z, sigmasum;
    double dt, w, theta, phi;
    double c[3], rr[3] = {0, 0, 0};
    double costheta = 0.94, widt = 10;
    double sigma;
    int Flag = 0;
    for( int itel = 0; itel < 15; itel++)
    {
        for( int i = 0; i < 3; i++)
        {
            for( int j = 0; j < 3; j++)
            {
                a[i][j] = 0;
            }
            b[i] = 0;

        }
        nus = 0; sigmasum = 0;
        for( int i = 0; i< np ; i++)
        {
            tHit = (LHHit*)((tHits)[i]);
            if (tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
            {
                continue;
            }
            if( geom->Getxyz(tHit->GetId(), x, y, z, 1, style) < 0)
            {
                continue;
            }
            if( itel == 0)
            {
                w = 1;
            }
            else 
            {
                dt = tHit->GetTime() - (rr[0] * x + rr[1] * y  + rr[2] - costheta * z)/0.2998;
                w  = exp(-0.5 * pow(dt/widt, 2));
                sigmasum += dt * dt;
                if( itel > 8)
                {
                    if( fabs(dt) > 3 * sigma && fabs(dt) > 30 )
                    {
                        tHit->SetStatus(1);
                    }
                }
            }
            nus++;
            c[0] = x * w;
            c[1] = y * w;
            c[2] = 1 * w;
            for( int k = 0; k < 3; k++)
            {
                for( int m = 0; m < 3; m++)
                {
                    a[m][k] += c[m] * c[k];
                }
                b[k] += ((tHit->GetTime()) * 0.2998 + costheta * z) * w * c[k];
            }

        }

        if( nus < 3)
        {
            break;
        }

        TDecompSVD svd(a);
        Bool_t ok;
        r = svd.Solve(b, ok);
        if( !ok )
        {
            break;
        }
        if( (r[0] * r[0] + r[1] * r[1]) > 1)
        {
            break;
        }
        if( ok )
        {
            sigma = sqrt(sigmasum/(nus -2)); 
            Flag  = 1;
            costheta = sqrt(1 - (r[0] * r[0] + r[1] * r[1]));
            rr[0] = r[0];
            rr[1] = r[1];
            rr[2] = r[2];
        }
    }
    dt = rr[0] * rr[0] + rr[1] * rr[1];
    if( Flag > 0.5 && dt <= 1. && dt > 0.)
    {
        theta = asin(sqrt(dt));
        phi   = atan2(rr[1], rr[0]);
        if( phi < 0)
        {
            phi += 2 * TMath::Pi();
        }
    }
    if( strcmp(style,"ED") == 0 )
    {
        trec->rec_Etheta_p = theta;
        trec->rec_Ephi_p   = phi;
        trec->rec_Ec0_p    = rr[2];
        trec->rec_Esigma_p = sigma;
    }
    if( theta > -1)
    {
        return  0;
    }
    else 
    {
        return -1;
    }

}


int G4KM2A_Reconstruction::eventrecline(LHEvent *tevent, KM2ARecEvent *trec)
{
    TClonesArray *HitsE, *HitsM;
    HitsE = tevent->GetHitsE();
    HitsM = tevent->GetHitsM();
    Init(trec);

    int twind = 400;
    int rwind = 100;
    int ned = geom->GetNED();
    trec->NtrigE = trigger(*HitsE, tevent->GetNhitE(), twind, "ED", 15);
    trec->NhitE  = tevent->GetNhitE();
    trec->NhitM  = tevent->GetNhitM();
    trec->NtrigW = tevent->GetNtrigE();      //no meaningful 

    trec->dt     = tevent->GetDt();
    trec->id     = tevent->GetId();
    trec->E      = tevent->GetE();
    trec->phi    = tevent->GetPhi();
    trec->theta  = tevent->GetTheta();
    trec->corex  = tevent->GetCorex();
    trec->corey  = tevent->GetCorey();
    
    trec->pNpE  =tevent->GetNpE();         // Particle number before entering the ED 
    trec->pNuM  =tevent->GetNuM();
    trec->pNuW  =tevent->GetNuW();
    trec->pNuW2 =tevent->GetNuW2();

    setnparticle(*HitsE, trec->NhitE, 40.);
    setnparticle(*HitsM, trec->NhitM, 75.);

    trec->NfiltE2 = spacetimefilter(*HitsE, trec->NhitE, twind, rwind, "ED", 0);
    if( trec->NfiltE2 >=3 && planarfit(*HitsE, trec->NhitE, trec, "ED") == 0)
    {
        trec->rec_theta = trec->rec_Etheta_p;
        trec->rec_phi   = trec->rec_Ephi_p;
        trec->rec_c0    = trec->rec_Ec0_p;
        trec->rec_sigma = trec->rec_Esigma_p;
        trec->rec_a     = 0;
    
    
        core_centre2(*HitsE, trec->NhitE, rwind, trec, "ED");
        double dr=trec->rec_Ez*tan(trec->rec_theta);
		trec->rec_x=trec->rec_Ex+dr*cos(trec->rec_phi);
		trec->rec_y=trec->rec_Ey+dr*sin(trec->rec_phi);
        if( conicalfit(*HitsE, trec->NhitE, trec, 0.035, "ED") == 0)
        {
            double alpha = 0.035;
            trec->rec_theta = trec->rec_Etheta_c;
            trec->rec_phi   = trec->rec_Ephi_c;
            trec->rec_a     = trec->rec_Ea;
            trec->rec_c0    = trec->rec_Ec0_c;
            trec->rec_sigma = trec->rec_Esigma_c;

            trec->NfiltE = noisefilter(*HitsE, trec->NhitE, -50, 100, rwind + 100, trec, "ED"); // Time : -50 ~ 100 ns Rwind 100
            if( trec->NfiltE > 10)
            {
                if( trec->NfiltE > 100 && trec->NfiltE < 200)
                {
                    noisefilter(*HitsE, trec->NhitE, -50, 100, rwind + 200, trec, "ED");
                }
                if( trec->NfiltE > 200)
                {
                    alpha = -1;
                    noisefilter(*HitsE, trec->NhitE, -50, 100, rwind + 300, trec, "ED");
                }

                core_likelihood(*HitsE, trec->NhitE, trec, "ED", 0);
                trec->rec_x  = trec->rec_Ex;
                trec->rec_y  = trec->rec_Ey;
            }

            if( conicalfit(*HitsE, trec->NhitE, trec, alpha, "ED") == 0)
            {
                trec->rec_theta = trec->rec_Etheta_c;
                trec->rec_phi   = trec->rec_Ephi_c;
                trec->rec_a     = trec->rec_Ea;
                trec->rec_c0    = trec->rec_Ec0_c;
                trec->rec_sigma = trec->rec_Esigma_c;
            }
            noisefilter(*HitsE, trec->NhitE, -30, 50, rwind, trec, "ED");

            /*
                NpE1 :  t :-30ns ~ 50ns r: 0m  ~ 100m
                NpE2 :  t :-30ns ~ 50ns r: 40m ~ 100m
                NpE3 :  t :-30ns ~ 50ns r：40m ~ 200m
                NpE4 :  t :-30ns ~ 50ns r: 0m  ~ 400m
                NpE5 :  t :-30ns ~ 50ns r: 0m  ~ 600m
                NpE6 :  t :-30ns ~ 50ns r: 0m  ~ 800m
            */
            trec->NpE1 = getnp(*HitsE, trec->NhitE, 0, trec, "ED");
            trec->NpE2 = getnp(*HitsE, trec->NhitE, 40, trec, "ED");

            noisefilter(*HitsE, trec->NhitE, -30, 50, 200, trec, "ED");
            trec->NpE3 = getnp(*HitsE, trec->NhitE, 40, trec, "ED");

            noisefilter(*HitsE, trec->NhitE, -30, 50, 400, trec, "ED");
            trec->NpE4 = getnp(*HitsE, trec->NhitE, 0, trec, "ED");

            noisefilter(*HitsE, trec->NhitE, -30, 50, 600, trec, "ED");
            trec->NpE5 = getnp(*HitsE, trec->NhitE, 0, trec, "ED");

            noisefilter(*HitsE, trec->NhitE, -30, 50, 800, trec, "ED");
            trec->NpE5 = getnp(*HitsE, trec->NhitE, 0, trec, "ED");


            noisefilter(*HitsE, trec->NhitE, -30, 50, 300, trec, "ED");
            
            /*
                NuM1 :  t :-30ns ~ 50ns r: 15m ~ 200m
                NuM2 :  t :-30ns ~ 50ns r: 40m ~ 200m
                NuM3 :  t :-30ns ~ 50ns r: 15m ~ 400m
                NuM4 :  t :-30ns ~ 50ns r: 40m ~ 400m
                NuM5 :  t :-30ns ~ 50ns r: 15m ~ 600m
                NuM6 :  t :-30ns ~ 50ns r: 15m ~ 800m
            */
            
            noisefilter(*HitsM, trec->NhitM, -30, 50, 200, trec, "MD");
            trec->NuM1 = getmuM(*HitsM, trec->NhitM, 15, trec);
            trec->NuM2 = getmuM(*HitsM, trec->NhitM, 40, trec);

            noisefilter(*HitsM, trec->NhitM, -30, 50, 400, trec, "MD");
            trec->NuM3 = getmuM(*HitsM, trec->NhitM, 15, trec);
            trec->NuM4 = getmuM(*HitsM, trec->NhitM, 40, trec);
            
            noisefilter(*HitsM, trec->NhitM, -30, 50, 600, trec, "MD");
            trec->NuM5 = getmuM(*HitsM, trec->NhitM, 15, trec);

            noisefilter(*HitsM, trec->NhitM, -30, 50, 800, trec, "MD");
            trec->NuM6 = getmuM(*HitsM, trec->NhitM, 15, trec);

            trec->Rec_Erho = recer50new3(trec->rec_Eage, trec->rec_Esize, trec->rec_theta);
            trec->Direction_Error = TMath::RadToDeg() * angle_between(trec->phi, TMath::Pi()/2 - trec->theta, trec->rec_phi, TMath::Pi()/2 - trec->rec_theta);


        }
    }
    return 0;
}

//get the sum number of particles for ED array
float G4KM2A_Reconstruction::getnp(TClonesArray &tHits,int np, int rcut, KM2ARecEvent *trec, const char *style){
    if( np < 4)
        return 0.;
    LHHit *tHit;
    int i;
    float dn=0.;
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus() > 2.5)
        {
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)>0)
            {
                dx = x-trec->rec_x;
                dy = y-trec->rec_y;
                dz = z;
                dr = sqrt(dx*dx + dy*dy + dz*dz-pow(dx*dirl + dy*dirm - dz*dirn,2.));
                if(dr > rcut) 
                    dn += tHit->GetPe();
            }
        }
   }
   return dn;
}
//get the particle density in certain region 
float G4KM2A_Reconstruction::getdensity(TClonesArray &tHits,int np, int rcut1, int rcut2, KM2ARecEvent *trec, const char *style){
    if(np<4)return 0.;
    LHHit *tHit;
    int i;
    float ndetector=0.,dn=0.;
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()>2.5){
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)>0){
              dx = x-trec->rec_x;
              dy = y-trec->rec_y;
              dz = z;
              dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
              if(dr>rcut1&&dr<rcut2){
                   dn += tHit->GetPe();
                   ndetector++;
              }
            }
        }
   }
   if(ndetector>2.5)  return dn/ndetector;
   else return 0;
}

float G4KM2A_Reconstruction::getdensityE(TClonesArray &tHits, int np, int rcut1, int rcut2, KM2ARecEvent *trec, const char *style)
{
    if( np < 4)
    {
        return 0;
    }
    LHHit *tHit;
    int nall;
    if( strcmp(style, "ED") == 0)
    {
        nall = geom->GetNED();
    }
    else if (strcmp(style, "MD") == 0) 
    {
        nall = geom->GetNMD();
    }
    double dirl, dirm, dirn;
    double x, y, z;
    double dx, dy, dz, dr;
    double dn = 0;
    double ndetector = 0;
    dirl = sin(trec->rec_theta) * cos( trec->rec_phi);
    dirm = sin(trec->rec_theta) * sin(trec->rec_phi);
    dirn = cos(trec->rec_theta);

    for( int i = 0; i < np; i++)
    {
        tHit = (LHHit *) ((tHits)[i]);
        if(tHit->GetStatus() > 2.5)
        {
            if(geom->Getxyz(tHit->GetId(), x, y, z, 1, style) > 0)
            {
                dx = x - trec->rec_x;
                dy = y - trec->rec_y;
                dz = z;
                dr = sqrt(dx*dx + dy*dy + dz*dz - pow(dx*dirl + dy*dirm + dz*dirn,2));
                if( dr > rcut1 && dr < rcut2)
                {
                    dn += tHit->GetPe();
                }

            }
        }
    } 
    for( int i = 0; i < nall; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if( geom->Getxyz(i, x, y, z, 0, style) < 0)
        {
            continue;
        }
        dx = x - trec->rec_x;
        dy = y - trec->rec_y;
        dz = z;
        dr = sqrt(dx*dx + dy*dy + dz*dz - pow(dx*dirl + dy*dirm + dz*dirn, 2));
        if( dr > rcut1 && dr < rcut2)
        {
            ndetector++; 
        }

    }
    if( ndetector > 2.5)
    {
        if(strcmp(style, "ED") == 0)
        {
            return dn/ndetector;
        }
        else if(strcmp(style, "MD") == 0)
        {
            return dn/ndetector/36.;
        }
        else 
        {
            return -1;
        }

    }
    else 
    {
        return 0.;
    }
}
//get the muon number of MD detector
float G4KM2A_Reconstruction::getmuM(TClonesArray &tHits,int np,int rcut,KM2ARecEvent *trec){
    LHHit *tHit;
    int i;
    float dn = 0;
    double dirl,dirm,dirn,dx,dy,dz,dr,x,y,z;
    dirl=sin(trec->rec_theta) * cos(trec->rec_phi);
    dirm=sin(trec->rec_theta) * sin(trec->rec_phi);
    dirn=cos(trec->rec_theta); 
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus() < 2.5 || tHit->GetPe() < 1.e-3)
            continue;
        if(geom->GetMDxyz(tHit->GetId(), x, y, z,1)<0)
            continue;
        dx = x - trec->rec_x;
        dy = y - trec->rec_y;
        dz = z;
        dr  = sqrt(dx*dx + dy*dy + dz*dz - pow(dx*dirl + dy*dirm - dz*dirn,2.));
        if(dr > rcut) 
            dn += tHit->GetPe();
   }
   return dn;
}


//get the muon number of WCDA detector
float G4KM2A_Reconstruction::getmuW(TClonesArray &tHits,int np,int rcut,KM2ARecEvent *trec){
    LHHit *tHit;
    int i;
    double dirl,dirm,dirn,dx,dy,dr,x,y,z,dn=0;

    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);

    dirl=sin(trec->rec_theta)*cos(trec->rec_phi);
    dirm=sin(trec->rec_theta)*sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(i=0; i<np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<2.5)
            continue;
        if(tHit->GetPe() > 16)
        {
            if(geom->GetWCDAxyz(tHit->GetId(),x,y,z,1)<0)continue;
            dx = x-trec->rec_x;
            dy = y-trec->rec_y;
            
            dr  = sqrt(dx*dx+dy*dy+z*z-pow(dx*dirl+dy*dirm-z*dirn,2.));
            if(dr>rcut) dn += tHit->GetPe()/45.;//from MC fit
        }  
   }
   return dn;
}
int G4KM2A_Reconstruction::conicalfit(TClonesArray &tHits,int np,KM2ARecEvent *trec,float alpha,const char *style){
    if(np<4)return -1;
    int i,j,k,nus,AA,Flag=0;
    double dt,w,the,phi,dirl,dirm,dirn,dr,c[4];
    double x,y,z,dx,dy,dz,sigma;
    if(alpha > 0)
        AA=1;
    else 
        AA=0; //to fit the alpha
    TMatrixD a(4,4);
    TVectorD b(4),r(4);
    LHHit *tHit;
    dirl=sin(trec->rec_theta) * cos(trec->rec_phi);
    dirm=sin(trec->rec_theta) * sin(trec->rec_phi);
    dirn=cos(trec->rec_theta);
    for(int iter = 0; iter < 30; iter++)
    {
        for(i = 0; i < 4; i++)
        {
            for(j = 0; j < 4; j++)
                a[i][j]=0;
            b[i]=0;
        }
        nus=0; sigma=0;
        for(i = 0; i < np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()<2.5 || tHit->GetPe()<1.e-3)
                continue; //not use the noise hit filter by planar fit
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style) < 0)
                continue;
            dx = x-trec->rec_x; 
            dy = y-trec->rec_y;
            dz = z;
            dr  = (dx*dx + dy*dy + dz*dz - pow(dx*dirl + dy*dirm - dz*dirn,2.));
            if(dr < 0)
                dr = 0;
            dr = sqrt(dr);

            if(iter == 0)
                w=1.;
            else
            {
              //  dt = tHit->GetTime()-(r[0]*x+r[1]*y-dirn*z+(r[2]+alpha*AA)*dr+r[3])/0.2998; //ns
              //  w  = exp(-0.5*pow(dt/10.,2));
                dt = tHit->GetTime()- (r[0]*x + r[1]*y - dirn*z + r[3])/0.2998; //ns please use this dt 
                if(dt > 0 || dt < -20) 
                    w  = exp(-0.5*pow(dt/10.,2)); //sig=10ns is better than 5,15,20 
                else 
                    w=1;
                dt = tHit->GetTime() - (r[0]*x + r[1]*y - dirn*z + (r[2] + alpha * AA) * dr + r[3])/0.2998;
                sigma += dt*dt;
            } 
            nus++;
            //w=sqrt(w); has a little worse for this 
            //w*=sqrt(tHit->GetNp()); not good
            c[0] = x * w;
            c[1] = y * w;
            c[2] = dr * (1 - AA) *w;
            c[3] = 1 * w;
            for(k = 0; k < 4; k++)
            {
                for(j = 0; j < 4; j++) 
                    a[j][k] += c[j]*c[k];
                b[k] += ((tHit->GetTime()) * 0.2998 - AA * alpha * dr + dirn*z) * w * c[k];
            }
            //if(DEBUG)printf("%d %lf %lf %lf %lf,%lf\n",i,c[0],c[1],tHit->GetTime(),dr,dirn*z);
        }
        //solve the equation
        if(nus < 4)
            break;
        TDecompSVD svd(a);
        Bool_t ok;
        r = svd.Solve(b,ok);
        printf(" ConicalFit : %d %d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
        if(!ok) break;
        if((r[0]*r[0] + r[1]*r[1]) > 1)
        {
             //printf(" ConicalFit >1: iter=%d Nhit=%d %lf %lf %lf\n",iter,nus,r[0],r[1],r[3]);
             break;
        }
        if(ok)
        {    
            Flag  = 1;
            sigma = sqrt(sigma/(nus-3.));
            dirl  = r[0];
            dirm  = r[1];
            dirn  = sqrt(1-(r[0] * r[0] + r[1] * r[1])); 
        }
    }
    //get the direction
    phi = -10.;  the = -10.;
    dt = r[0]*r[0] + r[1]*r[1];
    if(Flag > 0.5 && dt <= 1. && dt >= 0.)
    {
        the = asin(sqrt(dt));
        phi = atan2(r[1],r[0]);
        if(phi < 0.) 
            phi += PI*2.;
    }
    //set the result
    if(strcmp(style, "ED") == 0)
    {
        trec->rec_Etheta_c = the;
        trec->rec_Ephi_c   = phi;
        trec->rec_Ea       = r[2]*(1-AA)+AA*alpha;
        trec->rec_Ec0_c    = r[3];
        trec->rec_Esigma_c = sigma;
    }
    else if(strcmp(style,"WCDA") == 0){
        trec->rec_Wtheta_c = the;
        trec->rec_Wphi_c   = phi;
        trec->rec_Wa       = r[2]*(1-AA)+AA*alpha;
        trec->rec_Wc0_c    = r[3]; 
        trec->rec_Wsigma_c = sigma;
    }
    if(the>=0) return 0;
    else return -1;
}
//reconstruct the core of ED or WCDA array using centroid method 
void G4KM2A_Reconstruction::core_centre(TClonesArray &tHits,int np,KM2ARecEvent *trec,const char *style){        
    int i,j,nus=0;
    double x,y,z,w,mx=0,my=0,mz=0,total=0;
    LHHit *tHit;
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(i=0; i<np; i++){
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
        nus++;
        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style)<0)continue;
            // w = pow(tHit->GetNp(),2); //resolution is better than tHit->GetNp()
            //w = pow(tHit->GetPe(),2);
            w = tHit->GetPe(); 
            if(sqrt(pow(x,2)+pow(y,2))>CORE_R) w *=4.; 
            mx    += x*w;
            my    += y*w;
            mz    += z*w;
            total += w;    
    }
    //set result 
    if(strcmp(style, "ED") == 0)
    {   
        trec->rec_Ex = mx/total;
        trec->rec_Ey = my/total; 
        trec->rec_Ez = mz/total;
    }
    else if(strcmp(style, "WCDA") == 0)
    {
        trec->rec_Wx = mx/total;
        trec->rec_Wy = my/total;
        trec->rec_Wz = mz/total;
    }
}
//reconstruct the core of ED  array using optimized centroid method 
int G4KM2A_Reconstruction::core_centre2(TClonesArray &tHits,int np,int nr,KM2ARecEvent *trec,const char *style){
    if(np < 4)
        return -1;
    int i,j;
    double x,y,z,w,hx,hy,hz,dr,hx2,hy2;
    double Ax, Ay, Az, Aw;
    LHHit *tHit;
    for(int iter = 0; iter < 51; iter++)
    {
        Ax = 0; Ay = 0; Aw = 0;Az = 0; 
        for(i=0; i<np; i++)
        {
            tHit = (LHHit *)((tHits)[i]);
            if(tHit->GetStatus()<2.5 || tHit->GetPe()<1.e-3)
                continue;
            if(geom->Getxyz(tHit->GetId(),x,y,z,1,style) < 0)
                continue;
            if(iter == 0)
                w=1;
            else
            {
                dr = sqrt(pow(x-hx,2)+pow(y-hy,2));
                w  = exp(-0.5*pow(dr/15.,2)); //why use 15m?? can we optimized it???
            }
            //w *=(tHit->GetNp()); //Np is better than Np**2
            w *= (tHit->GetPe()); 
            //enhance the guard ring weight according the ED density 
            if(sqrt(pow(x,2)+pow(y,2))>CORE_R) 
                w *= 4.;
            Ax += x * w;
            Ay += y * w;
            Az += z * w;
            Aw += w;
        }
        hx = Ax/Aw;
        hy = Ay/Aw;
        hz = Az/Aw;
        dr = sqrt(pow(hx-hx2,2)+pow(hy-hy2,2));
        if( dr < 0.01 && iter > 10)
            break;
        hx2 = hx;  
        hy2 = hy;
    }
    //set the result 
    if(strcmp(style, "ED") == 0)
    {
        trec->rec_Ex = hx;
        trec->rec_Ey = hy;
        trec->rec_Ez = hz;
    }
    else if(strcmp(style,"WCDA") == 0)
    {
        trec->rec_Wx = hx;
        trec->rec_Wy = hy;
        trec->rec_Wz = hz;
    }
    return 0;
}
/*
    Par [0] : corex (m)
        [1] : corey (m)
        [2] : Size  (particle number)
        [3] : age s ()

*/

//NKG function1 for likelihood1 taken into account the detector with no hit 
void G4KM2A_Reconstruction::NKGfunction1(int &npar,double *gin,double &f,double *par,int iflag){
    int i,ned,nwc;
    //parameters use FengYouliang 2020-01-20 
   double x, y, z,r,u,dx,dy,dz,cs,A,rm=130.,sum=0; //the molliere radius is fixed to be 113m
    cs = TMath::Gamma(4.5-par[3])/(TMath::Gamma(par[3] - 0.5)*TMath::Gamma(5.0-2.*par[3]));
    cs = cs/(2 * PI * rm * rm);
    ned = geom->GetNED();
    nwc = geom->GetNWCDA();
    if(strcmp(_style,"ED")==0)
    {
        A = 1.0*_dirn;
        for(i = 0;i < ned; i++)
        {
            if(nED[i] > -0.5)
            {
                geom->Getxyz(i,x,y,z,0,_style);
                dx = x-par[0];
                dy = y-par[1];
                dz = z;
                r  = sqrt(dx*dx+ dy*dy + dz*dz - pow(dx*_dirl + dy*_dirm - dz*_dirn,2.));
                if( r < 0.3)
                    r=0.3;
                u  = A*cs * par[2]*pow(r/rm,par[3]-2)* pow(1+r/rm,par[3]-4.5);  
                if( u > 0)
                    sum += (nED[i])*log(u)-u; 
            }
        }
    }
    else if(strcmp(_style,"WCDA")==0){
        A=25.0*_dirn;
        for(i=0;i<nwc;i++){
            if(nWCDA[i]>-0.5){
                geom->Getxyz(i,x,y,z,0,_style);
                dx = x-par[0];
                dy = y-par[1];
                dz = z;
                r  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*_dirl+dy*_dirm-dz*_dirn,2.));
                if(r<0.1)r=0.1;
                u  = A*cs*par[2]*pow(r/rm,par[3]-2)*pow(1+r/rm,par[3]-4.5);
                if(u>0)sum += nWCDA[i]*log(u)-u;
            }
        }
    }
    f = -sum;
}
//NKG function2 for likelihood2 only taken into account the detector with  hit 
void G4KM2A_Reconstruction::NKGfunction2(int &npar,double *gin,double &f,double *par,int iflag){
    int i; 
    double x,y,z,r,u,dx,dy,dz,cs,A,rm = 130.,sum = 0;
    if(strcmp(_style,"ED")==0)
        A = 1.0*_dirn;
    else if(strcmp(_style,"WCDA")==0)
        A = 25.0*_dirn;
    cs = TMath::Gamma(4.5 - par[3])/(TMath::Gamma(par[3] - 0.5)*TMath::Gamma(5.0 - 2.*par[3]));
    cs = cs/(2*PI*rm*rm);
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
        
    TClonesArray* tHits=(TClonesArray*)minuit->GetObjectFit();
    LHHit *tHit;
    for(i=0; i<_np; i++){
        tHit = (LHHit *)((*tHits)[i]);
        if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)
            continue;
        geom->Getxyz(tHit->GetId(),x,y,z,1,_style);
        dx = x-par[0];
        dy = y-par[1];
        dz = z;
        r  = sqrt(dx*dx + dy*dy + dz*dz - pow(dx*_dirl + dy*_dirm - dz*_dirn,2.));
        //r  = sqrt(dx*dx+dy*dy-pow(dx*_dirl+dy*_dirm,2.));
        if(r < 0.2)
            r = 0.2;
        u  = A * cs*par[2] * pow(r/rm,par[3]-2) * pow(1+r/rm,par[3]-4.5);
        if(u>0)
            sum += (tHit->GetPe())*log(u)-u;
    }   
    f = -sum; 
}



/*
    Function to Reconstruct the core position use the NKG Fitting
*/
// likelihood using NKG function1 taken into account the detector with and without hit
// using NKG function2 only taken into account the detector with  hit
int G4KM2A_Reconstruction::core_likelihood(TClonesArray &tHits,int np,KM2ARecEvent *trec,const char *style,int method){
    if( np < 5)
        return -1;
    int  i,flag,ned,nwc;
    double size,x,y,sigma,esize,ex,ey,esigma,arglist[10];
    x = -1.e4;
    y = -1.e4;
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    ned = geom->GetNED();
    nwc = geom->GetNWCDA();
    if(np>4)
    {
        minuit->mncler();
        //init the value will be used for NKGfunction
        _np = np;
        size = 0;
        strcpy(_style,style);
        _dirl = sin(trec->rec_theta)*cos(trec->rec_phi);
        _dirm = sin(trec->rec_theta)*sin(trec->rec_phi);
        _dirn = cos(trec->rec_theta);
        //init all the detector hit value 
        if(strcmp(_style,"ED")==0){
            for(i = 0;i < ned; i++)
            { 
                if(nED[i] > -0.5)
                    nED[i] = 0;
            }
            LHHit *tHit;
            for(i=0; i<_np; i++)
            {
                tHit = (LHHit *)((tHits)[i]);
                if(tHit->GetStatus()<2.5 || tHit->GetPe()<1.e-3) 
                    continue;
                if(nED[geom->GetEDId2(tHit->GetId())] > -0.5)
                    nED[geom->GetEDId2(tHit->GetId())] += tHit->GetPe();
                size += tHit->GetPe();
            }
        }
        else if(strcmp(_style,"WCDA")==0){
            for(i=0;i<nwc;i++)    { if(nWCDA[i]>-0.5)nWCDA[i]=0;}
            LHHit *tHit;
            for(i=0; i<_np; i++){
                tHit = (LHHit *)((tHits)[i]);
                if(tHit->GetStatus()<2.5||tHit->GetPe()<1.e-3)continue;
                if(nWCDA[geom->GetWCDAId2(tHit->GetId())]>-0.5)nWCDA[geom->GetWCDAId2(tHit->GetId())] +=tHit->GetPe();
                size += tHit->GetPe();
            }
        }
        if(size < 4)
            return -1;
        arglist[0]=-1;
        minuit->mnexcm("SET PRINT",arglist,1,flag);//1:standard; 2: try to improve minuim (slower)
        minuit->mnparm(0,"corex",(double)trec->rec_x,20.,-1000.,1000.,flag);
        minuit->mnparm(1,"corey",(double)trec->rec_y,20.,-1000.,1000.,flag);
        minuit->mnparm(2,"size", size*160, 100,size,size*1.e7,flag);
        minuit->mnparm(3,"age",1.2,0.2,0.5,2.5,flag); //for NKG
        //flag=0 if no problems,  >0 if MNPARM unable to implement definition
        //minuit->FixParameter(2);
        minuit->FixParameter(3);
        if(method == 0)
            minuit->SetFCN(NKGfunction1);
        else 
            minuit->SetFCN(NKGfunction2);
        minuit->SetObjectFit(&tHits); 
        minuit->SetPrintLevel(-1); // -1  quiet (also suppresse all warnings),0  normal,1  verbose

        arglist[0] = 500; //loop
        arglist[1] = 0.1; //tolerance
        flag = 0;
        minuit->mnexcm("SIMPLEX",arglist,0,flag);
        minuit->FixParameter(2);
        minuit->FixParameter(3);
        arglist[0] = 1000; //loop
        arglist[1] = 0.01; //tolerance
        minuit->mnexcm("MIGRAD",arglist,0,flag);
        minuit->GetParameter(0,x,ex);
        minuit->GetParameter(1,y,ey);

        minuit->FixParameter(0);
        minuit->FixParameter(1);
        minuit->Release(2);
        minuit->Release(3);
        minuit->mnexcm("MIGRAD",arglist,0,flag);
        minuit->GetParameter(0,x,ex);
        minuit->GetParameter(1,y,ey);
        minuit->GetParameter(2,size,esize);
        minuit->GetParameter(3,sigma,esigma);
    }
    if(strcmp(style,"ED")==0)
    {
        trec->rec_Esize =float(size);
        trec->rec_Ex = float(x);
        trec->rec_Ey = float(y);
        trec->rec_Eage = float(sigma);
    }
    else if(strcmp(style,"WCDA")==0)
    {
        trec->rec_Wsize =float(size);
        trec->rec_Wx = float(x);
        trec->rec_Wy = float(y);
        trec->rec_Wage = float(sigma);
    }
    return 0;
}



int G4KM2A_Reconstruction::noisefilter(TClonesArray &tHits,int np,int twind,int twinu,int rwin,KM2ARecEvent *trec,const char *style){
    int i,nus=0;
    double dt,dirl,dirm,dirn,dr,tt;
    double x,y,z,dx,dy,dz;
    LHHit *tHit;
    if((strcmp(style, "WCDA")) == 0)
        tt = trec->dt;
    else 
        tt=0;
    //the direction from ED array is used here 
    dirl=sin(trec->rec_Etheta_c)*cos(trec->rec_Ephi_c);
    dirm=sin(trec->rec_Etheta_c)*sin(trec->rec_Ephi_c);
    dirn=cos(trec->rec_Etheta_c);
    //G4KM2A_Geometry *geom = G4KM2A_Geometry::GetInstance(ARRAY);
    for(i = 0; i < np; i++)
    {
        tHit = (LHHit *)((tHits)[i]);
        if(tHit->GetStatus()<0 || tHit->GetPe() < 1.e-3)
            continue; //not use the bad detector
        if(geom->Getxyz(tHit->GetId(),x,y,z,1,style) < 0)
            continue;
        dx = x-trec->rec_x;
        dy = y-trec->rec_y;
        dz = z;
        dr  = sqrt(dx*dx+dy*dy+dz*dz-pow(dx*dirl+dy*dirm-dz*dirn,2.));
        dt = tHit->GetTime()-tt-(dirl*x+dirm*y-dirn*z+trec->rec_Ea*dr+trec->rec_Ec0_c)/0.2998; //ns
        if(dt > twind && dt < twinu)
        {
            if(dr < rwin)
            {
                tHit->SetStatus(5);
                nus ++;
            }
            else 
                tHit->SetStatus(2); 
        } 
        else
        {
            if(dr < rwin)
                tHit->SetStatus(1);
            else 
                tHit->SetStatus(0);
        }
   }
   return nus;
}
void G4KM2A_Reconstruction::Draw(LHEvent *tevent, int num_id, KM2ARecEvent* trec)
{   
    eventrecline(tevent, trec);

    TCanvas* km2a_image = new TCanvas(Form("event id %d km2a image", num_id), "km2a_Image", 1600, 800);

    int ned = geom->GetNED();
    int nmd = geom->GetNMD();
    TGraph* g_ed = new TGraph();
    TGraph* g_md = new TGraph();

    TGraph* core_pos = new TGraph();
    core_pos->SetPoint(0, -trec->rec_y, trec->rec_x);
    core_pos->SetMarkerStyle(29);
    core_pos->SetMarkerColor(2);
    core_pos->SetMarkerSize(2);


    double x, y, z;

    for( int i = 0; i < ned; i ++)
    {
        geom->Getxyz(i, x, y, z, 1, "ED");
        z = x;
        x = -y;
        y = z;
        g_ed->SetPoint(i, x, y);
    }

    for( int i = 0; i < nmd; i++)
    {
        geom->Getxyz(i, x, y, z, 1, "MD");
        z = x;
        x = -y;
        y = z;
        g_md->SetPoint(i, x, y);
    }
    g_ed->SetMarkerStyle(97);
    g_ed->SetMarkerColor(18);
    g_ed->GetXaxis()->SetRangeUser(-750,650);
    g_ed->GetYaxis()->SetRangeUser(-750,650);
    g_ed->GetXaxis()->SetTitle("X (m)");
    g_ed->GetXaxis()->CenterTitle();
    g_ed->GetXaxis()->SetLabelSize(0.05);
    g_ed->GetXaxis()->SetTitleSize(0.05);
    g_ed->GetYaxis()->SetTitle("Y (m)");
    g_ed->GetYaxis()->CenterTitle();
    g_ed->GetYaxis()->SetLabelSize(0.05);
    g_ed->GetYaxis()->SetTitleSize(0.05);
    g_ed->SetMarkerSize(1);

    g_md->SetMarkerStyle(4);
    g_md->SetMarkerColor(16);
    g_md->SetMarkerSize(1);
    
    km2a_image->Divide(2,1);
    km2a_image->cd(1);
    gPad->SetMargin(0.2,0.2,0.2,0.2);
    gPad->SetTickx();
    gPad->SetTicky();
    g_ed->Draw("ap");
    g_md->Draw("psame");
    km2a_image->cd(2);
    gPad->SetMargin(0.2,0.2,0.2,0.2);
    gPad->SetTickx();
    gPad->SetTicky();
    g_ed->Draw("ap");
    g_md->Draw("psame");

    TH2Poly* ed_h2p = new TH2Poly();
    TH2Poly* md_h2p = new TH2Poly();

    km2a_image->cd(1);

    auto ed = tevent->GetHitsE();
    for( int i = 0; i < tevent->GetNhitE(); i++)
    {
        LHHit* tmp_ed = (LHHit* ) (*ed)[i];
        int id = tmp_ed->GetId();
        double x, y, z;
        double content = log10(tmp_ed->GetPe());
        if( content < 0)
        {
            continue;
        }
        geom->Getxyz(id, x, y, z, 1, "ED");
        z = x;
        x = -y;
        y = z;
        double x_ed[4] = {0};
        double y_ed[4] = {0};

        x_ed[0] = x - 5, y_ed[0] = y - 5;
        x_ed[1] = x + 5, y_ed[1] = y - 5;
        x_ed[2] = x + 5, y_ed[2] = y + 5;
        x_ed[3] = x - 5, y_ed[3] = y + 5;

        ed_h2p->AddBin(4, x_ed, y_ed);
        ed_h2p->Fill(x, y , content);
    }
    ed_h2p->Draw("colz same");
    TLatex *edtitle = new TLatex(560,600,"log_{10}(Ne)");
    edtitle->SetTextSize(0.04);
    edtitle->Draw("same");
    core_pos->Draw("psame");

    km2a_image->cd(2);
    auto md = tevent->GetHitsM();
    for( int i = 0; i < tevent->GetNhitM(); i++)
    {
        LHHit* tmp_md = (LHHit*) (*md)[i];
        int id = tmp_md->GetId();
        double content = log10(tmp_md->GetPe());
        double x, y, z;
        if(content <0 || tmp_md->GetStatus() <2.5)
        {
            continue;
        }

        geom->Getxyz(id, x, y, z, 1, "ED");
        z = x;
        x = -y;
        y = x;
        double x_md[4] = {0};
        double y_md[4] = {0};
        x_md[0] = x - 5, y_md[0] = y - 5;
        x_md[1] = x + 5, y_md[1] = y - 5;
        x_md[2] = x + 5, y_md[2] = y + 5;
        x_md[3] = x - 5, y_md[3] = y + 5;
        
        md_h2p->AddBin(4, x_md, y_md);
        md_h2p->Fill(x, y , content);
    }
    md_h2p->Draw("colz same");
    TLatex *mdtitle = new TLatex(500, 600, "log_{10}(N#mu)");
    mdtitle->SetTextSize(0.04);
    mdtitle->Draw("same");
    core_pos->Draw("psame");

    km2a_image->cd();
    TPaveText *pave_text =  new TPaveText(0.15,0.88,0.85,0.95);
    pave_text->SetFillStyle(0);
    pave_text->AddText(Form("E: %.2lfTeV, rec_theta %lf rec_phi %lf  direction_error : %lf", trec->E/1000, trec->rec_theta, trec->rec_phi, trec->Direction_Error));
    pave_text->Draw("same");

    km2a_image->SaveAs(Form("Event%d_km2a.png", num_id));
}