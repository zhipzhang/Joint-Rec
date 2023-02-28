#include <cmath>
#include <stdio.h>
#define PI 3.14159265358979312
double receNfit(double nhit){
    if(nhit>5)return 0.446739*log10(nhit)*log10(nhit)-0.159119*log10(nhit)+0.44573 ; 
    else return 0;
}
int CutKM2AErecNfit(double erec,double ratio){ //erec=log10(E/TeV)//ratio=log10(Nu+0.1/Ne)
    int i,n,m;
    double ll;
    //zenith<50 all >60% 
    //pgamma 0.81, 0.77, 0.73, 0.65, 0.64, 0.74, 0.76, 0.70, 0.87, 0.89, 
    //pback 1.30e-02, 6.90e-03, 1.74e-03, 4.48e-04, 2.33e-04, 1.56e-04, 1.30e-04, 6.41e-05, 1.42e-04, 3.23e-05, 
    //Q 7.16, 9.31, 17.40, 30.53, 41.96, 59.06, 66.98, 87.77, 72.76, 155.76, 
    double cut[10]={-2.20, -2.50, -2.60, -2.80, -2.60, -2.50, -2.50, -2.60, -2.40, -2.40};  
    //zenith<50 all >60% exclude sqrt((corex-20)**2+(corey+510)**2)<120
    //pgamma 0.81, 0.77, 0.72, 0.64, 0.69, 0.73, 0.76, 0.70, 0.87, 0.88, 
    //pback 1.19e-02, 6.35e-03, 1.61e-03, 4.28e-04, 2.60e-04, 1.49e-04, 1.29e-04, 4.95e-05, 1.37e-04, 3.38e-05, 
    //Q=7.45, 9.66, 18.01, 31.04, 42.85, 60.06, 66.76, 99.31, 73.81, 151.93, 
    //double cut[10]={-2.20, -2.50, -2.60, -2.80, -2.50, -2.50, -2.50, -2.60, -2.40, -2.40}; 
    if(erec<1)return -1;
     n=int((erec-1)/0.2);
     if(n>9)n=9;
     if(ratio<cut[n])return n;
     else return -1;
}
double rece2(double ne,double theta=-1){ //NpE2+theta+core rec energy,return energy log10(E) (TeV)
    double ene=0,ezen=0,er=0;
    if(ne>5&&theta>-0.1){ //NpE1/NpE2>2 dr<40m
        return 0.94789566*log10(ne)+0.60047041/pow(cos(theta),2)-0.91006520/cos(theta)+0.32395902;
    }
    else if(ne>5&&theta<0){ //all event
        return 0.955211*log10(ne)+0.13467;
    }
    else return -10;
};
int CutKM2AErecNe2(double erec,double ratio){ //erec=log10(E/TeV)//ratio=log10(Nu+0.1/Ne)
    int i,n,m;
    double ll;
    double cut[10]={-2.30, -2.70, -2.90, -3.00, -2.40, -2.40, -2.40, -2.30, -2.50, -2.30,};     
     if(erec<1)return -1;
     n=int((erec-1)/0.2);
     if(n>9)n=9;
     if(ratio<cut[n])return n;
     else return -1;
}
int  corefilter(double x, double y, int flag=0){
    double dr,phi;
    if(flag==0){
       dr=sqrt((x+355)*(x+355)+(y-195)*(y-195)); 
       if(dr<160)return 1;
       else return 0;
    }
    else if(flag==1) { 
        dr=sqrt(x*x+(y+20)*(y+20));
        phi=atan2(x,-(y+20))*180/PI; 
        if(phi<0)phi += 360;    
        if(dr>250&&dr<575&&phi>208&&phi<273)return 1;
        else return 0;   
    }
    else{
       printf("corefilter() flag error\n");
       return 0;
   }
};
double GetSpaceLC(double zen1, double phi1, double zen2, double phi2)
{
    if(zen1<0||zen2<0)return -1;

    double dr,dl0,dm0,dz0,dl1,dm1,dz1;

       dl0= sin(zen1)*cos(phi1);
       dm0= sin(zen1)*sin(phi1);
       dz0= cos(zen1);

       dl1 = sin(zen2)*cos(phi2);
       dm1 = sin(zen2)*sin(phi2);
       dz1 = cos(zen2);
     dr = acos(dl0*dl1+dm0*dm1+dz0*dz1)*57.295780;
    return dr;
}
