// -----------------------------------------------------------------------
//                                G4KM2A
//                       --- G4KM2A_Geometry ---
//                        class implementation file
// -----------------------------------------------------------------------
//  Description:
//
//    Managers the LHAASO detector position
//
// -----------------------------------------------------------------------
// History:
//
//    chen songzhan      -- 2017.3 created
//           
//========================================================================
//

#include "G4KM2A_Geometry.h"

#include <iostream>
#include <fstream>
//#include <map>
#include <algorithm>   
#include <iterator>
#include "math.h"
#include <string.h>

using namespace std;
int G4KM2A_Geometry::Flag=2;

G4KM2A_Geometry * G4KM2A_Geometry::m_myself = 0;
G4KM2A_Geometry::G4KM2A_Geometry(int tFlag)
{
    Flag=tFlag;
    fNED=0;
    fNMD=0;
    fNWCDA=0;
    fMinx=1.e10;  fMiny=1.e10;  fMinz=1.e10; 
    fMaxx=-1.e10; fMaxy=-1.e10; fMaxz=-1.e10;

    //read WCDA geometry file
    //ReadFileWCDA("./config/WCDA_pos_900.txt");

    if(Flag==0){//calibration mode
        ReadFileED("./config/ED_pos_unit.txt");
        ReadFileMD("./config/MD_pos_unit.txt");
    }
    else if(Flag==1){ //KM2A@ybj 
        ReadFileED("./config/ED_pos_ybj.txt");
        ReadFileMD("./config/MD_pos_ybj.txt");
    }
    else if(Flag==2){ //KM2A-33ED 
        ReadFileED("./config/ED_pos_33.txt");
        ReadFileMD("./config/MD_pos_0.txt");
    }
    else if(Flag==3){ //KM2A-71ED+10MD 
        ReadFileED("./config/ED_pos_71.txt");
        ReadFileMD("./config/MD_pos_10.txt");
    }
    else if(Flag==4){ //KM2A-quarter
        ReadFileED("./config/ED_pos_quarter.txt");
        ReadFileMD("./config/MD_pos_quarter.txt");
    }
    else if(Flag==5){ //KM2A-half 
        ReadFileED("./config/ED_pos_half.txt");
        ReadFileMD("./config/MD_pos_half.txt");
        //ReadFileMD("./config/MD_pos_halfnew.txt"); for MC 
    }
    else if(Flag==6){ //KM2A_all
        ReadFileED("./config/ED_pos_all.txt");
        ReadFileMD("./config/MD_pos_all.txt");
    }
    else {
        cout<<" No G4KM2A_Geometry !!!"<<endl;
        exit(22);
    }

    //adjust the range for detector range
    double dz = fMaxz-fMinz;
    fMinx -=dz;
    fMiny -=dz;
    fMaxx +=dz;
    fMaxy +=dz;
   if(fNED>0) SetArrayRange();
}
G4KM2A_Geometry::~G4KM2A_Geometry()
{

}
void G4KM2A_Geometry::SetArrayRange()
{
   int i,j,k,a,b,c;
   double dx,dy,dr=-30/57.295779;
   for(i=0;i<140;i++){
       fXX[i][0]=1000;
       fXX[i][1]=-1000;
       fYY[i][0]=1000;
       fYY[i][1]=-1000;
   }
   //get bin x range y, bin y range x
   for(i=0;i<fNED;i++){
         //exclude guard ring 
         if(sqrt(fEDx[i]*fEDx[i]+fEDy[i]*fEDy[i])>575)continue;

         dx=fEDx[i]*cos(dr)-fEDy[i]*sin(dr);
         dy=fEDx[i]*sin(dr)+fEDy[i]*cos(dr);
         j=int((dx+700)/10.);
         if(j>=0&&j<140){
            k=int(dy);
            if(k<fXX[j][0])fXX[j][0]=k;
            if(k>fXX[j][1])fXX[j][1]=k; 
         }
         j=int((dy+700)/10.);
         if(j>=0&&j<140){
            k=int(dx);
            if(k<fYY[j][0])fYY[j][0]=k;
            if(k>fYY[j][1])fYY[j][1]=k;
         }
   }
   //smooth y min
   a=fXX[0][0];
   for(i=0;i<140;i++){
       b=fXX[i][0];
       if(i<140)c=fXX[i+1][0];
       else c=fXX[i][0];
       if(fXX[i][0]>a)fXX[i][0]=a;
       if(fXX[i][0]>c)fXX[i][0]=c;
       a=b;
   }
   //smooth x min
   a=fYY[0][0];
   for(i=0;i<140;i++){
       b=fYY[i][0];
       if(i<140)c=fYY[i+1][0];
       else c=fYY[i][0];
       if(fYY[i][0]>a)fYY[i][0]=a;
       if(fYY[i][0]>c)fYY[i][0]=c;
       a=b;
   }
   //smooth y max
   a=fXX[0][1];
   for(i=0;i<140;i++){
       b=fXX[i][1];
       if(i<140)c=fXX[i+1][1];
       else c=fXX[i][1];
       if(fXX[i][1]<a)fXX[i][1]=a;
       if(fXX[i][1]<c)fXX[i][1]=c;
       a=b;
   }     
   //smooth x max
   a=fYY[0][1];
   for(i=0;i<140;i++){
       b=fYY[i][1];
       if(i<140)c=fYY[i+1][1];
       else c=fYY[i][1];
       if(fYY[i][1]<a)fYY[i][1]=a;
       if(fYY[i][1]<c)fYY[i][1]=c;
       a=b;
   }
   //for(i=0;i<140;i++)printf("x: %d %d %d\n ",i*10-700+5,fXX[i][0],fXX[i][1]);
   //for(i=0;i<140;i++)printf("y: %d %d %d\n ",i*10-700+5,fYY[i][0],fYY[i][1]);
}
int G4KM2A_Geometry::CheckCore(double dx,double dy, int nd)
{
    int i,j,k,nx,ny;
    double x,y,dr=-30/57.295779;
    x=dx*cos(dr)-dy*sin(dr);
    y=dx*sin(dr)+dy*cos(dr);
    //exclude ceter square region
    if(fabs(x-10)<140+nd&&fabs(y)<150+nd)return 0;

    //exclude the edge region
    nx=int((x+700.)/10.+2)-2;
    ny=int((y+700.)/10.+2)-2;
    if(nx<0||nx>139||ny<0||ny>139)return 0;
    if(y<fXX[nx][0]+nd||y>fXX[nx][1]-nd) return 0;
    if(x<fYY[ny][0]+nd||x>fYY[ny][1]-nd) return 0;
    return 1;
}
double G4KM2A_Geometry::GetCoreEdge(double dx,double dy)
{
    int i,j,k,nx,ny;
    double min,x,y,dr=-30/57.295779;
    x=dx*cos(dr)-dy*sin(dr);
    y=dx*sin(dr)+dy*cos(dr);
    //calculate the center region
    min=(fabs(x-10)-140);
    if(fabs(y+5)-165>min){
         min=fabs(y+5)-165;
    }

    nx=int((x+700.)/10.+2)-2;
    ny=int((y+700.)/10.+2)-2;
    if(nx<0||nx>139||ny<0||ny>139)min=-100;
    if(fXX[nx][0]<990){
      if(min>y-fXX[nx][0])min=y-fXX[nx][0];
      if(min>fXX[nx][1]-y)min=fXX[nx][1]-y;
    }
    if(fYY[ny][0]<990){
       if(min>x-fYY[ny][0])min=x-fYY[ny][0];
       if(min>fYY[ny][1]-x)min=fYY[ny][1]-x;
    }
    if(min<-101||min>300)min=-101; 
    return min;
}
void G4KM2A_Geometry::ReadFileED(const char filename[100] ) 
{
    ifstream file1;
    char sTemp[200];
    int i,id,n=0;
    double x,y,z,zeroZ;
    fEDIdmax=0;
    //file1.open("./config/ED_pos_ybj.txt");
   file1.open(filename);
    cout<< "Open ED file "<< filename<<endl;
    fRotationED=0;
    if(file1.good()){
        file1.getline(sTemp,200);
        i=sscanf(sTemp,"%*s %lf %*s %*s %lf",&x,&zeroZ);
        fRotationED=x;
    }
    while(file1.good()){
        file1.getline(sTemp,200);
        //i=sscanf(sTemp,"%d %lf %lf %lf",&id,&x,&y,&z);
        i=sscanf(sTemp,"%d %lf %lf %lf",&id,&y,&x,&z);
        if(i>1){
            y=-y; //change the detector coord to corsika coord
            z =z-zeroZ;
           // printf("%d %lf %lf %lf\n",id,y,x,z);
            if(id>fEDIdmax)fEDIdmax=id;
            fEDId.push_back(id);
            fEDx.push_back(x); 
            fEDy.push_back(y);
            fEDz.push_back(z);
            if(x-2<fMinx)fMinx=x-2;
            if(y-2<fMiny)fMiny=y-2;
            if(z-1<fMinz)fMinz=z-1;

            if(x+2>fMaxx)fMaxx=x+2;
            if(y+2>fMaxy)fMaxy=y+2;
            if(z+1>fMaxz)fMaxz=z+1;  
            n++;
        }
    }
    file1.close();
    cout<<"   "<<n<<" ED unit position"<<endl;
    fNED=n;

    for(i=0;i<fEDIdmax+1;i++){
         fEDId2.push_back(-1);
    }
    for(i=0;i<fNED;i++){
         id=fEDId[i];
         fEDId2[id]=i;
    }
} 

void G4KM2A_Geometry::ReadFileMD(const char filename[100])
{
    ifstream file1;
    char sTemp[200];
    int i,id,n=0;
    double x,y,z,zeroZ;
    fMDIdmax=0;
    file1.open(filename);
    cout<< "Open MD file "<<filename<<endl;
    if(file1.good()){
        file1.getline(sTemp,200);
        i=sscanf(sTemp,"%*s %lf",&zeroZ);
    }
    while(file1.good()){
        file1.getline(sTemp,200);
        i=sscanf(sTemp,"%d %lf %lf %lf",&id,&y,&x,&z);
        if(i>1){
            y=-y; //change the detector coord to corsika coord
            z =z-zeroZ;
            if(Flag>3) z=z-2.0; //now shift -5 for KM2A all
            //printf("%d %lf %lf %lf\n",id,y,x,z);
            if(id>fMDIdmax)fMDIdmax=id;
            fMDId.push_back(id);
            fMDx.push_back(x);
            fMDy.push_back(y);
            fMDz.push_back(z);
            if(x-8<fMinx)fMinx=x-8;
            if(y-8<fMiny)fMiny=y-8;
            if(z-2<fMinz)fMinz=z-2;

            if(x+8>fMaxx)fMaxx=x+8;
            if(y+8>fMaxy)fMaxy=y+8;
            if(z+2>fMaxz)fMaxz=z+2;
            n++;
        }
    }
    file1.close();
    cout<<"   "<<n<<" MD unit position"<<endl;
    fNMD=n;

    for(i=0;i<fMDIdmax+1;i++){
         fMDId2.push_back(-1);
    }
    for(i=0;i<fNMD;i++){
         id=fMDId[i];
         fMDId2[id]=i;
    }
}
void G4KM2A_Geometry::ReadFileWCDA(const char filename[100])
{
    ifstream file1;
    char sTemp[200];
    int i,id,n=0;
    double x,y,z,zeroZ;
    fWCDAIdmax=0;
    file1.open(filename);
    cout<< "Open WCDA file "<<filename<<endl;
    if(file1.good()){
        file1.getline(sTemp,200);
        i=sscanf(sTemp,"%*s %lf %*s %*s %lf",&x,&zeroZ);
        //fRotationWCDA=x;
    }
    while(file1.good()){
        file1.getline(sTemp,200);
        i=sscanf(sTemp,"%d %lf %lf %lf",&id,&y,&x,&z);
          y=-y; //change the detector coord to corsika coord
          z =z-zeroZ;
        if(i>1){
            if(id>fWCDAIdmax)fWCDAIdmax=id;
            fWCDAId.push_back(id);
            fWCDAx.push_back(x);
            fWCDAy.push_back(y);
            fWCDAz.push_back(z);
            if(x-8<fMinx)fMinx=x-8;
            if(y-8<fMiny)fMiny=y-8;
            if(z-2<fMinz)fMinz=z-2;

            if(x+8>fMaxx)fMaxx=x+8;
            if(y+8>fMaxy)fMaxy=y+8;
            if(z+2>fMaxz)fMaxz=z+2;
            n++;
        }
    }
    file1.close();
    cout<<"   "<<n<<" MD unit position"<<endl;
    fNWCDA=n;

    for(i=0;i<fWCDAIdmax+1;i++){
         fWCDAId2.push_back(-1);
    }
    for(i=0;i<fNWCDA;i++){
         id=fWCDAId[i];
         fWCDAId2[id]=i;
    }
}

//from real id to the id of array
int G4KM2A_Geometry::GetEDId2(int n){
      if(n<fEDIdmax+1)
          return fEDId2[n];
      else {
         // cout<<" G4KM2A_Geometry::GetEDId2 Error!!!"<<endl;
          return -1;
      }
}

//from real id to the id of array
int G4KM2A_Geometry::GetMDId2(int n){
      if(n<fMDIdmax+1)
          return fMDId2[n];
      else {
        //  cout<<" G4KM2A_Geometry::GetMDId2 Error!!!"<<endl;
          return -1;
      }
}
int G4KM2A_Geometry::GetWCDAId2(int n){
      if(n<fWCDAIdmax+1)
          return fWCDAId2[n];
      else {
          cout<<" G4KM2A_Geometry::GetWCDAId2 Error!!!"<<endl;
          return -1;
      }
}
int G4KM2A_Geometry::GetEDxyz(int id, double &x,double &y, double &z,int flag)
{
    int ii;
    if(flag==0) ii=id;
    else{
        //from real id to the id of array
        ii = GetEDId2(id);
    }
    if(ii>=fNED||ii<0){ 
      //     cout<<"Error: ED id is not right! flag="<<flag<<", ED="<<id<<", array="<<ii<<endl;
         return -1;
    }
    x=fEDx[ii];
    y=fEDy[ii];
    z=fEDz[ii];
    return 1;
}
int G4KM2A_Geometry::GetMDxyz(int id, double &x,double &y, double &z, int flag)
{
    int ii;
    if(flag==0) ii=id;
    else{
        //from real id to the id of array
        ii = GetMDId2(id);
    }
    if(ii>=fNMD||ii<0){
        // cout<<"Error: MD id is not right! MD="<<id<<", array="<<ii<<endl;
         return -1;
    }
    x=fMDx[ii];
    y=fMDy[ii];
    z=fMDz[ii];
    return 1;
}

int G4KM2A_Geometry::GetWCDAxyz(int id, double &x,double &y, double &z, int flag)
{
    int ii;
    if(flag==0) ii=id;
    else{
        ii = GetWCDAId2(id);
    }
    if(ii>=fNWCDA||ii<0){ 
        cout<<"Error: WCDA id is not right! "<<endl; 
        return -1;
    }
    x=fWCDAx[ii];
    y=fWCDAy[ii];
    z=fWCDAz[ii];
    return 1;
}

int G4KM2A_Geometry::Getxyz(int id, double &x,double &y, double &z, int flag,const char style[20])
{
     if(strcmp(style,"ED")==0){ return GetEDxyz(id,x,y,z,flag); }
     else if(strcmp(style,"MD")==0){ return GetMDxyz(id,x,y,z,flag); }
     else if(strcmp(style,"WCDA")==0){ return GetWCDAxyz(id,x,y,z,flag); }
}

double G4KM2A_Geometry::Get2EDr(int id1,int id2)
{
    if(id1>=fNED||id2>=fNED){ cout<<"Error: ED id is not right! "<<endl; exit(11);}
    return 1;
}
double G4KM2A_Geometry::Get2MDr(int id1,int id2)
{
    if(id1>=fNMD||id2>=fNMD){ cout<<"Error: MD id is not right! "<<endl; exit(11);}
    return 1;
}
void G4KM2A_Geometry::GetRangeX(double &minx,double &maxx)
{
    minx=fMinx;
    maxx=fMaxx;  
}
void G4KM2A_Geometry::GetRangeY(double &miny,double &maxy)
{
    miny=fMiny;
    maxy=fMaxy;
}
void G4KM2A_Geometry::GetRangeZ(double &minz,double &maxz)
{
    minz=fMinz;
    maxz=fMaxz;
}

double G4KM2A_Geometry::GetRangeX() 
{
    if(fabs(fMaxx)>fabs(fMinx))return fabs(fMaxx);
    else return fabs(fMinx);
}
double G4KM2A_Geometry::GetRangeY()                             
{
    if(fabs(fMaxy)>fabs(fMiny))return fabs(fMaxy);
    else return fabs(fMiny);
}
double G4KM2A_Geometry::GetRangeZ()                             
{
    if(fabs(fMaxz)>fabs(fMinz))return fabs(fMaxz);
    else return fabs(fMinz);
}
