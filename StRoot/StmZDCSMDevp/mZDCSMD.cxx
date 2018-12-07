#define mZDCSMD_cxx
#include "mZDCSMD.h"
#include <TH2.h>
#include "TH1.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#ifdef __ROOT__
ClassImp(mZDCSMD);
#endif

using namespace std;

mZDCSMD::mZDCSMD()
{
	//InitHist();
	mFileDirectory="./";
}

mZDCSMD::~mZDCSMD()
{

}

////////////////////////////////////////
void mZDCSMD::SetZDCSMDrefm(Float_t *nZDC,Float_t refmult,int runN)
{
	// input EastHor 0-7; EastVer, 8-15; WestHor 16-23; WestVer, 24-31;
	// East 0, West 1; horizontal 1, vertical 0
	for (int strip=1;strip<9;strip++) {              
		mZDCSMDadcR[0][1][strip-1] = nZDC[strip-1];      //  EastHor 
		mZDCSMDadcR[0][0][strip-1] = nZDC[strip-1+8];    //  EastVer 
		mZDCSMDadcR[1][1][strip-1] = nZDC[strip-1+16];   //  WestHor
		mZDCSMDadcR[1][0][strip-1] = nZDC[strip-1+24];   //  WestVer
		mZDCSMDadcC[0][1][strip-1] = nZDC[strip-1];      //  EastHor 
		mZDCSMDadcC[0][0][strip-1] = nZDC[strip-1+8];    //  EastVer 
		mZDCSMDadcC[1][1][strip-1] = nZDC[strip-1+16];   //  WestHor
		mZDCSMDadcC[1][0][strip-1] = nZDC[strip-1+24];   //  WestVer

		//cout<<mZDCSMD[0][0][strip-1] <<endl; 
	}

	mrefmult=refmult;
	mcent=FindCentrality(mrefmult);
	mrun=runN;
	runPointer(mrun);
}	

////////////////////////////////////////
void mZDCSMD::SetZDCSMDcent(Float_t *nZDC,int centrality,int runN)
{
	// input EastHor 0-7; EastVer, 8-15; WestHor 16-23; WestVer, 24-31;
	// East 0, West 1; horizontal 1, vertical 0
	for (int strip=1;strip<9;strip++) {              
		mZDCSMDadcR[0][1][strip-1] = nZDC[strip-1];      //  EastHor 
		mZDCSMDadcR[0][0][strip-1] = nZDC[strip-1+8];    //  EastVer 
		mZDCSMDadcR[1][1][strip-1] = nZDC[strip-1+16];   //  WestHor
		mZDCSMDadcR[1][0][strip-1] = nZDC[strip-1+24];   //  WestVer
		mZDCSMDadcC[0][1][strip-1] = nZDC[strip-1];      //  EastHor 
		mZDCSMDadcC[0][0][strip-1] = nZDC[strip-1+8];    //  EastVer 
		mZDCSMDadcC[1][1][strip-1] = nZDC[strip-1+16];   //  WestHor
		mZDCSMDadcC[1][0][strip-1] = nZDC[strip-1+24];   //  WestVer

		//cout<<mZDCSMD[0][0][strip-1] <<endl; 
	}

	mcent=centrality;
	// mrun=runN;
	mrun=17109018;
	runPointer(mrun);
}	

//---------------------------------------------
void mZDCSMD::SetMod(int mod){
	mMod=mod;

	cout<<endl;
	if(mMod==1)cout<<"###### step1: for get ZDCsmd pedstal and gain correction parameters ######"<<endl;
	if(mMod==2)cout<<"###### step2: for get ZDCsmd beam center correction parameters ######"<<endl;
	if(mMod==3)cout<<"###### step3: for get ZDCsmd psi shift (subevent) correction parameters ######"<<endl;
	if(mMod==4)cout<<"###### step4: for get ZDCsmd psi shift (fullevent) correction parameters ######"<<endl;
	if(mMod==5)cout<<"###### ZDCsmd ready (with QA Hists) ######"<<endl;
	if(mMod==6)cout<<"###### ZDCsmd ready (without Hists) ######"<<endl;

	if(mMod<1||mMod>6)cout<<"###### ZDCsmd mod wrong  ######"<<endl;
	cout<<endl;
}

//---------------------------------------------
bool mZDCSMD::InitRun(int mod=1){
	SetMod(mod);

	if(mMod<1||mMod>6)return false;

	InitHist();

	// raw
	mZDCSMD_rawQFull.Set(0.,0.);
	mZDCSMD_rawQEast.Set(0.,0.);
	mZDCSMD_rawQWest.Set(0.,0.);

	// gain corrected 
	mZDCSMD_corQFull.Set(0.,0.);
	mZDCSMD_corQEast.Set(0.,0.);
	mZDCSMD_corQWest.Set(0.,0.);

	// recenter 
	mZDCSMD_QFull.Set(0.,0.);
	mZDCSMD_QEast.Set(0.,0.);
	mZDCSMD_QWest.Set(0.,0.);

	// shifted
	mZDCSMD_QFullS.Set(0.,0.);
	mZDCSMD_QFulSS.Set(0.,0.);
	mZDCSMD_QEastS.Set(0.,0.);
	mZDCSMD_QWestS.Set(0.,0.);

	mZDCSMD_PsiEast=0.;
	mZDCSMD_PsiWest=0.;
	mZDCSMD_PsiFull=0.;
	mZDCSMD_PsiEastS=0.;
	mZDCSMD_PsiWestS=0.;
	mZDCSMD_PsiFullS=0.;
	mZDCSMD_PsiFulSS=0.;

	for (int strip=1;strip<9;strip++) {              
		//SetZDCSMDpeds(0,1,strip,mZDCSMDped0[0][1][strip-1]);  
		//SetZDCSMDpeds(0,0,strip,mZDCSMDped0[0][0][strip-1]);  
		//SetZDCSMDpeds(1,1,strip,mZDCSMDped0[1][1][strip-1]);  
		//SetZDCSMDpeds(1,0,strip,mZDCSMDped0[1][0][strip-1]);  
		SetZDCSMDpeds(0,1,strip,0.);  
		SetZDCSMDpeds(0,0,strip,0.);  
		SetZDCSMDpeds(1,1,strip,0.);  
		SetZDCSMDpeds(1,0,strip,0.);  
	}

	for(Int_t icent=0;icent<10;icent++){
		for(Int_t iDay=0; iDay<nDays; iDay++){
			for(Int_t ith=0;ith<nShis;ith++) {
				mZDCSMD_shiftEast_cos[icent][iDay][ith] = 0.;
				mZDCSMD_shiftEast_sin[icent][iDay][ith] = 0.;
				mZDCSMD_shiftWest_cos[icent][iDay][ith] = 0.;
				mZDCSMD_shiftWest_sin[icent][iDay][ith] = 0.;
				mZDCSMD_shiftFull_cos[icent][iDay][ith] = 0.;
				mZDCSMD_shiftFull_sin[icent][iDay][ith] = 0.;
				mZDCSMD_shiftFulS_cos[icent][iDay][ith] = 0.;
				mZDCSMD_shiftFulS_sin[icent][iDay][ith] = 0.;
			}
		}//j=0~nDays
	}//cent
	for(int iRun=0;iRun<totalRunNumber;iRun++){
			mZDCSMDCenterEx[iRun] = zdcsmd_ex0;
			mZDCSMDCenterEy[iRun] = zdcsmd_ey0;
			mZDCSMDCenterWx[iRun] = zdcsmd_wx0;
			mZDCSMDCenterWy[iRun] = zdcsmd_wy0;

	}//for
	// fill raw adc for ped gain	
	if(mod==1)return true;	

	//  
	TString mFileName;
	if(mod==2){mFileName=mFileDirectory+"/zdcsmd_ped.root";    openFile(mFileName);} // with pedstal  
	if(mod==3){mFileName=mFileDirectory+"/zdcsmd_center.root"; openFile(mFileName);} // with pedstal + beamcenter  
	if(mod==4){mFileName=mFileDirectory+"/zdcsmd_evpsub.root"; openFile(mFileName);} // with pedstal + beamcenter + shift factor (subevent) 
	if(mod>=5){mFileName=mFileDirectory+"/zdcsmd_evpfull.root";openFile(mFileName);} // with pedstal + beamcenter + shift factor (fullevent) 
	return true;	

}

//---------------------------------------------
bool mZDCSMD::InitRun(){
	return InitRun(1);
}

//--------------------------------------------------------------
Int_t mZDCSMD::openFile(TString name) {
	// Read the ZDCSMD constants root file
	TFile* pZDCSMDFile = new TFile(name, "READ");
	if (!pZDCSMDFile->IsOpen())cout<<"#####  No ZDCSMD related file. Will use default."<<endl;


	if(mMod>1)ReadZDCSMDFile(pZDCSMDFile); 
	if(mMod>3)ReadZDCSMDshiftFile(pZDCSMDFile); 


	// Close ZDCSMD constants file
	if (pZDCSMDFile->IsOpen())pZDCSMDFile->Close("R");
	return 1;	
}

////////////////////////////////////////
void mZDCSMD::InitEvent(){
	for (int strip=1;strip<9;strip++) {
		SetZDCSMDadcR(0,1,strip,0.);
		SetZDCSMDadcR(0,0,strip,0.);
		SetZDCSMDadcR(1,1,strip,0.);
		SetZDCSMDadcR(1,0,strip,0.);
		SetZDCSMDadcC(0,1,strip,0.);
		SetZDCSMDadcC(0,0,strip,0.);
		SetZDCSMDadcC(1,1,strip,0.);
		SetZDCSMDadcC(1,0,strip,0.);
	}
	mcent =0;
	mrun  =0;
	dayPtr=0;
	runPtr=0;
}

//---------------------------------------------
bool mZDCSMD::ZDCSMDbadrun()
{
	bool mBadZDC=false;
	for(int i=0;i<nZDCRuns;i++){
		//if(badRunzdc[i]==runPtr)mBadZDC=true;
		if(badRunNzdc[i]==mrun)mBadZDC=true;
		//if(badRunzdc[i]==runPtr)cout<<mrun<<endl;
		//if(badRunNzdc[i]==mrun)cout<<mrun<<endl;
	}

	return mBadZDC;
}

//-----------------------------------------------------------------------
Float_t mZDCSMD::ZDCSMD_GetPosition(int eastwest,int verthori,int strip) {
	//get position of each slat;strip starts from 1

	Float_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
	Float_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

	Float_t  mZDCSMDCenterex = mZDCSMDCenterEx[runPtr];//-0.07;
	Float_t  mZDCSMDCenterey = mZDCSMDCenterEy[runPtr];//-0.01;
	Float_t  mZDCSMDCenterwx = mZDCSMDCenterWx[runPtr];
	Float_t  mZDCSMDCenterwy = mZDCSMDCenterWy[runPtr];

	if(mMod>2){  // with beam center corrected
		if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1]-mZDCSMDCenterex;
		if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x[strip-1];
		if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterey;
		if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterwy;
	}
	else{        // raw beam center returned.
		if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1];
		if(eastwest==1 && verthori==0) return -zdcsmd_x[strip-1];
		if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.);
		if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.);
	}	

	return 0;
}

//----------------------------------------------------------------------
void mZDCSMD::InitHist(){

	gRandom=new TRandom3();
	gRandom->SetSeed(0);

	int nXs=(mMod==6)?1:1000;
	int nYs=(mMod==6)?1:nDays;

	char title[256];
	char yname[256];
	for (int i = 0; i < 8; i++) {
		sprintf(title,"ZDCSMDR_east_ver%d",i); mHistZDCSMDEastVR[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDR_east_hor%d",i); mHistZDCSMDEastHR[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDR_west_ver%d",i); mHistZDCSMDWestVR[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDR_west_hor%d",i); mHistZDCSMDWestHR[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDC_east_ver%d",i); mHistZDCSMDEastVC[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDC_east_hor%d",i); mHistZDCSMDEastHC[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDC_west_ver%d",i); mHistZDCSMDWestVC[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
		sprintf(title,"ZDCSMDC_west_hor%d",i); mHistZDCSMDWestHC[i] = new TH2F(title,title,nXs,-5.,995.,nYs, 0, Float_t(nDays));
	}

	for(int i=0;i<10;i++){
		//raw psi
		sprintf(title,"Flow_ZDCSMDPsi_R_East_cen%d" ,i); mHistZDCSMDPsi_R_East[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
		sprintf(title,"Flow_ZDCSMDPsi_R_West_cen%d" ,i); mHistZDCSMDPsi_R_West[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
		sprintf(title,"Flow_ZDCSMDPsi_R_Full_cen%d" ,i); mHistZDCSMDPsi_R_Full[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));

		//for shift file
		sprintf(title,"Flow_ZDCSMDshiftEast_c_cen%d",i); mHistZDCSMDshiftEast_c[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftEast_s_cen%d",i); mHistZDCSMDshiftEast_s[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftWest_c_cen%d",i); mHistZDCSMDshiftWest_c[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftWest_s_cen%d",i); mHistZDCSMDshiftWest_s[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftFulS_c_cen%d",i); mHistZDCSMDshiftFulS_c[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftFulS_s_cen%d",i); mHistZDCSMDshiftFulS_s[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftFull_c_cen%d",i); mHistZDCSMDshiftFull_c[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		sprintf(title,"Flow_ZDCSMDshiftFull_s_cen%d",i); mHistZDCSMDshiftFull_s[i] = new TProfile2D(title,title,nShis,0.5,nShis+0.5,nYs, 0, Float_t(nDays),-1,1,"");
		//ZDCSMD psi after shift
		sprintf(title,"Flow_ZDCSMDPsi_S_FulS_cen%d" ,i); mHistZDCSMDPsi_S_FulS[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
		sprintf(title,"Flow_ZDCSMDPsi_S_East_cen%d" ,i); mHistZDCSMDPsi_S_East[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
		sprintf(title,"Flow_ZDCSMDPsi_S_West_cen%d" ,i); mHistZDCSMDPsi_S_West[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
		sprintf(title,"Flow_ZDCSMDPsi_S_Full_cen%d" ,i); mHistZDCSMDPsi_S_Full[i]  = new TH2F(title,title,nPsiBins,0,twopi,nYs, 0, Float_t(nDays));
	}

	mHistCosZ = new TProfile("mHistCosZ","mHistCosZ",10,0,10, -1., 1.);  
	mHistSinZ = new TProfile("mHistSinZ","mHistSinZ",10,0,10, -1., 1.);  
	mZDCSMDBeamCenter = new TProfile2D("ZDCSMDBeamCenter","ZDCSMDBeamCenter",4,0.5,4.5,totalRunNumber,0,totalRunNumber,-20,20,"");

}


void mZDCSMD::WriteHist(){

	if(mMod>5)return;

	// Create mHist output file
	if(mHistOutputFileName == "") {
		cout<< "Please input output hist file" <<endl;
	}
	else {
		mhist_output = new TFile(mHistOutputFileName, "recreate");
		// Name was set previously in calling macro 
	}

	mHistCosZ ->Write();
	mHistSinZ ->Write();
	mZDCSMDBeamCenter->Write();  

	for(int i=0;i<8;i++) {
		mHistZDCSMDEastVR[i] ->Write();  mHistZDCSMDEastVC[i] ->Write();
		mHistZDCSMDEastHR[i] ->Write();  mHistZDCSMDEastHC[i] ->Write();
		mHistZDCSMDWestVR[i] ->Write();  mHistZDCSMDWestVC[i] ->Write();
		mHistZDCSMDWestHR[i] ->Write();  mHistZDCSMDWestHC[i] ->Write();
	}
	for(int i=0;i<10;i++){
		mHistZDCSMDPsi_R_East[i]->Write();
		mHistZDCSMDPsi_R_West[i]->Write();
		mHistZDCSMDPsi_R_Full[i]->Write();

		mHistZDCSMDshiftEast_c[i]->Write();
		mHistZDCSMDshiftEast_s[i]->Write();
		mHistZDCSMDshiftWest_c[i]->Write();
		mHistZDCSMDshiftWest_s[i]->Write();
		mHistZDCSMDshiftFulS_c[i]->Write();
		mHistZDCSMDshiftFulS_s[i]->Write();
		mHistZDCSMDshiftFull_c[i]->Write();
		mHistZDCSMDshiftFull_s[i]->Write();

		mHistZDCSMDPsi_S_FulS[i]->Write();
		mHistZDCSMDPsi_S_East[i] ->Write();
		mHistZDCSMDPsi_S_West[i] ->Write();
		mHistZDCSMDPsi_S_Full[i] ->Write();
	}

	 if(mhist_output!=NULL) mhist_output -> Write() ;
}

/////////////////
void mZDCSMD::Fill(){
	//mhist_output->cd(); 
	//mhist_output->Close(); 
	
	if(mMod>5)return;
	
	for(int i=0;i<8;i++) {
		mHistZDCSMDEastVR[i]->Fill(GetZDCSMDadcR(0,0,i+1),dayPtr);  mHistZDCSMDEastVC[i]->Fill(GetZDCSMDadcC(0,0,i+1),dayPtr);
		mHistZDCSMDEastHR[i]->Fill(GetZDCSMDadcR(0,1,i+1),dayPtr);  mHistZDCSMDEastHC[i]->Fill(GetZDCSMDadcC(0,1,i+1),dayPtr);
		mHistZDCSMDWestVR[i]->Fill(GetZDCSMDadcR(1,0,i+1),dayPtr);  mHistZDCSMDWestVC[i]->Fill(GetZDCSMDadcC(1,0,i+1),dayPtr);
		mHistZDCSMDWestHR[i]->Fill(GetZDCSMDadcR(1,1,i+1),dayPtr);  mHistZDCSMDWestHC[i]->Fill(GetZDCSMDadcC(1,1,i+1),dayPtr);
	}

	Float_t  mPsi_ZDCSMD[3];     //! without shift
	Float_t  mPsi_S_ZDCSMD[5];   //! after shift

	mZDCSMDBeamCenter->Fill(1,runPtr,GetZDCSMD_QEast().X());  
	mZDCSMDBeamCenter->Fill(2,runPtr,GetZDCSMD_QEast().Y()); 
	mZDCSMDBeamCenter->Fill(3,runPtr,GetZDCSMD_QWest().X()); 
	mZDCSMDBeamCenter->Fill(4,runPtr,GetZDCSMD_QWest().Y()); 

	//raw
	mPsi_ZDCSMD[0]   = GetZDCSMD_PsiFull();//[0,2Pi]
	mPsi_ZDCSMD[1]   = GetZDCSMD_PsiEast();//[0,2Pi]
	mPsi_ZDCSMD[2]   = GetZDCSMD_PsiWest();//[0,2Pi]

	//shift
	mPsi_S_ZDCSMD[0]   = GetZDCSMD_PsiFulSS();//[0,2Pi]
	mPsi_S_ZDCSMD[1]   = GetZDCSMD_PsiEastS();//[0,2Pi]
	mPsi_S_ZDCSMD[2]   = GetZDCSMD_PsiWestS();//[0,2Pi]
	mPsi_S_ZDCSMD[3]   = GetZDCSMD_PsiFullS();//[0,2Pi]

	//hist
	mHistZDCSMDPsi_R_East[mcent]->Fill(mPsi_ZDCSMD[1],dayPtr);
	mHistZDCSMDPsi_R_West[mcent]->Fill(mPsi_ZDCSMD[2],dayPtr);
	mHistZDCSMDPsi_R_Full[mcent]->Fill(mPsi_ZDCSMD[0],dayPtr);

	//correlation between them
	Float_t order  = 1.0;
	mHistCosZ->Fill(mcent, (Float_t)cos(order *(mPsi_S_ZDCSMD[1]-mPsi_S_ZDCSMD[2]+twopi/2.)));
	mHistSinZ->Fill(mcent, (Float_t)sin(order *(mPsi_S_ZDCSMD[1]-mPsi_S_ZDCSMD[2]+twopi/2.)));
	mHistCosZ->Fill(   0., (Float_t)cos(order *(mPsi_S_ZDCSMD[1]-mPsi_S_ZDCSMD[2]+twopi/2.)));
	mHistSinZ->Fill(   0., (Float_t)sin(order *(mPsi_S_ZDCSMD[1]-mPsi_S_ZDCSMD[2]+twopi/2.)));

	//ZDCSMD psi after shift
	mHistZDCSMDPsi_S_FulS[mcent]->Fill(mPsi_S_ZDCSMD[0],dayPtr);
	mHistZDCSMDPsi_S_East[mcent]->Fill(mPsi_S_ZDCSMD[1],dayPtr);
	mHistZDCSMDPsi_S_West[mcent]->Fill(mPsi_S_ZDCSMD[2],dayPtr);
	mHistZDCSMDPsi_S_Full[mcent]->Fill(mPsi_S_ZDCSMD[3],dayPtr);

	//shift An & Bn
	for(int i=1;i<nShis+1;i++) {
		mHistZDCSMDshiftEast_c[mcent]->Fill(i,dayPtr,cos((Float_t)i*mPsi_ZDCSMD[1]));
		mHistZDCSMDshiftEast_s[mcent]->Fill(i,dayPtr,sin((Float_t)i*mPsi_ZDCSMD[1]));
		mHistZDCSMDshiftWest_c[mcent]->Fill(i,dayPtr,cos((Float_t)i*mPsi_ZDCSMD[2]));
		mHistZDCSMDshiftWest_s[mcent]->Fill(i,dayPtr,sin((Float_t)i*mPsi_ZDCSMD[2]));
		mHistZDCSMDshiftFull_c[mcent]->Fill(i,dayPtr,cos((Float_t)i*mPsi_ZDCSMD[0]));
		mHistZDCSMDshiftFull_s[mcent]->Fill(i,dayPtr,sin((Float_t)i*mPsi_ZDCSMD[0]));
		mHistZDCSMDshiftFulS_c[mcent]->Fill(i,dayPtr,cos((Float_t)i*mPsi_S_ZDCSMD[3]));
		mHistZDCSMDshiftFulS_s[mcent]->Fill(i,dayPtr,sin((Float_t)i*mPsi_S_ZDCSMD[3]));
	}

	return;
}


void mZDCSMD::calibrateZDCSMDevp(){

	if(mMod>1)pedGain();

	//Q vector  from east ZDCSMD 
	Float_t eXsumR=0., eYsumR=0., eXWgtR=0., eYWgtR=0.;
	Float_t eXsumC=0., eYsumC=0., eXWgtC=0., eYWgtC=0.;
	for(Int_t strip = 1; strip < 8; strip++) {
		Float_t mADCpos=ZDCSMD_GetPosition(0,0,strip);
		eXsumR += mADCpos*GetZDCSMDadcR(0,0,strip);
		eXsumC += mADCpos*GetZDCSMDadcC(0,0,strip);
		eXWgtR += GetZDCSMDadcR(0,0,strip);
		eXWgtC += GetZDCSMDadcC(0,0,strip);
	}

	for(Int_t strip = 1; strip < 9; strip++) {
		Float_t mADCpos=ZDCSMD_GetPosition(0,1,strip);
		eYsumR += mADCpos*GetZDCSMDadcR(0,1,strip);
		eYsumC += mADCpos*GetZDCSMDadcC(0,1,strip);
		eYWgtR += GetZDCSMDadcR(0,1,strip);
		eYWgtC += GetZDCSMDadcC(0,1,strip);
	}

	//Q vector  from west ZDCSMD 
	Float_t wXsumR=0., wYsumR=0., wXWgtR=0., wYWgtR=0.;
	Float_t wXsumC=0., wYsumC=0., wXWgtC=0., wYWgtC=0.;
	for(Int_t strip = 1; strip < 8; strip++) {
		Float_t mADCpos=ZDCSMD_GetPosition(1,0,strip);
		wXsumR += mADCpos*GetZDCSMDadcR(1,0,strip);
		wXsumC += mADCpos*GetZDCSMDadcC(1,0,strip);
		wXWgtR += GetZDCSMDadcR(1,0,strip);
		wXWgtC += GetZDCSMDadcC(1,0,strip);
	}

	for(Int_t strip = 1; strip < 9; strip++) {
		Float_t mADCpos=ZDCSMD_GetPosition(1,1,strip);
		wYsumR += mADCpos*GetZDCSMDadcR(1,1,strip);
		wYsumC += mADCpos*GetZDCSMDadcC(1,1,strip);
		wYWgtR += GetZDCSMDadcR(1,1,strip);
		wYWgtC += GetZDCSMDadcC(1,1,strip);
	}

	mZDCSMD_rawQEast.Set((eXWgtR>0.) ? eXsumR/eXWgtR:0.,(eYWgtR>0.) ? eYsumR/eYWgtR:0.);
	mZDCSMD_rawQWest.Set((wXWgtR>0.) ? wXsumR/wXWgtR:0.,(wYWgtR>0.) ? wYsumR/wYWgtR:0.);
	mZDCSMD_corQEast.Set((eXWgtC>0.) ? eXsumC/eXWgtC:0.,(eYWgtC>0.) ? eYsumC/eYWgtC:0.);
	mZDCSMD_corQWest.Set((wXWgtC>0.) ? wXsumC/wXWgtC:0.,(wYWgtC>0.) ? wYsumC/wYWgtC:0.);
	mZDCSMD_QEast.Set(   (eXWgtC>0.) ? eXsumC/eXWgtC:0.,(eYWgtC>0.) ? eYsumC/eYWgtC:0.);
	mZDCSMD_QWest.Set(   (wXWgtC>0.) ? wXsumC/wXWgtC:0.,(wYWgtC>0.) ? wYsumC/wYWgtC:0.);

	mZDCSMD_rawQFull = mZDCSMD_rawQEast-mZDCSMD_rawQWest;
	mZDCSMD_corQFull = mZDCSMD_corQEast-mZDCSMD_corQWest;
	mZDCSMD_QFull    = mZDCSMD_QEast-mZDCSMD_QWest;

	mZDCSMD_PsiEast =mZDCSMD_QEast.Phi(); if(mZDCSMD_PsiEast <0.)mZDCSMD_PsiEast +=twopi;
	mZDCSMD_PsiWest =mZDCSMD_QWest.Phi(); if(mZDCSMD_PsiWest <0.)mZDCSMD_PsiWest +=twopi;
	mZDCSMD_PsiFull =mZDCSMD_QFull.Phi(); if(mZDCSMD_PsiFull <0.)mZDCSMD_PsiFull +=twopi;

	mZDCSMD_PsiEastS=mZDCSMD_QEast.Phi(); if(mZDCSMD_PsiEastS<0.)mZDCSMD_PsiEastS+=twopi;
	mZDCSMD_PsiWestS=mZDCSMD_QWest.Phi(); if(mZDCSMD_PsiWestS<0.)mZDCSMD_PsiWestS+=twopi;
	mZDCSMD_PsiFullS=mZDCSMD_QFull.Phi(); if(mZDCSMD_PsiFullS<0.)mZDCSMD_PsiFullS+=twopi;
	mZDCSMD_PsiFulSS=mZDCSMD_QFull.Phi(); if(mZDCSMD_PsiFulSS<0.)mZDCSMD_PsiFulSS+=twopi;


	////////////////////
	Float_t psi_e0 = mZDCSMD_PsiEast;
	Float_t psi_e  = psi_e0;
	for(Int_t i=1;i<nShis+1;i++) psi_e += 2*(-mZDCSMD_shiftEast_sin[mcent][dayPtr][i-1]*cos(i*psi_e0)+mZDCSMD_shiftEast_cos[mcent][dayPtr][i-1]*sin(i*psi_e0))/(Float_t)i;
	psi_e = atan2(sin(psi_e),cos(psi_e));
	if(psi_e<0)psi_e+=twopi;
	mZDCSMD_PsiEastS=psi_e;

	////////////////////
	Float_t psi_w0 = mZDCSMD_PsiWest;
	Float_t psi_w  = psi_w0;
	for(Int_t i=1;i<nShis+1;i++) psi_w += 2*(-mZDCSMD_shiftWest_sin[mcent][dayPtr][i-1]*cos(i*psi_w0)+mZDCSMD_shiftWest_cos[mcent][dayPtr][i-1]*sin(i*psi_w0))/(Float_t)i;
	psi_w = atan2(sin(psi_w),cos(psi_w));
	if(psi_w<0)psi_w+=twopi;
	mZDCSMD_PsiWestS=psi_w;

	// QfullS 
	Float_t mEQx=0., mEQy=0.;
	Float_t mWQx=0., mWQy=0.;
	Float_t mFQx=0., mFQy=0.;
	Float_t mQe  = mZDCSMD_QEast.Mod();
	Float_t mQw  = mZDCSMD_QWest.Mod();
	mEQx = mQe*cos(mZDCSMD_PsiEastS);
	mEQy = mQe*sin(mZDCSMD_PsiEastS);
	mWQx = mQw*cos(mZDCSMD_PsiWestS);
	mWQy = mQw*sin(mZDCSMD_PsiWestS);
	mFQx = mEQx - mWQx;
	mFQy = mEQy - mWQy;
	mZDCSMD_QEastS.Set(mEQx, mEQy);
	mZDCSMD_QWestS.Set(mWQx, mWQy);
	mZDCSMD_QFullS.Set(mFQx, mFQy);
	mZDCSMD_PsiFullS =mZDCSMD_QFullS.Phi();
	if(mZDCSMD_PsiFullS<0.)mZDCSMD_PsiFullS+=twopi;

	// shifted QfullSS 
	Float_t psi_f0 = mZDCSMD_PsiFullS;
	Float_t psi_f  = psi_f0;
	for(Int_t i=1;i<nShis+1;i++) psi_f += 2*(-mZDCSMD_shiftFulS_sin[mcent][dayPtr][i-1]*cos(i*psi_f0)+mZDCSMD_shiftFulS_cos[mcent][dayPtr][i-1]*sin(i*psi_f0))/(Float_t)i;
	psi_f = atan2(sin(psi_f),cos(psi_f));
	if(psi_f<0)psi_f+=twopi;
	mFQx = mZDCSMD_QFullS.Mod()*cos(psi_f);
	mFQy = mZDCSMD_QFullS.Mod()*sin(psi_f);
	mZDCSMD_QFulSS.Set(mFQx, mFQy);
	mZDCSMD_PsiFulSS=psi_f;


	// Fill histogram
	///////////////////////////
	if(ZDCSMDgoodEvent())Fill();


	return;
}

//-------------------------------------------------------------
bool mZDCSMD::ZDCSMDgoodEvent(){
	bool mGoodZDCSMD=true;

	//Float_t mZDCSMDadcSumEast=0.;
	//Float_t mZDCSMDadcSumWest=0.;
	//for(int iTile=0;iTile<16;iTile++){
	//	if(GetBBCadcR(0,iTile+1)>3995.  ) mGoodBBC=false; 
	//	if(GetBBCadcR(1,iTile+1)>3995.  ) mGoodBBC=false; 
	//	mBBCadcSumEast+=GetBBCadcR(0,iTile+1);
	//	mBBCadcSumWest+=GetBBCadcR(1,iTile+1);
	//}
	//if(mBBCadcSumEast<100. ) mGoodBBC=false;
	//if(mBBCadcSumWest<100. ) mGoodBBC=false;

	if(fabs(GetZDCSMD_rawQEast().X())<1e-6||fabs(GetZDCSMD_rawQEast().Y())<1e-6)mGoodZDCSMD=false;
	if(fabs(GetZDCSMD_rawQWest().X())<1e-6||fabs(GetZDCSMD_rawQWest().Y())<1e-6)mGoodZDCSMD=false;

	if(mMod>1 &&(fabs(GetZDCSMD_corQEast().X())<1e-6||fabs(GetZDCSMD_corQEast().Y())<1e-6))mGoodZDCSMD=false;
	if(mMod>1 &&(fabs(GetZDCSMD_corQWest().X())<1e-6||fabs(GetZDCSMD_corQWest().Y())<1e-6))mGoodZDCSMD=false;
	if(mMod>2 &&(fabs(GetZDCSMD_QEast().X())<1e-6||fabs(GetZDCSMD_QEast().Y())<1e-6))mGoodZDCSMD=false;
	if(mMod>2 &&(fabs(GetZDCSMD_QWest().X())<1e-6||fabs(GetZDCSMD_QWest().Y())<1e-6))mGoodZDCSMD=false;

	return mGoodZDCSMD;
}

//-------------------------------------------------------------
void mZDCSMD::pedGain(){
	//get ZDCSMD pedstal-subtracted and gain-corrected
	Float_t zdcsmdEastHorizontal = -1.;
	Float_t zdcsmdEastVertical   = -1.;
	Float_t zdcsmdWestHorizontal = -1.;
	Float_t zdcsmdWestVertical   = -1.;

	for (int strip=1;strip<9;strip++) {
		if (GetZDCSMDadcR(0,1,strip)) {
			zdcsmdEastHorizontal = (GetZDCSMDadcR(0,1,strip)-GetZDCSMDpeds(0,1,strip))/mZDCSMDGain[0][1][strip-1];
			SetZDCSMDadcC(0,1,strip,zdcsmdEastHorizontal);
		}
		if (GetZDCSMDadcR(0,0,strip)) {
			zdcsmdEastVertical   = (GetZDCSMDadcR(0,0,strip)-GetZDCSMDpeds(0,0,strip))/mZDCSMDGain[0][0][strip-1];
			SetZDCSMDadcC(0,0,strip,zdcsmdEastVertical);
		}
		if (GetZDCSMDadcR(1,1,strip)) {
			zdcsmdWestHorizontal = (GetZDCSMDadcR(1,1,strip)-GetZDCSMDpeds(1,1,strip))/mZDCSMDGain[1][1][strip-1];
			SetZDCSMDadcC(1,1,strip,zdcsmdWestHorizontal);
		}
		if (GetZDCSMDadcR(1,0,strip)) {
			zdcsmdWestVertical   = (GetZDCSMDadcR(1,0,strip)-GetZDCSMDpeds(1,0,strip))/mZDCSMDGain[1][0][strip-1];
			SetZDCSMDadcC(1,0,strip,zdcsmdWestVertical);
		}
	}
}

////--------------------------------------------------------------
Int_t mZDCSMD::ReadZDCSMDshiftFile(TFile* pPsiShiftFile){
	if (!pPsiShiftFile->IsOpen()) {
		cout<<"##### FlowMaker: No psi shift file. Will set shifts = 0."<<endl;
	}else cout<<"Read in ZDCSMD shift File"<<endl;

	char title[256];
	if (pPsiShiftFile->IsOpen()){
		for(int icent=0;icent<10;icent++){
			//ZDCSMD psi shift
			sprintf(title,"Flow_ZDCSMDshiftWest_c_cen%d",icent);	TProfile2D* HmZDCSMDshiftWest_c =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftWest_s_cen%d",icent);	TProfile2D* HmZDCSMDshiftWest_s =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftEast_c_cen%d",icent);	TProfile2D* HmZDCSMDshiftEast_c =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftEast_s_cen%d",icent);	TProfile2D* HmZDCSMDshiftEast_s =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftFull_c_cen%d",icent);	TProfile2D* HmZDCSMDshiftFull_c =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftFull_s_cen%d",icent);	TProfile2D* HmZDCSMDshiftFull_s =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftFulS_c_cen%d",icent);	TProfile2D* HmZDCSMDshiftFulS_c =(TProfile2D*) pPsiShiftFile->Get(title);
			sprintf(title,"Flow_ZDCSMDshiftFulS_s_cen%d",icent);	TProfile2D* HmZDCSMDshiftFulS_s =(TProfile2D*) pPsiShiftFile->Get(title);

			if(HmZDCSMDshiftEast_c && HmZDCSMDshiftWest_c) { //shift histograms exist
				if(icent==1){
					cout<<"get mZDCSMD_shiftEast_c && mZDCSMD_shiftWest_c"<<endl;
				}
				for(int iDay=0; iDay<nDays; iDay++){
					for(int ith=0;ith<nShis;ith++) {
						mZDCSMD_shiftEast_cos[icent][iDay][ith] = HmZDCSMDshiftEast_c->GetBinContent(ith+1,iDay+1);
						mZDCSMD_shiftEast_sin[icent][iDay][ith] = HmZDCSMDshiftEast_s->GetBinContent(ith+1,iDay+1);
						mZDCSMD_shiftWest_cos[icent][iDay][ith] = HmZDCSMDshiftWest_c->GetBinContent(ith+1,iDay+1);
						mZDCSMD_shiftWest_sin[icent][iDay][ith] = HmZDCSMDshiftWest_s->GetBinContent(ith+1,iDay+1);
						mZDCSMD_shiftFull_cos[icent][iDay][ith] = HmZDCSMDshiftFull_c->GetBinContent(ith+1,iDay+1);
						mZDCSMD_shiftFull_sin[icent][iDay][ith] = HmZDCSMDshiftFull_s->GetBinContent(ith+1,iDay+1);
					}
					if(HmZDCSMDshiftFull_c) {
						for(int ith=0;ith<nShis;ith++) {
							mZDCSMD_shiftFulS_cos[icent][iDay][ith] = HmZDCSMDshiftFulS_c->GetBinContent(ith+1,iDay+1);
							mZDCSMD_shiftFulS_sin[icent][iDay][ith] = HmZDCSMDshiftFulS_s->GetBinContent(ith+1,iDay+1);
						}
					} else{//mZDCSMDshiftFull doesn't exist
						for(int ith=0;ith<nShis;ith++) {
							mZDCSMD_shiftFulS_cos[icent][iDay][ith] = 0.;
							mZDCSMD_shiftFulS_sin[icent][iDay][ith] = 0.;
						}
					}
				}//j=0~nDays
			}else { //shift histograms don't exist
				for(int iDay=0; iDay<nDays; iDay++){
					for(int ith=0;ith<nShis;ith++) {
						mZDCSMD_shiftEast_cos[icent][iDay][ith] = 0.;
						mZDCSMD_shiftEast_sin[icent][iDay][ith] = 0.;
						mZDCSMD_shiftWest_cos[icent][iDay][ith] = 0.;
						mZDCSMD_shiftWest_sin[icent][iDay][ith] = 0.;
						mZDCSMD_shiftFull_cos[icent][iDay][ith] = 0.;
						mZDCSMD_shiftFull_sin[icent][iDay][ith] = 0.;
						mZDCSMD_shiftFulS_cos[icent][iDay][ith] = 0.;
						mZDCSMD_shiftFulS_sin[icent][iDay][ith] = 0.;
					}
				}//j=0~nDays
			}//else 
		}//cent
	}else {//file not exisite
		for(int icent =0;icent<10;icent++){
			for(int iDay=0; iDay<nDays; iDay++){
				for(int ith=0;ith<nShis;ith++) {
					mZDCSMD_shiftEast_cos[icent][iDay][ith] = 0.;
					mZDCSMD_shiftEast_sin[icent][iDay][ith] = 0.;
					mZDCSMD_shiftWest_cos[icent][iDay][ith] = 0.;
					mZDCSMD_shiftWest_sin[icent][iDay][ith] = 0.;
					mZDCSMD_shiftFull_cos[icent][iDay][ith] = 0.;
					mZDCSMD_shiftFull_sin[icent][iDay][ith] = 0.;
					mZDCSMD_shiftFulS_cos[icent][iDay][ith] = 0.;
					mZDCSMD_shiftFulS_sin[icent][iDay][ith] = 0.;
				}
			}//j=0~nDays
		}//cent
	}

	return 0;
}

//--------------------------------------------------------------
Int_t mZDCSMD::ReadZDCSMDFile(TFile* pZDCSMDConstFile) {
	// Read the ZDCSMD constants root file
	if (!pZDCSMDConstFile->IsOpen()) {
		cout<<"#####  No ZDCSMD constant file. Will use default."<<endl;
	} else {
		cout<<"##### ZDCSMD constant file read."<<endl; 
	}

	//for BeamCenter 
	if (pZDCSMDConstFile->IsOpen()){
		TProfile2D* mZDCSMDBeamCenter2D = (TProfile2D*) pZDCSMDConstFile->Get("ZDCSMDBeamCenter");
		for(int iRun=0;iRun<totalRunNumber;iRun++){
			//mZDCSMDCenterEx[iRun] = mZDCSMDBeamCenter2D->GetBinContent(1,iRun+1);
			//mZDCSMDCenterEy[iRun] = mZDCSMDBeamCenter2D->GetBinContent(2,iRun+1);
			//mZDCSMDCenterWx[iRun] = -1.*mZDCSMDBeamCenter2D->GetBinContent(3,iRun+1);
			////mZDCSMDCenterWx[iRun] = mZDCSMDBeamCenter2D->GetBinContent(3,iRun+1);
			//mZDCSMDCenterWy[iRun] = mZDCSMDBeamCenter2D->GetBinContent(4,iRun+1);

			//if(fabs(mZDCSMDCenterEx[iRun]-zdcsmd_ex0)>3.0)mZDCSMDCenterEx[iRun] = zdcsmd_ex0;;
			//if(fabs(mZDCSMDCenterEy[iRun]-zdcsmd_ey0)>3.0)mZDCSMDCenterEy[iRun] = zdcsmd_ey0;;
			//if(fabs(mZDCSMDCenterWx[iRun]-zdcsmd_wx0)>3.0)mZDCSMDCenterWx[iRun] = zdcsmd_wx0;;
			//if(fabs(mZDCSMDCenterWy[iRun]-zdcsmd_wy0)>3.0)mZDCSMDCenterWy[iRun] = zdcsmd_wy0;;
		}
	} else {
		for(int iRun=0;iRun<totalRunNumber;iRun++){
			mZDCSMDCenterEx[iRun] = zdcsmd_ex0;
			mZDCSMDCenterEy[iRun] = zdcsmd_ey0;
			mZDCSMDCenterWx[iRun] = zdcsmd_wx0;
			mZDCSMDCenterWy[iRun] = zdcsmd_wy0;
		}
	}

	return 0;
}

//-----------------------------------------------------------------------
void mZDCSMD::runPointer(int runId){

	dayPtr = (int)((runId)/1000%1000);
	if(nDays==1)dayPtr=0;

	int pointer=-999;
	for(int i=0; i<totalRunNumber; i++){
		if(runId==Runnumber[i]){
			pointer=i;
		}
	}
	if(pointer==-999)cout<<"Run number are not found! "<<runId<<endl;
	if(pointer==-999)pointer=0;

	runPtr =pointer;

	return ;
}

//-----------------------------------------------------------------------
int mZDCSMD::FindCentrality(Float_t refmult)
{
	int cent=0;
	for(int i=0;i<9;i++){
		if(refmult>=centZDCSMD[i])cent++;
	}

	return cent;
}



