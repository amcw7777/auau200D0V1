#ifndef mZDCSMD_h
#define mZDCSMD_h

#include "math.h"
#include "string.h"

#include "TString.h"
#include "TObject.h"
#include <vector>
#include <algorithm>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TDirectory.h"
#include "TVector3.h"

#include "runNumber.h"

#include <iostream>
using namespace std;


const double centZDCSMD[9] = {10,22,43,76,125,193,281,396,466}; 
//const Int_t nDays=//(Runnumber[totalRunNumber-1]/1000%1000 - Runnumber[0]/1000%1000);

const Double_t mZDCSMDGain[2][2][8] = {
{{1,1.20799,1.10295,1.10412,0.906926,0.914594,1.03453,1},
{1,1.27544,1.33041,1.32183,1.20831,1.47612,1.51068,1.5171}},
{{1,1.32221,1.55757,1.27181,1.06307,1.03059,1.10316,1},
{1,1.01041,1.13421,0.886185,0.960837,1.30616,1.11254,1.0531}}
};
//const Double_t zdcsmd_ex0 = 0.; 
//const Double_t zdcsmd_ey0 = 0.; 
//const Double_t zdcsmd_wx0 = 0.; 
//const Double_t zdcsmd_wy0 = 0.; 

const Double_t zdcsmd_ex0 = 4.96503; 
const Double_t zdcsmd_ey0 = 6.07288; 
const Double_t zdcsmd_wx0 = 4.83591; 
const Double_t zdcsmd_wy0 = 5.49622; 


const int nZDCRuns=34;
const int badRunzdc[nZDCRuns]={79, 80, 81, 82, 110, 113, 121, 122, 123, 124, 125, 126, 128, 129, 145, 259, 315, 316, 478, 479, 480, 481, 482, 483, 484, 485, 754, 782, 786, 1086, 1257, 1258, 1301, 1303, };
const int badRunNzdc[nZDCRuns]={17065016, 17065017, 17065018, 17065019, 17066004, 17066007, 17066021, 17066022, 17066023, 17066024, 17066025, 17066026, 17066028, 17066029, 17067001, 17071056, 17073022, 17073023, 17098037, 17098038, 17098039, 17098040, 17098043, 17098045, 17098046, 17098047, 17107042, 17108034, 17108038, 17116065, 17123021, 17123022, 17124039, 17124041, };

class mZDCSMD {

	public:

		mZDCSMD() ;       //  Constructor
		virtual          ~mZDCSMD() ;       //  Destructor

		//void SetmHistFileName(TString name) {mHistOutputFileName = name;}
		void  SetFileDirectory(TString name); 
		void  SetmHistFileName(TString name); 
		void  calibrateZDCSMDevp(); 
		void  WriteHist(); 

		void  SetZDCSMDrefm(Float_t *nZDCSMD,Float_t refmult,Int_t runN);  
		void  SetZDCSMDcent(Float_t *nZDCSMD,Int_t centrality,Int_t runN);  
		bool  ZDCSMDbadrun();
		bool  ZDCSMDgoodEvent();
		void  InitEvent(); 
		bool  InitRun();          
		bool  InitRun(Int_t mod);   // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6  
		void  SetMod(Int_t mod);    // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6  
		Int_t GetMod(){return mMod;} 
		// phiweight and shift are two independent method here!
		// pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,   ready+wQAHists 5,  ready+woQAHists 6  
		
		TVector2  GetZDCSMD_QFulSS();
		TVector2  GetZDCSMD_QFullS();
		TVector2  GetZDCSMD_QEastS();
		TVector2  GetZDCSMD_QWestS();

		// corrected + recentered Q
		TVector2  GetZDCSMD_QFull();
		TVector2  GetZDCSMD_QEast();
		TVector2  GetZDCSMD_QWest();
		// raw Q
		TVector2  GetZDCSMD_rawQFull();
		TVector2  GetZDCSMD_rawQEast();
		TVector2  GetZDCSMD_rawQWest();
		// ped/gain corrected Q
		TVector2  GetZDCSMD_corQFull();
		TVector2  GetZDCSMD_corQEast();
		TVector2  GetZDCSMD_corQWest();

		Float_t   GetZDCSMD_PsiEast();
		Float_t   GetZDCSMD_PsiWest();
		Float_t   GetZDCSMD_PsiFull();
		Float_t   GetZDCSMD_PsiEastS();
		Float_t   GetZDCSMD_PsiWestS();
		Float_t   GetZDCSMD_PsiFulSS();
		Float_t   GetZDCSMD_PsiFullS();
		Float_t   GetZDCSMDadcR(Int_t eastwest,Int_t verthori,Int_t strip) const;    // raw adc
		Float_t   GetZDCSMDadcC(Int_t eastwest,Int_t verthori,Int_t strip) const;    // ped/gain corrected
		Float_t   GetZDCSMDpeds(Int_t eastwest,Int_t verthori,Int_t strip) const;    // pedestal 
		Float_t   GetZDCSMD_shiftFull_sin(Int_t ith);
		Float_t   GetZDCSMD_shiftFull_cos(Int_t ith);
		Float_t   GetZDCSMD_shiftEast_sin(Int_t ith);
		Float_t   GetZDCSMD_shiftEast_cos(Int_t ith);
		Float_t   GetZDCSMD_shiftWest_sin(Int_t ith);
		Float_t   GetZDCSMD_shiftWest_cos(Int_t ith);
		Float_t   GetZDCSMD_shiftFulS_sin(Int_t ith);
		Float_t   GetZDCSMD_shiftFulS_cos(Int_t ith);
		Float_t   GetZDCSMD_EastX0();
		Float_t   GetZDCSMD_EastY0();
		Float_t   GetZDCSMD_WestX0();
		Float_t   GetZDCSMD_WestY0();
		Float_t   GetZDCSMD_moreshiftFull_sin(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftFull_cos(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftEast_sin(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftEast_cos(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftWest_sin(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftWest_cos(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftFulS_sin(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreshiftFulS_cos(Int_t cent,Int_t iDay,Int_t ith);
		Float_t   GetZDCSMD_moreEastX0(Int_t iRun);
		Float_t   GetZDCSMD_moreEastY0(Int_t iRun);
		Float_t   GetZDCSMD_moreWestX0(Int_t iRun);
		Float_t   GetZDCSMD_moreWestY0(Int_t iRun);


	protected:


	private:

		TRandom3 *gRandom;

		TString   mFileDirectory;
		TString   mHistOutputFileName ;
		TFile*    mhist_output;

		// phiweight and shift are two independent method here!
		Int_t     mMod;   // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,   ready 5,

		const Float_t twopi = 2.*TMath::Pi();
		static const Int_t nPsiBins=64;
		static const Int_t nDays = 365;

		static const Int_t nShis = 12;   // orders of shift 

		// raw
		TVector2  mZDCSMD_rawQFull;
		TVector2  mZDCSMD_rawQEast;
		TVector2  mZDCSMD_rawQWest;

		// gain corrected 
		TVector2  mZDCSMD_corQFull;
		TVector2  mZDCSMD_corQEast;
		TVector2  mZDCSMD_corQWest;

		// gain + recentered 
		TVector2  mZDCSMD_QFull;
		TVector2  mZDCSMD_QEast;
		TVector2  mZDCSMD_QWest;

		// shifted
		TVector2  mZDCSMD_QFulSS;
		TVector2  mZDCSMD_QFullS;
		TVector2  mZDCSMD_QEastS;
		TVector2  mZDCSMD_QWestS;

		Float_t   mZDCSMD_PsiEast;
		Float_t   mZDCSMD_PsiWest;
		Float_t   mZDCSMD_PsiFull;
		Float_t   mZDCSMD_PsiEastS;
		Float_t   mZDCSMD_PsiWestS;
		Float_t   mZDCSMD_PsiFulSS;
		Float_t   mZDCSMD_PsiFullS;


		//////////////////////////
		Float_t  mrefmult;
		Int_t    mcent;    // 1-9
		Int_t    mrun;
		Int_t    dayPtr; 
		Int_t    runPtr; 
		Float_t  mZDCSMDadcR[2][2][8];                       // ZDCSMD
		Float_t  mZDCSMDadcC[2][2][8];                       // ZDCSMD
		Float_t  mZDCSMDpeds[2][2][8];                       // ZDCSMDPed
		Float_t  mZDCSMDCenterEx[totalRunNumber];            //! ZDCSMD Beam Center 
		Float_t  mZDCSMDCenterEy[totalRunNumber];
		Float_t  mZDCSMDCenterWx[totalRunNumber];
		Float_t  mZDCSMDCenterWy[totalRunNumber];

		Float_t  mZDCSMD_shiftEast_cos[10][nDays][nShis];        //!ZDCSMD east shift
		Float_t  mZDCSMD_shiftEast_sin[10][nDays][nShis];        //!ZDCSMD east shift
		Float_t  mZDCSMD_shiftWest_cos[10][nDays][nShis];        //!ZDCSMD west shift
		Float_t  mZDCSMD_shiftWest_sin[10][nDays][nShis];        //!ZDCSMD west shift
		Float_t  mZDCSMD_shiftFull_sin[10][nDays][nShis];        //!ZDCSMD full shift
		Float_t  mZDCSMD_shiftFull_cos[10][nDays][nShis];        //!ZDCSMD full shift
		Float_t  mZDCSMD_shiftFulS_cos[10][nDays][nShis];        //!ZDCSMD fulS shift
		Float_t  mZDCSMD_shiftFulS_sin[10][nDays][nShis];        //!ZDCSMD fulS shift

		//ZDC beam center
		TProfile2D *mZDCSMDBeamCenter;
	
		// raw                        // corrected 	
		TH2F *mHistZDCSMDEastVR[8];   TH2F *mHistZDCSMDEastVC[8];
		TH2F *mHistZDCSMDEastHR[8];   TH2F *mHistZDCSMDEastHC[8];
		TH2F *mHistZDCSMDWestVR[8];   TH2F *mHistZDCSMDWestVC[8];
		TH2F *mHistZDCSMDWestHR[8];   TH2F *mHistZDCSMDWestHC[8];

		//raw psi
		TH2F *mHistZDCSMDPsi_R_East[10];
		TH2F *mHistZDCSMDPsi_R_West[10];
		TH2F *mHistZDCSMDPsi_R_Full[10];
		
		//for shift file
		TProfile2D *mHistZDCSMDshiftEast_c[10];
		TProfile2D *mHistZDCSMDshiftEast_s[10];
		TProfile2D *mHistZDCSMDshiftWest_c[10];
		TProfile2D *mHistZDCSMDshiftWest_s[10];
		TProfile2D *mHistZDCSMDshiftFull_c[10];
		TProfile2D *mHistZDCSMDshiftFull_s[10];
		TProfile2D *mHistZDCSMDshiftFulS_c[10];
		TProfile2D *mHistZDCSMDshiftFulS_s[10];
		
		//ZDCSMD psi after shift
		TH2F *mHistZDCSMDPsi_S_FulS[10];
		TH2F *mHistZDCSMDPsi_S_East[10];
		TH2F *mHistZDCSMDPsi_S_West[10];
		TH2F *mHistZDCSMDPsi_S_Full[10];
		
		//correlation between ZDCSMD
		// cos(n*delta_Psi) ZDCSMD
		TProfile *mHistCosZ;
		TProfile *mHistSinZ;


		Int_t  GetZDCSMD_PsiBins(Float_t psi){
			Int_t n=0;
			Int_t n0=fabs(psi)/twopi;
			if(psi<0.)psi+=twopi*(n0+1);
			if(psi>0.)psi-=twopi*n0;
			n = psi*nPsiBins/twopi;

			return n;
		}

		Int_t  FindCentrality(Float_t refmult);
		void   runPointer(Int_t runN);
		void   Fill();
		void   pedGain();

		void   InitHist();
		Int_t  openFile(TString name);  //open file 
		Int_t  ReadZDCSMDFile(TFile *name);                     // get the ZDCSMD constants
		Int_t  ReadZDCSMDshiftFile(TFile *name);                // get the ZDCSMD shift

		//for ZDCSMD file
		Float_t  ZDCSMD_GetPosition(Int_t eastwest,Int_t verthori,Int_t strip);
		void     SetZDCSMDadcR(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmd);
		void     SetZDCSMDadcC(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmd);
		void     SetZDCSMDpeds(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmd);

		
#ifdef __ROOT__
		ClassDef(mZDCSMD, 0)
#endif


};

#endif

inline void    mZDCSMD::SetZDCSMDadcR(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmdadcR) {mZDCSMDadcR[eastwest][verthori][strip-1] =(zdcsmdadcR>0.)? zdcsmdadcR:0.;}
inline void    mZDCSMD::SetZDCSMDadcC(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmdadcC) {mZDCSMDadcC[eastwest][verthori][strip-1] =(zdcsmdadcC>0.)? zdcsmdadcC:0.;}
inline void    mZDCSMD::SetZDCSMDpeds(Int_t eastwest,Int_t verthori,Int_t strip,const Float_t zdcsmdpeds) {mZDCSMDpeds[eastwest][verthori][strip-1] =(zdcsmdpeds>0.)? zdcsmdpeds:0.;}
inline Float_t mZDCSMD::GetZDCSMDadcR(Int_t eastwest,Int_t verthori,Int_t strip)const {return mZDCSMDadcR[eastwest][verthori][strip-1];}
inline Float_t mZDCSMD::GetZDCSMDadcC(Int_t eastwest,Int_t verthori,Int_t strip)const {return mZDCSMDadcC[eastwest][verthori][strip-1];}
inline Float_t mZDCSMD::GetZDCSMDpeds(Int_t eastwest,Int_t verthori,Int_t strip)const {return mZDCSMDpeds[eastwest][verthori][strip-1];}

inline TVector2  mZDCSMD::GetZDCSMD_QFulSS() {return mZDCSMD_QFulSS;}
inline TVector2  mZDCSMD::GetZDCSMD_QFullS() {return mZDCSMD_QFullS;}
inline TVector2  mZDCSMD::GetZDCSMD_QEastS() {return mZDCSMD_QEastS;}
inline TVector2  mZDCSMD::GetZDCSMD_QWestS() {return mZDCSMD_QWestS;}

// corrected + recentered Q
inline TVector2  mZDCSMD::GetZDCSMD_QFull() {return mZDCSMD_QFull;}
inline TVector2  mZDCSMD::GetZDCSMD_QEast() {return mZDCSMD_QEast;}
inline TVector2  mZDCSMD::GetZDCSMD_QWest() {return mZDCSMD_QWest;}

// raw Q
inline TVector2  mZDCSMD::GetZDCSMD_rawQFull() {return mZDCSMD_rawQFull;}
inline TVector2  mZDCSMD::GetZDCSMD_rawQEast() {return mZDCSMD_rawQEast;}
inline TVector2  mZDCSMD::GetZDCSMD_rawQWest() {return mZDCSMD_rawQWest;}

// ped/gain corrected Q
inline TVector2  mZDCSMD::GetZDCSMD_corQFull() {return mZDCSMD_corQFull;}
inline TVector2  mZDCSMD::GetZDCSMD_corQEast() {return mZDCSMD_corQEast;}
inline TVector2  mZDCSMD::GetZDCSMD_corQWest() {return mZDCSMD_corQWest;}

inline Float_t   mZDCSMD::GetZDCSMD_PsiEast()  {return mZDCSMD_PsiEast ;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiWest()  {return mZDCSMD_PsiWest ;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiFull()  {return mZDCSMD_PsiFull ;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiEastS() {return mZDCSMD_PsiEastS;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiWestS() {return mZDCSMD_PsiWestS;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiFulSS() {return mZDCSMD_PsiFulSS;}  
inline Float_t   mZDCSMD::GetZDCSMD_PsiFullS() {return mZDCSMD_PsiFullS;}  

inline Float_t   mZDCSMD::GetZDCSMD_shiftFull_sin(Int_t ith){return mZDCSMD_shiftFull_sin[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftFull_cos(Int_t ith){return mZDCSMD_shiftFull_cos[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftEast_sin(Int_t ith){return mZDCSMD_shiftEast_sin[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftEast_cos(Int_t ith){return mZDCSMD_shiftEast_cos[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftWest_sin(Int_t ith){return mZDCSMD_shiftWest_sin[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftWest_cos(Int_t ith){return mZDCSMD_shiftWest_cos[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftFulS_sin(Int_t ith){return mZDCSMD_shiftFulS_sin[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_shiftFulS_cos(Int_t ith){return mZDCSMD_shiftFulS_cos[mcent][dayPtr][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_EastX0(){return mZDCSMDCenterEx[runPtr];}
inline Float_t   mZDCSMD::GetZDCSMD_EastY0(){return mZDCSMDCenterEy[runPtr];}
inline Float_t   mZDCSMD::GetZDCSMD_WestX0(){return mZDCSMDCenterWx[runPtr];}
inline Float_t   mZDCSMD::GetZDCSMD_WestY0(){return mZDCSMDCenterWy[runPtr];}

inline Float_t   mZDCSMD::GetZDCSMD_moreshiftFull_sin(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftFull_sin[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftFull_cos(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftFull_cos[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftEast_sin(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftEast_sin[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftEast_cos(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftEast_cos[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftWest_sin(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftWest_sin[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftWest_cos(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftWest_cos[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftFulS_sin(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftFulS_sin[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreshiftFulS_cos(Int_t cent,Int_t iDay,Int_t ith) {return mZDCSMD_shiftFulS_cos[cent][iDay][ith];}  
inline Float_t   mZDCSMD::GetZDCSMD_moreEastX0(Int_t iRun){return mZDCSMDCenterEx[iRun];}
inline Float_t   mZDCSMD::GetZDCSMD_moreEastY0(Int_t iRun){return mZDCSMDCenterEy[iRun];}
inline Float_t   mZDCSMD::GetZDCSMD_moreWestX0(Int_t iRun){return mZDCSMDCenterWx[iRun];}
inline Float_t   mZDCSMD::GetZDCSMD_moreWestY0(Int_t iRun){return mZDCSMDCenterWy[iRun];}

//-------------------------------------------------------------
inline void  mZDCSMD::SetFileDirectory(TString name){ mFileDirectory = name; }
inline void  mZDCSMD::SetmHistFileName(TString name){
	if(mMod==1)mHistOutputFileName= name+"_ZDCSMD_ped.root"    ; 
	if(mMod==2)mHistOutputFileName= name+"_ZDCSMD_center.root"   ; 
	if(mMod==3)mHistOutputFileName= name+"_ZDCSMD_evpsub.root" ; 
	if(mMod==4)mHistOutputFileName= name+"_ZDCSMD_evpfull.root"; 
	if(mMod>=5)mHistOutputFileName= name+"_ZDCSMD_ready.root"  ; 
}	



