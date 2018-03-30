#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/***********************************************************************************
 **
 ** D0CorrelationV2Analyser
 **
 ** Author: Leon He
 ************************************************************************************
 **
 ** Description: 
 **
 ************************************************************************************
 **
 ** Log:
 **
 ********************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
//
#include "StThreeVectorF.hh"
#include "TLorentzVector.h"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
#include "StRoot/StmZDCSMDevp/mZDCSMD.h"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
//// 

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StPicoEvent;
class StKaonPion;
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StHFCuts;
class StPicoPrescales;
class StRefMultCorr;



class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);    
    ofstream fout1;

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();
    ofstream fout;

    bool  isGoodHadron(StPicoTrack const*) const;
    bool  isGoodEvent(StPicoEvent const*);
    int isD0Pair(StKaonPion const*) const;
    int isD0Pair50(StKaonPion const*) const;
    int isD0Pair150(StKaonPion const*) const;
    bool  isGoodTrack(StPicoTrack const*, StPicoEvent const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    bool isTofPion(StPicoTrack const* const, float beta) const;
    bool isTofKaon(StPicoTrack const* const, float beta, StPicoEvent const*) const;
    bool isTofPion(StPicoTrack const* const, float beta, StPicoEvent const*) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    float getD0Efficiency(TLorentzVector &);

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
    //StPicoDstMaker *
    StPicoDst *picoDst;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    //TFile* mPhi;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    //d0 v1 calculation
    mZDCSMD *miniZDCSMD;
    int mPrevRunId;
    double     PI;
    double     twoPI;
    std::vector<Int_t> mBadList;
		std::vector<Int_t> mRunList;
    int        mRun;            
    int mTempRunId;
    int mTempEvtId;
    double     mEnergy;            
    double     mWght;            
    TString    mListDir;            
    TH3D *d0MassPt;
    TH3D *d0BarMassPt;
    TH3D *d0MassPhiEta;
    TH3D *d0BarMassPhiEta;

    TH3D *d0MassPhiEta_50;
    TH3D *d0BarMassPhiEta_50;
    TH3D *d0MassPhiEta_150;
    TH3D *d0BarMassPhiEta_150;
    TH3D *d0MassPhiEta_5bin;
    TH3D *d0BarMassPhiEta_5bin;

    Int_t CheckrunNumber(int runnumber) const;            
    TH1D *zdcPsi;
    TH1D *zdcPsi_corr;
    TH1D *testEvent;
    TProfile *pionV1Plus;
    TProfile *pionV1Minus;
    TProfile *zdcResolution;
    TH1D *testTrack;
    TH2D *trackPhiEta;
    TH1D *testDPhi;
    TH2D *trackPhiEtaHFT;
    TH1D *hHitsDedx;
    TH2D *hSigmaPiBeta;
    // TNTuple *d0_canditates_tuple;
    enum 
    {
      NCENTRAPEFF=3,
      NRAPEFF=20,
    };
    TH3F *h3PtCentY_clone;
    TH3F *h3PtCentYcut_clone;

    double efficiency[4][6];
    ClassDef(StPicoD0AnaMaker, 1)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
