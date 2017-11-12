#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StAnaCuts.h"
#include "StThreeVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoTrack;
class TH1D;
class TH2D;
class TH2F;


class StPicoD0AnaMaker : public StMaker
{
  public:
    // StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
    StPicoD0AnaMaker(char const * name, 
        char const * outName,StPicoDstMaker* picoDstMaker);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();
    virtual void Clear(Option_t *opt="");

    int getEntries() const;


  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const* const ) const;
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    // bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const& , StThreeVectorF const& , double ) const;
    // int StPicoD0AnaMaker::getCentBin(float );
    // bool StPicoD0AnaMaker::isGoodTrack_basic(StPicoTrack const& , StThreeVectorF const&, double ) const;
    // bool StPicoD0AnaMaker::isGoodTrack_basic(StPicoTrack const& , StThreeVectorF const& , double ) const;
    // bool StPicoD0AnaMaker::isGoodTrack_hft(StPicoTrack const& , StThreeVectorF const& , double ) const
    //
      StPicoDstMaker* mPicoDstMaker;
    // StPicoD0Event* mPicoD0Event;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    TH2D *hPVZ;
    TH1D *deltaVZ;

    TH1D *hEta;
    TH1D *hTPCEta;
    TH1D *hHFTEta;
    TH1D *hPhi;
    TH1D *hTPCPhi;
    TH1D *hHFTPhi;
    TH1D *hPt;
    TH1D *hTPCPt;
    TH1D *hHFTPt;
    TH2D *hDcaXY;
    TH2D *hTPCDcaXY;
    TH2D *hHFTDcaXY;
    TH2D *hTPCDcaXY_wide;
    TH2D *hHFTDcaXY_wide;
    TH2D *hDcaZ;
    TH2D *hTPCDcaZ;
    TH2D *hHFTDcaZ;
    TH2D *hTPCDcaZ_wide;
    TH2D *hHFTDcaZ_wide;

    TH2D *hD0MassLike;
    TH2D *hD0MassUnlike;
    TH1D *hDecayLength;
    TH1D *hKaonDca;
    TH1D *hPionDca;
    TH1D *hDcaDaughters;
    TH1D *hPointing;
    TH1D *trackTest1;
    TH1D *trackTest2;
    TNtuple *tuplePair;
  
    //d0 background ntuples
    TH2F*      hVzVpdVz;
    TNtuple* ntpD0BackgroundSameSignPt01;
    TNtuple* ntpD0BackgroundSideBandPt01;
    TNtuple* ntpD0BackgroundSameSignPt12;
    TNtuple* ntpD0BackgroundSideBandPt12;
    TNtuple* ntpD0BackgroundSameSignPt23;
    TNtuple* ntpD0BackgroundSideBandPt23;
    TNtuple* ntpD0BackgroundSameSignPt35;
    TNtuple* ntpD0BackgroundSideBandPt35;
    TNtuple* ntpD0BackgroundSameSignPt510;
    TNtuple* ntpD0BackgroundSideBandPt510;


    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate

    ClassDef(StPicoD0AnaMaker, 0)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}


#endif
