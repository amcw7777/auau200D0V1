#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "TNtuple.h"


//StKFVertexMaker includes
//#include "TObjArray.h"
#include "StRoot/StmZDCSMDevp/mZDCSMD.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StEnumerations.h"
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
//
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
#include "TLorentzVector.h"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
class KFParticle; 
class StKFVerticesCollection; 
class TSpectrum; 
class StAnneling;
//// 

class StPicoDst;
class StPicoDstMaker;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH2D;
class TProfile;

class TMinuit;
class StKaonPion;
class StPicoPrescales;
class StRefMultCorr;


#include "TChain.h"
#include "StMaker.h"
// #include "StAnaCuts.h"
#include "StThreeVectorF.hh"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoTrack;
class StPicoEvent;
class TH1D;
class TH2D;
class TMinuit;

const int nCent = 10;
const int nVz=2;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
		void    setRunEnergyAndListDir(int run,double energy,char ListDir[256]);            

    void    DeclareHistograms();
    void    WriteHistograms();
    size_t  popcount(size_t) const;
    /*
    */
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;

    static StRefMultCorr* mRefMultCorr;
    TString    mOutName;

    TNtuple*   mD0Tuple;
    mZDCSMD *miniZDCSMD;
    bool isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isGoodTof(StPicoTrack const*) const;
    bool  isGoodEvent(StPicoEvent const*);
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    bool isTofPion(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    // int isD0Pair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const* const ) const;
    int isD0Pair50(StKaonPion const* const ) const;
    int isD0Pair150(StKaonPion const* const ) const;
    float getD0Efficiency(TLorentzVector &,bool);

    StRefMultCorr* mGRefMultCorrUtil; 
    //d0 v2 calculation
    bool  isGoodHadron(StPicoTrack const*) const;
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
    TH3D *d0MassPtY;
    TH3D *d0BarMassPtY;
    TH3D *d0MassPhiEta;
    TH3D *d0BarMassPhiEta;
    TH3D *d0MassPhiEta_noweight;
    TH3D *d0BarMassPhiEta_noweight;
    TH3D *d0MassPhiEta_5bin;
    TH3D *d0BarMassPhiEta_5bin;
    TH3D *d0MassPhiEta_50;
    TH3D *d0BarMassPhiEta_50;
    TH3D *d0MassPhiEta_150;
    TH3D *d0BarMassPhiEta_150;
    bool  readBadList();            
		bool  readRunList();            
    bool  removeBadID(int runnumber) const;            
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

    /////////for D0 efficiency // last is for all centrality bins
    TH3F *hD0Mc[10];
    TH3F *hD0Rc[10];
    TH3F *hD0BarMc[10];
    TH3F *hD0BarRc[10];

    //
    //	//StKFVertexMaker private
    ClassDef(StMyAnalysisMaker, 1)
};
#endif
