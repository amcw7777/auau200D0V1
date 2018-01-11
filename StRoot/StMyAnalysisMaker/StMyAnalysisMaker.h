#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "TNtuple.h"


//StKFVertexMaker includes
//#include "TObjArray.h"
#include "StRoot/StmZDCSMDevp/mZDCSMD.h"
// #include "StRoot/StRefMultCorr/StRefMultCorr.h"
// #include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StEnumerations.h"
#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
//
class StEvent;
//// 

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH2D;
class TProfile;

class StPicoPrescales;
class StRefMultCorr;



class TString;
class TFile;
class TNtuple;
class TH1D;
class TH2D;

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
    bool  isGoodEvent(StPicoEvent const*);

    /*
    */
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;

    static StRefMultCorr* mRefMultCorr;
    TString    mOutName;

    TNtuple*   mD0Tuple;
    mZDCSMD *miniZDCSMD;

    StRefMultCorr* mGRefMultCorrUtil; 
    //
    //	//StKFVertexMaker private
    ClassDef(StMyAnalysisMaker, 1)
};
#endif
