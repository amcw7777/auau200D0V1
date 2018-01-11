#include "StMyAnalysisMaker.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "phys_constants.h"
#include "StCuts.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"


#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

//
//

//
ClassImp(StMyAnalysisMaker)
  StRefMultCorr* StMyAnalysisMaker::mRefMultCorr = NULL;
  //-----------------------------------------------------------------------------
  StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
: StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {

  miniZDCSMD=new mZDCSMD();
  // phiweight and shift are two independent method here!
  // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6 
  // miniZDCSMD->SetFileDirectory("/star/u/amcw7777/d0V1AuAu2016/ZDCSMDFile");
  miniZDCSMD->SetFileDirectory("/global/homes/a/amcw7777/auau200GeVD0V1/ZDCSMDFile");
  if(!(miniZDCSMD->InitRun(5)))return kStFatal;  
  miniZDCSMD->SetmHistFileName(mOutName);
  if(!mRefMultCorr){
    mGRefMultCorrUtil = new StRefMultCorr("grefmult");
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    miniZDCSMD->WriteHist();
  }

  return kStOK;
}


//----------------------------------------------------------------------------- 
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Make() {
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  // const int  runID    = event->runId();
  const int  runID    = 17109018;//event->runId();
  // const int  evtID    = event->eventId();
  // const int refMult  = event->grefMult();
  if(!(isGoodEvent(event)))
  { 
    // LOG_WARN << " Not Min Bias! Skip! " << endm;
    return kStWarn;
  }

  // cout<<"this is a good event!"<<endl;
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  StThreeVectorF vtx = event->primaryVertex();
  float b = event->bField();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  if(event) {
    pVtx = event->primaryVertex();
  }
  ////////////  ZDCSMD ////////////////
  mGRefMultCorrUtil->init(mPicoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(mPicoDst->event()->grefMult(),pVtx.z(),mPicoDst->event()->ZDCx()) ;
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();

  // int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  // Bin       Centrality (16)   Centrality (9)
  //     -1           80-100%           80-100% // this one should be rejected in your centrality related analysis
  //     0            75-80%            70-80%
  //     1            70-75%            60-70%
  //     2            65-70%            50-60%
  //     3            60-65%            40-50%
  //     4            55-60%            30-40%
  //     5            50-55%            20-30%
  //     6            45-50%            10-20%
  //     7            40-45%             5-10%
  //     8            35-40%             0- 5%
  int cent = centrality+1;  

  //
  // ////////////  ZDCSMD ////////////////
  float ZDCSMDadc[32];
  for(int i=0;i<8;i++){
    ZDCSMDadc[i]   = 1.* event->ZdcSmdEastHorizontal(i);   // picoDst function i 0-7
    ZDCSMDadc[i+8] = 1.* event->ZdcSmdEastVertical(i);
    ZDCSMDadc[i+16]= 1.* event->ZdcSmdWestHorizontal(i);
    ZDCSMDadc[i+24]= 1.* event->ZdcSmdWestVertical(i);
  }


  miniZDCSMD->InitEvent();
  miniZDCSMD->SetZDCSMDcent(ZDCSMDadc,cent,runID);
  if(miniZDCSMD->ZDCSMDbadrun())return kStOK;
  miniZDCSMD->calibrateZDCSMDevp();
  if(!(miniZDCSMD->ZDCSMDgoodEvent()))return kStOK;
  int mMod=miniZDCSMD->GetMod();

  return kStOK;
}
//------------------------------------------------------------------------------
bool StMyAnalysisMaker::isGoodEvent(StPicoEvent const*mEvent)
{	

  const int  runID    = mEvent->runId();
  const int  evtID    = mEvent->eventId();
  const int  refMult  = mEvent->grefMult();
  // const int grefMult  = mEvent->grefMult();
  // const int  ranking  = mEvent->ranking();

  if(!mEvent->isTrigger(450050) &&
      !mEvent->isTrigger(450060)&&
      !mEvent->isTrigger(450005)&&
      !mEvent->isTrigger(450015)&&
      !mEvent->isTrigger(450025))
    return false; 

  StThreeVectorF Vertex3D=mEvent->primaryVertex();
  const double VertexX = Vertex3D.x(); 
  const double VertexY = Vertex3D.y(); 
  const double VertexZ = Vertex3D.z(); 
  const double vpdVz   = mEvent->vzVpd();

  //event cut
  if(refMult <=2 || refMult > 1000) return false;
  if(fabs(VertexZ) > 100) return false; 

  // hVertex2D ->Fill(VertexZ,vpdVz);
  // hDiffVz   ->Fill(VertexZ-vpdVz); 

  if(fabs(VertexZ) > 6) return false; 
  if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0)return false; 
  if(fabs(VertexZ-vpdVz)>3.)return false;       // no vpd cut in low energy?

  return true;

}

