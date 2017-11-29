#include "StMyAnalysisMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TLorentzVector.h"
// #include "StPicoDstMaker/StPicoV0.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include <vector>
#include "phys_constants.h"
#include "StCuts.h"

#include "StBTofUtil/tofPathLength.hh"
#include "TRandom3.h"
#include "StPhysicalHelixD.hh"

#include "TRMatrix.h"
#include "TRSymMatrix.h"
#include "TMinuit.h"

#include "StBTofUtil/tofPathLength.hh"
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
  mRun    = 16;
  mEnergy =200;
  mListDir="./";
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();
  PI   = 3.14159;
  twoPI= 2.*PI;
  mPrevRunId=-999;
  mTempRunId=-999;
  mTempEvtId=-999;

  mBadList.clear();
  mRunList.clear();

  if(!readRunList())return kStFatal;
  if(!readBadList())return kStFatal;

  miniZDCSMD=new mZDCSMD();
  // phiweight and shift are two independent method here!
  // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6 
  miniZDCSMD->SetFileDirectory("/star/u/amcw7777/d0V1AuAu2016/ZDCSMDFile");
  if(!(miniZDCSMD->InitRun(6)))return kStFatal;  
  miniZDCSMD->SetmHistFileName(mOutName);
  if(!mRefMultCorr){
    mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id() ;
    mRefMultCorr->setVzForWeight(6, -6.0, 6.0);
    mRefMultCorr->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt");
    /////////////////////////////////////
    // mGRefMultCorrUtil = new StRefMultCorr("grefmult");
  }

  d0MassPhiEta = new TH3D("d0MassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-1,1);
  d0MassPt = new TH2D("d0MassPt",";D^{0} mass (GeV/c^{2});p_{T} (GeV/c)",50,1.6,2.1,100,0,10);
  d0BarMassPhiEta = new TH3D("d0BarMassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-1,1);
  testDPhi = new TH1D("testDPhi","testDPhi",400,-1*twoPI,twoPI);
  testEvent = new TH1D("testEvent","",10,-0.5,9.5);
  zdcPsi = new TH1D("zdcPsi",";#psi_{ZDC}",1000,0,twoPI);
  zdcPsi_corr = new TH1D("zdcPsi_corr",";#psi_{ZDC}",1000,0,twoPI);
  pionV1Plus = new TProfile("pionV1Plus","",48,-1.2,1.2);
  pionV1Minus = new TProfile("pionV1Minus","",48,-1.2,1.2);
  zdcResolution = new TProfile("zdcResolution","",10,0,10);
  testTrack = new TH1D("testTrack","",10,0,10);
  trackPhiEta = new TH2D("trackPhiEta",";#phi;#eta",100,-1.*PI,PI,10,-1,1);
  trackPhiEtaHFT = new TH2D("trackPhiEtaHFT",";#phi;#eta",100,-1.*PI,PI,10,-1,1);
  pionV1Plus->Sumw2();
  pionV1Minus->Sumw2();
  zdcResolution->Sumw2();
  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Write();
    fout->Close();
  }

  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {

}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  d0MassPhiEta->Write();
  d0BarMassPhiEta->Write();
  testDPhi->Write();
  testEvent->Write();
  d0MassPt->Write();
  zdcPsi->Write();
  zdcPsi_corr->Write();
  pionV1Plus->Write();
  pionV1Minus->Write();
  zdcResolution->Write();
  testTrack->Write();
  trackPhiEta->Write();
  trackPhiEtaHFT->Write();
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
  const int  runID    = event->runId();
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
  // double vpdvz = mPicoDst->event()->vzVpd();
  // double vz = vtx.z();
  // if( fabs(vz)>6 && fabs(vz-vpdvz)>3 ) 
  //   return kStWarn;
  //StRefMultCorr
  // if ( runID != mPrevRunId ) {
  // 	mRefMultCorr->init(runID);
  // 	mPrevRunId = runID;
  // 	cout << "reset mPrevRunId = " << mPrevRunId << endl;
  // }
  // //mRefMultCorr->initEvent(refMult,VertexZ);
  // mRefMultCorr->initEvent(refMult,vz,event->ZDCx());
  //
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  // double mWght    = mRefMultCorr->getWeight();
  //
  // if( centrality<0||centrality>=(nCent-1)) return kStOK;
  // if( centrality<2||centrality>6) return kStOK; // 10 - 60 % centrality
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
  //cout<<refMult<<" "<<cent<<" "<<mRefMultCorr->getCentralityBin16()<<endl;


  ////////////  ZDCSMD ////////////////
  float ZDCSMDadc[32];
  for(int i=0;i<8;i++){
    ZDCSMDadc[i]   = 1.* event->ZdcSmdEastHorizontal(i);   // picoDst function i 0-7
    ZDCSMDadc[i+8] = 1.* event->ZdcSmdEastVertical(i);
    ZDCSMDadc[i+16]= 1.* event->ZdcSmdWestHorizontal(i);
    ZDCSMDadc[i+24]= 1.* event->ZdcSmdWestVertical(i);
  }


  miniZDCSMD->InitEvent();
  //miniZDCSMD->SetZDCSMDrefm(ZDCSMDadc,mult_corr,runID);;
  miniZDCSMD->SetZDCSMDcent(ZDCSMDadc,cent,runID);
  if(miniZDCSMD->ZDCSMDbadrun())return kStOK;
  testEvent->Fill(8);
  miniZDCSMD->calibrateZDCSMDevp();
  if(!(miniZDCSMD->ZDCSMDgoodEvent()))return kStOK;
  testEvent->Fill(9);
  //int mMod=miniZDCSMD->GetMod();

  float mZDC1Event_PsiF = miniZDCSMD->GetZDCSMD_PsiFulSS(); 
  float mZDC1Event_PsiOrigin = miniZDCSMD->GetZDCSMD_PsiFull(); 
  float mZDC1Event_PsiW = miniZDCSMD->GetZDCSMD_PsiWestS();
  float mZDC1Event_PsiE = miniZDCSMD->GetZDCSMD_PsiEastS();

  if(mZDC1Event_PsiF<0.)    mZDC1Event_PsiF+=twoPI;
  if(mZDC1Event_PsiW<0.)    mZDC1Event_PsiW+=twoPI;
  if(mZDC1Event_PsiE<0.)    mZDC1Event_PsiE+=twoPI;
  if(mZDC1Event_PsiF>twoPI) mZDC1Event_PsiF-=twoPI;
  if(mZDC1Event_PsiW>twoPI) mZDC1Event_PsiW-=twoPI;
  if(mZDC1Event_PsiE>twoPI) mZDC1Event_PsiE-=twoPI;
  /////////////////////////////////

  zdcPsi->Fill(mZDC1Event_PsiOrigin,mWght);
  zdcPsi_corr->Fill(mZDC1Event_PsiF,mWght);
	zdcResolution->Fill(cent,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 
	zdcResolution->Fill(0.,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 

  vector<int> kaonIndex;
  vector<int> pionIndex;
  vector<int> testPionIndex;
  kaonIndex.clear();
  pionIndex.clear();
  for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++)
  {
    StPicoTrack const* itrk = mPicoDst->track(i);
    if(!isGoodTrack(itrk))  continue;
    StThreeVectorF mom= itrk->gMom();
    double eta  = mom.pseudoRapidity();
    double phi  = mom.phi();
    double dca = (vtx - itrk->dcaPoint()).mag();
    trackPhiEta->Fill(phi,eta);
    if (isTpcPion(itrk) && fabs(dca)<1.5) 
      testPionIndex.push_back(i);
    if(!(itrk->isHft())) continue;
    trackPhiEtaHFT->Fill(phi,eta);
    testTrack->Fill(1);
    // if(!isGoodTof(itrk))  continue;
    testTrack->Fill(2);

    if (isTpcPion(itrk)) 
      pionIndex.push_back(i);
    bool tpcKaon = isTpcKaon(itrk,&vtx);
    float kBeta = getTofBeta(itrk,&vtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(itrk,kBeta);

    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(goodKaon) 
      kaonIndex.push_back(i);
  }
  for(unsigned int i=0;i<testPionIndex.size();i++)
  {
    if(centrality!=4) continue;
    StPicoTrack const* itrk = mPicoDst->track(testPionIndex[i]);
    StThreeVectorF mom= itrk->gMom();
    double eta  = mom.pseudoRapidity();
    TLorentzVector pionLorentz;
    pionLorentz.SetPtEtaPhiM(mom.perp(),eta,mom.phi(),0.139);
    double pionY = pionLorentz.Rapidity();
    double phi  = mom.phi();
    if(phi<0.0) phi += twoPI;
    double mcos1=cos(phi);
    double msin1=sin(phi);
    double v1ZDC1F = mcos1*cos(1.*mZDC1Event_PsiF) + msin1*sin(1.*mZDC1Event_PsiF);
    int charge = itrk->charge();
    if(charge > 0)
      pionV1Plus->Fill(pionY,v1ZDC1F,mWght);
    else
      pionV1Minus->Fill(pionY,v1ZDC1F,mWght);
  }


  for(unsigned int i=0;i<pionIndex.size();i++)
  {
    StPicoTrack const* pion= mPicoDst->track(pionIndex[i]);
    for(unsigned int j=0;j<kaonIndex.size();j++)
    {
      StPicoTrack const* kaon= mPicoDst->track(kaonIndex[j]);
      int charge = pion->charge() * kaon->charge();
      StKaonPion *kp = new StKaonPion(kaon,pion,j,i,vtx,b);
      if(charge>0) continue;//only unlike-sign pairs
      if(kp->pt()<1.) continue;// require D0 pT > 1.5 GeV/c
      // StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
      // StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());

      int isD0= isD0Pair(kp);
      if(isD0==0)
        continue;
      // cout<<"find D0 !"<<endl;
      d0MassPt->Fill(kp->m(),kp->pt());
      if(kp->pt()<1.5) continue;// require D0 pT > 1.5 GeV/c
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      double kpY = d0Lorentz.Rapidity();
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;
      deltaPhi = fabs(deltaPhi);
      // if(deltaPhi < 0)
      // {
      //   cout<<"kp phi = "<<kp->phi()<<"\tevet psi = "<<mZDC1Event_PsiF<<endl;
      //   cout<<"original dphi = "<<deltaPhi0<<endl;
      //   cout<<"folded dphi = "<<deltaPhi<<endl;
      // }
      testDPhi->Fill(deltaPhi);
      // cout<<"deltaPhi = "<<deltaPhi<<endl;
      if(kaon->charge() < 0)
        d0MassPhiEta->Fill(kp->m(),deltaPhi,kpY,mWght);
      if(kaon->charge() > 0)
        d0BarMassPhiEta->Fill(kp->m(),deltaPhi,kpY,mWght);
      delete kp;
    }
  }


  return kStOK;
}




bool StMyAnalysisMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    // double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gPtot();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}
float StMyAnalysisMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        // StPhysicalHelixD helix = trk->helix();
        StPhysicalHelixD helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > 0.6 && trk->nHitsFit() >= 20 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 && trk->nHitsDedx() >= 16;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
bool StMyAnalysisMaker::isGoodTof(StPicoTrack const * const trk) const
{
  int index2tof = trk->bTofPidTraitsIndex();
  if(index2tof < 0) return false;
  StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
  if(tofPid->btofMatchFlag()>0)
  {
    testTrack->Fill(3);
    if(fabs(tofPid->btofYLocal()) < 1.8)
    {
      testTrack->Fill(4);
      if(tofPid->btofBeta() > 0)
      {
        testTrack->Fill(5);
      }
    }
  }

  return tofPid->btofMatchFlag() > 0 && fabs(tofPid->btofYLocal()) < 1.8 && tofPid->btofBeta() > 0;
}


int StMyAnalysisMaker::isD0Pair(StKaonPion const* const kp) const
{

  // StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  // StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  // StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  TLorentzVector d0Lorentz;
  d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
  if(fabs(d0Lorentz.Rapidity())>1.) return 0;
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }
  // int charge = kaon->charge() * pion->charge();
  // if(charge>0)
  //   charge = kaon->charge()>0 ? 1:2;
  if(pairCuts)
    // return -1 * kaon->charge();
    return 1;
  else
    return 0;
}
//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= 15 &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::removeBadID(int runnumber) const
{
  // for (std::vector<int>::iterator it=mBadList.begin(); it!=mBadList.end(); ++it) 
  for (auto it=mBadList.begin(); it!=mBadList.end(); ++it) 
  { 
    if(runnumber==*it){
      cout<<runnumber<<"is a bad run!!!"<<endl; 
      return kTRUE;
    }
  }

  return kFALSE;
}  

//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::readBadList()
{ 
  if(0) return kTRUE;
  TString inf=mListDir + "/badList/";
  inf += Form("badrun%dList%.1f.list",mRun,mEnergy);
  ifstream inrun; 
  inrun.open(inf);
  if ( inrun.fail() ) {
    cout<< "cannot open " << inf.Data() << endl;
    return kFALSE;
  }
  Int_t runid;
  while ( inrun >> runid ) { mBadList.push_back(runid); }
  inrun.close();
  sort(mBadList.begin(),mBadList.end());

  vector<int>::iterator it;
  it = std::unique (mBadList.begin(), mBadList.end());
  mBadList.resize( std::distance(mBadList.begin(),it) );

  cout <<"badrun list :" <<inf.Data() << " loaded." << endl;
  cout <<"Total       :" <<mBadList.size()<< " bad runs. "<< endl;
  return kTRUE;
}
//------------------------------------------------------------------------------
bool StMyAnalysisMaker::isGoodEvent(StPicoEvent const*mEvent)
{	

  const int  runID    = mEvent->runId();
  const int  evtID    = mEvent->eventId();
  const int  refMult  = mEvent->grefMult();
  // const int grefMult  = mEvent->grefMult();
  // const int  ranking  = mEvent->ranking();

  testEvent->Fill(0);
  if((!mEvent->isTrigger(520001))
      &&(!mEvent->isTrigger(520011))
      &&(!mEvent->isTrigger(520021))
      &&(!mEvent->isTrigger(520031))
      &&(!mEvent->isTrigger(520041))
      &&(!mEvent->isTrigger(520051))
    )
    return false; 

  testEvent->Fill(1);

  // remove duplicate events
  ////////////////////////////////////
  if(mTempRunId==runID&&mTempEvtId==evtID)return false;
  mTempRunId=runID;  mTempEvtId=evtID;
  ////////////////////////////////////

  StThreeVectorF Vertex3D=mEvent->primaryVertex();
  const double VertexX = Vertex3D.x(); 
  const double VertexY = Vertex3D.y(); 
  const double VertexZ = Vertex3D.z(); 
  const double vpdVz   = mEvent->vzVpd();

  //event cut
  if(refMult <=2 || refMult > 1000) return false;
  testEvent->Fill(2);
  if(removeBadID(runID))return false;            
  testEvent->Fill(3);
  if(fabs(VertexZ) > 100) return false; 

  // hVertex2D ->Fill(VertexZ,vpdVz);
  // hDiffVz   ->Fill(VertexZ-vpdVz); 

  if(fabs(VertexZ) > 6) return false; 
  if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0)return false; 
  if(fabs(VertexZ-vpdVz)>3.)return false;       // no vpd cut in low energy?
  testEvent->Fill(4);

  //if(fabs(VertexZ) > mTreeCut::mVzMaxMap[mEnergy]) return kStOK; 
  //if(sqrt(pow(VertexX-mTreeCut::mVxMap[mEnergy],2.)+pow(VertexY-mTreeCut::mVyMap[mEnergy],2.))>mTreeCut::mVrMaxMap[mEnergy])return kStOK; 
  //if(mTreeCut::mVPDMap[mEnergy]&&fabs(VertexZ-vpdVz)>3.)return kStOK;       // no vpd cut in low energy?

  //check run number
  int runnumberPointer = -999;
  runnumberPointer=CheckrunNumber(runID);
  if(runnumberPointer == -999)return false;

  int dayPointer = (int)((runID)/1000%1000);
  int mRunL=mRunList.at(0);
  int mDayL=(int) ((mRunL)/1000%1000);
  dayPointer -= mDayL;
  // int timePointer = dayPointer/mStps;


  //StRefMultCorr
  if ( runID != mPrevRunId ) 
  {
    mRefMultCorr->init(runID);
    mPrevRunId = runID;
    cout << "reset mPrevRunId = " << mPrevRunId << endl;
  }
  // //mRefMultCorr->initEvent(refMult,VertexZ);
  mRefMultCorr->initEvent(refMult,VertexZ,mEvent->ZDCx());
  // mWght    = mRefMultCorr->getWeight();
  // double mult_corr= mRefMultCorr->getRefMultCorr() ;
  //
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  //
  if( centrality<0||centrality>=(nCent-1)) return false;
  testEvent->Fill(5);
  // int cent = centrality+1;  
  // //cout<<refMult<<" "<<cent<<" "<<mRefMultCorr->getCentralityBin16()<<endl;
  //
  // // careful cent 1-9    -nan  70-80, 60-70, 50-60, 40-50, 30-40, 20-30, 10-20, 5-10, 0-5 
  // double wCentSC[nCent]={-999, 1.,    1.,    1.,    1.,    1.,    1.,    0.4,   0.1,  0.1,};
  // //double wCentSC[nCent]={-999, 1.,    1.,    1.,    1.,    1.,    1.,    1.,   1.,  1.,};
  // //
  // if(wCentSC[cent]<1.){
  //   double mRand=gRandom->Rndm();
  //   if(mRand>wCentSC[cent]) return false;
  // }
  // testEvent->Fill(6);

  // mWght/=wCentSC[cent];
  mWght = 1;

  double mVz=6;
  double wVz=2.0*mVz/nVz;
  int    iVz=(VertexZ+mVz)/wVz;
  //cout<<wVz<<" "<<iVz<<endl;
  if(iVz<0||iVz>=nVz) return false;
  testEvent->Fill(7);

  return true;

}
bool StMyAnalysisMaker::readRunList()
{ 
  if(0) return kTRUE;
  // cout<<"Excuting readRunList function!"<<endl;
  TString inf=mListDir + "runList/";
  inf += Form("run%dList%.1f.list",mRun,mEnergy);
  ifstream inrun; 
  cout<<"file is located at: "<<inf<<endl;
  inrun.open(inf);
  if ( inrun.fail() ) {
    cout<< "cannot open " << inf.Data() << endl;
    return kFALSE;
  }
  Int_t runid;
  while ( inrun >> runid ) { mRunList.push_back(runid); }
  inrun.close();
  sort(mRunList.begin(),mRunList.end());

  vector<int>::iterator it;
  it = std::unique (mRunList.begin(), mRunList.end());
  mRunList.resize( std::distance(mRunList.begin(),it) );

  cout <<"Run list :" <<inf.Data() << " loaded. "<< endl;
  cout <<"Total    :" <<mRunList.size()<< " runs. "<< endl;

  if(mRunList.size()<1){cout<<"no run number found!!!"<<endl; return kFALSE;}

  return kTRUE;
}
//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::CheckrunNumber(int runnumber) const
{    
  int pointer=-999; 
  int id=0; 
  for (auto it=mRunList.begin(); it!=mRunList.end(); ++it) 
  { 
    if(runnumber==*it)pointer=id;
    id++;
  }

  if(pointer==-999)cout<<"Run number are not found! "<<runnumber<<endl;
  return pointer;
} 
//-----------------------------------------------------------------------------
// void StMyAnalysisMaker::setRunEnergyAndListDir(int run, double energy,char ListDir[256])
// {
// 	mRun    =run;
// 	mEnergy =energy;
// 	mListDir=ListDir;
// }

