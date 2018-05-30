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
#include "StCuts.cxx"

#include "StBTofUtil/tofPathLength.hh"
#include "TRandom3.h"
#include "StPhysicalHelixD.hh"

#include "TRMatrix.h"
#include "TRSymMatrix.h"
#include "TMinuit.h"

#include "TMath.h"
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
  cout<<"load picoDst"<<endl;
  cout<<"load picoEvent"<<endl;
  mOutName = outName;
  mRun    = 16;
  mEnergy =200;
  mListDir="./";
  cout<<"load others"<<endl;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();
  PI   = TMath::Pi();
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

  d0MassPtY = new TH3D("d0MassPtY",";D^{0} mass (GeV/c^{2});p_{T} (GeV/c)",50,1.6,2.1,100,0,10,4,-0.8,0.8);
  d0BarMassPtY = new TH3D("d0BarMassPtY",";D^{0} mass (GeV/c^{2});p_{T} (GeV/c)",50,1.6,2.1,100,0,10,4,-0.8,0.8);
  d0MassPhiEta = new TH3D("d0MassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta = new TH3D("d0BarMassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_noweight = new TH3D("d0MassPhiEta_noweight",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_noweight = new TH3D("d0BarMassPhiEta_noweight",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_50 = new TH3D("d0MassPhiEta_50",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_50 = new TH3D("d0BarMassPhiEta_50",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_150 = new TH3D("d0MassPhiEta_150",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_150 = new TH3D("d0BarMassPhiEta_150",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_pt2 = new TH3D("d0MassPhiEta_pt2",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_pt2 = new TH3D("d0BarMassPhiEta_pt2",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_pt25 = new TH3D("d0MassPhiEta_pt25",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_pt25 = new TH3D("d0BarMassPhiEta_pt25",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_pt3 = new TH3D("d0MassPhiEta_pt3",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_pt3 = new TH3D("d0BarMassPhiEta_pt3",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_5bin = new TH3D("d0MassPhiEta_5bin",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_5bin = new TH3D("d0BarMassPhiEta_5bin",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
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
  d0MassPhiEta->Sumw2();
  d0BarMassPhiEta->Sumw2();
  d0MassPhiEta_noweight->Sumw2();
  d0BarMassPhiEta_noweight->Sumw2();
  d0MassPhiEta_50->Sumw2();
  d0BarMassPhiEta_50->Sumw2();
  d0MassPhiEta_150->Sumw2();
  d0BarMassPhiEta_150->Sumw2();
  d0MassPhiEta_pt2->Sumw2();
  d0BarMassPhiEta_pt2->Sumw2();
  d0MassPhiEta_pt25->Sumw2();
  d0BarMassPhiEta_pt25->Sumw2();
  d0MassPhiEta_pt3->Sumw2();
  d0BarMassPhiEta_pt3->Sumw2();
  d0MassPhiEta_5bin->Sumw2();
  d0BarMassPhiEta_5bin->Sumw2();
  TFile *f_eff = new TFile("eff_simu.root");

  //X: pt (0,10,100), Y: eta(-1,1,20), Z: phi,(0,2pi,10);
  hD0Mc[9] = (TH3F *)f_eff->Get(Form("h3McD0PtYPhi_%i",0))->Clone(Form("hD0Mc_%i",9));
  hD0Rc[9] = (TH3F *)f_eff->Get(Form("h3RcD0PtYPhi_%i",0))->Clone(Form("hD0Rc_%i",9));
  hD0BarMc[9] = (TH3F *)f_eff->Get(Form("h3McD0BarPtYPhi_%i",0))->Clone(Form("hD0BarMc_%i",9));
  hD0BarRc[9] = (TH3F *)f_eff->Get(Form("h3RcD0BarPtYPhi_%i",0))->Clone(Form("hD0BarRc_%i",9));

  for(int i =0; i< 9;i++)
  {
    hD0Mc[i] = (TH3F *)f_eff->Get(Form("h3McD0PtYPhi_%i",i))->Clone(Form("hD0Mc_%i",i));
    hD0Rc[i] = (TH3F *)f_eff->Get(Form("h3RcD0PtYPhi_%i",i))->Clone(Form("hD0Rc_%i",i));
    hD0BarMc[i] = (TH3F *)f_eff->Get(Form("h3McD0BarPtYPhi_%i",i))->Clone(Form("hD0BarMc_%i",i));
    hD0BarRc[i] = (TH3F *)f_eff->Get(Form("h3RcD0BarPtYPhi_%i",i))->Clone(Form("hD0BarRc_%i",i));

    // X: pt , Y: Y
    hD0McPtY[i] = (TH2F *)hD0Mc[i]->Project3D("yx")->Clone(Form("hD0McPtY_%i",i));
    hD0McPtY[i]->RebinY(2);
    hD0BarMcPtY[i] = (TH2F *)hD0BarMc[i]->Project3D("yx")->Clone(Form("hD0BarMcPtY_%i",i));
    hD0BarMcPtY[i]->RebinY(2);
    hD0RcPtY[i] = (TH2F *)hD0Rc[i]->Project3D("yx")->Clone(Form("hD0RcPtY_%i",i));
    hD0RcPtY[i]->RebinY(2);
    hD0BarRcPtY[i] = (TH2F *)hD0BarRc[i]->Project3D("yx")->Clone(Form("hD0BarRcPtY_%i",i));
    hD0BarRcPtY[i]->RebinY(2);

    if(i!=0)
    {
      hD0Mc[9]->Add(hD0Mc[i]);
      hD0Rc[9]->Add(hD0Rc[i]);
      hD0BarMc[9]->Add(hD0BarMc[i]);
      hD0BarRc[9]->Add(hD0BarRc[i]);
    }
    D0Mc[i][0] = (TH1F *)hD0Mc[i]->ProjectionX(Form("D0McPt_%i",i))->Clone(Form("D0McPt_%i",i));
    D0Mc[i][1] = (TH1F *)hD0Mc[i]->ProjectionY(Form("D0McY_%i",i))->Clone(Form("D0McY_%i",i));
    D0Mc[i][2] = (TH1F *)hD0Mc[i]->ProjectionZ(Form("D0McPhi_%i",i))->Clone(Form("D0McPhi_%i",i));
    D0BarMc[i][0] = (TH1F *)hD0BarMc[i]->ProjectionX(Form("D0BarMcPt_%i",i))->Clone(Form("D0BarMcPt_%i",i));
    D0BarMc[i][1] = (TH1F *)hD0BarMc[i]->ProjectionY(Form("D0BarMcY_%i",i))->Clone(Form("D0BarMcY_%i",i));
    D0BarMc[i][2] = (TH1F *)hD0BarMc[i]->ProjectionZ(Form("D0BarMcPhi_%i",i))->Clone(Form("D0BarMcPhi_%i",i));
    D0Rc[i][0] = (TH1F *)hD0Rc[i]->ProjectionX(Form("D0RcPt_%i",i))->Clone(Form("D0RcPt_%i",i));
    D0Rc[i][1] = (TH1F *)hD0Rc[i]->ProjectionY(Form("D0RcY_%i",i))->Clone(Form("D0RcY_%i",i));
    D0Rc[i][2] = (TH1F *)hD0Rc[i]->ProjectionZ(Form("D0RcPhi_%i",i))->Clone(Form("D0RcPhi_%i",i));
    D0BarRc[i][0] = (TH1F *)hD0BarRc[i]->ProjectionX(Form("D0BarRcPt_%i",i))->Clone(Form("D0BarRcPt_%i",i));
    D0BarRc[i][1] = (TH1F *)hD0BarRc[i]->ProjectionY(Form("D0BarRcY_%i",i))->Clone(Form("D0BarRcY_%i",i));
    D0BarRc[i][2] = (TH1F *)hD0BarRc[i]->ProjectionZ(Form("D0BarRcPhi_%i",i))->Clone(Form("D0BarRcPhi_%i",i));
    for(int j = 0;j<3;j++)
    {
      D0Rc[i][j]->Divide(D0Mc[i][j]);
      D0BarRc[i][j]->Divide(D0Mc[i][j]);
    }
  }
  cout<<"Initialization done"<<endl;
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
  d0MassPhiEta_noweight->Write();
  d0BarMassPhiEta_noweight->Write();
  d0MassPhiEta_50->Write();
  d0BarMassPhiEta_50->Write();
  d0MassPhiEta_150->Write();
  d0BarMassPhiEta_150->Write();
  d0MassPhiEta_pt2->Write();
  d0BarMassPhiEta_pt2->Write();
  d0MassPhiEta_pt25->Write();
  d0BarMassPhiEta_pt25->Write();
  d0MassPhiEta_pt3->Write();
  d0BarMassPhiEta_pt3->Write();
  d0MassPhiEta_5bin->Write();
  d0BarMassPhiEta_5bin->Write();
  testDPhi->Write();
  testEvent->Write();
  d0MassPtY->Write();
  d0BarMassPtY->Write();
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
  mPicoEvent = (StPicoEvent *)mPicoDst->event();
  if(!mPicoEvent)
  {
    cout<<"no picoEvent at all"<<endl;
    LOG_WARN<<"no picoEvent at all"<<endl;
    return kStWarn;
  }
  // StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  mVtx = mPicoEvent->primaryVertex();
  mBField = mPicoEvent->bField();
  const int  runID    = mPicoEvent->runId();
  // const int  evtID    = event->eventId();
  // const int refMult  = event->grefMult();
  if(!(isGoodEvent()))
  { 
    return kStWarn;
  }

  // cout<<"this is a good event!"<<endl;
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  // StThreeVectorF vtx = mPicoEvent->primaryVertex();
  // float b = mPicoEvent->bField();
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  if(centrality == 7 || centrality == 8)
  {
    LOG_WARN << "centrality in 0-10%"<<endl;
    return kStWarn;
  }
  double reweight= mRefMultCorr->getWeight();
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
    ZDCSMDadc[i]   = 1.* mPicoEvent->ZdcSmdEastHorizontal(i);   // picoDst function i 0-7
    ZDCSMDadc[i+8] = 1.* mPicoEvent->ZdcSmdEastVertical(i);
    ZDCSMDadc[i+16]= 1.* mPicoEvent->ZdcSmdWestHorizontal(i);
    ZDCSMDadc[i+24]= 1.* mPicoEvent->ZdcSmdWestVertical(i);
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

  zdcPsi->Fill(mZDC1Event_PsiOrigin,reweight);
  zdcPsi_corr->Fill(mZDC1Event_PsiF,reweight);
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
    StThreeVectorF mom= itrk->gMom(mVtx,mBField);
    double eta  = mom.pseudoRapidity();
    double phi  = mom.phi();
    double dca = (mVtx - itrk->dcaPoint()).mag();
    trackPhiEta->Fill(phi,eta);
    // if (isTpcPion(itrk) && fabs(dca)<1.5) 
    //   testPionIndex.push_back(i);
    trackPhiEtaHFT->Fill(phi,eta);

    float kBeta = getTofBeta(itrk);
    bool tofAvailable = kBeta>0;

    bool tpcKaon = isTpcKaon(itrk);
    bool tofKaon = tofAvailable && isTofKaon(itrk,kBeta);
    bool goodKaon = tofAvailable? tofKaon&&tpcKaon : tpcKaon;
    if(goodKaon) 
      kaonIndex.push_back(i);

    bool tpcPion = isTpcPion(itrk);
    bool tofPion= tofAvailable && isTofPion(itrk,kBeta);
    bool goodPion = tofAvailable? tofPion&&tpcPion : tpcPion;
    if (goodPion) 
      pionIndex.push_back(i);
  }
  for(unsigned int i=0;i<pionIndex.size();i++)
  {
    StPicoTrack const* pion= mPicoDst->track(pionIndex[i]);
    for(unsigned int j=0;j<kaonIndex.size();j++)
    {
      StPicoTrack const* kaon= mPicoDst->track(kaonIndex[j]);

      int charge = pion->charge() * kaon->charge();
      if(charge>0) continue;//only unlike-sign pairs
      StKaonPion *kp = new StKaonPion(kaon,pion,j,i,mVtx,mBField);
      if(kp->pt()<1.5) continue;// require D0 pT > 1.5 GeV/c
      // if(kp->pt()<1.) continue;// require D0 pT > 1.5 GeV/c
      // StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
      // StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());

      int isD0 = isD0Pair(kp);
      int isD050 = isD0Pair50(kp);
      int isD0150 = isD0Pair150(kp);
      if(isD0==1)
      {
        // cout<<"find D0 !"<<endl;
        TLorentzVector d0Lorentz;
        bool isD0Meson = kaon->charge() < 0;
        d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
        float d0_eff = getD0Efficiency(d0Lorentz,isD0Meson);
        if(d0_eff<=0)
        {
          cout<<"efficiency = 0"<<endl;
          continue;
        }
        cout<<"efficiency = "<<d0_eff<<endl;
        double kpY = d0Lorentz.Rapidity();
        double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);

        if(deltaPhi > PI)
          deltaPhi = twoPI - deltaPhi;
        deltaPhi = fabs(deltaPhi);
        testDPhi->Fill(deltaPhi);
        if(kaon->charge() < 0)
        {
          d0MassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          if(kp->pt() >= 2)
            d0MassPhiEta_pt2->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          if(kp->pt() >= 2.5)
            d0MassPhiEta_pt25->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          if(kp->pt() >= 3)
            d0MassPhiEta_pt3->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);

          d0MassPhiEta_noweight->Fill(kp->m(),deltaPhi,kpY);
          d0MassPhiEta_5bin->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          d0MassPtY->Fill(kp->m(),kp->pt(),kpY,1.0/d0_eff*reweight);
        }
        if(kaon->charge() > 0)
        {
          if(kp->pt() >= 2)
            d0BarMassPhiEta_pt2->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          if(kp->pt() >= 2.5)
            d0BarMassPhiEta_pt25->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          if(kp->pt() >= 3)
            d0BarMassPhiEta_pt3->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          d0BarMassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
          d0BarMassPhiEta_noweight->Fill(kp->m(),deltaPhi,kpY);
          d0BarMassPtY->Fill(kp->m(),kp->pt(),kpY,1.0/d0_eff*reweight);
          d0BarMassPhiEta_5bin->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
        }
      }
      if(isD050==1)
      {
        // cout<<"find D0 !"<<endl;
        TLorentzVector d0Lorentz;
        bool isD0Meson = kaon->charge() < 0;
        d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
        float d0_eff = getD0Efficiency(d0Lorentz,isD0Meson);
        double kpY = d0Lorentz.Rapidity();
        double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);

        if(deltaPhi > PI)
          deltaPhi = twoPI - deltaPhi;
        deltaPhi = fabs(deltaPhi);
        testDPhi->Fill(deltaPhi);
        if(kaon->charge() < 0)
          d0MassPhiEta_50->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
        if(kaon->charge() > 0)
          d0BarMassPhiEta_50->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
      }
      if(isD0150==1)
      {
        // cout<<"find D0 !"<<endl;
        TLorentzVector d0Lorentz;
        bool isD0Meson = kaon->charge() < 0;
        d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
        float d0_eff = getD0Efficiency(d0Lorentz,isD0Meson);
        double kpY = d0Lorentz.Rapidity();
        double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);

        if(deltaPhi > PI)
          deltaPhi = twoPI - deltaPhi;
        deltaPhi = fabs(deltaPhi);
        testDPhi->Fill(deltaPhi);
        if(kaon->charge() < 0)
          d0MassPhiEta_150->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
        if(kaon->charge() > 0)
          d0BarMassPhiEta_150->Fill(kp->m(),deltaPhi,kpY,1.0/d0_eff*reweight);
      }
      delete kp;
    }//finish D0 loop
  }


  return kStOK;
}


bool StMyAnalysisMaker::isTofPion(StPicoTrack const * const trk, float beta) const
{
  bool tofPion = false;
  // StPicoEvent *mPicoEvent = (StPicoEvent *)mPicoDst->event();
  // float b = mPicoEvent->bField();
  // StThreeVectorF vtx = mPicoEvent->primaryVertex();
  if(beta>0)
  {
    // double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gMom(mVtx,mBField).mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
    tofPion = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofPion;
}



bool StMyAnalysisMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;
  // StPicoEvent *mPicoEvent = (StPicoEvent *)mPicoDst->event();
  // float b = mPicoEvent->bField();
  // StThreeVectorF vtx = mPicoEvent->primaryVertex();
  if(beta>0)
  {
    // double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gMom(mVtx,mBField).mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}
float StMyAnalysisMaker::getTofBeta(StPicoTrack const * const trk) const
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

        float L = tofPathLength(&mVtx, &btofHitPos, helix.curvature());
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
bool StMyAnalysisMaker::isTpcKaon(StPicoTrack const * const trk) const
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
  // return trk->gPt() > 0.6 && trk->nHitsFit() >= 20 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 && trk->nHitsDedx() >= 16;
  if (!trk->isHft())
    return false;
  // StPicoEvent *mPicoEvent = (StPicoEvent *)mPicoDst->event();
  // StThreeVectorF vtx = mPicoEvent->primaryVertex();
  // float b = mPicoEvent->bField();
  StThreeVectorF gmom = trk->gMom(mVtx,mBField);
  StPhysicalHelixD helix = trk->helix(mBField);
  return gmom.perp() > 0.6 && trk->nHitsFit() >= 20 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 && fabs(gmom.pseudoRapidity()) < 1 && fabs(helix.geometricSignedDistance(mVtx)) > 0.005;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
////////////////////////////////////////////////////////////////////////////////////
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
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  const int indexCent = getCentIndex(centrality);
  const int indexPt = getPtIndex(kp->pt());
  bool pairCuts = false;
  if(indexCent<0) return false;
  if(indexPt<0) return false;
  pairCuts = ( cos(kp->pointingAngle())> mycuts::cosTheta
      && kp->kaonDca() > mycuts::kDca[indexCent][indexPt]
      && kp->pionDca() > mycuts::pDca[indexCent][indexPt]
      && kp->dcaDaughters() < mycuts::dcaDaughters[indexCent][indexPt]
      && sin(kp->pointingAngle())*kp->decayLength() < mycuts::dcaV0ToPv[indexCent][indexPt]
      && kp->decayLength() > mycuts::decayLength[indexCent][indexPt]
      );

  
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
int StMyAnalysisMaker::isD0Pair50(StKaonPion const* const kp) const
{
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  const int indexCent = getCentIndex(centrality);
  const int indexPt = getPtIndex(kp->pt());
  bool pairCuts = false;
  if(indexCent<0) return false;
  if(indexPt<0) return false;
  pairCuts = ( cos(kp->pointingAngle())> mycuts::cosTheta
      && kp->kaonDca() > mycuts::kDca1[indexCent][indexPt]
      && kp->pionDca() > mycuts::pDca1[indexCent][indexPt]
      && kp->dcaDaughters() < mycuts::dcaDaughters1[indexCent][indexPt]
      && sin(kp->pointingAngle())*kp->decayLength() < mycuts::dcaV0ToPv1[indexCent][indexPt]
      && kp->decayLength() > mycuts::decayLength1[indexCent][indexPt]
      );

  
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
int StMyAnalysisMaker::isD0Pair150(StKaonPion const* const kp) const
{
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  const int indexCent = getCentIndex(centrality);
  const int indexPt = getPtIndex(kp->pt());
  bool pairCuts = false;
  if(indexCent<0) return false;
  if(indexPt<0) return false;
  pairCuts = ( cos(kp->pointingAngle())> mycuts::cosTheta
      && kp->kaonDca() > mycuts::kDca2[indexCent][indexPt]
      && kp->pionDca() > mycuts::pDca2[indexCent][indexPt]
      && kp->dcaDaughters() < mycuts::dcaDaughters2[indexCent][indexPt]
      && sin(kp->pointingAngle())*kp->decayLength() < mycuts::dcaV0ToPv2[indexCent][indexPt]
      && kp->decayLength() > mycuts::decayLength2[indexCent][indexPt]
      );

  
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
  // TString inf=mListDir + "/badList/";
  // inf += Form("badrun%dList%.1f.list",mRun,mEnergy);
  TString inf = mListDir + "goodruns2016_wpxl_grefhftrcut.txt";
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
bool StMyAnalysisMaker::isGoodEvent()
{	

  const int  runID    = mPicoEvent->runId();
  const int  evtID    = mPicoEvent->eventId();
  const int  refMult  = mPicoEvent->grefMult();
  // const int grefMult  = mPicoEvent->grefMult();
  // const int  ranking  = mPicoEvent->ranking();

  testEvent->Fill(0);
  if((!mPicoEvent->isTrigger(520001))
      &&(!mPicoEvent->isTrigger(520011))
      &&(!mPicoEvent->isTrigger(520021))
      &&(!mPicoEvent->isTrigger(520031))
      &&(!mPicoEvent->isTrigger(520041))
      &&(!mPicoEvent->isTrigger(520051))
      &&(!mPicoEvent->isTrigger(570001))//VPDMB-5-sst for production2
    )
    return false; 

  testEvent->Fill(1);

  // remove duplicate events
  ////////////////////////////////////
  if(mTempRunId==runID&&mTempEvtId==evtID)return false;
  mTempRunId=runID;  mTempEvtId=evtID;
  ////////////////////////////////////

  StThreeVectorF Vertex3D=mPicoEvent->primaryVertex();
  const double VertexX = Vertex3D.x(); 
  const double VertexY = Vertex3D.y(); 
  const double VertexZ = Vertex3D.z(); 
  const double vpdVz   = mPicoEvent->vzVpd();

  //event cut
  if(refMult <=2 || refMult > 1000) return false;
  testEvent->Fill(2);
  auto it = find (mBadList.begin(),mBadList.end(), runID);
  if (it == mBadList.end())
  {
      return false;
  }
  // if(removeBadID(runID))return false;            
  testEvent->Fill(3);
  if(fabs(VertexZ) > 100) return false; 

  // hVertex2D ->Fill(VertexZ,vpdVz);
  // hDiffVz   ->Fill(VertexZ-vpdVz); 

  if(fabs(VertexZ) > 6) return false; 
  if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0)return false; 
  if(fabs(VertexZ-vpdVz)>3.)return false;       // no vpd cut in low energy?
  if(fabs(VertexX)<1e-5 && fabs(VertexY)<1e-5 && fabs(VertexZ)<1e-5) return false;
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
  mRefMultCorr->initEvent(refMult,VertexZ,mPicoEvent->ZDCx());
  // reweight    = mRefMultCorr->getWeight();
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

  // reweight/=wCentSC[cent];

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

float StMyAnalysisMaker::getD0Efficiency(TLorentzVector &d0Lorentz,bool isD0)
{
  int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  int ptBin = (int)(d0Lorentz.Pt()/0.1)+1;
  int ptBinRef = (int)(5./0.1)+1;
  if(d0Lorentz.Pt()>10)  ptBin = 100;
  int yBin = (int)((d0Lorentz.Rapidity()+1)/0.1)+1;
  int yBin_new = (int)((d0Lorentz.Rapidity()+1)/0.2)+1;
  int phiBin = (int)((d0Lorentz.Phi()+TMath::Pi())/(0.2*TMath::Pi()))+1;
  float mcYield,rcYield;
  float mcYieldRef,rcYieldRef;
  float mcYieldRef_err,rcYieldRef_err;
  float mcYield_err,rcYield_err;
  double eff_pt,eff_phi,eff_y;
  double eff_pt_err;
  if (isD0)
  {
    mcYield = hD0Mc[centrality]->GetBinContent(ptBin,yBin,phiBin);
    rcYield = hD0Rc[centrality]->GetBinContent(ptBin,yBin,phiBin);
    mcYield_err = hD0Mc[centrality]->GetBinError(ptBin,yBin,phiBin);
    rcYield_err = hD0Rc[centrality]->GetBinError(ptBin,yBin,phiBin);

    mcYieldRef = hD0McPtY[centrality]->GetBinContent(ptBin,yBin_new);
    rcYieldRef = hD0RcPtY[centrality]->GetBinContent(ptBin,yBin_new);
    mcYieldRef_err = hD0McPtY[centrality]->GetBinError(ptBin,yBin_new);
    rcYieldRef_err = hD0RcPtY[centrality]->GetBinError(ptBin,yBin_new);

    eff_pt = D0Rc[centrality][0]->GetBinContent(ptBin);
    eff_pt_err = D0Rc[centrality][0]->GetBinError(ptBin);
    eff_y = D0Rc[centrality][1]->GetBinContent(yBin);
    eff_phi = D0Rc[centrality][2]->GetBinContent(phiBin);
  }
  else
  {
    mcYield = hD0BarMc[centrality]->GetBinContent(ptBin,yBin,phiBin);
    rcYield = hD0BarRc[centrality]->GetBinContent(ptBin,yBin,phiBin);
    mcYield_err = hD0BarMc[centrality]->GetBinError(ptBin,yBin,phiBin);
    rcYield_err = hD0BarRc[centrality]->GetBinError(ptBin,yBin,phiBin);

    mcYieldRef = hD0BarMcPtY[centrality]->GetBinContent(ptBin,yBin_new);
    rcYieldRef = hD0BarRcPtY[centrality]->GetBinContent(ptBin,yBin_new);
    mcYieldRef_err = hD0BarMcPtY[centrality]->GetBinError(ptBin,yBin_new);
    rcYieldRef_err = hD0BarRcPtY[centrality]->GetBinError(ptBin,yBin_new);

    eff_pt = D0BarRc[centrality][0]->GetBinContent(ptBin);
    eff_pt_err = D0BarRc[centrality][0]->GetBinError(ptBin);
    eff_y = D0BarRc[centrality][1]->GetBinContent(yBin);
    eff_phi = D0BarRc[centrality][2]->GetBinContent(phiBin);
  }
  // cout<<"pt bin = "<<ptBin<<"\ty bin = "<<yBin<<"\tphi bin = "<<phiBin<<endl;
  // cout<<"mc = "<<mcYield<<"\trc = "<<rcYield<<endl;
  double eff = rcYield/mcYield;
  double eff_err = sqrt(pow(rcYield_err/mcYield,2)+pow(rcYield*mcYield_err/mcYield/mcYield,2));
  double effRef = rcYieldRef/mcYieldRef;
  double effRef_err = sqrt(pow(rcYieldRef_err/mcYieldRef,2)+pow(rcYieldRef*mcYieldRef_err/mcYieldRef/mcYieldRef,2));
  // cout<<"############ efficiency = "<<eff<<" +/- "<<eff_err<<endl;
  // cout<<"############ relative error = "<<eff_err/eff<<endl;
  // cout<<"############ pt efficiency = "<<eff_pt<<" +/- "<<eff_pt_err<<endl;
  // cout<<"############ relative error = "<<eff_pt_err/eff_pt<<endl;
  // cout<<"############ ref efficiency = "<<effRef<<" +/- "<<effRef_err<<endl;
  // cout<<"############ relative error = "<<effRef_err/effRef<<endl;
  // double eff = eff_pt*eff_y*eff_phi/0.005/0.005/0.05;
  return effRef;
}
// get pt bin
int StMyAnalysisMaker::getPtIndex(float const pt) const 
{
    int bin = -1;
    for (int i = 0; i < mycuts::nPtBins; i++)
    {
        if ((pt >= mycuts::PtEdge[i]) && (pt < mycuts::PtEdge[i + 1]))
            bin = i;
    }
    return bin;
}

// get cent bin
int StMyAnalysisMaker::getCentIndex(int const cent) const 
{
    int bin = -1;
    for (int i = 0; i < mycuts::nCent; i++)
    {
        if ((cent > mycuts::CentEdge[i]) && (cent < mycuts::CentEdge[i + 1]))
            bin = i;
    }
    return bin;
}
