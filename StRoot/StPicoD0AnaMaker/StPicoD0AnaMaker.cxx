#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
#include <stdio.h>
#include <time.h>
#include <algorithm>

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mOutputFile->cd();


  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");
  miniZDCSMD=new mZDCSMD();
  miniZDCSMD->SetFileDirectory("/global/homes/a/amcw7777/auau200GeVD0V1/ZDCSMDFile");
  // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6 
  if(!(miniZDCSMD->InitRun(5)))return kStFatal;
  miniZDCSMD->SetmHistFileName(mOutFileName+"zdc");

  PI   = 3.14159;
  twoPI= 2.*PI;
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
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  fout1.close();
  mOutputFile->cd();
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
  miniZDCSMD->WriteHist();

  // save user variables here
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();


  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  // if(!(isGoodEvent()) || !event->isMinBias())//minBias trigger requires
  if(!(isGoodEvent(event)))//minBias trigger requires
  {
    LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  if(event) {
    pVtx = event->primaryVertex();
  }
  ////////////  ZDCSMD ////////////////
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int cent = centrality+1;
  if(centrality<0) {
    LOG_WARN << "not minBias sample!" << endl;
    return kStWarn;
  }
  // if(!(centrality==4||centrality==5||centrality==6))//minBias trigger requires
  // {
  //   LOG_WARN << " Not Good Event! Skip! " << endm;
  //   return kStWarn;
  // }
  const int  runID    = 17109018;//event->runId();
  float ZDCSMDadc[32];
  for(int i=0;i<8;i++){
    ZDCSMDadc[i]   =  event->ZdcSmdEastHorizontal(i);   // picoDst function i 0-7
    ZDCSMDadc[i+8] =  event->ZdcSmdEastVertical(i);
    ZDCSMDadc[i+16]=  event->ZdcSmdWestHorizontal(i);
    ZDCSMDadc[i+24]=  event->ZdcSmdWestVertical(i);
  }
  // cout<<"zdcsmd init."<<endl;
  miniZDCSMD->InitEvent();
  // cout<<"zdcsmd set."<<endl;
  miniZDCSMD->SetZDCSMDcent(ZDCSMDadc,cent,runID);
  // miniZDCSMD->SetZDCSMDcent(ZDCSMDadc,4,17109018);
  // cout<<"zdcsmd calibrate."<<endl;
  miniZDCSMD->calibrateZDCSMDevp();
  if(!(miniZDCSMD->ZDCSMDgoodEvent()))return kStOK;

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

  zdcPsi->Fill(mZDC1Event_PsiOrigin);
  zdcPsi_corr->Fill(mZDC1Event_PsiF);
	zdcResolution->Fill(cent,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 
	zdcResolution->Fill(0.,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 
////////////////////////////////////////////////////
  // cout<<"zdcsmd done."<<endl;
  double reweight = mGRefMultCorrUtil->getWeight();
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
    if (!kaon->isHFTTrack() || !pion->isHFTTrack()) continue;
    if (!isTpcPion(pion)) continue;
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    int charge=0;
    float d0Pt = kp->pt();
    if(d0Pt < 1.5)  continue;
    // double dMass = kp->m();
    // double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    // double reweight_eff = 1;//= (efficiency[0][fitindex]/efficiency[centBin][fitindex]);

    if((charge=isD0Pair(kp))!=0 )
    {
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      d0MassPt->Fill(kp->m(),kp->pt());
      double kpY = d0Lorentz.Rapidity();
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;
      deltaPhi = fabs(deltaPhi);
      testDPhi->Fill(deltaPhi);
      // cout<<"deltaPhi = "<<deltaPhi<<endl;
      if(kaon->charge() < 0)
        d0MassPhiEta->Fill(kp->m(),deltaPhi,kpY);
      if(kaon->charge() > 0)
        d0BarMassPhiEta->Fill(kp->m(),deltaPhi,kpY);

    }//D loop

  }

  return kStOK;
}
//-----------------------------------------------------------------------------

int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
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

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/

int StPicoD0AnaMaker::isD0Pair50(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0044 &&
      kp->pionDca() > 0.0120 && kp->kaonDca() > 0.0119 &&
      kp->dcaDaughters() < 0.0069 && kp->decayLength()>0.0144;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0036 &&
      kp->pionDca() > 0.0102 && kp->kaonDca() > 0.0110 &&
      kp->dcaDaughters() < 0.0048 && kp->decayLength()>0.0204;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0031 &&
      kp->pionDca() > 0.0118 && kp->kaonDca() > 0.0109 &&
      kp->dcaDaughters() < 0.0044 && kp->decayLength()>0.0242;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0026 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0106 &&
      kp->dcaDaughters() < 0.0049 && kp->decayLength()>0.0245;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0032 &&
      kp->pionDca() > 0.0096 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0047 && kp->decayLength()>0.0300;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

int StPicoD0AnaMaker::isD0Pair150(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0072 &&
      kp->pionDca() > 0.0092 && kp->kaonDca() > 0.0105 &&
      kp->dcaDaughters() < 0.0077 && kp->decayLength()>0.0110;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0053 &&
      kp->pionDca() > 0.0078 && kp->kaonDca() > 0.0068 &&
      kp->dcaDaughters() < 0.0078 && kp->decayLength()>0.0168;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0074 && kp->decayLength()>0.0187;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.0065 && kp->kaonDca() > 0.0066 &&
      kp->dcaDaughters() < 0.0068 && kp->decayLength()>0.0199;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0047 && kp->kaonDca() > 0.0041 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0180;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
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
        StPhysicalHelixD helix = trk->helix();

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
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    double ptot = trk->dcaGeometry().momentum().mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}

//------------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodEvent(StPicoEvent const*mEvent)
{	

  const int  runID    = mEvent->runId();
  const int  evtID    = mEvent->eventId();
  const int  refMult  = mEvent->grefMult();
  // const int grefMult  = mEvent->grefMult();
  // const int  ranking  = mEvent->ranking();

  testEvent->Fill(0);
  if(!mEvent->isTrigger(450050) &&
      !mEvent->isTrigger(450060)&&
      !mEvent->isTrigger(450005)&&
      !mEvent->isTrigger(450015)&&
      !mEvent->isTrigger(450025))
    return false; 

  testEvent->Fill(1);


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
  testEvent->Fill(4);
  // mRefMultCorr->initEvent(refMult,VertexZ,mEvent->ZDCx());
  // //
  // int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
  // //
  // // if( centrality<0||centrality>=(nCent-1)) return false;
  // testEvent->Fill(5);
  // mWght = 1;
  //
  // double mVz=6;
  // double wVz=2.0*mVz/nVz;
  // int    iVz=(VertexZ+mVz)/wVz;
  // //cout<<wVz<<" "<<iVz<<endl;
  // if(iVz<0||iVz>=nVz) return false;
  // testEvent->Fill(7);

  return true;

}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= 15 &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}

bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > 0.6 && trk->nHitsFit() >= 20 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 && trk->nHitsDedx() >= 16;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}