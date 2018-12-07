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

  ifstream ifs("efficiency.txt");
  for(int i=0; i<6; i++)
    for(int j=0; j<4; j++)
      ifs>>efficiency[j][i];

  // -------------- USER VARIABLES -------------------------
  // mGRefMultCorrUtil = new StRefMultCorr("grefmult");
  miniZDCSMD=new mZDCSMD();
  miniZDCSMD->SetFileDirectory("/global/homes/a/amcw7777/auau200GeVD0V1/ZDCSMDFile");
  // pedgain 1,  beam center 2,  phiweight + shift (subevent) 3,  Fullevent 4,  ready+wQAHists 5,  ready+woQAHists 6 
  if(!(miniZDCSMD->InitRun(5)))return kStFatal;
  miniZDCSMD->SetmHistFileName(mOutFileName+"zdc");

  PI   = 3.14159;
  twoPI= 2.*PI;
  d0MassPt = new TH3D("d0MassPt",";D^{0} mass (GeV/c^{2});p_{T} (GeV/c);y",50,1.6,2.1,100,0,10,4,-0.8,0.8);
  d0BarMassPt = new TH3D("d0BarMassPt",";D^{0} mass (GeV/c^{2});p_{T} (GeV/c);y",50,1.6,2.1,100,0,10,4,-0.8,0.8);
  d0MassPhiEta = new TH3D("d0MassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta = new TH3D("d0BarMassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_pt3 = new TH3D("d0MassPhiEta_pt3",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_pt3 = new TH3D("d0BarMassPhiEta_pt3",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_50 = new TH3D("d0MassPhiEta_50",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_50 = new TH3D("d0BarMassPhiEta_50",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_150 = new TH3D("d0MassPhiEta_150",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_150 = new TH3D("d0BarMassPhiEta_150",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_5bin = new TH3D("d0MassPhiEta_5bin",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_5bin = new TH3D("d0BarMassPhiEta_5bin",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
  d0MassPhiEta_noeff = new TH3D("d0MassPhiEta_noweight",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_noeff = new TH3D("d0BarMassPhiEta_noweight",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,5,0,PI,4,-0.8,0.8);
  d0MassPhiEta_nHits25 = new TH3D("d0MassPhiEta_nHits25",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_nHits25 = new TH3D("d0BarMassPhiEta_nHits25",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0MassPhiEta_noHitsRatio = new TH3D("d0MassPhiEta_noHitsRatio",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  d0BarMassPhiEta_noHitsRatio = new TH3D("d0BarMassPhiEta_noHitsRatio",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",50,1.6,2.1,4,0,PI,4,-0.8,0.8);
  testDPhi = new TH1D("testDPhi","testDPhi",400,-2.*twoPI,2.*twoPI);
  testEvent = new TH1D("testEvent","",10,-0.5,9.5);
  zdcPsi = new TH1D("zdcPsi",";#psi_{ZDC}",1000,0,twoPI);
  zdcPsi_corr = new TH1D("zdcPsi_corr",";#psi_{ZDC}",1000,0,twoPI);
  zdcPsi_east = new TH1D("zdcPsi_east",";#psi_{ZDC}",1000,0,twoPI);
  zdcPsi_west = new TH1D("zdcPsi_west",";#psi_{ZDC}",1000,0,twoPI);
  centra_weight = new TH2D("centra_weight", "",10,0,10,1000,-10,10);
  pionV1Plus = new TProfile("pionV1Plus","",48,-1.2,1.2);
  pionV1Minus = new TProfile("pionV1Minus","",48,-1.2,1.2);
  zdcResolution = new TProfile("zdcResolution","",10,0,10);
  testTrack = new TH1D("testTrack","",10,0,10);
  trackPhiEta = new TH2D("trackPhiEta",";#phi;#eta",100,-1.*PI,PI,10,-1,1);
  trackPhiEtaHFT = new TH2D("trackPhiEtaHFT",";#phi;#eta",100,-1.*PI,PI,10,-1,1);
  hHitsDedx = new TH1D("hitsDedx","",100,0,100);
  hSigmaPiBeta = new TH2D("sigmaPiBeta","",100,0,5,50,0,0.05);
  hCentrality = new TH1D("hCentrality","",10,0,10);
  hCentralityWeight = new TH1D("hCentralityWeight","",10,0,10);
  h2RefMult = new TH2D("h2RefMult",";gRefMult;refMult",1500,0,1500,1500,0,1500);
  h2RefMult_0_80 = new TH2D("h2RefMult_0_80",";gRefMult;refMult",1500,0,1500,1500,0,1500);
  h2RefMult_10_80 = new TH2D("h2RefMult_10_80",";gRefMult;refMult",1500,0,1500,1500,0,1500);
  pionV1Plus->Sumw2();
  pionV1Minus->Sumw2();
  zdcResolution->Sumw2();
  d0MassPhiEta->Sumw2();
  d0BarMassPhiEta->Sumw2();
  d0MassPhiEta_50->Sumw2();
  d0BarMassPhiEta_50->Sumw2();
  d0MassPhiEta_pt3->Sumw2();
  d0BarMassPhiEta_pt3->Sumw2();
  d0MassPhiEta_150->Sumw2();
  d0BarMassPhiEta_150->Sumw2();
  d0MassPhiEta_5bin->Sumw2();
  d0BarMassPhiEta_5bin->Sumw2();
  d0MassPhiEta_noeff->Sumw2();
  d0BarMassPhiEta_noeff->Sumw2();
  d0MassPhiEta_nHits25->Sumw2();
  d0BarMassPhiEta_nHits25->Sumw2();
  d0MassPhiEta_noHitsRatio->Sumw2();
  d0BarMassPhiEta_noHitsRatio->Sumw2();
  fout1.open("xcheck_topo_d0_candidates.csv");
  fout1<<"D0phi,pipt,DL,kpt,evId,kdca,dca12i,D0rap,keta,pidca,D0mass,nD0,nD0bar,D0pt,D0dca,pieta,cosTheta,kaonIdx,pionIdx,from"<<endl;
  // d0_candidates_tuple = new TNtuple("d0_candidates_tuple","","D0phih:pipt:DL:kpt:evId:kdca:dca12i:D0rap:keta:pidca:D0mass:nD0:nD0bar:D0pt:D0dca:pieta:cosTheta");
  
///////////////Load efficiency files/////////////
  TFile* fRapEff = new TFile("D0_eff_Run14.root");
  h3PtCentY_clone = (TH3F *)fRapEff->Get("h3PtCentY")->Clone("h3PtCentY_clone");
  h3PtCentYcut_clone = (TH3F *)fRapEff->Get("h3PtCentYcut")->Clone("h3PtCentYcut_clone");
  h3PtY_clone = (TH2F *)h3PtCentY_clone->Project3D("zx")->Clone("h3PtY_clone");
  h3PtYcut_clone = (TH2F *)h3PtCentYcut_clone->Project3D("zx")->Clone("h3PtYcut_clone");

  TFile* fPionPlus= new TFile("Eff_PionPlus_embedding.root");
  for(int i = 0;i<9;i++)
  {
    hPionPlus[i] = (TH1F *)fPionPlus->Get(Form("h1Ratiocent_%d",i))->Clone(Form("hPionPlus_%d",i));
    hPionPlus_noRatio[i] = (TH1F *)fPionPlus->Get(Form("h1Ratiocent_noRatio_%d",i))->Clone(Form("hPionPlus_noRatio_%d",i));
    hPionPlus_nHits25[i] = (TH1F *)fPionPlus->Get(Form("h1Ratiocent_nHits25_%d",i))->Clone(Form("hPionPlus_nHits25_%d",i));
  }
  TFile* fPionMinus= new TFile("Eff_PionMinus_embedding.root");
  for(int i = 0;i<9;i++)
  {
    hPionMinus[i] = (TH1F *)fPionMinus->Get(Form("h1Ratiocent_%d",i))->Clone(Form("hPionMinus_%d",i));
    hPionMinus_noRatio[i] = (TH1F *)fPionMinus->Get(Form("h1Ratiocent_noRatio_%d",i))->Clone(Form("hPionMinus_noRatio_%d",i));
    hPionMinus_nHits25[i] = (TH1F *)fPionMinus->Get(Form("h1Ratiocent_nHits25_%d",i))->Clone(Form("hPionMinus_nHits25_%d",i));
  }

  TFile* fKaonPlus= new TFile("Eff_KaonPlus_embedding.root");
  for(int i = 0;i<9;i++)
  {
    hKaonPlus[i] = (TH1F *)fKaonPlus->Get(Form("h1Ratiocent_%d",i))->Clone(Form("hKaonPlus_%d",i));
    hKaonPlus_noRatio[i] = (TH1F *)fKaonPlus->Get(Form("h1Ratiocent_noRatio_%d",i))->Clone(Form("hKaonPlus_noRatio_%d",i));
    hKaonPlus_nHits25[i] = (TH1F *)fKaonPlus->Get(Form("h1Ratiocent_nHits25_%d",i))->Clone(Form("hKaonPlus_nHits25_%d",i));
  }
  TFile* fKaonMinus= new TFile("Eff_KaonMinus_embedding.root");
  for(int i = 0;i<9;i++)
  {
    hKaonMinus[i] = (TH1F *)fKaonMinus->Get(Form("h1Ratiocent_%d",i))->Clone(Form("hKaonMinus_%d",i));
    hKaonMinus_noRatio[i] = (TH1F *)fKaonMinus->Get(Form("h1Ratiocent_noRatio_%d",i))->Clone(Form("hKaonMinus_noRatio_%d",i));
    hKaonMinus_nHits25[i] = (TH1F *)fKaonMinus->Get(Form("h1Ratiocent_nHits25_%d",i))->Clone(Form("hKaonMinus_nHits25_%d",i));
  }


  pionDca = new TH2D("pionDca","",10,0,10,1000,0,0.1);
  kaonDca = new TH2D("kaonDca","",10,0,10,1000,0,0.1);
  decayLength = new TH2D("decayLength","",10,0,10,1000,0,0.1);
  dcaDaughter = new TH2D("dcaDaughter","",10,0,10,1000,0,0.1);
  pointingAngle = new TH2D("pointingAngle","",10,0,10,1000,0,0.1);
  pionPtEta = new TH2D("pionPtEta","",200,-1,1,100,0,10);
  kaonPtEta = new TH2D("kaonPtEta","",200,-1,1,100,0,10);

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
  // fRapEff->Close();
  mOutputFile->cd();
  pionPtEta->Write();
  kaonPtEta->Write();
  pionDca->Write();
  kaonDca->Write();
  decayLength->Write();
  dcaDaughter->Write();
  pointingAngle->Write();
  d0MassPhiEta->Write();
  d0BarMassPhiEta->Write();
  // d0MassPhiEta_pt3->Write();
  // d0BarMassPhiEta_pt3->Write();
  // d0MassPhiEta_50->Write();
  // d0BarMassPhiEta_50->Write();
  // d0MassPhiEta_150->Write();
  // d0BarMassPhiEta_150->Write();
  // d0MassPhiEta_5bin->Write();
  // d0BarMassPhiEta_5bin->Write();
  d0MassPhiEta_noeff->Write();
  d0BarMassPhiEta_noeff->Write();
  d0MassPhiEta_noHitsRatio->Write();
  d0BarMassPhiEta_noHitsRatio->Write();
  d0MassPhiEta_nHits25->Write();
  d0BarMassPhiEta_nHits25->Write();
  testDPhi->Write();
  testEvent->Write();
  d0MassPt->Write();
  d0BarMassPt->Write();
  zdcPsi->Write();
  zdcPsi_corr->Write();
  zdcPsi_east->Write();
  zdcPsi_west->Write();
  centra_weight->Write();
  pionV1Plus->Write();
  pionV1Minus->Write();
  zdcResolution->Write();
  testTrack->Write();
  trackPhiEta->Write();
  trackPhiEtaHFT->Write();
  hHitsDedx->Write();
  hSigmaPiBeta->Write();
  hCentrality->Write();
  hCentralityWeight->Write();
  h2RefMult->Write();
  h2RefMult_0_80->Write();
  h2RefMult_10_80->Write();
  miniZDCSMD->WriteHist();

  cout<<"writing new histo"<<endl;
  // save user variables here
  delete mPrescales;
  mOutputFile->Close();

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
  // mGRefMultCorrUtil->print();
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int cent = centrality+1;
  h2RefMult->Fill(event->grefMult(), event->refMult());
  if(picoDst->event()->runId() < 15107000 || mGRefMultCorrUtil->isBadRun(picoDst->event()->runId())) 
  {
    LOG_WARN << "bad run!" << endl;
    return kStWarn;
  }
  if(centrality<0) 
  {
    LOG_WARN << "not minBias sample!" << endl;
    return kStWarn;
  }
  testEvent->Fill(3);
  h2RefMult_0_80->Fill(event->grefMult(), event->refMult());
  if(centrality >= 7) {
  // if(centrality<0 ) {
    LOG_WARN << "0-10\% centrality" << endl;
    return kStWarn;
  }
  testEvent->Fill(4);
  h2RefMult_10_80->Fill(event->grefMult(), event->refMult());
  double reweight = mGRefMultCorrUtil->getWeight();
  hCentralityWeight->Fill(centrality, reweight);
  centra_weight->Fill(centrality, reweight);
  hCentrality->Fill(centrality);

  // cout<<"centrality = "<<centrality<<endl;
  for(unsigned int i=0;i<picoDst->numberOfTracks();i++)
  {
    StPicoTrack const* itrk = picoDst->track(i);
    testTrack->Fill(0);
    if(!itrk->isHFTTrack()) continue;
    if(!isGoodTrack(itrk,event))  continue;
    testTrack->Fill(1);

    StThreeVectorF Vertex3D=event->primaryVertex();
    float b = event->bField();
    StThreeVectorF gmom = itrk->gMom(Vertex3D,b);

    bool tpcPion = isTpcPion(itrk);
    bool tpcKaon = isTpcKaon(itrk,&pVtx);
    float kBeta = getTofBeta(itrk,&pVtx);
    float pBeta = getTofBeta(itrk,&pVtx);
    bool kTofAvailable = kBeta>0;
    bool pTofAvailable = pBeta>0;
    bool tofKaon = kTofAvailable && isTofKaon(itrk,kBeta);
    bool tofPion = pTofAvailable && isTofPion(itrk,pBeta);
    bool goodKaon = kTofAvailable ? tofKaon&&tpcKaon : tpcKaon;
    bool goodPion = pTofAvailable ? tofPion&&tpcPion : tpcPion;
    if(goodKaon)
    {
      testTrack->Fill(2);
      kaonPtEta->Fill(gmom.pseudoRapidity(),gmom.perp());
    }
    if(goodPion)
    {
      testTrack->Fill(3);
      pionPtEta->Fill(gmom.pseudoRapidity(),gmom.perp());
    }
  }
  // if(!(centrality==4||centrality==5||centrality==6))//minBias trigger requires
  // {
  //   LOG_WARN << " Not Good Event! Skip! " << endm;
  //   return kStWarn;
  // }
  // const int  runID    = 17109018;//event->runId();
  const int  runID    = event->runId();
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
  zdcPsi_east->Fill(mZDC1Event_PsiE);
  zdcPsi_west->Fill(mZDC1Event_PsiW);
	zdcResolution->Fill(cent,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 
	zdcResolution->Fill(0.,cos(1.*(mZDC1Event_PsiE-mZDC1Event_PsiW+PI))); 
////////////////////////////////////////////////////
  // cout<<"zdcsmd done."<<endl;
  // if (reweight < 0.5)
    // cout<<"event reweight =  "<<reweight<<endl;
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    // int kaonId = kp->kaonIdx();
    // int pionId = kp->pionIdx();
    float d0Pt = kp->pt();
    if(d0Pt < 1.5)  continue;
    // if(d0Pt < 2.5)  continue;
    if (cos(kp->pointingAngle()) < 0.99) continue;
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());
    StThreeVectorF pimom = pion->gMom(pVtx,event->bField());
    StThreeVectorF kmom = kaon->gMom(pVtx,event->bField());

    if (!kaon->isHFTTrack() || !pion->isHFTTrack()) continue;

    // if (!isGoodTrack(kaon,event) || !isGoodTrack(pion,event)) continue;
    if (!isGoodTrack(kaon,event)) continue;
    if (!isGoodTrack(pion,event)) continue;

    // if (!isTpcPion(pion)) continue;
    bool tpcPion = isTpcPion(pion);
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    float pBeta = getTofBeta(pion,&pVtx);
    bool kTofAvailable = kBeta>0;
    bool pTofAvailable = pBeta>0;
    bool tofKaon = kTofAvailable && isTofKaon(kaon,kBeta,event);
    bool tofPion = pTofAvailable && isTofPion(pion,pBeta,event);
    bool goodKaon = kTofAvailable ? tofKaon&&tpcKaon : tpcKaon;
    bool goodPion = pTofAvailable ? tofPion&&tpcPion : tpcPion;
    // bool goodKaon = tpcKaon && tofKaon;
    // bool goodPion = tpcPion && tofPion;

    // double diffBeta = 0;
    // if(pBeta>0)
    // {
    //   double ptot = pion->dcaGeometry().momentum().mag();
    //   float beta_p = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
    //   diffBeta = fabs(1/pBeta - 1/beta_p);
    // }
    // hHitsDedx->Fill(pion->nHitsDedx());
    // hSigmaPiBeta->Fill(pion->nSigmaPion(),diffBeta);
    if(!goodKaon) continue;
    if(!goodPion) continue;
    int charge=0;
    // double dMass = kp->m();
    // int centBin = 0;
    // if(centrality>=7) centBin=1;
    // else if(centrality>=4)  centBin=2;
    // else centBin=3;
    // int fitindex = 5;
    // if(d0Pt<5)
    //   fitindex = static_cast<int>(d0Pt);
    // double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    // double reweight_eff = 1;//= (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    charge=isD0Pair(kp);
    int count_d0 = 0;
    int count_d0bar = 0;
    if(charge<0)
    {
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      double kpY = d0Lorentz.Rapidity();
      if(fabs(kpY) > 0.8)
        continue;
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      kaonDca->Fill(kp->pt(),kp->kaonDca());
      pionDca->Fill(kp->pt(),kp->pionDca());
      dcaDaughter->Fill(kp->pt(),kp->dcaDaughters());
      decayLength->Fill(kp->pt(),kp->decayLength());
      pointingAngle->Fill(kp->pt(),sin(kp->pointingAngle())*kp->decayLength());
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;
      deltaPhi = fabs(deltaPhi);

      double deltaPhi1 = fabs(kp->phi()-mZDC1Event_PsiF);
      if (deltaPhi1 > twoPI)
        deltaPhi1 = 2.*twoPI - deltaPhi1;
      if (deltaPhi1 > PI)
        deltaPhi1 = twoPI - deltaPhi1;
      cout<<deltaPhi1<<","<<deltaPhi<<endl;

      float d0Eff = getD0Efficiency(d0Lorentz);
      // cout<<"D0 efficiency = "<<d0Eff<<endl;

      testDPhi->Fill(deltaPhi);
      // cout<<"deltaPhi = "<<deltaPhi<<endl;
      if(kaon->charge() < 0)
      {
        count_d0++;
        testTrack->Fill(5);
        d0MassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        if(d0Pt >= 3)
          d0MassPhiEta_pt3->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        d0MassPhiEta_5bin->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        d0MassPhiEta_noeff->Fill(kp->m(),deltaPhi,kpY,1.*reweight);
        d0MassPt->Fill(kp->m(),kp->pt(),kpY,1.);
      }
      if(kaon->charge() > 0)
      {
        count_d0bar++;
        testTrack->Fill(6);
        d0BarMassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        if(d0Pt >= 3)
          d0BarMassPhiEta_pt3->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        d0BarMassPhiEta_5bin->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        d0BarMassPhiEta_noeff->Fill(kp->m(),deltaPhi,kpY,1.*reweight);
        d0BarMassPt->Fill(kp->m(),kp->pt(),kpY,1.);
      }
      fout1<<kp->phi()<<","<<pimom.perp()<<","<<kp->decayLength()<<","<<kmom.perp()<<","<<mPicoD0Event->eventId()<<",";
      fout1<<kp->kaonDca()<<","<<kp->dcaDaughters()<<","<<kpY<<","<<kmom.pseudoRapidity()<<","<<kp->pionDca()<<","<<kp->m()<<","<<count_d0<<",";
      fout1<<count_d0bar<<","<<kp->pt()<<","<<sin(kp->pointingAngle())*kp->decayLength()<<","<<pimom.pseudoRapidity()<<","<<cos(kp->pointingAngle())<<",";
      fout1<<kp->kaonIdx()<<","<<kp->pionIdx()<<",Leon"<<","<<centrality<<endl;

    }//D loop
    /*
    if(isD0Pair50(kp)<0)
    {
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      double kpY = d0Lorentz.Rapidity();
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;

      float d0Eff = getD0Efficiency(d0Lorentz);
      deltaPhi = fabs(deltaPhi);
      if(kaon->charge() < 0)
      {
        d0MassPhiEta_50->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
      }
      if(kaon->charge() > 0)
      {
        d0BarMassPhiEta_50->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
      }
    }// D0 50 
    if(isD0Pair150(kp)<0)
    {
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      double kpY = d0Lorentz.Rapidity();
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;

      float d0Eff = getD0Efficiency(d0Lorentz);
      deltaPhi = fabs(deltaPhi);
      if(kaon->charge() < 0)
      {
        d0MassPhiEta_150->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
      }
      if(kaon->charge() > 0)
      {
        d0BarMassPhiEta_150->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
      }
    }// D0 50 
    */

  }
  /*
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    int kaonId = kp->kaonIdx();
    int pionId = kp->pionIdx();
    float d0Pt = kp->pt();
    if(d0Pt < 1.5)  continue;
    if (cos(kp->pointingAngle()) < 0.99) continue;
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());
    StKaonPion *kp_test = new StKaonPion(kaon,pion,kaonId,pionId,pVtx,event->bField());
    StThreeVectorF pimom = pion->gMom(pVtx,event->bField());
    StThreeVectorF kmom = kaon->gMom(pVtx,event->bField());

    if (!kaon->isHFTTrack() || !pion->isHFTTrack()) continue;

    if (!isLooseTrack(kaon,event)) continue;
    if (!isLooseTrack(pion,event)) continue;

    bool tpcPion = isTpcPion(pion);
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    float pBeta = getTofBeta(pion,&pVtx);
    bool kTofAvailable = kBeta>0;
    bool pTofAvailable = pBeta>0;
    bool tofKaon = kTofAvailable && isTofKaon(kaon,kBeta,event);
    bool tofPion = pTofAvailable && isTofPion(pion,pBeta,event);
    bool goodKaon = kTofAvailable ? tofKaon&&tpcKaon : tpcKaon;
    bool goodPion = pTofAvailable ? tofPion&&tpcPion : tpcPion;
    if(!goodKaon) continue;
    if(!goodPion) continue;
    int charge=0;
    int centBin = 0;
    if(centrality>=7) centBin=1;
    else if(centrality>=4)  centBin=2;
    else centBin=3;
    int fitindex = 5;
    if(d0Pt<5)
      fitindex = static_cast<int>(d0Pt);
    double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    // double reweight_eff = 1;//= (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    charge=isD0Pair(kp);
    int count_d0 = 0;
    int count_d0bar = 0;
    if(charge<0)
    {
    cout<<"x,y,z = "<<pVtx.x()<<","<<pVtx.y()<<","<<pVtx.z()<<",field = "<<event->bField()<<endl;
    cout<<"kp: pionDca = "<<kp->pionDca()<<"\ttestkp = "<<kp_test->pionDca()<<endl;
    cout<<"kp: kaonDca = "<<kp->kaonDca()<<"\ttestkp = "<<kp_test->kaonDca()<<endl;
    cout<<"kp: pointingAngle = "<<kp->pointingAngle()<<"\ttestkp = "<<kp_test->pointingAngle()<<endl;
    cout<<"kp: dcaDaughters = "<<kp->dcaDaughters()<<"\ttestkp = "<<kp_test->dcaDaughters()<<endl;
    cout<<"kp: decayLength = "<<kp->decayLength()<<"\ttestkp = "<<kp_test->decayLength()<<endl;
      TLorentzVector d0Lorentz;
      d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
      double kpY = d0Lorentz.Rapidity();
      double deltaPhi = fabs(kp->phi()-mZDC1Event_PsiF);
      if(deltaPhi > PI)
        deltaPhi = twoPI - deltaPhi;
      deltaPhi = fabs(deltaPhi);
      float d0Eff = getD0Efficiency(d0Lorentz);
      int ptBin = (int)(d0Pt/0.1+1);
      double tpcEff = hPionMinus[centrality]->GetBinContent(ptBin) * hKaonPlus[centrality]->GetBinContent(ptBin);
      double tpcEff_noRatio = hPionMinus_noRatio[centrality]->GetBinContent(ptBin) * hKaonPlus_noRatio[centrality]->GetBinContent(ptBin);
      double tpcEff_nHits25 = hPionMinus_nHits25[centrality]->GetBinContent(ptBin) * hKaonPlus_nHits25[centrality]->GetBinContent(ptBin);
      cout<<"tpcEff = "<<tpcEff<<","<<tpcEff_noRatio<<","<<tpcEff_nHits25<<endl;
      deltaPhi = fabs(deltaPhi);
      if(kaon->charge() < 0)
      {
        if(pion->nHitsFit() >= 20 && kaon->nHitsFit() >= 20)
          d0MassPhiEta_noHitsRatio->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        if(pion->nHitsFit() >= 25 && kaon->nHitsFit() >= 25
            && (1.0*pion->nHitsFit()/pion->nHitsMax())>=0.52 && (1.0*kaon->nHitsFit()/kaon->nHitsMax())>=0.52)
          d0MassPhiEta_nHits25->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff*tpcEff_noRatio/tpcEff_nHits25);
        if(pion->nHitsFit() >= 20 && kaon->nHitsFit() >= 25
            && (1.0*pion->nHitsFit()/pion->nHitsMax())>=0.52 && (1.0*kaon->nHitsFit()/kaon->nHitsMax())>=0.52)
          d0MassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff*tpcEff/tpcEff_nHits25);
      }
      if(kaon->charge() > 0)
      {
        if(pion->nHitsFit() >= 20 && kaon->nHitsFit() >= 20)
          d0BarMassPhiEta_noHitsRatio->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff);
        if(pion->nHitsFit() >= 25 && kaon->nHitsFit() >= 25
            && (1.0*pion->nHitsFit()/pion->nHitsMax())>=0.52 && (1.0*kaon->nHitsFit()/kaon->nHitsMax())>=0.52)
          d0BarMassPhiEta_nHits25->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff*tpcEff_nHits25/tpcEff_nHits25);
        if(pion->nHitsFit() >= 25 && kaon->nHitsFit() >= 25
            && (1.0*pion->nHitsFit()/pion->nHitsMax())>=0.52 && (1.0*kaon->nHitsFit()/kaon->nHitsMax())>=0.52)
          d0BarMassPhiEta->Fill(kp->m(),deltaPhi,kpY,1.*reweight/d0Eff*tpcEff/tpcEff_nHits25);
      }

    }//D loop
  }
  */
  return kStOK;
}
//-----------------------------------------------------------------------------

int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  // TLorentzVector d0Lorentz;
  // d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
  // if(fabs(d0Lorentz.Rapidity())>1.) return 0;
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
  // if(charge>0)
  //   charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
    // return -1;
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
  // if(charge>0)
  //   charge = kaon->charge()>0 ? 1:2;


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
  // if(charge>0)
  //   charge = kaon->charge()>0 ? 1:2;


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

//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofPion(StPicoTrack const * const trk, float beta) const
{
  bool tofPion= false;

  if(beta>0)
  {
    double ptot = trk->dcaGeometry().momentum().mag();
    float beta_p = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
    tofPion = fabs(1/beta - 1/beta_p) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofPion;
}
///-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta,StPicoEvent const*mEvent) const
{
  bool tofKaon = false;

  StThreeVectorF vtx=mEvent->primaryVertex();
  float b = mEvent->bField();
  if(beta>0)
  {
    double ptot = trk->gMom(vtx,b).mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}

//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofPion(StPicoTrack const * const trk, float beta,StPicoEvent const*mEvent) const
{
  bool tofPion= false;

  StThreeVectorF vtx=mEvent->primaryVertex();
  float b = mEvent->bField();
  if(beta>0)
  {
    // double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gMom(vtx,b).mag();
    float beta_p = ptot/sqrt(ptot*ptot+M_PION_PLUS*M_PION_PLUS);
    tofPion = fabs(1/beta - 1/beta_p) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofPion;
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
  // if(refMult <=2 || refMult > 1000) return false;
  // if(fabs(VertexZ) > 100) return false; 

  // hVertex2D ->Fill(VertexZ,vpdVz);
  // hDiffVz   ->Fill(VertexZ-vpdVz); 

  if(fabs(VertexZ) > 6) return false; 
  if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0)return false; 
  if(fabs(VertexZ-vpdVz)>3.)return false;       // no vpd cut in low energy?
  if(fabs(VertexX)<1e-5 && fabs(VertexY)<1e-5 && fabs(VertexZ)<1e-5 ) return false;
  testEvent->Fill(2);
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

bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk,StPicoEvent const*mEvent) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  // return trk->gPt() > 0.6 && trk->nHitsFit() >= 20 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 && trk->nHitsDedx() >= 16;
  StThreeVectorF Vertex3D=mEvent->primaryVertex();
  float b = mEvent->bField();
  StThreeVectorF gmom = trk->gMom(Vertex3D,b);
  StPhysicalHelixD helix = trk->helix();
  return gmom.perp() > 0.6 && trk->nHitsFit() >= 20 && fabs(gmom.pseudoRapidity()) < 1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>=0.52 && fabs(helix.geometricSignedDistance(Vertex3D)) > 0.005;
  // return trk->gPt() > 0.6 && trk->nHitsFit() >= 20  && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52 ;
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}

bool StPicoD0AnaMaker::isLooseTrack(StPicoTrack const * const trk,StPicoEvent const*mEvent) const
{
  StThreeVectorF Vertex3D=mEvent->primaryVertex();
  float b = mEvent->bField();
  StThreeVectorF gmom = trk->gMom(Vertex3D,b);
  StPhysicalHelixD helix = trk->helix();
  return gmom.perp() > 0.6 && fabs(gmom.pseudoRapidity()) < 1.&& fabs(helix.geometricSignedDistance(Vertex3D)) > 0.005;
}
float StPicoD0AnaMaker::getD0Efficiency(TLorentzVector &d0)
{
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  float d0Rap = d0.Rapidity();
  float d0Pt = d0.Pt();
  if (d0Pt > 7.9)
    d0Pt = 7.9;
  int yBin = (int)((d0Rap+1.0)/0.1); 
  int ptBin = (int)(d0Pt/0.1+1);
  // int ptBinRef = (int)(5./0.1+1);
  int centBin = centrality+1;
  double eff = 1.0 * (h3PtCentYcut_clone->GetBinContent(ptBin,centBin,yBin)) / (h3PtCentY_clone->GetBinContent(ptBin,centBin,yBin));
  double eff_ref = 1.0 * (h3PtYcut_clone->GetBinContent(ptBin,yBin)) / (h3PtY_clone->GetBinContent(ptBin,yBin));
 
  // double eff_5 = 1.0 * (h3PtCentYcut_clone->GetBinContent(ptBinRef,centBin,yBin)) / (h3PtCentY_clone->GetBinContent(ptBinRef,centBin,yBin));
  if(eff_ref > 0)
    return eff;
  else
    return 0;
}



