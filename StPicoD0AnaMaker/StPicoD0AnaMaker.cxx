#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"

#include "phys_constants.h"
#include "StCuts.h"
#include "StBTofUtil/tofPathLength.hh"
ClassImp(StPicoD0AnaMaker)

// StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, 
    char const * outName,StPicoDstMaker* picoDstMaker): 
  StMaker(name),mPicoDstMaker(picoDstMaker), mOutFileName(outName), 
  mOutputFile(NULL),  mEventCounter(0)
{}

Int_t StPicoD0AnaMaker::Init()
{
   // mPicoD0Event = new StPicoD0Event();

   // mChain = new TChain("T");
   // std::ifstream listOfFiles(mInputFileList.Data());
   // if (listOfFiles.is_open())
   // {
   //    std::string file;
   //    while (getline(listOfFiles, file))
   //    {
   //       LOG_INFO << "StPicoD0AnaMaker - Adding :" << file << endm;
   //       mChain->Add(file.c_str());
   //    }
   // }
   // else
   // {
   //    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
   //    return kStErr;
   // }
   //
   // mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   // mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mOutputFile->cd();

   hPVZ = new TH2D("pvz",";tpc vz;vpd vz",500,-100,100,500,-100,100);
   deltaVZ = new TH1D("deltaVZ","",500,-10,10);

   hEta = new TH1D("eta","",200,-1.5,1.5);
   hPhi = new TH1D("phi","",200,-4,4);
   hPt = new TH1D("pt","",200,0,5);
   hDcaXY = new TH2D("dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hDcaZ = new TH2D("dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);

   hTPCEta = new TH1D("tpc_eta","",200,-1.5,1.5);
   hTPCPhi = new TH1D("tpc_phi","",200,-4,4);
   hTPCPt = new TH1D("tpc_pt","",200,0,5);
   hTPCDcaXY = new TH2D("tpc_dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hTPCDcaZ = new TH2D("tpc_dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);
   hTPCDcaXY_wide = new TH2D("tpc_dcaXY_wide",";pt;dcaXY",200,0,5,300,-1.5,1.5);
   hTPCDcaZ_wide = new TH2D("tpc_dcaZ_wide",";pt;dcaZ",200,0,5,300,-1.5,1.5);

   hHFTEta = new TH1D("hft_eta","",200,-1.5,1.5);
   hHFTPhi = new TH1D("hft_phi","",200,-4,4);
   hHFTPt = new TH1D("hft_pt","",200,0,5);
   hHFTDcaXY = new TH2D("hft_dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hHFTDcaZ = new TH2D("hft_dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);
   hHFTDcaXY_wide = new TH2D("hft_dcaXY_wide",";pt;dcaXY",200,0,5,300,-1.5,1.5);
   hHFTDcaZ_wide = new TH2D("hft_dcaZ_wide",";pt;dcaZ",200,0,5,300,-1.5,1.5);

   hD0MassLike = new TH2D("d0MassLike","",100,0,10,500,1.6,2.1);
   hD0MassUnlike = new TH2D("d0MassUnlike","",100,0,10,500,1.6,2.1);
   hDecayLength = new TH1D("decayLength","",500,0,0.05);
   hKaonDca = new TH1D("kaonDca","",200,0,0.02);
   hPionDca = new TH1D("pionDca","",200,0,0.02);
   hDcaDaughters = new TH1D("dcaDaughters","",500,0,0.01);
   hPointing = new TH1D("pointing","",500,0,0.01);
   trackTest1 = new TH1D("trackTest1","",10,0,10);
   trackTest2 = new TH1D("trackTest2","",10,0,10);

   tuplePair = new TNtuple("tuplePair","","mass:pt:kDca:pDca:kpDca:theta:decayLenght:charge:dpt");

   hVzVpdVz = new TH2F("hVzVpdVz","hVzVpdVz",200,-100,100,200,-100,100);
   ntpD0BackgroundSameSignPt01 = new TNtuple("ntpD0BackgroundSameSignPt01", "ntpD0BackgroundSameSignPt01", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt01 = new TNtuple("ntpD0BackgroundSideBandPt01", "ntpD0BackgroundSideBandPt01", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt12 = new TNtuple("ntpD0BackgroundSameSignPt12", "ntpD0BackgroundSameSignPt12", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt12 = new TNtuple("ntpD0BackgroundSideBandPt12", "ntpD0BackgroundSideBandPt12", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt23 = new TNtuple("ntpD0BackgroundSameSignPt23", "ntpD0BackgroundSameSignPt23", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt23 = new TNtuple("ntpD0BackgroundSideBandPt23", "ntpD0BackgroundSideBandPt23", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt35 = new TNtuple("ntpD0BackgroundSameSignPt35", "ntpD0BackgroundSameSignPt35", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt35 = new TNtuple("ntpD0BackgroundSideBandPt35", "ntpD0BackgroundSideBandPt35", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt510 = new TNtuple("ntpD0BackgroundSameSignPt510", "ntpD0BackgroundSameSignPt510", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt510 = new TNtuple("ntpD0BackgroundSideBandPt510", "ntpD0BackgroundSideBandPt510", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");


   // -------------- USER VARIABLES -------------------------
   
   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   hD0MassLike->Write();
   hD0MassUnlike->Write();
   hDecayLength->Write();
   hKaonDca->Write();
   hPionDca->Write();
   hDcaDaughters->Write();
   hPointing->Write();

   hPVZ->Write();
   hVzVpdVz->Write();
   deltaVZ->Write();
   hEta->Write();
   hTPCEta->Write();
   hHFTEta->Write();
   hPhi->Write();
   hTPCPhi->Write();
   hHFTPhi->Write();
   hPt->Write();
   hTPCPt->Write();
   hHFTPt->Write();
   hDcaZ->Write();
   hTPCDcaZ->Write();
   hHFTDcaZ->Write();
   hTPCDcaZ_wide->Write();
   hHFTDcaZ_wide->Write();
   hDcaXY->Write();
   hTPCDcaXY->Write();
   hHFTDcaXY->Write();
   hTPCDcaXY_wide->Write();
   hHFTDcaXY_wide->Write();

   // tuplePair->Write();
   trackTest1->Write();
   trackTest2->Write();

   mOutputFile->Write();
   mOutputFile->Close();

   return kStOK;
}
//-----------------------------------------------------------------------------
void StPicoD0AnaMaker::Clear(Option_t *opt)
{
}
Int_t StPicoD0AnaMaker::Make()
{
   readNextEvent();
   // cout<<"new event coming"<<endl;

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   // if(mPicoD0Event->runId() != picoDst->event()->runId() ||
   //     mPicoD0Event->eventId() != picoDst->event()->eventId())
   // {
   //   LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
   //   LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
   //   exit(1);
   // }

   // -------------- USER ANALYSIS -------------------------
   
   // check if good event (including bad run)
   // if(!mHFCuts->isGoodEvent(const_cast<const StPicoDst*>(picoDst), NULL)) 
   //   return kStOk;

   // Loop all tracks
  bool isGoodTrigger = false;
   std::array<unsigned int, 5> const triggers = { 530003,    // vpdmb-5-p-nobsmd-hlt 
                                                  430903,    // vpdmb-5-p-nobsmd-hlt 
                                                  430904,    // vpdmb-5-p-nobsmd-hlt 
                                                  430906,    // vpdmb-5-p-nobsmd-hlt 
                                                  430907    // vpdmb-5-p-nobsmd-hlt 
                                                };
  for(auto trg: triggers)
  {
    if(picoDst->event()->isTrigger(trg)) isGoodTrigger = true;
  }
   StThreeVectorF vtx = picoDst->event()->primaryVertex();
   float b = picoDst->event()->bField();
   double vpdvz = picoDst->event()->vzVpd();
   double vz = vtx.z();
  // if(vz==0)
  //    exit(1);
   hPVZ->Fill(vz,vpdvz);
   hVzVpdVz->Fill(vz,vpdvz);
   deltaVZ->Fill(vz-vpdvz);

   // if(!isGoodTrigger)
   // {
   //   // cout<<"not good trigger"<<endl;
   //   // exit(1);
   //   return kStWarn;
   // }
   if( fabs(vz)>6 && fabs(vz-vpdvz)>3 ) 
     return kStWarn;

   
   // cout<<"number of tracks = "<<picoDst->numberOfTracks()<<endl;
   for(int i=0;i<picoDst->numberOfTracks();i++)
   {
     StPicoTrack const* trk = picoDst->track(i);
     double eta = trk->gMom(vtx, b).pseudoRapidity();
     double phi = trk->gMom(vtx, b).phi() ;
     double pt = trk->gPt();
     StPhysicalHelixD helix = trk->helix(b);
     double dcaZ = (trk->dcaPoint().z() - vtx.z());
     double dcaXY = helix.geometricSignedDistance(vtx.x(),vtx.y());
     double dca = (vtx-trk->dcaPoint()).mag();
     // TOF 
     int tofIndex = trk->bTofPidTraitsIndex();
     bool TofMatch = kFALSE;
     StPicoBTofPidTraits* tofPidTraits;
     if (tofIndex >= 0)  tofPidTraits = picoDst->btofPidTraits(tofIndex); //GNX
     if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
     if(!TofMatch) continue;
     if(dca>3)  continue;

     if(pt<0.2)  continue;
     if(fabs(eta)>1) continue;
     hEta->Fill(eta);
     hPhi->Fill(phi);
     hPt->Fill(pt);
     hDcaXY->Fill(pt,dcaXY);
     hDcaZ->Fill(pt,dcaZ);
     if(trk->nHitsFit()<20) continue;
     hTPCEta->Fill(eta);
     hTPCPhi->Fill(phi);
     hTPCPt->Fill(pt);
     hTPCDcaXY->Fill(pt,dcaXY);
     hTPCDcaZ->Fill(pt,dcaZ);
     hTPCDcaXY_wide->Fill(pt,dcaXY);
     hTPCDcaZ_wide->Fill(pt,dcaZ);
     if(!(trk->isHFTTrack())) continue;
     hHFTEta->Fill(eta);
     hHFTPhi->Fill(phi);
     hHFTPt->Fill(pt);
     hHFTDcaXY->Fill(pt,dcaXY);
     hHFTDcaZ->Fill(pt,dcaZ);
     hHFTDcaXY_wide->Fill(pt,dcaXY);
     hHFTDcaZ_wide->Fill(pt,dcaZ);
   }


   for(int i=0;i<picoDst->numberOfTracks();i++)
   {
     StPicoTrack const* itrk = picoDst->track(i);
     trackTest1->Fill(1);
     if(!isGoodTrack(itrk))  continue;
     trackTest1->Fill(2);
     if (!isTpcPion(itrk)) continue;
     trackTest1->Fill(3);
     if(!(itrk->isHft())) continue;
     trackTest1->Fill(4);
     for(int j=0;j<picoDst->numberOfTracks();j++)
     {
       StPicoTrack const* jtrk = picoDst->track(j);
       trackTest2->Fill(1);
       if(!isGoodTrack(jtrk))  continue;
       trackTest2->Fill(2);
       if(!(jtrk->isHFTTrack())) continue;
       trackTest2->Fill(3);
       bool tpcKaon = isTpcKaon(jtrk,&vtx);
       float kBeta = getTofBeta(jtrk,&vtx);
       bool tofAvailable = kBeta>0;
       bool tofKaon = tofAvailable && isTofKaon(jtrk,kBeta);
       bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
       if(!goodKaon) continue;
       trackTest2->Fill(4);
       int charge = itrk->charge() * jtrk->charge();
       StKaonPion *kp = new StKaonPion(jtrk,itrk,j,i,vtx,b);
       hDecayLength->Fill(kp->decayLength());
       hKaonDca->Fill(kp->kaonDca());
       hPionDca->Fill(kp->pionDca());
       hDcaDaughters->Fill(kp->dcaDaughters());
       // cout<<"find pair"<<endl;
       hPointing->Fill(kp->pointingAngle() * kp->decayLength());
       StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
       StPicoTrack const* pion = picoDst->track(kp->pionIdx());
       if(kaon->charge()*pion->charge()>0 && kp->m()>1.813 && kp->m()<1.919)
       {
         // if(randomSample(kp->pt()))
         {
           if(kp->pt()>0 && kp->pt()<1)
             ntpD0BackgroundSameSignPt01->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>1 && kp->pt()<2)
             ntpD0BackgroundSameSignPt12->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>2 && kp->pt()<3)
             ntpD0BackgroundSameSignPt23->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>3 && kp->pt()<5)
             ntpD0BackgroundSameSignPt35->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>5 && kp->pt()<10)
             ntpD0BackgroundSameSignPt510->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
         }
       }
       if(kaon->charge()*pion->charge()<0 && ((kp->m()>1.761 && kp->m()<1.813) || (kp->m()>1.919 && kp->m()<1.971)))
       {
         // if(randomSample(kp->pt()))
         {
           if(kp->pt()>0 && kp->pt()<1)
             ntpD0BackgroundSideBandPt01->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>1 && kp->pt()<2)
             ntpD0BackgroundSideBandPt12->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>2 && kp->pt()<3)
             ntpD0BackgroundSideBandPt23->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>3 && kp->pt()<5)
             ntpD0BackgroundSideBandPt35->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>5 && kp->pt()<10)
             ntpD0BackgroundSideBandPt510->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
         }
       }
       // }
       tuplePair->Fill(kp->m(),kp->pt(),kp->kaonDca(),kp->pionDca(),kp->dcaDaughters(),kp->pointingAngle(), kp->decayLength(),charge,fabs(itrk->gPt()-jtrk->gPt()));
       if((charge=isD0Pair(kp))!=0 )
       {
         cout<<"find good pair"<<endl;
         if(charge>0)
           hD0MassLike->Fill(kp->pt(),kp->m());
         if(charge==-1)
           hD0MassUnlike->Fill(kp->pt(),kp->m());
       }
       delete kp;
   }
}
//
// TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
// cout<<"number of pairs = "<<aKaonPion->GetEntries()<<endl;
//
// for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
   // {
   //    // this is an example of how to get the kaonPion pairs and their corresponsing tracks
   //    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
   //    // if(!isGoodPair(kp)) continue;
   //
   //    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
   //    StPicoTrack const* pion = picoDst->track(kp->pionIdx());
   //
   //    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
   //    if (!isTpcPion(pion)) continue;
   //    bool tpcKaon = isTpcKaon(kaon,&vtx);
   //    float kBeta = getTofBeta(kaon,&vtx);
   //    bool tofAvailable = kBeta>0;
   //    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
   //    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
   //    if(!goodKaon) continue;
   //
   //    int charge = kaon->charge() * pion->charge();
   //    hDecayLength->Fill(kp->decayLength());
   //    hKaonDca->Fill(kp->kaonDca());
   //    hPionDca->Fill(kp->pionDca());
   //    hDcaDaughters->Fill(kp->dcaDaughters());
   //    hPointing->Fill(kp->pointingAngle() * kp->decayLength());
   //    if((charge=isD0Pair(kp))!=0 )
   //    {
   //      if(charge>0)
   //        hD0MassLike->Fill(kp->pt(),kp->m());
   //      if(charge==-1)
   //        hD0MassUnlike->Fill(kp->pt(),kp->m());
   //    }
   // }

   return kStOK;
}


int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  TLorentzVector d0Lorentz;
  d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
  if(fabs(d0Lorentz.Rapidity())>1.) return 0;
  bool pairCuts = false;
  // if(kp->pt()<1)
  // {
  //   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
  //     kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
  //     kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  // }
  // else if(kp->pt()<2)
  // {
  //   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
  //     kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
  //     kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  // }
  // else if(kp->pt()<3)
  // {
  //   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
  //     kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
  //     kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  // }
  // else if(kp->pt()<5)
  // {
  //   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
  //     kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
  //     kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  // }
  // else 
  // {
  //   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
  //     kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
  //     kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  // }
   pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.005 &&
      kp->pionDca() > 0.01 && kp->kaonDca() > 0.01 &&
      kp->dcaDaughters() < 0.01 && kp->decayLength()>0.015;  
  
  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  return trk->gPt() > 0.2&& trk->nHitsFit() >= 20;
}
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
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
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


  //return  trk->nHitsFit() >= mycuts::nHitsFit;
//-----------------------------------------------------------------------------

// bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const& trk, StThreeVectorF const& pVtx, double dca_sign) const {
//    int tofIndex = trk.bTofPidTraitsIndex();
//    bool TofMatch = kFALSE;
//
//    StPicoBTofPidTraits* tofPidTraits;
//    if (tofIndex >= 0)  tofPidTraits = mPicoDstMaker->picoDst()->btofPidTraits(tofIndex); //GNX
//    if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
//    bool vzSel = false; if (fabs(mPicoEvent->primaryVertex().z()) < 4.0) vzSel = true;
//
//    return trk.gPt()      >= anaCuts::minPt  &&
//       trk.nHitsFit() >= anaCuts::nHitsFit &&
//       fabs(trk.gMom(pVtx, mPicoEvent->bField()).pseudoRapidity()) <= anaCuts::eta
//       && vzSel;
// }
//
// bool StPicoD0AnaMaker::isGoodTrack_hft(StPicoTrack const& trk, StThreeVectorF const& pVtx, double dca_sign) const {
//    int tofIndex = trk.bTofPidTraitsIndex();
//    bool TofMatch = kFALSE;
//    StPicoBTofPidTraits* tofPidTraits;
//    if (tofIndex >= 0)  tofPidTraits = mPicoDstMaker->picoDst()->btofPidTraits(tofIndex); //GNX
//    if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
//    bool vzSel = false; if (fabs(mPicoEvent->primaryVertex().z()) < 4.0) vzSel = true;
//
//    return (!anaCuts::requireHFT || trk.isHFTTrack()) &&
//       trk.gPt()      >= anaCuts::minPt  &&
//       trk.nHitsFit() >= anaCuts::nHitsFit &&
//       fabs(trk.gMom(pVtx, mPicoEvent->bField()).pseudoRapidity()) <= anaCuts::eta
//       && TofMatch && vzSel;
// }
//
// bool StPicoD0AnaMaker::isGoodTrack_tpc(StPicoTrack const& trk, StThreeVectorF const& pVtx, double dca_sign) const {
//    int tofIndex = trk.bTofPidTraitsIndex();
//    bool TofMatch = kFALSE;
//    StPicoBTofPidTraits* tofPidTraits;
//    if (tofIndex >= 0)  tofPidTraits = mPicoDstMaker->picoDst()->btofPidTraits(tofIndex); //GNX 
//    if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
//    bool vzSel = false; if (fabs(mPicoEvent->primaryVertex().z()) < 4.0) vzSel = true; 
//
//    return trk.gPt()      >= anaCuts::minPt  &&
//       trk.nHitsFit() >= anaCuts::nHitsFit &&
//       fabs(trk.gMom(pVtx, mPicoEvent->bField()).pseudoRapidity()) <= anaCuts::eta
//       && TofMatch && vzSel;
// }
//
// bool StPicoD0AnaMaker::isGoodTrack_basic(StPicoTrack const& trk, StThreeVectorF const& pVtx, double dca_sign) const {
//    return trk.gPt()      >= anaCuts::minPt  &&
//       trk.nHitsFit() >= anaCuts::nHitsFit &&
//       fabs(trk.gMom(pVtx, mPicoEvent->bField()).pseudoRapidity()) <= anaCuts::eta;
// }
//
// int StPicoD0AnaMaker::getCentBin(float ntrk) {
//    int centb = -1;
//    for (int ic=0;ic<anaCuts::NCENT;ic++) {
//       if (ntrk > anaCuts::centCuts[ic]) {centb=ic; break;}
//    }
//         return centb;
// }
//
