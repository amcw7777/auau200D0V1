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
	mRun    =0;
	mEnergy =0;
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

	// badList.clear();
	// runList.clear();

	// if(!readRunList())return kStFatal;
	// if(!readBadList())return kStFatal;

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

  d0MassPhiEta = new TH3D("d0MassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",500,1.6,2.1,4,0,PI,4,-1,1);
  d0BarMassPhiEta = new TH3D("d0BarMassPhiEta",";D^{0} mass (GeV/c^{2});#phi_{D^{0}}-#psi_{ZDC};#eta",500,1.6,2.1,4,0,PI,4,-1,1);
  zdcPsi = new TH1D("zdcPsi",";#psi_{ZDC}",1000,0,twoPI);
  zdcPsi_corr = new TH1D("zdcPsi_corr",";#psi_{ZDC}",1000,0,twoPI);
  pionV1Plus = new TProfile("pionV1Plus","",48,-1.2,1.2);
  pionV1Minus = new TProfile("pionV1Minus","",48,-1.2,1.2);
  pionV1Plus->Sumw2();
  pionV1Minus->Sumw2();
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
  zdcPsi->Write();
  zdcPsi_corr->Write();
  pionV1Plus->Write();
  pionV1Minus->Write();
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
	const int  evtID    = event->eventId();
	const int refMult  = event->grefMult();
  if(!(isGoodEvent(event)))
  { 
    // LOG_WARN << " Not Min Bias! Skip! " << endm;
    return kStWarn;
  }

	mTempRunId = runID;  
  mTempEvtId = evtID;
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }

  StThreeVectorF vtx = event->primaryVertex();
  float b = event->bField();
  double vpdvz = mPicoDst->event()->vzVpd();
  double vz = vtx.z();
  // if( fabs(vz)>6 && fabs(vz-vpdvz)>3 ) 
  //   return kStWarn;
	//StRefMultCorr
	if ( runID != mPrevRunId ) {
		mRefMultCorr->init(runID);
		mPrevRunId = runID;
		cout << "reset mPrevRunId = " << mPrevRunId << endl;
	}
	//mRefMultCorr->initEvent(refMult,VertexZ);
	mRefMultCorr->initEvent(refMult,vz,event->ZDCx());

	int centrality = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
	double mWght    = mRefMultCorr->getWeight();

	if( centrality<0||centrality>=(nCent-1)) return kStOK;
	if( centrality<2||centrality>6) return kStOK; // 10 - 60 % centrality
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
	miniZDCSMD->calibrateZDCSMDevp();
	if(!(miniZDCSMD->ZDCSMDgoodEvent()))return kStOK;
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

  vector<int> kaonIndex;
  vector<int> pionIndex;
  kaonIndex.clear();
  pionIndex.clear();
  for(int i=0;i<mPicoDst->numberOfTracks();i++)
  {
    StPicoTrack const* itrk = mPicoDst->track(i);
    if(!isGoodTrack(itrk))  continue;
    if(!(itrk->isHft())) continue;
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
  for(int i=0;i<pionIndex.size();i++)
  {
    if(centrality!=4) continue;
    StPicoTrack const* itrk = mPicoDst->track(pionIndex[i]);
		StThreeVectorF mom= itrk->gMom();
		double p    = mom.mag();
		double pt   = mom.perp();
		double eta  = mom.pseudoRapidity();
		double phi  = mom.phi();
		if(phi<0.0) phi += twoPI;
		double mcos1=cos(phi);
		double msin1=sin(phi);
		double v1ZDC1F = mcos1*cos(1.*mZDC1Event_PsiF) + msin1*sin(1.*mZDC1Event_PsiF);
    int charge = itrk->charge();
    if(charge > 0)
      pionV1Plus->Fill(eta,v1ZDC1F,mWght);
    else
      pionV1Minus->Fill(eta,v1ZDC1F,mWght);
  }


  for(int i=0;i<pionIndex.size();i++)
  {
    StPicoTrack const* itrk = mPicoDst->track(pionIndex[i]);
    for(int j=0;j<kaonIndex.size();j++)
    {
      StPicoTrack const* jtrk = mPicoDst->track(kaonIndex[j]);
       int charge = itrk->charge() * jtrk->charge();
       StKaonPion *kp = new StKaonPion(jtrk,itrk,j,i,vtx,b);
       // if(kp->pt()<1.5) continue;// require D0 pT > 1.5 GeV/c
       if(charge>0) continue;//only unlike-sign pairs
       StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
       StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());

       int isBar = isD0Pair(kp);
       if(isBar != 0)
         cout<<"find D0 !"<<endl;
       if(isBar>0)
         d0MassPhiEta->Fill(kp->m(),kp->phi()-mZDC1Event_PsiF,kp->eta(),mWght);
       if(isBar<0)
         d0BarMassPhiEta->Fill(kp->m(),kp->phi()-mZDC1Event_PsiF,kp->eta(),mWght);
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
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}

int StMyAnalysisMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
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
  // int charge = kaon->charge() * pion->charge();
  // if(charge>0)
  //   charge = kaon->charge()>0 ? 1:2;
  if(pairCuts)
    return -1 * kaon->charge();
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

bool StMyAnalysisMaker::removeBadID(int runnumber)
{
	for (std::vector<int>::iterator it=badList.begin(); it!=badList.end(); ++it) 
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
	while ( inrun >> runid ) { badList.push_back(runid); }
	inrun.close();
	sort(badList.begin(),badList.end());

	vector<int>::iterator it;
	it = std::unique (badList.begin(), badList.end());
	badList.resize( std::distance(badList.begin(),it) );

	cout <<"badrun list :" <<inf.Data() << " loaded." << endl;
	cout <<"Total       :" <<badList.size()<< " bad runs. "<< endl;
	return kTRUE;
}
//------------------------------------------------------------------------------
bool StMyAnalysisMaker::isGoodEvent(StPicoEvent const*event) const
{	
  const int  runID    = event->runId();
	const int  evtID    = event->eventId();
  bool isBadRun = false;
  std::array<unsigned int, 87> const badRunN={17062049, 17063007, 17063008, 17063039, 17065037, 17065041, 17065042, 17065054, 17065055, 17065056, 17075018, 17075019, 17075024, 17075025, 17075026, 17075099, 17075101, 17075102, 17075105, 17075106, 17075107, 17075108, 17075110, 17075111, 17076001, 17076002, 17076003, 17076011, 17098017, 17098035, 17098036, 17098037, 17098038, 17098039, 17098040, 17098043, 17098045, 17098046, 17098047, 17099031, 17099035, 17099042, 17103041, 17103061, 17104001, 17104002, 17104003, 17104004, 17104005, 17104006, 17104007, 17104008, 17104009, 17104011, 17104012, 17104013, 17104014, 17104015, 17104016, 17104017, 17109051, 17109052, 17109055, 17109056, 17109057, 17109058, 17109059, 17109060, 17109061, 17109068, 17109069, 17109070, 17114024, 17114025, 17114032, 17114033, 17114034, 17115038, 17116065, 17117003, 17119016, 17120009, 17121002, 17121029, 17122039, 17124039, 17124041, };
  for(auto badRun: badRunN)
  { 
    if(badRun == runID)  isBadRun = true;
  }
  /////////////////////////////////////////////////////////////////////
  bool isGoodTrigger = false;
  std::array<unsigned int, 6> const triggers = { 520001,    // vpdmb-5-p-nobsmd-hlt 
    520011,    // vpdmb-5-p-nobsmd-hlt 
    520021,    // vpdmb-5-p-nobsmd-hlt 
    520031,    // vpdmb-5-p-nobsmd-hlt 
    520041,// vpdmb-5-p-nobsmd-hlt 
    520051// vpdmb-5-p-nobsmd-hlt 
  };
  for(auto trg: triggers)
  {
    if(event->isTrigger(trg)) isGoodTrigger = true;
  }
  /////////////////////////////////////////////////////////////////////

	StThreeVectorF Vertex3D=event->primaryVertex();
	const double VertexX = Vertex3D.x(); 
	const double VertexY = Vertex3D.y(); 
	const double VertexZ = Vertex3D.z(); 
	const double vpdVz   = event->vzVpd();
  bool isGoodVtx = (fabs(VertexZ) < 6) && (sqrt(pow(VertexX,2.)+pow(VertexY,2.))<2.0) && (fabs(VertexZ-vpdVz)<3.);

  /////////////////////////////////////////////////////////////////////
  bool isDuplicateEvent = false;
	if(mTempRunId==runID && mTempEvtId==evtID)
    isDuplicateEvent = true;

  return !isBadRun && isGoodTrigger && isGoodVtx && !isDuplicateEvent;
}

