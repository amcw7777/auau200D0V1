//
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
//
//
StChain *chain;
void readPicoDst(const Char_t *inputFile="test.list", const Char_t *outputFile="test.root")
{
  Int_t nEvents = 10000000;
  //Load all the System libraries
  string SL_version = "SL16d";
  string env_SL = getenv("STAR");
  if (env_SL.find(SL_version) == string::npos)
  {
    cout << "Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..." << endl;
    exit(1);
  }

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StBTofUtil");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StPicoPrescales");
  gSystem->Load("StPicoD0EventMaker");
  gSystem->Load("StMyAnalysisMaker");
  gSystem->Load("StmZDCSMDevp");
  chain = new StChain();

  cout<<inputFile<<endl;
  StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
   // StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDst");

  StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("ana",picoMaker,outputFile);
	// anaMaker->setRunEnergyAndListDir(run,energy,ListDir);            

  chain->Init();
  cout<<"chain->Init();"<<endl;
  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;
  for (Int_t i=0; i<nEvents; i++){
  // for (Int_t i=0; i<nEvents; i+=5){
    if(i%100==0)
      cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) { cout << "Bad return code!" << iret << endl; break;}

    total++;

  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  delete chain;


}
