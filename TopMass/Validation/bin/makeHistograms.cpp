//SYS
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <utility>
//ROOT
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TTree.h"
#include "TAxis.h"
#include "TPaveStats.h"
#include "TMath.h"
//#include "TRandom2.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TParameter.h"
//USER
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "interface/tdrStyle.h"

using namespace std;


void makeHistograms( string  myInFile  ){

  //Set TDR style
  TDRStyle tdrSTY; 
  tdrSTY.setTDRStyle();

  TFile * inFile = new TFile( myInFile.c_str() );
  TTree * psexpTree = (TTree*)inFile->Get("psexpTree");
  double mem_mass;
  double mem_mass_err, mem_jes, mem_jes_err, mem_corr, 
    mem_jes_intErr, mem_jes_err_intErr, mem_mass_intErr, mem_mass_err_intErr; 
  psexpTree->SetBranchAddress("mem_mass", &mem_mass);
  psexpTree->SetBranchAddress("mem_mass_err",&mem_mass_err);
  psexpTree->SetBranchAddress("mem_mass_err_intErr",&mem_mass_err_intErr);
  psexpTree->SetBranchAddress("mem_mass_intErr",&mem_mass_intErr);
  psexpTree->SetBranchAddress("mem_jes", &mem_jes);
  psexpTree->SetBranchAddress("mem_jes_err",&mem_jes_err);
  psexpTree->SetBranchAddress("mem_jes_err_intErr",&mem_jes_err_intErr);
  psexpTree->SetBranchAddress("mem_jes_intErr",&mem_jes_intErr);
  psexpTree->SetBranchAddress("mem_corr",&mem_corr);
  //  psexpTree->SetBranchAddress("mem_prob",&mem_prob);

  //Make tree to hold genral info on pseudo-exp. Important for resampling corrections.
  TTree * psexpInfoTree = (TTree*)inFile->Get("psexpInfoTree");
  int bkgnumEvents, signumEvents, psexp, signumEvtPool, bkgnumEvtPool;
  double trueJES, trueMass;
  psexpInfoTree->SetBranchAddress("signalEvents",&signumEvents);
  psexpInfoTree->SetBranchAddress("bkgEvents",&bkgnumEvents);
  psexpInfoTree->SetBranchAddress("numPsexp",&psexp);
  psexpInfoTree->SetBranchAddress("signalPool",&signumEvtPool);
  psexpInfoTree->SetBranchAddress("bkgPool",&bkgnumEvtPool);
  psexpInfoTree->SetBranchAddress("trueJES",&trueJES);
  psexpInfoTree->SetBranchAddress("trueMass",&trueMass);

  TParameter<Int_t>* mem_npes = (TParameter<Int_t>*)inFile->Get("mem_npes");
 
  if(psexpInfoTree->GetEntries() != 1) {return;}
  psexpInfoTree->GetEntry(0);
 
  // psexpTree->Draw("mem_mass");

  string myOutFile = myInFile.substr(0,myInFile.length()-5)+"_psExp.root"; 
  TFile * calF = new TFile(myOutFile.c_str() , "recreate");
  //Define useful histograms
  //NOTE TO USER: You must define the histogram boundaries to suit the results of your pseudoexperiments.
  //If your boundaries are too narrow, however, this macro will seg-fault.

  stringstream sstr;
  sstr.str("");
  sstr<<"Mass_"<<trueMass;
  TCanvas * cmass = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hmass = new TH1F("massH", "Mtop Measurement from MEM", 60, trueMass, trueMass + 6.0);
  (hmass->GetXaxis())->SetTitle("Mtop [GeV]");  (hmass->GetYaxis())->SetTitle("#"); (hmass->GetXaxis())->SetNdivisions(506);

 
  sstr.str("");
  sstr<<"MassErr_"<<trueMass;
  TCanvas * cmasserr = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hmasserr = new TH1F("masserrH", "Mtop stat error from MEM",  60, 0.1, 0.3);
  (hmasserr->GetXaxis())->SetTitle("Mtop Stat. error [GeV]");  (hmasserr->GetYaxis())->SetTitle("#"); (hmasserr->GetXaxis())->SetNdivisions(508);
  
  sstr.str("");
  sstr<<"MassPull_"<<trueMass;
  TCanvas * cmasspull = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hmasspull = new TH1F("masspullH", "Pull Distribution - Mtop form MEM", 200, 0, 40);
  (hmasspull->GetXaxis())->SetTitle("pull");  (hmasspull->GetYaxis())->SetTitle("#");
  
  sstr.str("");
  sstr<<"MassRelErr_"<<trueMass;
  TCanvas * cmassrelerr = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hmassrelerr = new TH1F("massrelerrH", "Mtop relative stat error from MEM", 60, 0.000, 0.009);
  (hmassrelerr->GetXaxis())->SetTitle("Mtop rel. stat. error");  (hmassrelerr->GetYaxis())->SetTitle("#");
    
  // sstr.str("");
  // sstr<<"MassProb_"<<trueMass;
  // TCanvas* cmassprob = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600,600);
  // TH1F* hmassprob = new TH1F("massprobH","Mtop prob from MEM",100,0,3);
  // (hmassprob->GetXaxis())->SetTitle("Mtop Probability"); (hmassprob->GetYaxis())->SetTitle("#");

  sstr.str("");
  sstr<<"MassErr2D_"<<trueMass;
  TCanvas* cmasserr2D = new TCanvas(sstr.str().c_str(),sstr.str().c_str(),600,600);
  TH2F* hmasserr2D = new TH2F("masserr2D","Mtop vs Uncert",60, trueMass - 1.0, trueMass + 2.0,60, 0.1, 0.3);
  (hmasserr2D->GetXaxis())->SetTitle("Mtop [GeV]"); (hmasserr2D->GetYaxis())->SetTitle("Mtop Stat. err [GeV]");

  sstr.str("");
  sstr<<"JES_"<<trueJES;
  TCanvas * cjes = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hjes = new TH1F("jesH", "JES Measurement from MEM", 341, 0.7995, 1.1405);
  (hjes->GetXaxis())->SetTitle("JES");  (hjes->GetYaxis())->SetTitle("#");
  
  sstr.str("");
  sstr<<"JESErr_"<<trueJES;
  TCanvas * cjeserr = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hjeserr = new TH1F("jeserrH", "JES stat error from MEM", 40, 0.00105, 0.00505);
  (hjeserr->GetXaxis())->SetTitle("JES Stat. error [GeV]");  (hjeserr->GetYaxis())->SetTitle("#");
  
  sstr.str("");
  sstr<<"JESPull_"<<trueJES;
  TCanvas * cjespull = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hjespull = new TH1F("jespullH", "Pull Distribution - JES form MEM", 21, -4.1, 4.1);
  (hjespull->GetXaxis())->SetTitle("pull"); (hjespull->GetYaxis())->SetTitle("#");
  
  sstr.str("");
  sstr<<"JESRelErr_"<<trueJES;
  TCanvas * cjesrelerr = new TCanvas(sstr.str().c_str(),sstr.str().c_str(), 600, 600);
  TH1F * hjesrelerr = new TH1F("jesrelerrH", "JES relative stat error from MEM", 40, 0.0020, 0.0060);
  (hjesrelerr->GetXaxis())->SetTitle("JES rel. stat. error");  (hjesrelerr->GetYaxis())->SetTitle("#");  

  TH2F * hmassjes = new TH2F("massJESH", "JES vs Mtop Measurement from MEM", 41, trueMass - 8.2, trueMass + 8.2, 341, 0.7995, 1.1405);

  double numInd= min ((double)signumEvtPool/(double)signumEvents, (double)bkgnumEvtPool/(double)bkgnumEvents);

  
  for (int i = 0 ; i< psexpTree->GetEntries(); i++ ){
   
    psexpTree->GetEntry(i);

    hmass->Fill(mem_mass);  
    double tmp =sqrt( pow(mem_mass_err,2) + pow(mem_mass_intErr,2) + pow(mem_mass_err_intErr,2));
    //double tmp = mem_mass_err;
    if(tmp != 0.0) {
      hmasspull->Fill((double)(mem_mass-trueMass)/tmp);
      hmasserr->Fill(tmp);  hmassrelerr->Fill(tmp/mem_mass);
      hmasserr2D->Fill(mem_mass,tmp);
    }
    //    hmassprob->Fill(mem_prob);

    hjes->Fill(mem_jes);  
    tmp =sqrt( pow(mem_jes_err,2) + pow(mem_jes_intErr,2) + pow(mem_jes_err_intErr,2));
    if(tmp != 0.0) {
      hjespull->Fill((double)(mem_jes-trueJES)/tmp);
      hjeserr->Fill(tmp);  hjesrelerr->Fill(tmp/mem_jes);
    }
  }
  
  stringstream mystr;  

  //############################################## MASS ##########################################################
  //Mass distribution ###################
  double hmassmax =  hmass->GetBinContent(hmass->GetMaximumBin());
  double hmassmean  = hmass->GetMean(); 
  double hmassrms  = hmass->GetRMS(); 
  
  TLine l(hmassmean, 0, hmassmean, hmassmax);
  l.SetLineColor(kRed);
  if(  hmass->GetEntries()>0){
    hmass->Fit("gaus"); 
  
  TF1 * gaus_hmass = hmass->GetFunction("gaus");
  gaus_hmass->SetLineColor(kBlue);
  }

  TPaveStats ptmass(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptmass.SetName("massMeas");
  ptmass.SetBorderSize(0);
  ptmass.SetFillColor(kWhite);   ptmass.SetLineColor(0);
  ptmass.SetTextAlign(12);  ptmass.SetTextColor(kRed);
  mystr.str("");   mystr<<hmassmean<<" #pm "<<(hmassrms/sqrt(numInd));

  ptmass.AddText( mystr.str().c_str());
  cout << "finished mass fit" << endl; 
  //Error distribution ###################
  double hmasserrmax =  hmasserr->GetBinContent(hmasserr->GetMaximumBin());
  if(  hmasserr->GetEntries()>0){
    hmasserr->Fit("gaus");
  TF1 * gaus_hmasserr = hmasserr->GetFunction("gaus");
  gaus_hmasserr->SetLineColor(kBlue);
  }

  double hmasserrmean  = hmasserr->GetMean();  
  double hmasserrrms  = hmasserr->GetRMS(); 

  TLine lerr(hmasserrmean, 0, hmasserrmean, hmasserrmax);
  lerr.SetLineColor(kRed);

  TPaveStats ptmasserr(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptmasserr.SetName("massErrMeas");
  ptmasserr.SetBorderSize(0);
  ptmasserr.SetFillColor(kWhite);ptmasserr.SetLineColor(0);
  ptmasserr.SetTextAlign(12); ptmasserr.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean "<<hmasserrmean<<" #pm "<<(hmasserrrms/sqrt(numInd));

  ptmasserr.AddText( mystr.str().c_str());
  cout << "finished mass err fit" << endl;
  //RelativeError distribution ###################
  double hmassrelerrmax =  hmassrelerr->GetBinContent(hmassrelerr->GetMaximumBin());
 if(  hmassrelerr->GetEntries()>0){
  hmassrelerr->Fit("gaus");
 
  TF1 * gaus_hmassrelerr = hmassrelerr->GetFunction("gaus");
  gaus_hmassrelerr->SetLineColor(kBlue);
 }
  double hmassrelerrmean  = hmassrelerr->GetMean();  
  double hmassrelerrrms  = hmassrelerr->GetRMS(); 

  TLine lrelerr(hmassrelerrmean, 0, hmassrelerrmean, hmassrelerrmax);
  lrelerr.SetLineColor(kRed);

  TPaveStats ptmassrelerr(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptmassrelerr.SetName("massErrMeas");
  ptmassrelerr.SetBorderSize(0);
  ptmassrelerr.SetFillColor(kWhite);ptmassrelerr.SetLineColor(0);
  ptmassrelerr.SetTextAlign(12); ptmassrelerr.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean "<<hmassrelerrmean<<" #pm "<<(hmassrelerrrms/sqrt(numInd));

  ptmassrelerr.AddText( mystr.str().c_str());
  cout << "finished relative mass error fit" << endl;
  //Pull distribution ###################
  double hmasspullmax =  hmasspull->GetBinContent(hmasspull->GetMaximumBin());
 if(  hmasspull->GetEntries()>0){
  hmasspull->Fit("gaus"); 
 
  TF1 * gaus_hmasspull = hmasspull->GetFunction("gaus");
  gaus_hmasspull->SetLineColor(kBlue);
 }
  double hmasspullmean  = hmasspull->GetMean(); 
  double hmasspullrms  =  hmasspull->GetRMS(); 
  //double hmasspullrmserror =  hmasspull->GetRMSError();

    TLine lpull(hmasspullmean, 0, hmasspullmean, hmasspullmax);
     lpull.SetLineColor(kRed);

  TPaveStats ptmasspull(0.5994803,0.5,0.9984164,0.6514842,"brNDC");
  ptmasspull.SetName("masspullMeas");
  ptmasspull.SetBorderSize(0);
  ptmasspull.SetFillColor(kWhite);ptmasspull.SetLineColor(0);
  ptmasspull.SetTextAlign(12); ptmasspull.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean pull: "<<hmasspullmean<<" #pm "<<(hmasspullrms/sqrt(numInd));
  ptmasspull.AddText( mystr.str().c_str());
  mystr.str("");  
  mystr<<"Pull width: "<<hmasspullrms<<" #pm "<<(hmasspullrms * sqrt(0.5*((1./((double)signumEvtPool+(double)bkgnumEvtPool))+( 1./((double)psexp-1.)))));
  ptmasspull.AddText( mystr.str().c_str());
  cout << "finished mass pull fit" << endl; 

  //############################################## JES ##########################################################


  //JES distribution #####################################

  double hjesmax =  hjes->GetBinContent(hjes->GetMaximumBin());
  double hjesmean  = hjes->GetMean(); 
  double hjesrms  = hjes->GetRMS(); 
 
  if(  hjesmean == 99.){
    hjes->Fit("gaus");   
 
    TF1 * gaus_hjes = hjes->GetFunction("gaus");
    gaus_hjes->SetLineColor(kBlue);
  }
 
  TLine lj(hjesmean, 0, hjesmean, hjesmax);
  lj.SetLineColor(kRed);

  TPaveStats ptjes(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptjes.SetName("jesMeas");
  ptjes.SetBorderSize(0);
  ptjes.SetFillColor(kWhite);   ptjes.SetLineColor(0);
  ptjes.SetTextAlign(12);  ptjes.SetTextColor(kRed);
  mystr.str("");   mystr<<hjesmean<<" #pm "<<(hjesrms/sqrt(numInd));

  ptjes.AddText( mystr.str().c_str());
 
  cout << "finished jes fit" << endl;
  //Error distribution ###################
  double hjeserrmax =  hjeserr->GetBinContent(hjeserr->GetMaximumBin());
  
  if(  hjesmean == 99.){
    hjeserr->Fit("gaus");
    TF1 * gaus_hjeserr = hjeserr->GetFunction("gaus");
    gaus_hjeserr->SetLineColor(kBlue);
  }
  double hjeserrmean  = hjeserr->GetMean();  
  double hjeserrrms  = hjeserr->GetRMS(); 
  
  TLine lerrj(hjeserrmean, 0, hjeserrmean, hjeserrmax);
  lerrj.SetLineColor(kRed);
  
  TPaveStats ptjeserr(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptjeserr.SetName("jesErrMeas");
  ptjeserr.SetBorderSize(0);
  ptjeserr.SetFillColor(kWhite);ptjeserr.SetLineColor(0);
  ptjeserr.SetTextAlign(12); ptjeserr.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean "<<hjeserrmean<<" #pm "<<(hjeserrrms/sqrt(numInd));
  
  ptjeserr.AddText( mystr.str().c_str());

  cout << "finished jes error fit" << endl; 
  //Relative Error distribution ###################
  double hjesrelerrmax =  hjesrelerr->GetBinContent(hjeserr->GetMaximumBin());
  if(  hjesmean == 99.){
    hjesrelerr->Fit("gaus");
  
    TF1 * gaus_hjesrelerr = hjesrelerr->GetFunction("gaus");
    gaus_hjesrelerr->SetLineColor(kBlue);
  }  
  
  double hjesrelerrmean  = hjesrelerr->GetMean();  
  double hjesrelerrrms  = hjesrelerr->GetRMS(); 
  
  TLine lrelerrj(hjeserrmean, 0, hjesrelerrmean, hjesrelerrmax);
  lrelerrj.SetLineColor(kRed);
  
  TPaveStats ptjesrelerr(0.5386566,0.6600824,0.8819273,0.7649317,"brNDC");
  ptjesrelerr.SetName("jesErrMeas");
  ptjesrelerr.SetBorderSize(0);
  ptjesrelerr.SetFillColor(kWhite);ptjesrelerr.SetLineColor(0);
  ptjesrelerr.SetTextAlign(12); ptjesrelerr.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean "<<hjesrelerrmean<<" #pm "<<(hjesrelerrrms/sqrt(numInd));

  ptjesrelerr.AddText( mystr.str().c_str());

  
  //pULL distribution ###################
  double hjespullmax =  hjespull->GetBinContent(hjespull->GetMaximumBin());
  double hjespullmean  = hjespull->GetMean();
  double hjespullrms  =  hjespull->GetRMS(); 

  TLine lpullj(hjespullmean, 0, hjespullmean, hjespullmax);
  lpullj.SetLineColor(kRed);
  
  if(  hjesmean == 99.){
    hjespull->Fit("gaus"); 
  
    TF1 * gaus_hjespull = hjespull->GetFunction("gaus");
    gaus_hjespull->SetLineColor(kBlue);
  }


  TPaveStats ptjespull(0.5994803,0.5,0.9984164,0.8014842,"brNDC");
  ptjespull.SetName("jespullMeas");
  ptjespull.SetBorderSize(0);
  ptjespull.SetFillColor(kWhite);ptjespull.SetLineColor(0);
  ptjespull.SetTextAlign(12); ptjespull.SetTextColor(kRed);
  mystr.str("");   mystr<<"Mean pull: "<<hjespullmean<<" #pm "<<(hjespullrms/sqrt(numInd));
  ptjespull.AddText( mystr.str().c_str());
  mystr.str("");   
  mystr<<"Pull width: "<<hjespullrms<<" #pm "<<(hjespullrms * sqrt(0.5*((1./((double)signumEvtPool+(double)bkgnumEvtPool))+( 1./((double)psexp-1.)))));  
  ptjespull.AddText( mystr.str().c_str());
  cout << "finished jes pull fit" << endl;   
  
  cjes->cd(); 
  hjes->Draw();lj.Draw("same"); 
  ptjes.Draw("same");
    
  cjeserr->cd(); 
  hjeserr->Draw();lerrj.Draw("same");  
  ptjeserr.Draw("same");
    
  cjespull->cd();     
  hjespull->Draw(); lpullj.Draw("same"); 
  ptjespull.Draw("same");
  
  cjesrelerr->cd(); 
  hjesrelerr->Draw(); lrelerrj.Draw("same");  
  ptjesrelerr.Draw("same");

  cmass->cd(); 
  hmass->Draw();l.Draw("same"); 
  ptmass.Draw("same");

  cmasserr->cd(); 
  hmasserr->Draw();lerr.Draw("same");  
  ptmasserr.Draw("same");

  cmasspull->cd();     
  hmasspull->Draw(); lpull.Draw("same"); 
  ptmasspull.Draw("same");

  cmassrelerr->cd();     
  hmassrelerr->Draw(); lrelerr.Draw("same"); 
  ptmassrelerr.Draw("same");

  //  cmassprob->cd();
  //  hmassprob->Draw();

  cmasserr2D->cd();
  hmasserr2D->Draw();

  cout << "finished drawing" << endl;

  calF->cd();
  cmass->Write();  cmasserr->Write();  cmasspull->Write();  cmassrelerr->Write(); 
  cjes->Write();  cjeserr->Write();  cjespull->Write();  cjesrelerr->Write(); 
  cmasserr2D->Write();

  hmassjes->Write();

  TTree* ctree = psexpInfoTree->CloneTree();
  ctree->Write();

  cout << "Number of Attempted PEs: " << mem_npes->GetVal() << endl;

  cout<<"Mass: "<<hmassmean<<" #pm "<<(hmassrms/sqrt(numInd))<<endl;
  cout<<"Error: "<<hmasserrmean<<" #pm "<<(hmasserrrms/sqrt(numInd))<<endl;
  cout<<"Rel. Error: "<<hmassrelerrmean<<" #pm "<<(hmassrelerrrms/sqrt(numInd))<<endl;
  cout<<"Mean pull: "<<hmasspullmean<<" #pm "<<(hmasspullrms/sqrt(numInd))<<endl;
  cout<<"Pull width: "<<hmasspullrms<<" #pm "<<(hmasspullrms * sqrt(0.5*((1./((double)signumEvtPool+(double)bkgnumEvtPool))+( 1./((double)psexp-1.)))))<<endl;
  

  cout<<"Jes: "<<hjesmean<<" #pm "<<(hjesrms/sqrt(numInd))<<endl;
  cout<<"Error: "<<hjeserrmean<<" #pm "<<(hjeserrrms/sqrt(numInd))<<endl;
  cout<<"Rel. Error: "<<hjesrelerrmean<<" #pm "<<(hjesrelerrrms/sqrt(numInd))<<endl;
  cout<<"Mean pull: "<<hjespullmean<<" #pm "<<(hjespullrms/sqrt(numInd))<<endl;
  cout<<"Pull width: "<<hjespullrms<<" #pm "<<(hjespullrms * sqrt(0.5*((1./((double)signumEvtPool+(double)bkgnumEvtPool))+( 1./((double)psexp-1.)))))<<endl;  

  //Clean up
 
  delete calF; 
   return;
}


int main(int argc,  char* argv[]){

  if(argc==2){
    //    std::string  myInFileFolder = (string)argv[1];
    makeHistograms((string)argv[1]);
    return 0;
  }else{
    cout<<"ERROR! \n \t Usage: \n \t signalOnly [Input File] "<<endl;
    return -1;
  }
  
  
}



