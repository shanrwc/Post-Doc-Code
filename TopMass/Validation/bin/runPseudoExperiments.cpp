//SYS
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <utility>
//ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMarker.h"
#include "TParameter.h"
//USER
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "Mass/CommonTools/acceptanceTools.h"
#include "Mass/CommonTools/normalizationTools.h"
//#include "Mass/CommonTools/valid_sampleSelector.h"
#include "Mass/CommonTools/calib_sampleSelector.h"
#include "interface/tdrStyle.h"

using namespace RooFit;
using namespace std;

Double_t Likelihood2D(Double_t *x, Double_t *par){

  Double_t xx =x[0];
  Double_t yy =x[1]; 

  Double_t a1 = par[0]; Double_t b1 = par[1];
  Double_t a2 = par[2]; Double_t b2 = par[3];
  Double_t c_c = par[4];
  Double_t d_d = par[5];

  Double_t f = (0.5/(1-pow(c_c,2))) * (pow((xx-a1)/b1,2) + pow((yy-a2)/b2,2) - (2.*c_c*(xx-a1)*(yy-a2)/(b1*b2)))  + d_d;

  return f;
}

//1D Parabola for mass Only fit
Double_t Likelihood1D(Double_t *x, Double_t *par){   
  Float_t xx =x[0];
  Float_t a1 = par[0]; Float_t b1 = par[1];  Float_t c1 = par[2];
  Double_t f =  pow((xx-a1)/b1,2)/2  + c1;
  return f;
}



Double_t LikeProjection(Double_t *x, Double_t *par){

  Double_t xx =x[0];
  Double_t yy =par[6]; 

  Double_t a1 = par[0]; Double_t b1 = par[1];
  Double_t a2 = par[2]; Double_t b2 = par[3];
  Double_t c_c = par[4];
  Double_t d_d = par[5];

  Double_t f = (0.5/(1-pow(c_c,2))) * (pow((xx-a1)/b1,2) + pow((yy-a2)/b2,2) - (2.*c_c*(xx-a1)*(yy-a2)/(b1*b2)))  + d_d;

  return f;
}


void runPseudotest( string  myInFile, int psexp  ){

  //Set TDR style
  TDRStyle tdrSTY; 
  tdrSTY.setTDRStyle();

  // //Open event signal likelihood file
  //  sampleSelector * sigSelector = new sampleSelector();
  stringstream sstr;
  sstr.str("");

  //The sampleSelector class reads the config file, accesses the input files,
  //and creates the ensemble of events used in a pseudoexperiment.
  sampleSelector selector;

  //The initialize function MUST be called before pseudoexperiments are performed;
  //it causes the sampleSelector to read in the input file and set up objects for
  //pseudoexperiments.  Each pair in the vector contains first the number of events
  //from a given sample in the pseudoexperiment and second the number of events in
  //the pool from that sample.  The signal sample is always in the zeroth slot;
  //additional background samples are in the order listed in the drive file.
  //It also prints lots of information to screen and returns an error code--
  //Pay attention to the output and do not run pseudoexperiments if the 
  //intialize function returns something besides 0!
  //1 - config file not found or not openable
  //2 - input file not found
  //3 - input file does not contain likelihood tree
  //4 - input signal file does not contain truth information
  //5 - input files' mass hypotheses don't match
  //6 - input files' JES hypothese don't match
  //7 - background sample fractions excede 1.0
  //8 - signal fraction != 1 - sum(background fractions)
  //9 - input samples using different event weights
  //10 - input file has no mass or jes hypotheses compatible with bounds
  //11 - input file name and trigger string are not an acceptable match
  vector<pair<int,int> > numEvtPools;
  int status =  selector.initialize(myInFile,numEvtPools);
  if (status != 0){
    cout << "[ERROR]: sampleSelector failed to initialize!" << endl;
    cout << "[ERROR]: Status = " << status << endl;
    return;
  }
  
  
  //Truth information is taken from the signal sample; this function must be called to 
  //properly set things up for pseudoexperiments.
  //You also get useful output from it to make sure you are running want you want to run.
  double trueMass=0.0;
  double trueJES= 0.0;
  if( selector.getTruthInfo(trueMass,trueJES) == -1) {
    cout<<"[ERROR]: Truth information not found! EXIT"<<endl;
    return;
  }
  cout << "True Mass = " << trueMass << endl;
  cout << "True JES = " << trueJES << endl;


  set<double> set_mass=  selector.getMassHypos();
  set<double> set_jes=  selector.getJESHypos();
  int numMassHyp = set_mass.size();
  int numJESHyp = set_jes.size();
  
  //Create output file
  sstr.str("");
  sstr<<"MassJES_Bkg_M_"<<trueMass<<"_JES_"<<trueJES<<"_"<<selector.getSystematicLabel()<<"_"<<selector.getSeed()<<".root";
   
  TFile * outputFile = new TFile(sstr.str().c_str(), "recreate");
  //Make tree to hold results from pseudo-exp.
  TTree * psexpTree = new TTree("psexpTree","Pseudo-experiment Tree");
  double mem_mass, mem_mass_err, mem_jes, mem_jes_err, mem_corr, mem_pmin;

  psexpTree->Branch("mem_mass", &mem_mass, "mem_mass/D");
  psexpTree->Branch("mem_mass_err",&mem_mass_err,"mem_mass_err/D");
  psexpTree->Branch("mem_jes",&mem_jes,"mem_jes/D");
  psexpTree->Branch("mem_jes_err",&mem_jes_err,"mem_jes_err/D");
  psexpTree->Branch("mem_corr",&mem_corr," mem_corr/D");
  psexpTree->Branch("mem_pmin",&mem_pmin,"mem_pmin/D");

  double mem_mass_intErr, mem_mass_err_intErr, mem_jes_intErr, mem_jes_err_intErr, mem_corr_intErr; 

  psexpTree->Branch("mem_mass_intErr", &mem_mass_intErr, "mem_mass_intErr/D");
  psexpTree->Branch("mem_mass_err_intErr",&mem_mass_err_intErr,"mem_mass_err_intErr/D");
  psexpTree->Branch("mem_jes_intErr",&mem_jes_intErr,"mem_jes_intErr/D");
  psexpTree->Branch("mem_jes_err_intErr",&mem_jes_err_intErr,"mem_jes_err_intErr/D");
  psexpTree->Branch("mem_corr_intErr",&mem_corr_intErr," mem_corr_intErr/D");

  //Make a couple histograms
  TH1D* edm_pass = new TH1D("edm_pass","EDM of Good Fits",100,0,0.0001);
  TH1D* edm_fail = new TH1D("edm_fail","EDM of Bad Fits",100,0,10); 

  //Make tree to hold genral info on pseudo-exp. Important for resampling corrections.
  TTree * psexpInfoTree = new TTree("psexpInfoTree","Pseudo-experiment info for Resampling");
  int bkgnumEvents,signumEvents, signumEvtPool,bkgnumEvtPool ;
  psexpInfoTree->Branch("signalEvents",&signumEvents,"signumEvents/I");
  psexpInfoTree->Branch("bkgEvents",&bkgnumEvents,"bkgnumEvents/I");
  psexpInfoTree->Branch("numPsexp",&psexp,"psexp/I");
  psexpInfoTree->Branch("signalPool",&signumEvtPool,"signumEvtPool/I");
  psexpInfoTree->Branch("bkgPool",&bkgnumEvtPool,"bkgnumEvtPool/I");
  psexpInfoTree->Branch("trueJES",&trueJES,"trueJES/D");
  psexpInfoTree->Branch("trueMass",&trueMass,"trueMass/D");

  signumEvents = numEvtPools[0].first;  
  signumEvtPool = numEvtPools[0].second;

  if( numEvtPools.size() >1 ){
    bkgnumEvents = numEvtPools[1].first;
    bkgnumEvtPool =  numEvtPools[1].second;
  }else{
    bkgnumEvtPool = 0.0;
    bkgnumEvents = 0.0;
  }
  
  psexpInfoTree->Fill();
  
  stringstream selectEvents;
  stringstream dataCut;

  //run pseudoexperiments
  int skip=0;

  //  Double_t pe_nll_start=10000000.;

  for (int p = 0 ; p < psexp+skip; p++){
    //Get sample
    RooDataSet* sample = selector.getEnsemble();
   
    cout<<"Starting pseudo-exp: "<<p-skip<<endl;
           
    Double_t pe_nll[500];     Double_t pe_nll_err[500];  
    Double_t pe_mass[500];    Double_t pe_jes[500];
    Double_t  ex[500]; 
     
    for( int y = 0 ; y<500; y++){
      pe_nll[y]=0.; pe_nll_err[y]=0.;
      pe_mass[y]=0.;pe_jes[y]=0.;
      ex[y]=0.;
    }
     
    selector.reduce2Arrays(sample, pe_mass, pe_jes, pe_nll, pe_nll_err );

    // //Rescalling for graphical reasons
    // for (int i = 0; i<(numMassHyp*numJESHyp); i++  ){
    //   pe_nll[i] -= pe_nll_start;   
    // }

    //Use Minuit2 --- it's better
    TVirtualFitter::SetDefaultFitter("Minuit2");

    //Define graphs for...
    //...2D visualization
    sstr.str(""); 
    sstr<<"canMassJES"<<p;  
    TCanvas * canMassJES = new TCanvas(sstr.str().c_str(),"", 600, 600);
    TGraph2DErrors * grMassJES = new TGraph2DErrors(numMassHyp*numJESHyp, &pe_mass[0], &pe_jes[0],  &pe_nll[0],  &ex[0],  &ex[0],  &pe_nll_err[0] );
    //A bit of rescaling
    double minLike = grMassJES->GetZmin();
    for (int k=0;k<grMassJES->GetN(); k++) { grMassJES->GetZ()[k] -= minLike;}

    //Define graphs for...
    //...1D visualization
    sstr.str(""); 
    sstr<<"canMass"<<p;  
    TCanvas * canMass = new TCanvas(sstr.str().c_str(),"MTop", 600); 
    TGraphErrors * grMass = new TGraphErrors(numMassHyp, &pe_mass[0], &pe_nll[0], &ex[0], &pe_nll_err[0] );
    //A bit of rescaling
    // minLike = grMass->GetMinimum();
    //cout<<minLike<<endl;
    for (int k=0;k<grMass->GetN(); k++) { grMass->GetY()[k] -= minLike;}

    //a bit of editing
    grMass->SetTitle("Negative Log likelihood - Mass Measurement");
    (grMass->GetXaxis())->SetTitle("Mtop [GeV]");
    (grMass->GetYaxis())->SetTitle("-ln(like) [Arb. Units]");

    //Define parabolas...
    sstr.str("");  
    sstr<<"parabola2D_"<<p;
    TF2 *  f2 = new TF2(sstr.str().c_str(), Likelihood2D ,  (*set_mass.begin())-1.0, (*set_mass.rbegin())+ 1.0, (*set_jes.begin())-0.14, (*set_jes.rbegin())+ 0.14, 6);
    f2->SetParNames("mass_min", "mass_err", "jes_min", "jes_err", "corr","min");
    f2->SetParameters(trueMass, 0.3, trueJES,  0.0035, -0.63, 30.0 ); 

    sstr.str("");  
    sstr<<"parabola1D_"<<p;
    TF1 *  f3 = new TF1(sstr.str().c_str(), Likelihood1D ,  (*set_mass.begin())-1.0, (*set_mass.rbegin())+ 1.0, 3);
    f3->SetParNames("mass_min", "mass_err", "min");
    f3->SetParameters(trueMass+2.0,  0.162 , 0.0); 
    f3->SetLineWidth(1);
    f3->SetLineColor(kBlue);

    TFitResultPtr  pr;


    //Fit parabola to likelihoods
    if(numJESHyp > 1){
      sstr.str("");  
      sstr<<"parabola2D_"<<p;   
      pr = grMassJES->Fit(sstr.str().c_str(),"S");
    }else{
      sstr.str("");  
      sstr<<"parabola1D_"<<p;   
      pr = grMass->Fit(sstr.str().c_str(),"RS");
    }

    TFitResult * rr =pr.Get();

    if (rr->Status() != 0)
    {
      if(numJESHyp > 1){
	sstr.str("");  
	sstr<<"parabola2D_"<<p;   
	pr = grMassJES->Fit(sstr.str().c_str(),"S");
      }else{
	cout << f3->GetParameter(1)<<" "<<f3->GetParameter(2)<<endl;
	sstr.str("");  
	sstr<<"parabola1D_"<<p;   
	pr = grMass->Fit(sstr.str().c_str(),"RS");
      }

      rr = pr.Get();
    }

    if(rr->Status() != 0.0){
      skip++;

      if (skip < 3)
      {

	edm_fail->Fill(rr->Edm());
	sstr.str("");
	sstr<<"fail"<<p<<"_"<<(rr->Status());
	canMass->SetName(sstr.str().c_str());
	grMass->SetName(sstr.str().c_str());
	canMass->cd();
	grMass->Draw("a p");
	canMass->Update();	
	canMass->Write();
      }

      delete  grMassJES; delete canMassJES;
      continue;
    }
    
    double absmin;
    if(numJESHyp > 1){
      mem_mass =f2->GetParameter(0);  mem_mass_err =fabs(f2->GetParameter(1));  
      mem_jes =f2->GetParameter(2) ;  mem_jes_err= fabs(f2->GetParameter(3));
      mem_corr = f2->GetParameter(4); absmin= f2->GetParameter(5);
      mem_pmin = f2->GetParameter(5); edm_pass->Fill(rr->Edm());

      mem_mass_intErr =f2->GetParError(0);  mem_mass_err_intErr =fabs(f2->GetParError(1));  
      mem_jes_intErr =f2->GetParError(2) ;  mem_jes_err_intErr= fabs(f2->GetParError(3));
      mem_corr_intErr = f2->GetParError(4); 

    }else{
      
      mem_mass =f3->GetParameter(0);  mem_mass_err =fabs(f3->GetParameter(1));  
      mem_jes = 99;  mem_jes_err = 99.;
      mem_corr = 99. ; absmin = f3->GetParameter(2);
      mem_pmin = f3->GetParameter(2); edm_pass->Fill(rr->Edm());

      mem_mass_intErr =f3->GetParError(0);  mem_mass_err_intErr =fabs(f3->GetParError(1));  
      mem_jes_intErr = 99. ;  mem_jes_err_intErr = 99.;
      mem_corr_intErr = 99.; 
    }

    psexpTree->Fill();
     
    //Safe only the first 15 fits --- space reasons
    if(p-skip < 5){

      if( numJESHyp > 1){
	//a bit of editing
	grMassJES->SetTitle("Negative Log likelihood - Mass Measurement");
	(grMassJES->GetXaxis())->SetTitle("Mtop [GeV]");
	(grMassJES->GetYaxis())->SetTitle("-ln(like) [Arb. Units]");
	(grMassJES->GetXaxis())->SetRangeUser( (*set_mass.begin())-1.0, (*set_mass.rbegin())+ 1.0);
	(grMassJES->GetYaxis())->SetRangeUser( (*set_jes.begin())-0.14, (*set_jes.rbegin())+ 0.14);
	
	TPaveText *myText = new TPaveText(0.6105528,0.236014,0.8894472,0.3863636,"brNDC");
	//NDC sets coords relative to pad
	myText->SetTextSize(0.04);
	myText->SetFillColor(0); //white background
	myText->SetTextAlign(12);
	
	outputFile->cd();
	sstr.str("");
	sstr<<"pseudoexp"<<p-skip<<endl;
	TDirectory *cdtof = outputFile->mkdir(sstr.str().c_str());
	cdtof->cd();  
	
	TF1 f1X("f1X", LikeProjection, (*set_mass.begin())-1.0, (*set_mass.rbegin())+ 1.0, 7);    
	TF1 f1Y("f1Y", LikeProjection, (*set_jes.begin())-0.14, (*set_jes.rbegin())+ 0.14, 7);
	
	Double_t var_array[500];
	Double_t like_array[500];       Double_t likeErr_array[500];
	
	for( set<double>::iterator j_it = set_jes.begin(); j_it!= set_jes.end(); j_it++){
	  selector.reduce2MassArray( sample ,(*j_it), var_array, like_array, likeErr_array);
	  for (int i = 0; i<numMassHyp; i++  ){
	    like_array[i] -= minLike;   
	  }
	  
	  sstr.str("");
	  sstr<<"MassSlice"<<p-skip<<"_JES"<<(*j_it);
	  TCanvas canSlice(sstr.str().c_str(),"", 600, 600); 
	  TGraphErrors grX(numMassHyp, &var_array[0],  &like_array[0],   &ex[0],  &likeErr_array[0] );
	  (grX.GetXaxis())->SetTitle("Mtop [GeV]");    (grX.GetYaxis())->SetTitle("-ln(like) [Arb. Units]");
	  f1X.SetParameters(mem_mass,  mem_mass_err, mem_jes,  mem_jes_err, mem_corr, absmin, (*j_it)); 
	  grX.Draw("a p"); f1X.Draw("same");
	  canSlice.Update();
	  canSlice.Write();
	}
	
	for( set<double>::iterator m_it = set_mass.begin(); m_it!= set_mass.end(); m_it++){

	  selector.reduce2JESArray( sample, (*m_it), var_array, like_array, likeErr_array);
	  
	  for (int i = 0; i<numJESHyp; i++  ){
	    like_array[i] -= minLike;   
	  }
	  
	  sstr.str("");
	  sstr<<"JESSlice"<<p-skip<<"_Mass"<<(*m_it);
	  TCanvas canSlice(sstr.str().c_str(), "", 600, 600); 
	  
	  TGraphErrors grY(numJESHyp, &var_array[0],  &like_array[0],   &ex[0],  &likeErr_array[0] );
	  (grY.GetXaxis())->SetTitle("JES");  (grY.GetYaxis())->SetTitle("-ln(like) [Arb. Units]");
	  
	  f1Y.SetParameters(mem_jes,  mem_jes_err, mem_mass,  mem_mass_err, mem_corr, absmin, (*m_it)); 
	  
	  grY.Draw("a p"); f1Y.Draw("same");
	  
	  canSlice.Update();
	  canSlice.Write();
	}
      

	
	sstr.str("");
	sstr<<"M_{top}: "<<fixed<<setprecision(2)<<mem_mass<<" #pm "
	    <<fixed<<setprecision(2)<<sqrt(pow(mem_mass_err,2) + pow(mem_mass_intErr,2)+pow(mem_mass_err_intErr,2));    
	myText->AddText(sstr.str().c_str());
	
	
	sstr.str("");
	sstr<<"JES: "<<fixed<<setprecision(3)<<mem_jes<<" #pm "
	    <<fixed<<setprecision(3)<<sqrt(pow(mem_jes_err,2) + pow(mem_jes_intErr,2)+pow(mem_jes_err_intErr,2));    
	myText->AddText(sstr.str().c_str());
     
	
	sstr.str(""); 
	sstr<<"contour"<<p;  
	TCanvas * canCont = new TCanvas(sstr.str().c_str(),"MTop", 600, 600);
	Double_t contours[3];
	
	contours[0]= absmin+0.5; contours[1]= absmin+1.0; contours[2]= absmin+1.5;
	
	f2->SetRange(trueMass-1.0, trueJES-0.01, trueMass+1.0, trueJES+0.01);
	f2->SetContour(3, contours);
	canCont->cd();
	f2->Draw("cont1");
	
	TMarker mem_marker;
	mem_marker.SetMarkerStyle(2); mem_marker.SetMarkerColor(kBlack);
	canCont->cd();  mem_marker.DrawMarker(mem_mass,mem_jes);
	
	TMarker true_marker;
	true_marker.SetMarkerStyle(29);true_marker.SetMarkerColor(kBlack);
	canCont->cd();  true_marker.DrawMarker(trueMass, trueJES);
	
	canCont->Update();   
	myText->Draw();
	canCont->Write(); 
	delete canCont;
      

	canMassJES->cd();
	grMassJES->Draw("a p");      
	canMassJES->Write();
	
      }else{
	
	TPaveText *myText = new TPaveText(0.3982412,0.6765734,0.5979899,0.8269231,"brNDC");
	//NDC sets coords relative to pad
	myText->SetName("massResult");
	myText->SetTextSize(0.04);
	myText->SetFillColor(0); //white background
	myText->SetTextAlign(12);
    
	sstr.str("");
	sstr<<fixed<<setprecision(2)<<mem_mass<<" #pm "
	    <<fixed<<setprecision(2)<<sqrt(pow(mem_mass_err,2) + pow(mem_mass_intErr,2)+pow(mem_mass_err_intErr,2)); 

	myText->AddText(sstr.str().c_str());
	
	canMass->cd();
	grMass->Draw("a p");
	canMass->Update();   
	myText->Draw();
	
	canMass->Write();
	
      }
    }//end loop to safe only 20 pseudo-experiments

    delete  grMass;
    delete canMass;

    delete  grMassJES;
    delete canMassJES;
  }//end loop on pseudo-experiments!!!

  TParameter<Int_t>* mem_npes = new TParameter<Int_t>("mem_npes",skip+psexp);
  
  outputFile->cd();
  edm_pass->Write();
  edm_fail->Write();
  mem_npes->Write();
  psexpTree->Write();
  psexpInfoTree->Write();

  //Clean up
  delete outputFile; 
  return;
}

int main(int argc,  char* argv[]){

  if(argc==3){
    runPseudotest((string)argv[1], atof(argv[2]));
    return 0;
  }else{
    cout<<"ERROR! \n \t Usage: \n \t signalBackgroundMassJES [Input File]  [Number of Pseudoexperiments]"<<endl;
    return -1;
  }
}



