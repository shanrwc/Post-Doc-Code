#include <math.h>
//#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string.h>
#include <sstream> 
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "TROOT.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "RooPolyVar.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TH2.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TMatrixD.h"
#include "TFitResultPtr.h"
#include "TProfile.h"
#include "interface/tdrStyle.h"

void calculateAccptanceError(float Mass, float JES, float err[], float corr[][8], float &error){

  float M = Mass-172.5;  float J = JES-1.00;
  float errorSquare = 0.;

  for (int i=0; i<4; i++){
    errorSquare += pow(err[i]*pow(J,i),2);
    errorSquare += pow(err[i+4]*pow(J,i)*M,2);
  }

  for (int i=0; i<8; i++ ){
    for(int j=0; j<8; j++){
      if(j<=i) continue;
      int JESPower = i+j;
      int MassPower = 0;
      if( i > 3 && j > 3){
	JESPower -= 8;
	MassPower = 2;
      }else if(i>3 || j>3){
	JESPower -=4;
	MassPower = 1;
      }
      errorSquare += 2*err[i]*err[j]*corr[i][j]*pow(J, JESPower)*pow(M, MassPower);
    }
  }
  error = sqrt( errorSquare);
}



const int numPara = 8;

Double_t fun(Double_t *x, Double_t *par) {
 
   Double_t xx =  x[0] - 172.5;
   Double_t yy =  x[1] - 1.00;
 
   Double_t a = par[0]+ par[1] * yy + par[2]*yy*yy + par[3]*yy*yy*yy ;
   Double_t b = par[4]+ par[5] * yy + par[6]*yy*yy + par[7]*yy*yy*yy ;

   Double_t result = a + b * xx;

  return result;
}


Double_t m_fun(Double_t *x, Double_t *par) {
 
  Double_t xx[2] = {x[0], par[numPara]};
  return  fun(xx, par);
}


Double_t j_fun(Double_t *x, Double_t *par) {
  Double_t xx[2] = {par[numPara], x[0]};
  return  fun(xx, par);
}




void acceptanceTerm(){

  //Set TDR style
  TDRStyle tdrSTY; 
  tdrSTY.setTDRStyle();

  cout<<"################################################################################################################"<<endl;
  cout<<"#                                                                                                              #"<<endl;
  cout<<"#  Acceptance Function used :  (p0 + p1*JES + p2*JES^2 + p3*JES^3) + (p4 + p5*JES + p6*JES^2 + p7*JES^3)*Mass  #"<<endl;
  cout<<"#                                                                                                              #"<<endl;
  cout<<"################################################################################################################"<<endl;

  const int numJES = 10;   const int numMass = 8;

  float JESHyp[numJES];
  float jj = 0.80; 
  for (int ii = 0; ii<numJES; ii++) {
    JESHyp[ii] = jj;
    cout<<jj<<" ";
    jj+=0.05;
  } 

  cout<<""<<endl;
  float MassHyp[numMass] ={ 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
  
  //Get data
  TFile f1("/afs/cern.ch/work/c/cimmino/public/Acceptance/selectedRejected.root") ;
 

  float mass, jes;
  int acc;
  TTree *dataTree = (TTree*)f1.Get("accTree");
  dataTree->SetBranchAddress("mass",&mass);
  dataTree->SetBranchAddress("jes",&jes);
  dataTree->SetBranchAddress("acc",&acc);
   
  //Define Output file
  TFile f2("accTest1.root", "recreate") ;
  //Define acceptance tree. Will hold info for likelihood acceptance correction
  float p[numPara];
  float perr[numPara];
  float pcorr[numPara][numPara];
  
  TTree * m_accTree = new TTree("accTermTree","Acceptance Tree");
  stringstream bname;
  stringstream ss;
  
  for (int i = 0 ; i<numPara; i++){
    bname.str("");
    bname<<"p"<<i;
    ss.str("");
    ss<<"p"<<i<<"/F";
    m_accTree->Branch(bname.str().c_str(),&p[i],ss.str().c_str());
  
    bname.str("");
    bname<<"perr"<<i;
    ss.str("");
    ss<<"perr"<<i<<"/F";
    m_accTree->Branch(bname.str().c_str(),&perr[i],ss.str().c_str());
  
    for (int j = 0 ; j<numPara; j++){
      bname.str("");
      bname<<"pcorr"<<i<<"_"<<j;
      ss.str("");
      ss<<"perr"<<i<<"_"<<j<<"/F";
      m_accTree->Branch(bname.str().c_str(),&pcorr[i][j],ss.str().c_str());
    }
  }
  
  TCanvas *c1 = new TCanvas("massCan","",600,600);
  TH1F * m_hn[numJES]; 
  TH1F * m_hd[numJES]; 
  TF1 * m_fitfun[numJES];
  TLegend *leg1 =  new TLegend(0.7952261,0.2202797,0.9798995,0.8723776,NULL,"brNDC");
  leg1->SetHeader("");
  leg1->SetFillColor(kWhite);
  TCanvas *c4 = new TCanvas("jesCan", "", 600, 600);
  TH1F * j_hn[numMass]; 
  TH1F * j_hd[numMass]; 
  TF1 * j_fitfun[numMass];
  TLegend *leg2 = new TLegend(0.1545226,0.1555944,0.3103015,0.5367133,NULL,"brNDC");
  leg2->SetFillColor(kWhite);
  leg2->SetHeader("");
  
  TCanvas *c2 = new TCanvas("AccCan", "", 600, 600);
  TH2F * Accn= new TH2F("Accn", "Acceptance", 300, 160.5, 190.5, 51, 0.745, 1.255);
  TH2F * Accd= new TH2F("Accd", "Acceptance", 300, 160.5, 190.5, 51, 0.745, 1.255);
  TH2F * Acc= new TH2F("Acc", "Acceptance", 300, 160.5, 190.5, 51, 0.745, 1.255);
  Acc->GetXaxis()->SetTitle("Mtop [GeV]");
  Acc->GetYaxis()->SetTitle("JES");  
  Acc->GetZaxis()->SetTitle("Fraction of Accepted Events");
 
  //  TH2F * pull2D= new TH2F("pull2D", "", 300, 160., 190., 40, 0.80, 1.20);

  TCanvas *c_pull = new TCanvas("pull", "", 600, 600);
  TH1F * pull= new TH1F("pull", "", 16, -0.004, 0.004);
  TCanvas *c_pull1 = new TCanvas("pull1", "", 600, 600);
  TH1F * pull_1= new TH1F("pull_1", "", 28, -3.5, 3.5);
  
  TCanvas *c3 = new TCanvas("FunCan", "",600,600);
  TF2* fitfun = new TF2("acceptanceTerm", fun, 160., 185., 0.7, 1.3,  numPara);
  fitfun->SetParameter(0, 4.04266e-01); 
  fitfun->SetParameter(1, -7.67641e-01);
  fitfun->SetParameter(2, 5.e-01);
  fitfun->SetParameter(3, 3.5e-01 );
  
  stringstream hName;
  
  //initialize histograms and 1D functions
  for (int j = 0 ; j<numJES; j++){
    hName.str("");
    hName<<"m_hn_"<<setfill('0')<<setw(3)<<100*JESHyp[j];
    m_hn[j] = new TH1F(hName.str().c_str(), "", 300, 160.5, 190.5);
    m_hn[j]->SetMarkerStyle(7);
    m_hn[j]->SetMarkerColor(1+j);
    if(JESHyp[j]>=1.25) {  m_hn[j]->SetLineColor(kOrange);m_hn[j]->SetMarkerColor(kOrange); }
    else  {m_hn[j]->SetLineColor(1+j);}
    m_hn[j]->Sumw2();
    m_hn[j]->GetXaxis()->SetTitle("Mtop [GeV]");
    m_hn[j]->GetYaxis()->SetTitle("Acc./Rej.");

    hName.str("");
    hName<<"m_hd_"<<setfill('0')<<setw(3)<<100*JESHyp[j];
    m_hd[j] = new TH1F(hName.str().c_str(), "", 300, 160.5, 190.5);
    m_hd[j]->Sumw2();
   
    hName.str("");
    hName<<"m_f_"<<j;
    m_fitfun[j] = new TF1(hName.str().c_str(),m_fun, 160, 185.,  numPara+1);
    if(JESHyp[j]>=1.) {  m_fitfun[j]->SetLineColor(1+j);}
      else   { m_fitfun[j]->SetLineColor(1+j);}
  }
  
  for (int j = 0 ; j<numMass; j++){
    hName.str("");
    hName<<"j_hn_"<<setfill('0')<<setw(4)<<10*MassHyp[j];
    j_hn[j] = new TH1F(hName.str().c_str(), "", 51, 0.745, 1.255);
    j_hn[j]->SetMarkerStyle(7);
    j_hn[j]->SetMarkerColor(1+j);
    j_hn[j]->SetLineColor(1+j);
    j_hn[j]->Sumw2();
    j_hn[j]->GetXaxis()->SetTitle("JES");
    j_hn[j]->GetYaxis()->SetTitle("Acc./Rej.");
    hName.str("");
    hName<<"j_hd_"<<setfill('0')<<setw(4)<<10*MassHyp[j];
    j_hd[j] = new TH1F(hName.str().c_str(), "", 51, 0.745, 1.255);
    j_hd[j]->Sumw2();
   
    hName.str("");
    hName<<"j_f_"<<j;
    j_fitfun[j] = new TF1(hName.str().c_str(),j_fun, 0.745, 1.255,  numPara+1);
    j_fitfun[j]->SetLineColor(1+j);
  }

  
  int inEntries = dataTree->GetEntries();
  //loop on inoput data
  
  for (int i = 0 ; i< inEntries; i++){
  
    dataTree->GetEntry(i);

    hName.str("");
    hName<<"m_hd_"<<setfill('0')<<setw(3)<<100* jes;
    ((TH1F*)gDirectory->Get(hName.str().c_str()))->Fill(mass);

    hName.str("");
    hName<<"j_hd_"<<setfill('0')<<setw(4)<<10* mass;
    ((TH1F*)gDirectory->Get(hName.str().c_str()))->Fill(jes);
    
    Accd->Fill(mass, jes);      
  
    if(acc==1){
      hName.str("");
      hName<<"m_hn_"<<setfill('0')<<setw(3)<<100*jes;
    
      ((TH1F*)gDirectory->Get(hName.str().c_str()))->Fill(mass);

      hName.str("");
      hName<<"j_hn_"<<setfill('0')<<setw(4)<<10*mass;
      ((TH1F*)gDirectory->Get(hName.str().c_str()))->Fill(jes);
      
      Accn->Fill(mass, jes);
    }
  }

  //divide to calculate acceptance
  for (int j = 0 ; j<numJES; j++){
    m_hn[j]->Divide( m_hn[j], m_hd[j], 1, 1, "B");
    hName.str("");
  }


  for (int j = 0 ; j<numMass; j++){
    j_hn[j]->Divide(j_hn[j], j_hd[j], 1, 1, "B");
  }


  Acc->Sumw2();
  Acc->Divide(Accn, Accd, 1, 1, "B");

  c2->cd();
  Acc->Draw("P"); //Draw BEFORE fitting!



    //Perform Fit and save results
  TFitResultPtr  pr =  Acc->Fit("acceptanceTerm", "S E");
  



  TFitResult * r =pr.Get();

  cout<<"Chi2/NDF = "<<(Acc->GetFunction("acceptanceTerm"))->GetChisquare()/(Acc->GetFunction("acceptanceTerm"))->GetNDF()<<endl;
  
  fitfun->Draw("same");
  c2->Write();

  TH1D * prMass = (TH1D * )Accn->ProjectionX();
  TH1D * prJES =  (TH1D * )Accn->ProjectionY();

  prMass->Write();
  prJES->Write();

  //Get Correlation Matrix
  TMatrixDSym mr = r->GetCorrelationMatrix();
  
  //Fill the tree with fit parameters
  for(int i = 0 ; i<numPara; i++) {
    
    //Tree parameters 
    p[i]=fitfun->GetParameter(i);
    perr[i]=fitfun->GetParError(i);
    cout<<"$p_"<<i<<"$ ";
    for(int j = 0 ; j<numPara; j++){
      //Tree parameters     
      pcorr[i][j]= mr(i, j);
       cout<<pcorr[i][j]<<"  ";

    }
        cout<<""<<endl;
  }
  
  //Write Tree to file
  m_accTree->Fill();
  m_accTree->Write();
 
  map< pair<float, float> , float > accMap;
  map< pair<float, float> , float > allMap;

 //  // Do stuff to calculate pull distributions
 //   for (int i = 0; i<numMass; i++){
 //     for(int j = 0 ; j<numJES; j++){

 //       pair <float, float> mypair;
 //       mypair = make_pair(MassHyp[i],JESHyp[j]);

 //       if( allMap.size() == 0 || allMap.count(mypair)== 0){
 // 	accMap[mypair]=1.;  allMap[mypair]=1.; 
 //       }
      
 //       for (int entry = 0 ; entry< inEntries; entry++){
 // 	dataTree->GetEntry(entry);
 // 	if (mass== MassHyp[i] && jes==JESHyp[j] ){
 // 	  allMap[mypair]+=1.;
 // 	  if(acc==1){accMap[mypair]+=1.;}
 // 	}
 //       }
 //     }
 //   }


  for(int entry = 0 ; entry <inEntries; entry++){
    
	dataTree->GetEntry(entry);
	pair <float, float> mypair;
	mypair = make_pair(mass,jes);
	
	if( allMap.size() == 0 || allMap.count(mypair)== 0){
	  accMap[mypair]=0.;  allMap[mypair]=0.; 
	}

	allMap[mypair]+=1.;

	if(acc==1){accMap[mypair]+=1.;}	 
  }






  map< pair<float, float> , float>::iterator itr;

  for (itr = allMap.begin(); itr!= allMap.end(); itr++){
    pair <float, float> mypair = (*itr).first;
    float mymass = mypair.first; float myjes = mypair.second;
    
   
    //  if(mymass == 172.5 ) continue; 

    float mytotal = (*itr).second;
 
    float myacc = accMap[mypair]/mytotal;
 
    float myError = sqrt( myacc*(1 - myacc)/(mytotal-1));
    
    float er;

    calculateAccptanceError(mymass, myjes,  perr,  pcorr,  er);


    float thValue = fitfun->Eval(mymass, myjes);
    //  float thError = fitfunError->Eval(mymass, myjes);
    //  cout<<mymass<<" "<<myjes<<" "<<myacc<<" "<<er<<" "<<thValue<<" "<<myError<<endl;
    pull->Fill( -myacc + thValue  );
    // cout<< myacc - thValue / mymass<<endl;
    pull_1->Fill( (thValue - myacc)/sqrt(pow(myError,2) - pow(er,2)) );
    //   cout<<(thValue - myacc)/er<<endl;
    //pull_1->Fill( (thValue - myacc)/myError);
  }
  
  pull->Fit("gaus");

  pull_1->Fit("gaus");
  // pull_2->Fit("gaus");
  
  // for (int ibin = 1; ibin<= pull2D->GetNbinsX(); ibin++){
  //     for (int jbin = 1; jbin<= pull2D->GetNbinsY(); jbin++){
  //       if ( pull2D->GetBinContent(ibin,jbin) != 0){
  // 	pull->Fill((Accn->GetBinContent(ibin,jbin) - pull2D->GetBinContent(ibin,jbin)));
  // 	pull_1->Fill((Accn->GetBinContent(ibin,jbin) - pull2D->GetBinContent(ibin,jbin))/Accn->GetBinError(ibin,jbin));
  //       }
  //     }
  //   }
  

  c_pull->cd(); pull->Draw();
  c_pull1->cd(); pull_1->Draw();
  
   c_pull1->Write();
  
  //Set parameters for 1D functions
  for (int i = 0 ; i<numJES; i++){
    
    for(int j = 0 ; j<numPara; j++){
      m_fitfun[i]->SetParameter(j, fitfun->GetParameter(j));
    }

    m_fitfun[i]->SetParameter(numPara, JESHyp[i]);
  }

  //Set parameters for 1D functions
  for (int i = 0 ; i<numMass; i++){
    
    for(int j = 0 ; j<numPara; j++){
      j_fitfun[i]->SetParameter(j, fitfun->GetParameter(j));
    }
    
    j_fitfun[i]->SetParameter(numPara, MassHyp[i]);
  }


  //Start Drawing and saving
  c1->cd(); 
  m_hn[0]->GetYaxis()->SetRangeUser(0.0, 0.5);
  m_hn[0]->Draw("e p");

  m_fitfun[0]->Draw("same");

  hName.str("");
  hName<<"JES="<<JESHyp[0];
  leg1->AddEntry(m_fitfun[0],  hName.str().c_str(), "l");
  for (int j = 1 ; j<numJES; j++){

    m_hn[j]->Draw("same e p");
    m_fitfun[j]->Draw("same");
    
    hName.str("");
    hName<<"JES="<<JESHyp[j];
    leg1->AddEntry(m_fitfun[j],  hName.str().c_str(), "l");    
  }

  leg1->Draw("same");
  c1->Write();
  
 


  c4->cd(); 
  j_hn[numMass-1]->Draw("e p");
  j_hn[numMass-1]->GetYaxis()->SetRangeUser(0.0, 0.5);
  j_fitfun[numMass-1]->Draw("same");

  hName.str("");
  hName<<"Mtop="<<MassHyp[numMass-1];
  leg2->AddEntry(j_fitfun[numMass-1],  hName.str().c_str(), "l");  
  for (int j = numMass-2 ; j>=0 ; j--){

    j_hn[j]->Draw("same e p");
    j_fitfun[j]->Draw("same");
    
    hName.str("");
    hName<<"Mtop="<<MassHyp[j];
    leg2->AddEntry(j_fitfun[j],  hName.str().c_str(), "l");
  }
  leg2->Draw("same");
  c4->Write();
  
  c3->cd();
  fitfun->Draw("SURF1");
  c3->Write();




  //A bit of cleaning  
  for (int j = 0 ; j<numJES; j++){
    delete m_fitfun[j];
    delete m_hd[j];
    delete m_hn[j];
  }
  
  for (int j = 0 ; j<numMass; j++){
    delete j_fitfun[j];
    delete j_hd[j];
    delete j_hn[j]; 
    
  }

  delete m_accTree;

  delete c1; delete c2; delete c3; delete c4;

  return;
}


int main ( int argc, const char* argv[] ) {
 

   gROOT->SetStyle("Plain");
  acceptanceTerm();
  
  return 0 ;
}

