#include<stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "TRandom.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TMath.h"
#include "TH2.h"
#include "TH1.h"
#include <vector>
#include "TLorentzVector.h"
using namespace std;

#define MAXEVENT -1

typedef struct {
  Float_t eta, phi, pt, jma, ntracks, btag, hadem, dummy1, dummy2;
  Int_t evtn, type, index;
} PARTICLE;


Double_t simpleTF(Double_t *v, Double_t *par){
  
  //  Double_t x = v[0];
  Double_t y = par[10];
  Double_t d = v[0];


  Double_t a1 = par[0];
  Double_t b1 = par[1];
  Double_t p1 = a1 + b1*y; 

  Double_t a2 = par[2];
  Double_t b2 = par[3];
  Double_t p2 = a2 + b2*y;

  Double_t a3 = par[4];
  Double_t b3 = par[5];
  Double_t p3 = a3 + b3*y;

  Double_t a4 = par[6];
  Double_t b4 = par[7];
  Double_t p4 = a4 + b4*y;

  Double_t a5 = par[8];
  Double_t b5 = par[9];
  Double_t p5 = a5 + b5*y;

  Double_t gaus1 = TMath::Exp(-0.5 * pow((d-p1)/ p2, 2));
  Double_t gaus2 = TMath::Exp(-0.5 * pow((d-p4)/ p5, 2));

  Double_t result = (1./((p2+p3*p5)*TMath::Sqrt(2.*TMath::Pi())))*(gaus1 + p3*gaus2);

  return result ;   
}

void gaussianSmear( string  INFILE_NAME ,  string JES_str){

  float JES = atof(JES_str.c_str());

  string  OUTFILE_NAME =  INFILE_NAME.substr(0, INFILE_NAME.find(".lhco")) + "_smeared_JES"+JES_str+".lhco";
  string  TreeFILE_NAME = INFILE_NAME.substr(0, INFILE_NAME.find(".lhco")) + "_smeared_JES"+JES_str+".root";
 
  // Define simple particle structures
  static PARTICLE parton;
  static PARTICLE smear;
  
  //create a file TreeFILE_NAME to hold event information;
  TFile treeFile(TreeFILE_NAME.c_str(),"recreate");
  //Define event tree 
  TTree partonTree("eventTree","event level info");
  //TTree smearedTree("smearedTree","event level info");
  
  partonTree.Branch("parton",&parton,"eta:phi:pt:jma:ntracks:btag:hadem:dummy1:dummy2/F:evtn:type:index/I");
  partonTree.Branch("smear",&smear,"eta:phi:pt:jma:ntracks:btag:hadem:dummy1:dummy2/F:evtn:type:index/I");

  //JetTree
  Float_t x, y, z;  
  Int_t btag;

  TTree jetTree("JetTree","Tree containing only jet info");
  jetTree.Branch("x",&x, "x/F"); jetTree.Branch("y",&y, "y/F"); 
  jetTree.Branch("z",&z, "z/F"); jetTree.Branch("btag",&btag, "btag/I");


  //Define smearing function
  TF1 *smearF = new TF1("smearFunc", simpleTF, -500, 500, 11);
  smearF->SetParNames("a1", "b1","a2", "b2","a3", "b3","a4", "b4","a5", "b5","parton");
  TLorentzVector jetV, muonV;
  
  //Define histograms
  TH2F EpVsEj("EpartonVsEjet", "Eparton vs Ejet", 1501, -0.5, 1500.5, 1501, -0.5, 1500.5);
  EpVsEj.SetXTitle("Ejet (GeV)");
  EpVsEj.SetYTitle("Eparton (GeV)");

  TH2F deltaE("deltaE", "",  1500, 0.0, 1500.0,  600, -300., 300.);
  deltaE.SetXTitle("Eparton (GeV)");
  deltaE.SetYTitle("Eparton-Ejet (GeV)");
  
  TH2F deltaE1("deltaE1", "",  600, -300., 300.,  1500, 0.0, 1500.0);

  deltaE1.SetYTitle("Eparton (GeV)");
  deltaE1.SetXTitle("Eparton-Ejet (GeV)");
   
    
   //Read input file
   int i = -999;
   if( MAXEVENT != -1) i = -1;

   ifstream myInfile(INFILE_NAME.c_str());
   ofstream myOutfile(OUTFILE_NAME.c_str()); // LHCO file with smeared events
   string s_type, s_index, s_eta, s_phi, s_pt, s_jma, s_ntracks, s_btag, s_hadem, s_dummy1, s_dummy2, s_evt, s_trigger;



   if ( myInfile.is_open() ) {
     cout<<"Opening file "<<INFILE_NAME<<" JES = "<<JES_str<<endl;
  
     while ( myInfile>>s_index && i < MAXEVENT){
  
       if(s_index == "#" ) { //label line
	 myInfile>>s_type;  myInfile>>s_eta;   myInfile>>s_phi;  myInfile>>s_pt;  myInfile>>s_jma;	 
	 myInfile>>s_ntracks;	 myInfile>>s_btag;  myInfile>>s_hadem; myInfile>>s_dummy1; myInfile>>s_dummy2;
	 i++;
	 //write to LHCO file
	 if (i < MAXEVENT )myOutfile<<s_index<<"\t"<<s_type<<"\t"<<s_eta<<"\t"<<s_phi<<"\t"<<s_pt<<"\t"<<s_jma<<"\t"<<s_ntracks<<"\t"<<s_btag<<"\t"<<s_hadem<<"\t"<<s_dummy1<<"\t"<<s_dummy2<<"\n";

       }else if(s_index == "0"){
	 //line with event info is shorter

	 myInfile>>s_evt; //just fill in two (temporary) two variables
	 myInfile>>s_trigger;
	 
	 //	 cout<<i<<" "<<s_evt<<endl;; 
	 parton.evtn = atoi(s_evt.c_str());
	 smear.evtn = atoi(s_evt.c_str());
	 //write to LHCO file
	 myOutfile<<s_index<<"\t"<<s_evt<<"\t"<<s_trigger<<"\n";

       }else{ //read complete line (remaining 10 columns)
	 // cout<<i<<endl;
	 myInfile>>s_type;  myInfile>>s_eta;   myInfile>>s_phi;  myInfile>>s_pt;  myInfile>>s_jma;	 
	 myInfile>>s_ntracks;	 myInfile>>s_btag;  myInfile>>s_hadem; myInfile>>s_dummy1;	 myInfile>>s_dummy2;

	 parton.index= atoi(s_index.c_str());
	 parton.type= atoi(s_type.c_str());
	 parton.eta = atof(s_eta.c_str());
	 parton.phi = atof(s_phi.c_str());
	 parton.pt  = atof(s_pt.c_str());
	 parton.jma= atof(s_jma.c_str());
	 parton.ntracks= atof(s_ntracks.c_str());
	 parton.btag= atof(s_btag.c_str());
	 parton.hadem= atof(s_hadem.c_str());
	 parton.dummy1= atof(s_dummy1.c_str());
	 parton.dummy2= atof(s_dummy2.c_str());
	 
	 smear.index= atoi(s_index.c_str()); 
	 smear.type= atoi(s_type.c_str());
	 smear.eta = atof(s_eta.c_str());
	 smear.phi = atof(s_phi.c_str());
	 smear.pt  = atof(s_pt.c_str());
	 smear.jma= atof(s_jma.c_str());
	 smear.ntracks= atof(s_ntracks.c_str());
	 smear.btag= atof(s_btag.c_str());
	 smear.hadem= atof(s_hadem.c_str());
	 smear.dummy1= atof(s_dummy1.c_str());
	 smear.dummy2= atof(s_dummy2.c_str());
	 	
	 if(s_type == "4"){ //select only jet
	   jetV.SetPtEtaPhiM(parton.pt , parton.eta, parton.phi, parton.jma);
	 
	   Double_t Eparton = jetV.E(); //get parton energy
           Double_t Ejet = 0.0;
           int t = 0 ;	
           
	   //	   cout<<"I found a jet"<<endl;
	   // cout<<"I am smearing it "<<endl;
	   if( parton.btag== 2){
	     smearF->SetRange( -500.0, Eparton );
	     smearF->SetParameters( -1.99091, -0.0243212,2.89767, 0.0886286, 0.0618142,0.00103996,-14.7397,0.00169109,5.1975,0.166501, Eparton);
	   }else{
	     smearF->SetRange( -500.0, Eparton );
	     smearF->SetParameters( -1.99091, -0.0243212, 2.89767, 0.0886286, 0.0618142,0.00103996,-14.7397,0.00169109,5.1975,0.166501, Eparton);
	   }

	   // smearF->Update();
	     Double_t d =smearF->GetRandom();
	     
	     Ejet = Eparton - d;

           while ( pow(Ejet,2) <  pow(parton.jma,2) && t<5 ){

	     //cout<<"I am smearing it "<<endl;
	     
	     if( parton.btag== 2){
	       smearF->SetRange( -500.0, Eparton );
	       smearF->SetParameters( -1.99091, -0.0243212, 2.89767, 0.0886286, 0.0618142,0.00103996,-14.7397,0.00169109,5.1975,0.166501, Eparton);
	     }else{
	       smearF->SetRange( -500.0, Eparton );
	       smearF->SetParameters( -1.99091, -0.0243212, 2.89767, 0.0886286, 0.0618142,0.00103996,-14.7397,0.00169109,5.1975,0.166501, Eparton);
	     }
	     
	     
	     // smearF->Update();
	     d =smearF->GetRandom();
	     
	     Ejet = Eparton - d;
	     
	   }
	   
	   Ejet = Ejet*JES;

	   EpVsEj.Fill(Ejet, Eparton);
	   deltaE.Fill(Eparton, Eparton-Ejet);
	   deltaE1.Fill(Eparton-Ejet, Eparton);
	   
	   Float_t smearedPt = TMath::Sqrt( ( pow(Ejet,2) - pow(smear.jma,2))/ (1+ pow(sinh(smear.eta),2)));



	   smear.pt = smearedPt;
	   
	   x = Ejet;
	   y = Eparton;
	   z = Eparton - Ejet;
	   btag = (parton.btag==2.0)? 1:0 ;

	   jetTree.Fill();
	   // if (Ejet > 1500 || Eparton > 1500)  cout<<Eparton<<" "<<Ejet<<endl;
	   
	 }else if(s_type == "2"){ //select muons and fix mass
	   muonV.SetPtEtaPhiM(parton.pt , parton.eta, parton.phi, parton.jma);
	   Double_t eMuon= muonV.E();
	   smear.jma = 0.1056583668; //set mass to correct value
	   //recalculate Pt
	   smear.pt =  TMath::Sqrt( ( pow(eMuon,2) - pow(smear.jma,2))/ (1+ pow(sinh(smear.eta),2)));
	 }
 
	 //Write to LHCO file
	 myOutfile<<smear.index<<"\t"<<smear.type<<"\t"<<smear.eta<<"\t"<<smear.phi<<"\t"<<smear.pt<<"\t"<<smear.jma<<"\t"
		  <<smear.ntracks<<"\t"<<smear.btag<<"\t"<<smear.hadem<<"\t"<<smear.dummy1<<"\t"<<smear.dummy2<<"\n";
	 
	 //Fill parton and smear trees
	 partonTree.Fill();
	 //	 smearedTree.Fill();
       
       }
       
       if(MAXEVENT == -1 ) i = -99; //consider all events
       
     }
     

     //Write to file     
     EpVsEj.Write();
     deltaE.Write();
     deltaE1.Write();
     
     jetTree.Write();
     partonTree.Write();
     //     smearedTree.Write();
     //Close files
     myInfile.close();  
     myOutfile.close();
//       myOutfile2.close();
      return;

  }    else {
    cout << "Unable to open file"; 
    return;}
}


#ifndef __CINT__

 int main ( int argc, const char* argv[] ) {
 
   if( argc == 3){  
     gaussianSmear ( (string)argv[1], (string)argv[2]) ;
     return 0;
   }else if(argc == 2){
     gaussianSmear ( (string)argv[1], "1");
   }

   return -1;
 }

#endif
