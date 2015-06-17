#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "TFile.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "TTree.h"
using namespace std;

bool applyleptonCuts(float type, float eta,  float pt){
  if(type == 2.){//muon
    if(eta <= 2.1 && pt>=26.0 ){
      return true;
    }else{
      return false;
    }
  }else{//electron
    if(eta <= 2.5 && pt>=45.0 ){
      return true;
    }else{
      return false;
    }

  }
}


bool applyjetCuts(float eta, float pt, float dR ,float ntracks){
 if(eta <= 2.4 && pt>=30.0 && ntracks >=1 && dR>= 0.3){
   return true;
  }else{
    return false;
  }
}


float deltaR(float eta1,  float phi1,  float eta2, float phi2) {

   double deta = eta1 - eta2;

   float result = phi1 - phi2;
   while (result > float(M_PI)) result -= float(2*M_PI);
   while (result <= -float(M_PI)) result += float(2*M_PI);
   float dphi = result;
   return std::sqrt(deta*deta + dphi*dphi);
}


void applycuts(string  INFILE_NAME, float m , float j ){


  cout<<"Applying cuts to file: "<<INFILE_NAME<<endl;
  
  float events=0; float muons=0; float electrons;
  float acceptedmuons=0;   float acceptedevents=0;
 
  //Define output files
  string  OUTFILE_NAME =  INFILE_NAME.substr(0, INFILE_NAME.find(".lhco")) + "_selEvt.lhco";
  ofstream myOutfile(OUTFILE_NAME.c_str());
 
  //Open input file
  ifstream myInfile(INFILE_NAME.c_str());
  string  s_index, s_type, s_eta, s_phi, s_pt, s_jma, s_ntracks, s_btag, s_hadem, s_dummy1, s_dummy2;
  string  e_index, e_evt, e_trigger;
  string  p_index[6], p_type[6], p_eta[6], p_phi[6], p_pt[6], p_jma[6], p_ntracks[6], p_btag[6], p_hadem[6], p_dummy1[6], p_dummy2[6];
  string  l_index, l_type, l_eta, l_phi, l_pt, l_jma;
  string  j_eta[6], j_phi[6], j_pt[6], j_jma[6],  j_ntracks[6];
 
 //create a file datasetFILE_NAME to hold event information;
  TFile f("/tmp/cimmino/mytest.root","update");

  //create a tree  for all events 
  TTree *t1;
 
  t1 = (TTree*)f.Get("accTree");
 Float_t sample, mass, jes;
  Int_t acc;
  if(t1==NULL){

    t1 = new TTree("accTree","Acceptance Tree");
    t1->Branch("sample",&sample,"sample/F");
    t1->Branch("mass",&mass,"mass/F");
    t1->Branch("jes",&jes,"jes/F");
    t1->Branch("acc",&acc,"acc/I");

  }else{
  
    t1->SetBranchAddress("sample",&sample);
    t1->SetBranchAddress("mass",&mass);
    t1->SetBranchAddress("jes",&jes);
    t1->SetBranchAddress("acc",&acc);
 
  }
 
  if ( myInfile.is_open() ) {

    while ( myInfile>>s_index){
      
      if(s_index == "#" ){
	int jet = 0;

	events++;
	
	//int flag ;
	//title row
	myInfile>>s_type;
	myInfile>>s_eta;
	myInfile>>s_phi;
	myInfile>>s_pt; 
	myInfile>>s_jma;	
	myInfile>>s_ntracks;
	myInfile>>s_btag;
	myInfile>>s_hadem;
	myInfile>>s_dummy1;
	myInfile>>s_dummy2;

	//event row
	myInfile>>e_index;  myInfile>>e_evt;  myInfile>>e_trigger;

	//loop on partons
	for (int p = 0 ; p<6; p++){

       	//lepton row
	myInfile>>p_index[p];  myInfile>>p_type[p];  myInfile>>p_eta[p];   myInfile>>p_phi[p];  myInfile>>p_pt[p];  myInfile>>p_jma[p];	 
	myInfile>>p_ntracks[p];	 myInfile>>p_btag[p];  myInfile>>p_hadem[p]; myInfile>>p_dummy1[p]; myInfile>>p_dummy2[p];


	//check process type
	if (atof(p_type[p].c_str())==2) {		

	  l_eta = p_eta[p]; l_phi = p_phi[p]; l_pt = p_pt[p]; l_jma=p_jma[p];  
	  muons++;
	}
	else if (atof(p_type[p].c_str())==1) {		

	  l_eta = p_eta[p]; l_phi = p_phi[p]; l_pt = p_pt[p]; l_jma=p_jma[p];
	  electrons++;

	}else if(atof(p_type[p].c_str())==4) {	

	  j_eta[jet] = p_eta[p]; j_phi[jet] = p_phi[p]; j_pt[jet] = p_pt[p]; j_jma[jet]=p_jma[p];  j_ntracks[jet]=p_ntracks[p];
	  jet++;
	}}


	//measure muon isolation
        float dR[4];
	for (int jetnum= 0 ; jetnum<4; jetnum++) {
	  dR[jetnum] =deltaR(atof(l_eta.c_str()), atof(l_phi.c_str()),  atof(j_eta[jetnum].c_str()),atof(j_phi[jetnum].c_str())) ;
	}

	if (jet <=4 ){
	  if (applyleptonCuts(atof(l_type.c_str()), atof(l_eta.c_str()),  atof(l_pt.c_str())) && 
	      applyjetCuts(atof(j_eta[0].c_str()), atof(j_pt[0].c_str()), dR[0], atof(j_ntracks[0].c_str())) && 
	      applyjetCuts(atof(j_eta[1].c_str()), atof(j_pt[1].c_str()), dR[1], atof(j_ntracks[1].c_str())) &&
	      applyjetCuts(atof(j_eta[2].c_str()), atof(j_pt[2].c_str()), dR[2], atof(j_ntracks[2].c_str())) && 
	      applyjetCuts(atof(j_eta[3].c_str()), atof(j_pt[3].c_str()), dR[3], atof(j_ntracks[3].c_str()))){
	    
	    acceptedmuons++;
	    acc = 1;//event accepted
	    
	    //Write to LHCO file
	    myOutfile<<s_index<<"\t"<<s_type<<"\t"<<s_eta<<"\t"<<s_phi<<"\t"<<s_pt<<"\t"<<s_jma<<"\t"
		     <<s_ntracks<<"\t"<<s_btag<<"\t"<<s_hadem<<"\t"<<s_dummy1<<"\t"<<s_dummy2<<"\n";
	    
	    myOutfile<<e_index<<"\t"<<acceptedevents<<"\t"<<e_trigger<<"\n";
	    
	    for (int p = 0 ; p<6; p++) {
	      
	      myOutfile<<p_index[p]<<"\t"<<p_type[p]<<"\t"<<p_eta[p]<<"\t"<<p_phi[p]<<"\t"<<p_pt[p]<<"\t"<<p_jma[p]<<"\t"
		       <<p_ntracks[p]<<"\t"<<p_btag[p]<<"\t"<<p_hadem[p]<<"\t"<<p_dummy1[p]<<"\t"<<p_dummy2[p]<<"\n";
	      
	    }
	  }else{ acc = 0;}
	  
	  
	  mass = m; jes=j;
	  t1->Fill();
	}else{events--;}
      }else{
	cout<<"[ERROR] First line of event does not start with :   # typ eta phi pt jmass  ntrk  btag  had/em   dummy  dummy"<<endl;
	cout<<s_index<<endl;
	cout<<"Exiting!"<<endl;
	return;

      }
    
    }
   

    //    cout<<"JES is "<<j<<endl;

    cout<<"Cuts appplied on "<<events<<" events."<<endl;
    cout<<"Accepted  "<<acceptedmuons<<"/"<<muons<<" muon events."<<endl;
    //   cout<<"Accepted  "<<acceptedelectrons<<"/"<<electrons<<" electron events."<<endl;
  
    myInfile.close();
    myOutfile.close();
    
    t1->Write();
 //     if(muon_old != NULL){
//        muon_old->append(*muon);
//        muon_old->Write();
//      }else{
//        muon->Write();
//      }

    //    delete muon;

  }
}


#ifndef __CINT__

 int main ( int argc, const char* argv[] ) {



   //  s_mass.replace(s_mass.find("_"), 1, ".");

   applycuts ( argv[1], atof(argv[2]),atof(argv[3]) ) ;
   return 0 ;
 }

#endif
