#ifndef acc_acceptanceTools_H
#define acc_acceptanceTools_H

#include <string>
#include <sstream>
//ROOFIT
#include <RooDataSet.h>
#include <RooRealVar.h>

using namespace std;
using namespace RooFit;

namespace acc{

  class acceptanceTools{

  public:
    explicit acceptanceTools(std::string type){
      std::cout<<"Using acceptance for "<<type<<std::endl;
      type_= type;

      if (type_=="validation" || type_=="Validation"){
	numberOfPar = 8;
      } else if (type_== "muon_calibration" ){
	numberOfPar = 2;
      }else if (type_== "muon_hitfitcut" ){
	numberOfPar = 2;
      }else if (type_ == "muon_hitfit_tf" ){
	numberOfPar = 2;
      }else if (type_== "electron_calibration" ){
	numberOfPar = 2;
      }else{
	numberOfPar = 0;
      }


      for (int i = 0 ; i< numberOfPar; i++){
	parameters_[i] = 0.; parErrors_[i] = 0.;
	for (int j=0; j<numberOfPar; j++) { 
	  corrMatrix_[i][j]=0.; 
	}
      }



    };

    ~acceptanceTools(){};
    
    int numberOfPar;
    
    
    int loadAcceptanceFile(){
      
      //Define variables for acceptance term
      RooRealVar *accVar[8];  RooRealVar *accErrVar[8]; RooRealVar * accCorrVar[64];
      RooArgSet accArgSet;
      stringstream ss;   int corrIndex=0;
      for (int i=0; i<numberOfPar; i++){
	ss.str("");
	ss<<"p"<<i;
	accVar[i]= new RooRealVar(ss.str().c_str() ,ss.str().c_str(), -1., 1.);
	accArgSet.add(*accVar[i]);
	ss.str("");
	ss<<"perr"<<i;
	accErrVar[i]= new RooRealVar(ss.str().c_str() ,ss.str().c_str(), -1., 1.);
	accArgSet.add(*accErrVar[i]);
	for (int j=0; j<numberOfPar; j++){
	  ss.str("");
	  ss<<"pcorr"<<i<<"_"<<j;
	  accCorrVar[corrIndex]= new RooRealVar(ss.str().c_str() ,ss.str().c_str(), -1., 1.);
	  accArgSet.add(*accCorrVar[corrIndex]);					    
	}
      }
      

      std::stringstream inFilePath;
      inFilePath.str("");  
      
      if (type_=="validation" || type_=="Validation"){
	inFilePath<<"data/acceptanceTerm.root"; 
      }else if(type_=="muon_calibration"){
	inFilePath<<"data/muon_calibration.root"; 
      }else if (type_=="muon_hitfitcut"){
	inFilePath<<"data/acceptanceTerm_wcut.root";
      }else if (type_ == "muon_hitfit_tf") {
	inFilePath<<"data/accHitFitcut_TFmatch.root";
      }else if(type_=="electron_calibration"){
	inFilePath<<"data/electron_calibration.root";
      }

      
      TFile * file_acc = new TFile(inFilePath.str().c_str()); 
      if(file_acc == NULL) {
	std::cout<<"[ERROR]: Acceptance file not "<<inFilePath.str()<<" found"<<std::endl;
	return 0;   
      }
      
      //Get Acceptance term tree
      TTree * accTree = (TTree*)file_acc->Get("accTermTree");
      if (accTree==NULL) {
	std::cout<<"[ERROR]: Acceptance tree not valid."<<std::endl;
	return 0;
      }
      std::cout<<"Getting acceptance file: "<<inFilePath.str()<<std::endl;

      RooDataSet * accData = new  RooDataSet("accData","Acceptance Data", accArgSet, Import(*accTree)); 
      if( accData == 0 ){ 
	cout<<"[ERROR]: Acceptance file not found!";
	return 0; 
      }else{
	const  RooArgSet * accRow = accData->get(0); 
	corrIndex=0; //resetIndex
	if(accRow != 0  ){
	  for (int i = 0 ; i<numberOfPar ; i++){//Read acceptance tree and fill terms for acceptance parameters
	    ss.str("");
	    ss<<"p"<<i;
	    parameters_[i]=((RooRealVar*)accRow->find(ss.str().c_str()))->getVal();
	    ss.str("");
	    ss<<"perr"<<i;
	    parErrors_[i]=((RooRealVar*)accRow->find(ss.str().c_str()))->getVal();
	    for (int j=0; j<numberOfPar; j++){
	      ss.str("");
	      ss<<"pcorr"<<i<<"_"<<j;
	      corrMatrix_[i][j] = ((RooRealVar*)accRow->find(ss.str().c_str()))->getVal();
	      corrIndex++;
	    }
	  }//End read acceptance tree and terms for acceptance parameters
	  return 1;
	}else{
	  cout<<"[ERROR]: Acceptance file found but EMPTY!";
	  return 0;
	}
      }
    }



    void calculateAcceptance(double myMass, double myJES, double &myAcceptance){
	myAcceptance = 0.;

      if (type_=="validation" || type_=="Validation"){
	double accMass = myMass-172.5;  double accJES =  myJES-1.00;
	double p0 = parameters_[0]+(parameters_[1]*accJES)+(parameters_[2]*pow(accJES,2))+(parameters_[3]*pow(accJES,3));
	double p1 = parameters_[4]+(parameters_[5]*accJES)+(parameters_[6]*pow(accJES,2))+(parameters_[7]*pow(accJES,3));
	
	myAcceptance = p0 + p1*accMass;
      
      }else if(type_=="muon_calibration"){
	double accMass = myMass-172.5;
	double p0 = parameters_[0];
	double p1 = parameters_[1];
	myAcceptance = p0 + p1*accMass;
      }else if(type_=="electron_calibration"){
	double accMass = myMass-172.5;
	double p0 =  1.149308e-02; //parameters_[0];
	double p1 = 1.12941e-04;//parameters_[1];
	myAcceptance = p0 + p1*accMass;
      }else if (type_=="muon_hitfitcut"){
	double accMass = myMass-172.5;
	double p0 = parameters_[0];
	double p1 = parameters_[1];
	myAcceptance = p0 + p1*accMass;
      }else if (type_ == "muon_hitfit_tf"){
	double accMass = myMass - 172.5;
	double p0 = parameters_[0];
	double p1 = parameters_[1];
	myAcceptance = p0 + p1*accMass;
      }
      else{
	
	myAcceptance = 0.;
     
      }
      
    }




    
    
    void calculateAcceptanceError(double Mass, double JES, double &error){
      
      double M = Mass - 172.5;
      double J = JES - 1.00;

      if (type_=="validation" || type_=="Validation"){
	
	double errorSquare = 0.;
	
	for (int i=0; i<4; i++){
	  errorSquare += pow(parErrors_[i]*pow(J,i),2);
	  errorSquare += pow(parErrors_[i+4]*pow(J,i)*M,2);
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
	    errorSquare += 2*parErrors_[i]*parErrors_[j]*corrMatrix_[i][j]*pow(J, JESPower)*pow(M, MassPower);
	  }
	}
	error = sqrt( errorSquare);
      }else if (type_ == "muon_calibration"){

	double p0_err = parErrors_[0];
	double p1_err = parErrors_[1];

	double corrTerm = corrMatrix_[0][1];
        
	error = sqrt(pow(p0_err,2)+pow(p1_err*M,2) + 2*M*corrTerm*p0_err*p1_err);
      }else if (type_ == "electron_calibration"){

	double p0_err = 2.01552e-05;//parErrors_[0];
	double p1_err = 6.28081e-06;//parErrors_[1];

	double corrTerm = 0.211601;//corrMatrix_[0][1];
        
	error = sqrt(pow(p0_err,2)+pow(p1_err*M,2) + 2*M*corrTerm*p0_err*p1_err);
      }else if (type_ == "muon_hitfitcut"){
	double p0_err = parErrors_[0];
	double p1_err = parErrors_[1];
	double corrTerm = corrMatrix_[0][1];
	error = sqrt(pow(p0_err,2)+pow(p1_err*M,2) + 2*M*corrTerm*p0_err*p1_err);
      }else if (type_ == "muon_hitfit_tf"){
	double p0_err = parErrors_[0];
	double p1_err = parErrors_[1];
	double corrTerm = corrMatrix_[0][1];
	error = sqrt(pow(p0_err,2)+pow(p1_err*M,2) + 2*M*corrTerm*p0_err*p1_err);
    }else{
	error = -1.;
      }
    }

  
  

  private:
    std::string type_;
    double parameters_[50] ;
    double parErrors_[50] ;
    double corrMatrix_[50][50] ;

    
  };





}


#endif
