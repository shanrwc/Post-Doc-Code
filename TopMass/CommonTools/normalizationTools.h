#ifndef xsec_normalizationTools_H
#define xsec_normalizationTools_H

#include <string>
#include <sstream>
//ROOFIT
#include <RooDataSet.h>
#include <RooRealVar.h>

using namespace std;
using namespace RooFit;

namespace xsec{

  class normalizationTools{

  public:
    explicit normalizationTools(std::string type){
      std::cout<<"Using cross section for "<<type<<std::endl;
      type_= type;

    };

    ~normalizationTools(){};
    
    
    RooDataSet* loadNormalizationFile( ){

      //Open normalizzation file
      stringstream inFilePath;

      if(type_=="muon" || type_=="MUON"){
	inFilePath.str(""); 
	inFilePath<<"data/crossSection.root";
      }else if (type_=="muon_low"){
	inFilePath.str("");
	inFilePath<<"data/extCrossSection.root";
      }else{
	inFilePath.str("");
      }

      TFile * f_cross = new TFile(inFilePath.str().c_str());
      if(f_cross == NULL) {
	std::cout<<"[ERROR]: Acceptance file "<<inFilePath.str()<<" not found."<<std::endl;
	return NULL;
      }
      //Access normalization tree
      TTree *tree_cross =(TTree*)f_cross->Get("crossSection_m");  
      if (tree_cross==NULL) {
	std::cout<<"[ERROR]: Acceptance tree not valid."<<std::endl;
	return NULL;
      }
      
      RooRealVar cross_mass("mass","mass", 150., 200.) ; 
      RooRealVar cross("cross","cross sectiom", 0., 100.) ; 
      RooRealVar cross_err("cross_err","cross sectiom", 0., 100.) ; 

      RooDataSet * myData = new RooDataSet("crossTree","Acceptance Data",RooArgSet(cross_mass, cross, cross_err), Import(*tree_cross));

      return myData;
    }

    
  



   

  private:
    std::string type_;
  
    
  };





}


#endif
