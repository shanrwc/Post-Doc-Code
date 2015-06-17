#ifndef SAMPLESELECTOR_H
#define SAMPLESELECTOR_H

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TRandom2.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

#include "Mass/CommonTools/normalizationTools.h"
#include "Mass/CommonTools/acceptanceTools.h"

using namespace std;
using namespace RooFit;

//This class handles reading in a sample or samples for use in pseudoexperiments.
//Eventually, it will be able to mix different samples and apply cuts,
//but that will require more implementation in the future.
//The main point is to give that analysis code a RooDataSet or RooDataSets to play with.

class sampleSelector{

 public:

   sampleSelector()
   { 
     normalizationType_ = "";
     acceptanceType_ = "";
     useResampling = true;
   }

   ~sampleSelector()
   {
     ;
   }
  
   int initialize(string cardfile, vector<pair<int,int> > &numEvtPools)
   {
     //Read in parameters from card file 
     ifstream drivecard (cardfile.c_str());
     
     string line; 

     double tSigFrac = 0.0;
     double totalFrac = 0.0;
 
     double mlow = 0;
     double mhigh = 0;
     double jlow = 0;
     double jhigh = 0;
     
     vector<string> filenames; 
     
     if (drivecard.is_open()) 
     { 
       while (drivecard.good()) 
       { 
	 
	 getline(drivecard,line); 

	 string first = line.substr(0,1); 
	 if (first.compare("#") != 0) 
	 { 
	   vector<string> pieces = splitString(line);
	   if (pieces.size()  == 0) continue;
	   if (pieces[0].find("N_EVENTS:") != string::npos && pieces.size() > 1) 
	   {
	     nEvents = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("LOW_MASS_HYPOTHESIS") != string::npos && pieces.size() > 1)
	   {
	     mlow = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("HIGH_MASS_HYPOTHESIS") != string::npos && pieces.size() > 1)
	   {
	     mhigh = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("LOW_JES_HYPOTHESIS") != string::npos && pieces.size() > 1)
	   {
	     jlow = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("HIGH_JES_HYPOTHESIS") != string::npos && pieces.size() > 1)
	   {
	     jhigh = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("NORMALIZATION_TYPE:") != string::npos && pieces.size() > 1)
	   {
	     normalizationType_ = pieces[1];
	   }
	   if (pieces[0].find("USE_RESAMPLING:") != string::npos && pieces.size() > 1)
	   {
	     if (pieces[1].find("no") != string::npos) useResampling = false;
	   }
	   if (pieces[0].find("APPLY_ACCEPTANCE:") != string::npos && pieces.size() > 1)
	   {
	     if (pieces[1].find("yes") != string::npos) applyAcceptance = true;
	     if (pieces[1].find("no") != string::npos) applyAcceptance = false;
	   }
	   if (pieces[0].find("ACCEPTANCE_TYPE:") != string::npos && pieces.size() > 1)
	   {
	     cout << pieces[1] << endl;
	     acceptanceType_ = pieces[1];
	   }
	   if (pieces[0].find("SIGNAL_FILE:") != string::npos && pieces.size() > 2) 
	   { 
	     TFile *sfile = new TFile(pieces[1].c_str());
	     if (sfile == NULL)
	     { 
	       cout << "[ERROR}: Input file " << pieces[1] <<" not found!"<<endl;
	       return (2);
	     }

	     RooRealVar evt("evt","evt",0,500000000);
	     RooRealVar mass("mass","mass",150.,200.);
	     RooRealVar jes("jes","jes",0.5,1.5);
	     RooRealVar like("like","like",1.0e-40,100);
	     RooRealVar like_err("like_err","like_err",-1000,1000);
 
	     TTree* stree = (TTree*)sfile->Get("likeTree");
	     if (stree == NULL)
	     {
	       cout << "[ERROR]: Likelihood tree not found!" << endl;
	       return (3);
	     }
 
	     sigData = new RooDataSet("sigData","ProcessLikelihoods",RooArgSet(evt,mass,jes,like,like_err),Import(*stree));
	     //Also, get tree with truth information
	     infoTree_ = (TTree*)sfile->Get("infoTree");
	     if (infoTree_ == NULL)
	     { 
	       cout << "[ERROR]: Truth information not found!" << endl;
	       return(3);
	     }
	     tSigFrac = makeDouble(pieces[2]);
	   }
	   if (pieces[0].find("BACKGROUND_FILE:") != string::npos && pieces.size() > 2)
	   {  
	     TFile* bfile = new TFile(pieces[1].c_str());
	     if (bfile == NULL)
	     {
	       cout << "[ERROR]: Input file " << pieces[1] << " not found!" << endl;
	       return (1);
	     }

	     RooRealVar evt("evt","evt",0,500000000);
	     RooRealVar mass("mass","mass",150.,200.);
	     RooRealVar jes("jes","jes",0.5,1.5);
	     RooRealVar like("like","like",1.0e-40,100);
	     RooRealVar like_err("like_err","like_err",-1000,1000);

	     TTree* btree = (TTree*)bfile->Get("likeTree");
	     if (btree == NULL)
	     {
	       cout << "[ERROR]: Likelihood tree not found!" << endl;
	       return (2);
	     }

	     RooDataSet* temp = new RooDataSet("bkgData","ProcessLikelihoods",RooArgSet(evt,mass,jes,like,like_err),Import(*btree));

	     bkgDatas.insert(pair<RooDataSet*,double>(temp,makeDouble(pieces[2])));
	     totalFrac += makeDouble(pieces[2]);
	   }
	 }//closes if-statement over drivecard's active lines
       }//closes while-loop over drivecard's goodness
     }//closes if-statement over drivecard's openness 
     else
     {
       cout << "[ERROR]: Config text file not found!" << endl;
       return(1);
     }

     //Now, check that hypothesis sets are consistent to signal
     //Loop on Signal sample
     for (int k = 0; k != sigData->numEntries(); ++k)
     {
       const RooArgSet* row = sigData->get(k);
       double mhypo = ((RooRealVar*)row->find("mass"))->getVal();
       double jhypo = ((RooRealVar*)row->find("jes"))->getVal();
       if(mhypo >= mlow && mhypo<= mhigh && jhypo >= jlow && jhypo<= jhigh ){
	 set_mass_.insert(mhypo);
	 set_jes_.insert(jhypo);
       }
     }

     if( set_mass_.size()==0 || set_jes_.size()==0 ){
       cout << "[ERROR]: Mass and/or JES hypothesis not compatible with request" << endl;
       return(1); 
     }

     //Loop on all background samples
     map<RooDataSet*, double>::iterator iter;
     set<double> bMHypos;
     set<double> bJHypos;
     for (iter = bkgDatas.begin(); iter != bkgDatas.end(); ++iter){    
       bMHypos.clear(); bJHypos.clear();
       for (int k = 0; k != (iter->first)->numEntries(); ++k)
       {  
	 const RooArgSet* row = (iter->first)->get(k);
	 double mhypo = ((RooRealVar*)row->find("mass"))->getVal();
	 double jhypo = ((RooRealVar*)row->find("jes"))->getVal();
	 if(  set_mass_.find(mhypo)!=set_mass_.end() &&   set_jes_.find(jhypo)!=set_jes_.end()){
	   bMHypos.insert(mhypo);
	   bJHypos.insert(jhypo);
	 }
       } 
       
       //Now final comparison      
       if(set_mass_ != bMHypos || set_jes_ != bJHypos ){ 
	 cout << "[ERROR]: Mismatch between signal and background mass hypotheses!" << endl; 
	 return(5); 
       } 
     } 
     
     //Set up normalization and acceptance tools
     xsec::normalizationTools normTool(normalizationType_);
     crossData = normTool.loadNormalizationFile();
     accTool = new acc::acceptanceTools(acceptanceType_);
     accTool->loadAcceptanceFile();

     //Finally, calculate signal fraction of ensembles
     if (totalFrac > 1.0)
     {
       cout << "[ERROR]: Fractions assigned to backgrounds exceed 1.0!" << endl;
       return(7);
     }
     if (tSigFrac != (1.0-totalFrac))
     {
       cout << "[ERROR]: Given signal fraction does not equal 1-sum(background fractions)!" << endl;
       return(8);
     }
     nSigEvents = nEvents*(1.0-totalFrac);

     //Now fill sets with event numbers and weights, and write out pool information for the user
     //Start by writing out signal sample info
     int likeEntries = sigData->numEntries();
     cout<<"############################"<<endl;
     cout<<"Entries in signal sample likelihood file: "<<likeEntries<<endl;
     set<int> set_evt;
     
     for (int i = 0; i< likeEntries; i++)
     {
       const  RooArgSet * row = sigData->get(i);
       set_evt.insert(((RooRealVar*)row->find("evt"))->getVal());
     }
     
     numEvtPools.push_back(pair<int,int>(nSigEvents,set_evt.size()));
     totalHyp_ =  set_mass_.size()*set_jes_.size();
     copy(set_evt.begin(),set_evt.end(),back_inserter(sig_evts));
     
     cout<<"# Found events: "<<sig_evts.size()<<endl;
     cout<<"# Found mass hyp: "<<set_mass_.size()<<endl;
     cout<<"# First mass hyp: "<<(double)(*set_mass_.begin())<<endl;
     cout<<"# Found JES hyp: "<<set_jes_.size()<<endl;
     cout<<"# First JES hyp: "<<(double)(*set_jes_.begin())<<endl;
     cout<<"############################"<<endl;
     
     //Now print out a limited set for the background files
     for (iter = bkgDatas.begin(); iter != bkgDatas.end(); ++iter)
     {
       RooDataSet* temp = iter->first;
       int nEntries = temp->numEntries();
       cout << "Entries in background sample likelihood file: " << nEntries << endl;
       set_evt.clear();
       for (int i = 0; i != nEntries; ++i)
       {
	 const RooArgSet* row = temp->get(i);
	 set_evt.insert(((RooRealVar*)row->find("evt"))->getVal());
       }
       cout << "# Found events: " << set_evt.size() << endl;
       numEvtPools.push_back(pair<int,int>(iter->second*nEvents,set_evt.size()));
       cout << "# Fraction of Sample: " << iter->second << endl;
       
       vector<int> temp2;
       copy(set_evt.begin(),set_evt.end(),back_inserter(temp2));
       bkg_evts.push_back(temp2);
     }
     cout << "############################"<<endl;
     
     cout << "sampleSelector successfully initialized." << endl;
     return(0);
   } 


   int getTruthInfo(double& trueMass, double& trueJES)
   {

     if(infoTree_){

     RooRealVar tJES("trueJES", "trueJES", 0.5, 1.5);
     RooRealVar tMass("trueMass", "trueMass", 160., 200.);

     RooDataSet * infoData = new  RooDataSet("infoData","True Mass and JES",RooArgSet(tMass, tJES), Import(*infoTree_)); 
     const  RooArgSet * infoRow = infoData->get(0);
     trueJES = ((RooRealVar*)infoRow->find("trueJES"))->getVal();
     trueMass = ((RooRealVar*)infoRow->find("trueMass"))->getVal();

     return 1;    
     }else{  return -1;}

   }
  

  void reduce2MassArray(RooDataSet * myData, double JESHyp, Double_t array[], Double_t like_array[], Double_t likeErr_array[]){

    int m =0;
    double min = 100000000;
    std::stringstream dataCut;

    for ( set<double>::iterator m_it = set_mass_.begin(); m_it!= set_mass_.end(); m_it++, m++){
      double M= (*m_it);
      dataCut.str("");
      dataCut<<"jes<"<<JESHyp+0.01<<" && "<<"jes>"<<JESHyp-0.01<<" && mass<"<<M+0.01<<" && mass>"<<M-0.01;
      RooDataSet* cutdata = (RooDataSet*)myData->reduce(dataCut.str().c_str());
      array[m]=M;
      const  RooArgSet * row = cutdata->get(0);
      like_array[m]=((RooRealVar*)row->find("norm_nl_like"))->getVal();
      if (like_array[m] < min) min = like_array[m];
      likeErr_array[m]=((RooRealVar*)row->find("norm_nl_err"))->getVal();
    }

    for (int i = 0; i != m; ++i)
    {
      like_array[i] = like_array[i] - min;
    }
  }

  void reduce2JESArray( RooDataSet * myData, double MassHyp, Double_t array[], Double_t like_array[], Double_t likeErr_array[]){

   int m =0;
   double min = 100000000;
   std::stringstream dataCut;
   for ( set<double>::iterator j_it = set_jes_.begin(); j_it!= set_jes_.end(); j_it++, m++){
     double J= (*j_it);
     dataCut.str("");
     dataCut<<"jes<"<<J+0.01<<" && "<<"jes>"<<J-0.01<<" && mass<"<<MassHyp+0.01<<" && mass>"<<MassHyp-0.01;
     RooDataSet* cutdata = (RooDataSet*)myData->reduce(dataCut.str().c_str());
     array[m]=J;

     const  RooArgSet * row = cutdata->get(0);
     like_array[m]=((RooRealVar*)row->find("norm_nl_like"))->getVal();
     if (like_array[m] < min) min = like_array[m];
     likeErr_array[m]=((RooRealVar*)row->find("norm_nl_err"))->getVal();
   }

   for (int i = 0; i != m; ++i)
   {
     like_array[i] = like_array[i] - min;
   }
   
 }
 
  RooDataSet* getEnsemble()
  {
    //This is the function that should be called externally to get one ensemble or pseudoexperiment
    //It needs to get selected events from each component sample via the selectEvents function
    int seed = 10;

    RooDataSet* ensemble = selectEvents(sigData,sig_evts,nSigEvents,seed);
    map<RooDataSet*,double>::iterator iter;
    int n = 0;
    for (iter = bkgDatas.begin(); iter != bkgDatas.end(); ++iter)
    {
      seed += 10;
      ensemble->append(*selectEvents(iter->first,bkg_evts[n],nEvents*(iter->second),seed));
    }

    //It needs to combine all the samples and then combine event likelihoods into process likelihoods
    //And normalize/add acceptance to everything
    RooDataSet* normEnsemble = mergeEventLLs(ensemble);

    //And then it needs to return the dataset
    return(normEnsemble);
  }
 
  RooDataSet* selectEvents( RooDataSet* dataset, vector<int> evtNums, int nEvts, int seed )
  {
    //Select the requested number of events with replacement and return the datset
    RooRealVar evt("evt","evt", 0, 500000000);
    RooRealVar mass("mass","mass", 150., 200.);
    RooRealVar jes("jes","jes", 0.50, 1.50);
    RooRealVar like("like", "like", 1.0e-40, 100);
    RooRealVar nl_like("nl_like", "nl_like", 0, 1000000);
    RooRealVar nl_err("nl_err", "nl_err", 0, 1000000);
    RooRealVar like_err("like_err", "like_err", 0, 100);
    
    std::stringstream cut2;  

    RooDataSet*  selData  = new RooDataSet("selData","Selected Data", RooArgSet(evt, mass, jes, like, like_err));
   
    ++count;
    TRandom2 myrnd(count + seed);

    while ((selData->numEntries()/ totalHyp_  < nEvts) && (evtNums.size() > 0))
    {
      int rEvts = min(190, nEvts-selData->numEntries()/ totalHyp_ );
      string cut = "(";
      
      for (int i = 0; i < rEvts; ++i)
      {
	if (evtNums.size() == 0) continue;
	else if (evtNums.size() == 1)
	{
	  cut += "evt==" + makeString(evtNums[0],0) + " || ";
	  if (!useResampling) evtNums.clear();
	}
	else
	{
	  int tempN = (int)myrnd.Uniform(evtNums.size());
	  cut += "evt==" + makeString(evtNums[tempN],0) + " || ";
	  
	  if (!useResampling)
	  {
	    evtNums.erase(evtNums.begin() + tempN);
	  }
	}
      }
      cut += "evt==-99) && mass>"+makeString(*set_mass_.begin()-0.1,3)+" && mass<"+makeString(*set_mass_.rbegin()+0.1,3)+" && jes>"+makeString(*set_jes_.begin()-0.1,3)+" && jes<"+makeString(*set_jes_.rbegin()+0.1,3);

      RooDataSet* cutdata = (RooDataSet*)dataset->reduce(cut.c_str());
      selData->append(*cutdata);
     }
  
    return(selData);
  }

  RooDataSet* mergeEventLLs(RooDataSet* dataset)
  {
    std::stringstream cut2;

    RooRealVar mass("mass","mass", 150., 200.);
    RooRealVar jes("jes","jes", 0.50, 1.50);
    RooRealVar norm_nl_like("norm_nl_like", "norm_nl_like", 0, 1000000);
    RooRealVar norm_nl_err("norm_nl_err", "norm_nl_err", 0, 100000);
    RooDataSet *  sampleData = new RooDataSet("sampleData","Sample likelihood", RooArgSet(mass, jes, norm_nl_like, norm_nl_err ) );
    
    RooRealVar cross ("cross","cross section",0.,100.);
    RooRealVar cross_err("cross_err","cross section",0.,100.);

    double accept = 1.0;
    double accept_err = 0.0;

    for( set<double>::iterator m_it = set_mass_.begin(); m_it!= set_mass_.end(); m_it++)
    {
      for( set<double>::iterator j_it = set_jes_.begin(); j_it!= set_jes_.end(); j_it++)
      {
	cut2.str("");
	cut2<<"mass>"<<*m_it-0.01<<" && mass<"<<*m_it+0.01<<" && jes<"<<*j_it+0.01<<" && jes>"<<*j_it-0.01;

	double nl_like;
	double nl_err;

	RooDataSet * massJesData = (RooDataSet*) dataset->reduce(cut2.str().c_str()) ;

     	double nll=0.0; double nll_err=0.0;
	mass =*m_it; jes =*j_it;
	
	//Calculate negative log likelihood
	for (int i = 0 ; i< massJesData->numEntries(); i++)
	{
	  const  RooArgSet * row = massJesData->get(i);
	  nll+=-log(((RooRealVar*)row->find("like"))->getVal());
	  //	  cout << nll << " " << ((RooRealVar*)row->find("like"))->getVal() << endl;
	  nll_err+=((RooRealVar*)row->find("like_err"))->getVal()/((RooRealVar*)row->find("like"))->getVal();
	}
	
	nl_like = nll; nl_err = nll_err;

	//Add the normalization and acceptance, if using
	cut2.str("");
	cut2<<"mass>"<<*m_it-0.01<<" && mass<"<<*m_it+0.01;
	RooDataSet* massCrossData = (RooDataSet*) crossData->reduce(RooArgSet(cross,cross_err),cut2.str().c_str());
	if (massCrossData->numEntries() > 1)
	{
	  cout << "[ERROR]: More than one normalization found for " << cut2.str() << endl;
	  return 0;
	}

	const RooArgSet* rowCross = massCrossData->get(0);
	double norm = ((RooRealVar*)rowCross->find("cross"))->getVal();
	double norm_err = ((RooRealVar*)rowCross->find("cross_err"))->getVal();

	if (applyAcceptance)
	{
	  accTool->calculateAcceptance(*m_it,*j_it,accept);
	  accTool->calculateAcceptanceError(*m_it,*j_it,accept_err);
	}

	norm_nl_like = nl_like + nEvents*log(norm) + nEvents*log(accept);
	norm_nl_err = nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept);

	sampleData -> add(RooArgSet(mass, jes, norm_nl_like, norm_nl_err) );
      }
    }
    return(sampleData);
  }

void reduce2Arrays( RooDataSet * myData, Double_t mass_array[],  Double_t jes_array[], Double_t like_array[], Double_t likeErr_array[]){

  int k =0;
  std::stringstream dataCut;
  for( set<double>::iterator m_it = set_mass_.begin(); m_it!= set_mass_.end(); m_it++){
      
    for( set<double>::iterator j_it = set_jes_.begin(); j_it!= set_jes_.end(); j_it++){
      dataCut.str("");
      dataCut<<"mass<"<<*m_it+0.01<<" && "<<"mass>"<<*m_it-0.01<<" && "<< "jes<"<<*j_it+0.01<<" && "<<"jes>"<<*j_it-0.01;     

      RooDataSet * massJesData = (RooDataSet*) myData->reduce(dataCut.str().c_str()) ;

      const  RooArgSet * row = massJesData->get(0);
      if (row==0)
      {
	cout << "reduce2Arrays function produced an empty dataset!" << endl;
	return;
      }
      like_array[k]=((RooRealVar*)row->find("norm_nl_like"))->getVal();
      likeErr_array[k]=((RooRealVar*)row->find("norm_nl_err"))->getVal();
      mass_array[k]=*m_it; jes_array[k]=*j_it;
            
      k++;
    }    
  }
}

 int getPseudoSize()
 {
   return(nEvents);
 }

 set<double> getMassHypos()
 {
   return(set_mass_);
 }

 set<double> getJESHypos()
 {
   return(set_jes_);
 }

 private:
  int count;
  int totalHyp_;
  TTree* infoTree_;

  //Parameters and objects for normalizationa and acceptance
  string normalizationType_;
  RooDataSet* crossData;
  bool applyAcceptance;
  string acceptanceType_;
  acc::acceptanceTools* accTool;

  //Sets to store mass and jes hypotheses
  set<double> set_mass_;
  set<double> set_jes_;    
  set<int> set_evt_;

  //Pseudoexperiment parameters and the data themselves
  bool useResampling;
  int nEvents;
  int nSigEvents;
  RooDataSet* sigData;
  map<RooDataSet*,double> bkgDatas;
    
  //vectors to store lists of event numbers (useful to have when making pseudoexperiments, I promise)
  vector<int> sig_evts;
  vector< vector<int> > bkg_evts;

  vector<string> splitString(string phrase)
  {
    stringstream ss(phrase);
    string buffer;
    vector<string> pieces;
    
    while(ss >> buffer)
    {
      pieces.push_back(buffer);
    }
    return(pieces);
  }

  template<typename T>
  std::string makeString(T const& value, int prec)
  {
    std::stringstream sstr;
    sstr << fixed;
    sstr << setprecision(prec);
    sstr << value;
    return sstr.str();
  }

  double makeDouble(string input)
  {
    istringstream stm;
    stm.str(input);
    double num;
    stm >> num;

    return(num);
  }
};

#endif
