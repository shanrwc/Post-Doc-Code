#ifndef SAMPLESELECTOR_H
#define SAMPLESELECTOR_H

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TMap.h"
#include "TParameter.h"
#include "TClonesArray.h"
#include "TIterator.h"


#include "RooDataSet.h"
#include "RooRealVar.h"

#include "Mass/CommonTools/normalizationTools.h"
#include "Mass/CommonTools/acceptanceTools.h"

using namespace std;
using namespace RooFit;

//This class handles reading in a sample or samples and returning datasets
//for use in pseudoexperiments.

class sampleSelector{

 public:
  
   sampleSelector()
   { 
     seed = 10;
     seed_int = 10;
     err_corr = 1.0;
     normalizationType_ = "";
     acceptanceType_ = "";
     useResampling = true;
     trueMass = 0;
     trueJES = 0;
     treeCut = "";
     systLabel = "";
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

     double totalFrac = 0.0;
     double tSigFrac = 0.0;

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
	   if (pieces[0].find("RANDOM_SEED") != string::npos && pieces.size() > 1)
	   {
	     seed = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("SEED_INTERVAL") != string::npos && pieces.size() > 1)
	   {
	     seed_int = makeDouble(pieces[1]);
	   }
	   if (pieces[0].find("PARABOLA_ERROR_CORRECTION") != string::npos && pieces.size() > 1)
	   {
	     err_corr = makeDouble(pieces[1]);
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
	   if (pieces[0].find("APPLY_CUT:") != string::npos && pieces.size() > 1)
	   {
	     treeCut = pieces[1];
	   }
	   if (pieces[0].find("SIGNAL_FILE:") != string::npos && pieces.size() > 3) 
	   { 
	     TFile *sfile = new TFile(pieces[1].c_str(),"READ");
	     if (sfile == NULL)
	     {
	       cout << "[ERROR}: Input file " << pieces[1] <<" not found!"<<endl;
	       return(2);
	     }
	     TTree* stree = (TTree*)(sfile->Get("LHCO/Particle_tree"));
	     if (stree == NULL)
	     {
	       cout << "[ERROR]: Signal particle tree not found!" << endl;
	       return(3);
	     }
	     TTree* sltree = (TTree*)(sfile->Get("Likelihood_tree"));
	     if (sltree == NULL)
	     {
	       cout << "[ERROR]: Signal likelihood tree not found!" << endl;
	       return(3);
	     }

	     //Check that sample and systematic label are consistant
	     bool sysMatch = checkSampleVSyst(pieces[1], pieces[3]);
	     if (!sysMatch)
	     {
	       cout << "[ERROR]: Mismatch between input file and systematic name!" << endl;
	       return(11);
	     }

	     sigData = convertTree(stree,sltree,"sigData",pieces[3]);
	     //Also, check truth information and store systematic label
	     systLabel = pieces[3];
	     if (trueMass + trueJES < 0.00001)
	     {
	       cout << "[ERROR]: Truth information not found!" << endl;
	       return(4);
	     }
	     tSigFrac = makeDouble(pieces[2]);
	   }
	   if (pieces[0].find("BACKGROUND_FILE:") != string::npos && pieces.size() > 3)
	   {
	     TFile* bfile = new TFile(pieces[1].c_str(),"READ");
	     if (bfile == NULL)
	     {
	       cout << "[ERROR]: Input file " << pieces[1] << " not found!" << endl;
	       return (2);
	     }

	     TTree* btree = (TTree*)bfile->Get("LHCO/Particle_tree");
	     if (btree == NULL)
	     {
	       cout << "[ERROR]: Background particle tree not found!" << endl;
	       return (3);
	     }
	     TTree* bltree = (TTree*)bfile->Get("Likelihood_tree");
	     if (bltree == NULL)
	     {
	       cout << "[ERROR]: Background likelihood tree not found!" << endl;
	       return(3);
	     }

	     RooDataSet* temp = convertTree(btree,bltree,"bkgData",pieces[3]);
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
       return(10); 
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
     if (fabs(1.0-(tSigFrac+totalFrac)) > 0.001)
     {
       cout << "[ERROR]: Given signal fraction does not equal 1-sum(background fractions)!" << endl;
       return(8);
     }
     nSigEvents = nEvents*(1.0-totalFrac);

     //Check that all samples are using the same weight type
     bool match = true;
     for (int i = 0; i != (int)weight_types.size(); ++i)
     {
       for (int j = i + 1; j != (int)weight_types.size(); ++j)
       {
	 if (weight_types[i].compare(weight_types[j]) == 0)
	 {
	   match = false;
	 }
       }
     }
     if (!match)
     {
       cout << "[ERROR]: Input samples are using different types of event weights!" << endl;
       return(9);
     }

     //Now fill sets with event numbers and weights, and write out pool information for the user
     //Start by writing out signal sample info
     int likeEntries = sigData->numEntries();
     cout<<"############################"<<endl;
     cout<<"Entries in signal sample likelihood file: "<<likeEntries<<endl;
     set<pair<int,double> > set_evt;
     double wTotal = 0;

     for (int i = 0; i< likeEntries; i++)
     {
       const  RooArgSet * row = sigData->get(i);
       set_evt.insert(pair<int,double>(((RooRealVar*)row->find("evt"))->getVal(),((RooRealVar*)row->find("weight"))->getVal()));
       wTotal += ((RooRealVar*)row->find("weight"))->getVal();
     }

     numEvtPools.push_back(pair<int,int>(nSigEvents,wTotal));
     totalHyp_ =  set_mass_.size()*set_jes_.size();
     copy(set_evt.begin(),set_evt.end(),back_inserter(sig_evts));    
     
     cout<<"# Found events: "<<sig_evts.size()<<endl;
     cout<<"# Found mass hyp: "<<set_mass_.size()<<endl;
     cout<<"# First mass hyp: "<<(double)(*set_mass_.begin())<<endl;
     cout<<"# Found JES hyp: "<<set_jes_.size()<<endl;
     cout<<"# First JES hyp: "<<(double)(*set_jes_.begin())<<endl;
     cout<<"############################"<<endl;
     
     //Now print out a limited set for the background files
     //     map<RooDataSet*, double>::iterator iter;
     for (iter = bkgDatas.begin(); iter != bkgDatas.end(); ++iter)
     {
       RooDataSet* temp = iter->first;
       int nEntries = temp->numEntries();
       cout << "Entries in background sample likelihood file: " << nEntries << endl;
       set_evt.clear();
       wTotal = 0;
       for (int i = 0; i != nEntries; ++i)
       {
	 const RooArgSet* row = temp->get(i);
	 set_evt.insert(pair<int,double>(((RooRealVar*)row->find("evt"))->getVal(),((RooRealVar*)row->find("weight"))->getVal()));
	 wTotal += ((RooRealVar*)row->find("weight"))->getVal();
       }
       cout << "# Found events: " << set_evt.size() << endl;
       numEvtPools.push_back(pair<int,int>(iter->second*nEvents,wTotal));
       cout << "# Fraction of Sample: " << iter->second << endl;
       
       vector<pair<int,double> > temp2;
       copy(set_evt.begin(),set_evt.end(),back_inserter(temp2));
       bkg_evts.push_back(temp2);
     }
     cout << "############################"<<endl;
   
     cout << "sampleSelector successfully initialized." << endl;
     return(0);
   } 

   RooDataSet* convertTree(TTree* etree, TTree* ltree, string name, string weightName)
   {
     RooRealVar evt("evt","evt",0,5000000000);
     RooRealVar run("run","run",0,5000000000);
     RooRealVar mass("mass","mass",150.,200.);
     RooRealVar jes("jes","jes",0.5,1.5);
     RooRealVar like("like","like",1.0e-40,100);
     RooRealVar like_err("like_err","like_err",0,1000);
     RooRealVar weight("weight","weight",-100,100);

     RooDataSet* sample = new RooDataSet(name.c_str(),"ProcessLikelihoods",RooArgSet(evt,run,mass,jes,like,like_err,weight));

     double eventNumber;
     double runNumber;
     double truthMass;
     double truthJES;
     double eventWeight;
     double evtW_num;
     double evtW_denom;
     vector<double> *pdf_weight = NULL;

     etree->SetBranchStatus("*",0);
     etree->SetBranchStatus("event_number",1);
     etree->SetBranchAddress("event_number",&eventNumber);
     etree->SetBranchStatus("run_number",1);
     etree->SetBranchAddress("run_number",&runNumber);
     etree->SetBranchStatus("trueMass",1);
     etree->SetBranchAddress("trueMass",&truthMass);
     etree->SetBranchStatus("trueJES",1);
     etree->SetBranchAddress("trueJES",&truthJES);
     etree->SetBranchStatus("event_weight",1);
     etree->SetBranchAddress("event_weight",&eventWeight);

     if (weightName.compare("NOMINAL") == 0)
     {
       evtW_num = 1.0;
       evtW_denom = 1.0;
     }
     else if (weightName.compare("TRIGGER_UP") == 0)
     {
       etree->SetBranchStatus("trigger_weight_up",1);
       etree->SetBranchAddress("trigger_weight_up",&evtW_num);
       etree->SetBranchStatus("trigger_weight",1);
       etree->SetBranchAddress("trigger_weight",&evtW_denom);
     }
     else if (weightName.compare("TRIGGER_DOWN") == 0)
     {
       etree->SetBranchStatus("trigger_weight_down",1);
       etree->SetBranchAddress("trigger_weight_down",&evtW_num);
       etree->SetBranchStatus("trigger_weight",1);
       etree->SetBranchAddress("trigger_weight",&evtW_denom);
     }
     else if (weightName.compare("LEP_ID_UP") == 0)
     {
       etree->SetBranchStatus("lepID_weight_up",1);
       etree->SetBranchAddress("lepID_weight_up",&evtW_num);
       etree->SetBranchStatus("lepID_weight",1);
       etree->SetBranchAddress("lepID_weight",&evtW_denom);
     }
     else if (weightName.compare("LEP_ID_DOWN") == 0)
     {
       etree->SetBranchStatus("lepID_weight_down",1);
       etree->SetBranchAddress("lepID_weight_down",&evtW_num);
       etree->SetBranchStatus("lepID_weight",1);
       etree->SetBranchAddress("lepID_weight",&evtW_denom);
     }
     else if (weightName.compare("BTAG_UP") == 0)
     {
       etree->SetBranchStatus("btag_weight_up",1);
       etree->SetBranchAddress("btag_weight_up",&evtW_num);
       etree->SetBranchStatus("btag_weight",1);
       etree->SetBranchAddress("btag_weight",&evtW_denom);
     }
     else if (weightName.compare("BTAG_DOWN") == 0)
     {
       etree->SetBranchStatus("btag_weight_down",1);
       etree->SetBranchAddress("btag_weight_down",&evtW_num);
       etree->SetBranchStatus("btag_weight",1);
       etree->SetBranchAddress("btag_weight",&evtW_denom);
     }
     else if (weightName.compare("MISTAG_UP") == 0)
     {
       etree->SetBranchStatus("mistag_weight_up",1);
       etree->SetBranchAddress("mistag_weight_up",&evtW_num);
       etree->SetBranchStatus("btag_weight",1);
       etree->SetBranchAddress("btag_weight",&evtW_denom);
     }
     else if (weightName.compare("MISTAG_DOWN") == 0)
     {
       etree->SetBranchStatus("mistag_weight_down",1);
       etree->SetBranchAddress("mistag_weight_down",&evtW_num);
       etree->SetBranchStatus("btag_weight",1);
       etree->SetBranchAddress("btag_weight",&evtW_denom);
     }
     else if (weightName.compare("TOPPT") == 0)
     {
       evtW_num = 1.0;
       etree->SetBranchStatus("topPt_weight",1);
       etree->SetBranchAddress("topPt_weight",&evtW_denom);
     }
     else if (weightName.compare("BFRAG") == 0)
     {
       evtW_num = 1.0;
       etree->SetBranchStatus("bfrag_weight",1);
       etree->SetBranchAddress("bfrag_weight",&evtW_denom);
     }
     else if (weightName.compare("PILE_UP_UP") == 0)
     {
       etree->SetBranchStatus("pu_weight_up",1);
       etree->SetBranchAddress("pu_weight_up",&evtW_num);
       etree->SetBranchStatus("pu_weight",1);
       etree->SetBranchAddress("pu_weight",&evtW_denom);
     }
     else if (weightName.compare("PILE_UP_DOWN") == 0)
     {
       etree->SetBranchStatus("pu_weight_down",1);
       etree->SetBranchAddress("pu_weight_down",&evtW_num);
       etree->SetBranchStatus("pu_weight",1);
       etree->SetBranchAddress("pu_weight",&evtW_denom);
     }
     else if (weightName.find("PDF_") != string::npos)
     {
       etree->SetBranchStatus("pdf_weight",1);
       etree->SetBranchAddress("pdf_weight",&pdf_weight);
       evtW_num = 1.0;
       evtW_denom = 1.0;
     }
     else if (weightName.find("CROSS_SECTION") != string::npos || weightName.find("ACCEPTANCE") != string::npos)
     {
       evtW_num = 1.0;
       evtW_denom = 1.0;
     }
     else
     {
       cout << "Weight type " << weightName << " unknown: using NOMINAL" << endl;
       evtW_num = 1.0;
       evtW_denom = 1.0;
     }
     double leventNumber;
     double lrunNumber;
     TH2D* llhisto = new TH2D();
     TH2D* luhisto = new TH2D();

     ltree->SetBranchStatus("*",1);
     ltree->SetBranchAddress("event_number",&leventNumber);
     ltree->SetBranchAddress("run_number",&lrunNumber);
     ltree->SetBranchAddress("likes",&llhisto);
     ltree->SetBranchAddress("luncerts",&luhisto);
     etree->AddFriend(ltree);

     if(treeCut.length() > 0)
     {
       TTree* tempTree = etree;
       etree = tempTree->CopyTree(treeCut.c_str());
     }

     int eventNr = 0;
     while(etree->GetEntry(eventNr++))
     {
       evt = eventNumber;
       run = runNumber;
       weight = eventWeight*evtW_num/evtW_denom;

       if (weightName.find("PDF_") != string::npos)
       {
	 int nPDF = int(makeDouble(weightName.substr(weightName.length()-1,1)));
	 weight = eventWeight*pdf_weight->at(nPDF);
       }

       if (name.find("sig") != string::npos)
       {
	 trueMass = truthMass;
	 //	 trueJES = truthJES;
	 trueJES = 1.0;

	 if (truthMass == -1 || truthJES == -1)
	 {
	   trueMass = 172.3;
	   trueJES = 1.0;
	 }
       }

       if (eventNumber == leventNumber && runNumber == lrunNumber)
       {
	 TAxis* xaxis = llhisto->GetXaxis();
	 TAxis* yaxis = llhisto->GetYaxis();
	 
	 int nX = xaxis->GetNbins();
	 int nY = yaxis->GetNbins();

	 for (int i = 1; i <= nX; ++i)
	 {
	   for (int j = 1; j <= nY; ++j)
	   {
	     mass = xaxis->GetBinCenter(i);
	     jes = yaxis->GetBinCenter(j);
	     like = llhisto->GetBinContent(i,j);
	     like_err = luhisto->GetBinContent(i,j);

	     sample->add(RooArgSet(evt,run,mass,jes,like,like_err,weight));
	   }
	 }
       }
     }

     return(sample);
   }

   int getTruthInfo(double& tMass, double& tJES)
   {

     if(trueMass != 0 && trueJES != 0)
     {

       tMass = trueMass;
       tJES = trueJES;

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
    RooDataSet* ensemble = selectEvents(sigData,sig_evts,nSigEvents,seed);

    map<RooDataSet*,double>::iterator iter;
    int n = 0;
    for (iter = bkgDatas.begin(); iter != bkgDatas.end(); ++iter)
    {
      seed += seed_int;
      ensemble->append(*selectEvents(iter->first,bkg_evts[n],nEvents*(iter->second),seed));
      ++n;
    }

    //It needs to combine all the samples and then combine event likelihoods into process likelihoods
    //And normalize/add acceptance to everything

    RooDataSet* normEnsemble = mergeEventLLs(ensemble);

    //And then it needs to return the dataset
    return(normEnsemble);
  }
 
  RooDataSet* selectEvents( RooDataSet* dataset, vector<pair<int,double> > evtNums, int nEvts, int seed )
  {
    //Select the requested number of events with replacement and return the datset
    RooRealVar evt("evt","evt", 0, 500000000);
    RooRealVar mass("mass","mass", 150., 200.);
    RooRealVar jes("jes","jes", 0.50, 1.50);
    RooRealVar like("like", "like", 1.0e-40, 100);
    RooRealVar nl_like("nl_like", "nl_like", 0, 1000000);
    RooRealVar nl_err("nl_err", "nl_err", 0, 1000000);
    RooRealVar like_err("like_err", "like_err", 0, 1000);
    RooRealVar run("run","run",0,5000000000);
    RooRealVar weight("weight","weight",-100,100);
    
    std::stringstream cut2;  

    RooDataSet*  selData  = new RooDataSet("selData","Selected Data", RooArgSet(evt, run, mass, jes, like, like_err));
   
    ++count;
    TRandom2 myrnd(count + seed);
    
    double maxWeight;
    double minWeight;
    getWeightLimits(evtNums,maxWeight,minWeight);

    while ((selData->numEntries()/ totalHyp_  < nEvts) && (evtNums.size() > 0))
    {
      int rEvts = min(100, nEvts-selData->numEntries()/ totalHyp_ );
      string cut = "(";
      
      for (int i = 0; i < rEvts; ++i)
      {
	if (evtNums.size() == 0) continue;
	else
	{
	  int tempN = (int)myrnd.Uniform(evtNums.size());
	  pair<int,double> entry = evtNums[tempN];
	  while (entry.second < (double)myrnd.Uniform(minWeight,maxWeight))
	  {
	    tempN = (int)myrnd.Uniform(evtNums.size());
	    entry = evtNums[tempN];
	  }
	  
	  cut += "evt==" + makeString(entry.first,0) + " || ";
	  
	  if (!useResampling)
	  {
	    evtNums.erase(evtNums.begin()+tempN);
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
    
    //    xsec::normalizationTools normTool(normalizationType_);
    //    RooDataSet* crossData = normTool.loadNormalizationFile();
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
	  nll_err+=((RooRealVar*)row->find("like_err"))->getVal()/((RooRealVar*)row->find("like"))->getVal();
	  //	  cout << nll_err << endl;
	}
	
	nl_like = nll; nl_err = nll_err;

	//Add the normalization and acceptance, if using
	cut2.str("");
	cut2<<"mass>"<<*m_it-0.01<<" && mass<"<<*m_it+0.01;
	RooDataSet* massCrossData = (RooDataSet*) crossData->reduce(RooArgSet(cross,cross_err),cut2.str().c_str());
	if (massCrossData->numEntries() > 1)
	{
	  cout << "[ERROR]: More than one normalization!" << endl;
	  return 0;
	}
	else if (massCrossData->numEntries() < 1)
	{
	  cout << "[ERROR]: No normalization found!" << endl;
	  return 0;
	}

	const RooArgSet* rowCross = massCrossData->get(0);
	double norm = ((RooRealVar*)rowCross->find("cross"))->getVal();
	double norm_err = ((RooRealVar*)rowCross->find("cross_err"))->getVal();

	if (applyAcceptance)
	{
	  //	  acc::acceptanceTools accTool(acceptanceType_);
	  //	  accTool.loadAcceptanceFile();

	  accTool->calculateAcceptance(*m_it,*j_it,accept);
	  accTool->calculateAcceptanceError(*m_it,*j_it,accept_err);
	}

	if (systLabel.compare("CROSS_SECTION_UP") == 0)
	{
	  norm_nl_like = nl_like + nEvents*log(norm+norm_err) + nEvents*log(accept);
	  norm_nl_err = (nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept))/err_corr;
	}
	else if (systLabel.compare("CROSS_SECTION_DOWN") == 0)
	{
	  norm_nl_like = nl_like + nEvents*log(norm-norm_err) + nEvents*log(accept);
	  norm_nl_err = (nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept))/err_corr;
	}
	else if (systLabel.compare("ACCEPTANCE_UP") == 0)
	{
	  norm_nl_like = nl_like + nEvents*log(norm) + nEvents*log(accept+accept_err);
	  norm_nl_err = (nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept))/err_corr;
	}
	else if (systLabel.compare("ACCEPTANCE_DOWN") == 0)
	{
	  norm_nl_like = nl_like + nEvents*log(norm) + nEvents*log(accept-accept_err);
	  norm_nl_err = (nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept))/err_corr;
	}
	else
	{
	  norm_nl_like = nl_like + nEvents*log(norm) + nEvents*log(accept);
	  norm_nl_err = (nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept))/err_corr;
	}
	//	cout << "Final Uncertainties: " << nl_err << " " << nEvents*fabs(norm_err/norm) << " " << nEvents*(accept_err/accept) << endl;
	//	cout << "Final Likelihood: " << mass.getVal() << " " << norm_nl_like.getVal() << " " << norm_nl_err.getVal() << endl;
	//	norm_nl_err = nl_err + nEvents*fabs(norm_err/norm) + nEvents*(accept_err/accept) - nEvents*fabs(norm_err/norm) - nEvents*(accept_err/accept);
	//	cout << nl_like << " " << nEvents*log(norm) << " " << nEvents*log(accept) << endl;

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
      //      dataCut<<"mass=="<<*m_it<<" && "<<"jes<"<<*j_it+0.01<<" && "<<"jes>"<<*j_it-0.01;
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

 void getWeightLimits(vector<pair<int,double> > evtInfo,double& max, double& min)
 {
   double tempMax = 0;
   double tempMin = 0;
   for (int i = 0; i != (int)evtInfo.size(); ++i)
   {
     if (evtInfo[i].second > tempMax) tempMax = evtInfo[i].second;
     if (evtInfo[i].second < tempMin) tempMin = evtInfo[i].second;
   }

   if (tempMin < 0)
   {
     cout << "[getWeightLimits]: WARNING-NEGATIVE EVENT WEIGHT FOUND!" << endl;
   }
   max = tempMax;
   min = tempMin;
 }

 bool checkSampleVSyst(string sample, string label)
 {
   if (sample.find("nominal") != string::npos)
   {
     if (label.find("NOMINAL") != string::npos || label.find("TRIGGER") != string::npos || label.find("LEP_ID") != string::npos || label.find("BTAG") != string::npos || label.find("MISTAG") != string::npos || label.find("CROSS_SECTION") != string::npos || label.find("ACCEPTANCE") != string::npos || label.find("TOPPT") != string::npos || label.find("PILE_UP") != string::npos || label.find("BFRAG") != string::npos)
     {
       return(true);
     }
     else
     {
       return(false);
     }
   }
   else
   {
     string snippet1 = label.substr(0,label.find("_"));
     string snippet2 = label.substr(label.rfind("_")+1,label.length());
     transform(snippet1.begin(),snippet1.end(),snippet1.begin(),::tolower);
     transform(snippet2.begin(),snippet2.end(),snippet2.begin(),::tolower);
     //     cout << snippet1 << " " << snippet2 << endl;
     if (sample.find(snippet1) != string::npos && sample.find(snippet2) != string::npos)
     {
       return (true);
     }
     else
     {
       return(false);
     }
   }
 }

 string getSystematicLabel()
 {
   return(systLabel);
 }

 double getSeed()
 {
   return(seed);
 }

 private:
  int count;
  int totalHyp_;
  double trueMass;
  double trueJES;

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
  double seed;
  double seed_int;
  double err_corr;
  bool useResampling;
  string treeCut;
  string systLabel;
  int nEvents;
  int nSigEvents;
  vector<string> weight_types;
  RooDataSet* sigData;
  map<RooDataSet*,double> bkgDatas;
    
  //vectors to store lists of event numbers and weights (useful to have when making pseudoexperiments, I promise)
  vector<pair<int,double> > sig_evts;
  vector< vector<pair<int,double> > > bkg_evts;

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
