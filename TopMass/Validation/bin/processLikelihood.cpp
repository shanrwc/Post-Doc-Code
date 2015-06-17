#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <set>

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TMath.h"
#include "TH2.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMap.h"
#include "TParameter.h"
#include "TClonesArray.h"
#include "TArrayD.h"

using namespace std;

template<typename T>
std::string makeString(T const& value)
{
  std::stringstream sstr;
  sstr << fixed;
  sstr << setprecision(3);
  sstr << value;
  return sstr.str();
}

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

double makeDouble(string input)
{
  istringstream stm;
  stm.str(input);
  double num;
  stm >> num;

  return(num);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//These lines must be edited by the user
const int numJES = 1;
double JESHyp[numJES] = {1.0};

const int numMass = 7;
//double MassHyp[numMass] = {154.0,158.0,162.0,166.0,170.0,174.0,178.0};//166
//double MassHyp[numMass] = {159.0,163.0,167.0,171.0,175.0,179.0,183.0,};//169
//double MassHyp[numMass] = {161.0,165.0,169.0,173.0,177.0,181.0,185.0};//171
//double MassHyp[numMass] = {162.0,166.0,170.0,174.0,178.0,182.0,186.0};//172
//double MassHyp[numMass] = {161.0,165.0,169.0,173.0,177.0,181.0,185.0};//173
double MassHyp[numMass] = {165.0,169.0,173.0,177.0,181.0,185.0,189.0};//175
//double MassHyp[numMass] = {166.0,170.0,174.0,178.0,182.0,186.0,190.0};//178

string inDir = "/afs/cern.ch/user/s/swalch/likelihoods/Electron/ttbar_signal/nominal/175_5";

//The order of input files here must match the order of JES hypotheses
//given above.
string fileList[numJES] = {inDir+"/JES100/mergedFile.txt"};
//This outfile should contain a petuple, and it will be updated to include
//a second TTree object with the likelihoods
string outFile = "LHCOTree_175.5_1.00_TTbar_electron_01_cutset_nominal.root";

void processLikelihood()
{
  //Create two maps to hold the likelihoods and their uncertainties
  map<string,TH2D*> likes;
  map<string,TH2D*> uncerts;

  double minL = 100;
  double maxL = 0;
  double minU = 100;
  double maxU = 0;

  for (int nJES = 0; nJES != numJES; ++nJES)
  {
    ifstream input (fileList[nJES].c_str());
    string line;

    if (input.is_open())
    {
      while (input.good())
      {
	getline(input,line);
	vector<string> pieces = splitString(line);
	if ((int)pieces.size() < 4) continue;
	//split up run and event number
	string runEvtNum = pieces[0];
	double tlike = makeDouble(pieces[3]);
	double tuncert = makeDouble(pieces[4]);
	if (pieces[1].find("format") == string::npos && pieces[1].find("Weight") == string::npos && (tlike == 0 || tuncert == 0))
	{
	  cout << line << endl;
	  //	  cout << "Event " << runEvtNum << " has zero likelihood/uncertainty!" << endl;
	}
	if (tlike > maxL) maxL = tlike;
	if (tlike < minL) minL = tlike;
	if (tuncert > maxU) maxU = tuncert;
	if (tuncert < minU) minU = tuncert;

	if (likes.find(runEvtNum) != likes.end() && uncerts.find(runEvtNum) != uncerts.end())
	{
	  //Update histograms with additional likelihoods and uncertainties
	  int mindex = makeDouble(pieces[1])-1;
	  likes[runEvtNum]->Fill(MassHyp[mindex],JESHyp[nJES],tlike);
	  uncerts[runEvtNum]->Fill(MassHyp[mindex],JESHyp[nJES],tuncert);
	}
	else if (likes.find(runEvtNum)  == likes.end() && uncerts.find(runEvtNum) == uncerts.end())
	{
	  //Add the histograms if they cannot be found
	  string llname = pieces[0] + "_likes";
	  double mSpace = 1;
	  if (numMass > 1) mSpace = MassHyp[1]-MassHyp[0];
	  double jesSpace = 0.1;
	  if (numJES > 1) jesSpace = JESHyp[1]-JESHyp[0];
	  TH2D* lltemp = new TH2D(llname.c_str(),llname.c_str(),numMass,MassHyp[0]-mSpace/2,MassHyp[numMass-1]+mSpace/2,numJES,JESHyp[0]-jesSpace/2,JESHyp[numJES-1]+jesSpace/2);
	  int mindex = makeDouble(pieces[1])-1;
	  lltemp->Fill(MassHyp[mindex],JESHyp[nJES],tlike);
	  likes.insert(pair<string,TH2D*>(runEvtNum,lltemp));

	  string luname = pieces[0] + "_uncerts";
	  TH2D* lutemp = new TH2D(luname.c_str(),luname.c_str(),numMass,MassHyp[0]-mSpace/2,MassHyp[numMass-1]+mSpace/2,numJES,JESHyp[0]-jesSpace/2,JESHyp[numJES-1]+jesSpace/2);
	  lutemp->Fill(MassHyp[mindex],JESHyp[nJES],tuncert);
	  uncerts.insert(pair<string,TH2D*>(runEvtNum,lutemp));
	}
	else
	{
	  //Spit out an error if only one of the two histograms can be found
	  cout << "[ERROR]: Macro is not loading likelihoods correctly!" << endl;
	  return;
	}
      }//closes while-loop over input file's goodness
    }//closes if-statment of input file's openness
  }//closes loop over JESs/input files

  //Now, with all the likelihoods in hand, make an output tree
  TFile* file = new TFile(outFile.c_str(),"UPDATE");
  if (file == NULL)
  {
    cout << "input file not found" << endl;
    return;
  }
  TTree* evtTree = (TTree*)(file->Get("LHCO/Particle_tree"));
  if (evtTree == NULL)
  {
    cout << "input event tree not found" << endl;
    return;
  }
  double run_Number;
  evtTree->SetBranchStatus("*",0);
  evtTree->SetBranchStatus("run_number",1);
  evtTree->SetBranchAddress("run_number",&run_Number);
  evtTree->GetEntry(0);

  //  file->cd("LHCO");
  double eventNumber;
  double runNumber;
  TH2D* llhisto = new TH2D();
  TH2D* luhisto = new TH2D();
  TTree* lltree = new TTree("Likelihood_tree","tree containing event likelihoods");
  lltree->Branch("event_number",&eventNumber);
  lltree->Branch("run_number",&runNumber);
  lltree->Branch("likes",&llhisto);
  lltree->Branch("luncerts",&luhisto);

  map<string,TH2D*>::iterator iter;
  for (iter = likes.begin(); iter != likes.end(); ++iter)
  {
    string temp = iter->first;
    size_t pos = temp.find(".");
    if (pos != string::npos)
    {
      runNumber = makeDouble(temp.substr(0,pos));
      eventNumber = makeDouble(temp.substr(pos+1,100));
    }
    else
    {
      runNumber = run_Number;
      eventNumber = makeDouble(temp);
    }
   
    llhisto = iter->second;
    luhisto = uncerts[iter->first];

    lltree->Fill();
  }

  //Finally, attach the two trees together and save everything
  lltree->BuildIndex("event_number","run_number");
  evtTree->AddFriend(lltree);

  //Write out some useful information
  cout << "Event tree contains " << evtTree->GetEntries() << " events." << endl;
  cout << "Likelihood tree contains " << lltree->GetEntries() << " events." << endl;
  if (evtTree->GetEntries() != lltree->GetEntries())
  {
    cout << "Warning: possible mismatch of samples" << endl;
  }
  cout << "Likelihood range: " << minL << " - " << maxL << endl;
  cout << "Uncertainty range: " << minU << " - " << maxU << endl;

  lltree->Write();
  file->Close();
}//closes processLikelihood function

#ifndef __CINT__

int main (int argc, const char* argv[])
{
  processLikelihood();

  return(0);
}

#endif
