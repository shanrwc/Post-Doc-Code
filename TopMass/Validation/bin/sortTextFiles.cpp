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

string updateHypNum(string input, int num)
{
  vector<string> pieces = splitString(input);
  if ((int)pieces.size() != 5)
  {
    cout << "Trying to rewrite the wrong line!" << endl;
    return(" ");
  }
  string newline = pieces[0]+" "+makeString(num)+" "+pieces[2]+" "+pieces[3]+" "+pieces[4];
  return(newline);
}

int indexHyp(double mass)
{
  double mhypos [33] = {154.0,157.0,158.0,159.0,160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0,171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0,181.0,182.0,183.0,184.0,185.0,186.0,187.0,190.0};

  int index = -1;
  for (int i = 0; i != 33; ++i)
  {
    if (mhypos[i] == mass) index = i;
  }

  return(index);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//This macro takes text files produces with several sets of mass hypotheses
//and sorts them into sets that match those of signal samples.  It is meant
//to be used on w+jets, z+jets, and similar samples.

void sortTextFiles(int argc, const char* argv[])
{
  //Build array of vectors that will hold lines of text
  const int nHYPOS = 33;
  double mhypos [nHYPOS] = {154.0,157.0,158.0,159.0,160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0,171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0,181.0,182.0,183.0,184.0,185.0,186.0,187.0,190.0};
  vector<string> lines [nHYPOS];
  for (int i = 0; i != nHYPOS; ++i)
  {
    vector<string> temp;
    lines[i] = temp;
  }

  double stMasses [7] = {166.0,169.0,171.0,172.0,173.0,175.0,178.0};

  //Loop over input files, sorting lines of text into appropriate vectors
  ifstream infile;
  string line;
  for (int i = 1; i < argc; ++i)
  {
    cout << "Accessing file" << argv[i] << endl;
    infile.open(argv[i]);

    if (!infile) {
      cout << "FILE NOT FOUND! Skipping File!" << endl;
      continue;
  }
    if (infile.is_open())
    {
      while (infile.good())
      {
	getline(infile,line);
	vector<string> pieces = splitString(line);
	if ((int)pieces.size() != 5) continue;
	int index = makeDouble(pieces[1]);
	(lines[index-1]).insert(line);
      }//closes while-loop over file's goodness
    }//closes if-statement over file's opennes

  //Now loop over output files, and write out the appropriate lines to each one
  //Note that hypothesis numbers will need to get REASSIGNED appropriately
    for (int nMass = 0; nMass != 7; ++nMass)
    {
      double tMass = stMasses[nMass];
      //open output file
      ofstream outfile;
      outfile.open(("mergedFile_"+makeString(tMass)+".txt").c_str());

      //now loop over mass hypotheses
      for (int nHyp = tMass-12; nHyp <=tMass+12; nHyp += 4)
      {
	//get vector of lines
	vector<string> collect = lines[indexHyp(nHyp)];
	for (int i = 0; i != (int)collect.size(); ++i)
	{
	  outfile << collect[i] << endl;
	}
      }

      outfile.close();
    }

}

#ifndef __CINT__

int main (int argc, const char* argv[] ) {

  sortTextFiles( argc, argv );

  return(0);
}

#endif
