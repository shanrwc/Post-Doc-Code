//SYS
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <utility>
#include <vector>
#include <map>
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
#include "TMatrixDSym.h"
//USER
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "Mass/CommonTools/acceptanceTools.h"
#include "Mass/CommonTools/normalizationTools.h"
#include "Mass/CommonTools/calib_sampleSelector.h"
#include "interface/tdrStyle.h"

using namespace std;

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

template<typename T>
string makeString(T const& value, int prec)
{
  stringstream sstr;
  sstr << fixed;
  sstr << setprecision(prec);
  sstr << value;
  return sstr.str();
}

string convertLabel(string stub)
{
  if (stub.compare("CALIBRATION") == 0)
  {
    return("Calibration");
  }
  else if (stub.compare("TRIGGER") == 0)
  {
    return("Trigger Efficiency");
  }
  else if (stub.compare("LEP_ID") == 0)
  {
    return("Lepton ID");
  }
  else if (stub.compare("BTAG") == 0)
  {
    return("$b$-tag Efficiency");
  }
  else if (stub.compare("MISTAG") == 0)
  {
    return("Mistag Rate");
  }
  else if (stub.compare("CROSS_SECTION") == 0)
  {
    return("Normalization");
  }
  else if (stub.compare("ACCEPTANCE") == 0)
  {
    return("Acceptance");
  }
  else if (stub.compare("SCALE") == 0)
  {
    return("Renormalization-Factorization");
  }
  else if (stub.compare("MATCHING") == 0)
  {
    return("ME-PS Matching");
  }
  else if (stub.compare("COLOR") == 0)
  {
    return("Color Reconnection");
  }
  else if (stub.compare("TOPPT") == 0)
  {
    return ("Top $p_{T}$~Reweighting");
  }
  else if (stub.compare("PILEUP") == 0)
  {
    return ("Pile-Up");
  }
  else if (stub.compare("UNDER_EVT") == 0)
  {
    return("Underlying Event");
  }
  else if (stub.compare("MPFBIAS") == 0)
  {
    return("In-Situ Group");
  }
  else if (stub.compare("RELFSR") == 0)
  {
    return("Inter-Calibration Group");
  }
  else
  {
    return (stub);
  }

}

int calcResults(string inFile)
{
  // TDRStyle tdrSTY;
  // tdrSTY.setTDRStyle();

  //Set the offset/central value used for plotting
  const double OFFSET = 172.5;

  //First, read the input file and sort out the filenames
  ifstream drivecard (inFile.c_str());
  string line;

  string outfilename;
  string outrootname;
  TFile* datafile = 0;
  bool calcUncert = false;
  vector<TFile*> calibfiles;
  map<string, vector<TFile*> > systfiles;
  map<string,double> upSysts;
  map<string,double> downSysts;

  if (drivecard.is_open())
  {
    while (drivecard.good())
    {
      getline(drivecard,line);

      string first = line.substr(0,1);
      //skip commented lines
      if (first.compare("#") != 0)
      {
	vector<string> pieces = splitString(line);
	if (pieces.size() == 0) continue;
	if (pieces[0].find("OUTPUT_TABLE") != string::npos && pieces.size() > 1)
	{
	  outfilename = pieces[1];
	}
	if (pieces[0].find("OUTPUT_HISTOGRAMS") != string::npos && pieces.size() > 1)
	{
	  outrootname = pieces[1];
	}
	if (pieces[0].find("DATA") != string::npos && pieces.size() > 1)
	{
	  datafile = new TFile(pieces[1].c_str(),"READ");
	  if (datafile == NULL)
	  {
	    cout << "[ERROR]: Data file not found!" << endl;
	    return(1);
	  }
	}
	if (pieces[0].find("CALIBRATION") != string::npos && pieces.size() > 1)
	{
	  TFile* temp = new TFile(pieces[1].c_str(),"READ");
	  if (temp == NULL)
	  {
	    cout << "[ERROR]: Calibration file not found!" << endl;
	    return(2);
	  }
	  calibfiles.push_back(temp);
	}
	if (pieces[0].find("CALIBRATE_UNCERTAINTY") != string::npos && pieces.size() > 1)
	{
	  if (pieces[1].find("yes") != string::npos) calcUncert = true;
	}
	if (pieces[0].find("SYSTEMATIC") != string::npos && pieces.size() > 2)
	{
	  TFile* temp = new TFile(pieces[1].c_str(),"READ");
	  if (temp == NULL)
	  {
	    cout <<"[ERROR]: Systematic file not found!" << endl;
	    return(3);
	  }

	  if (systfiles.find(pieces[2]) != systfiles.end())
	  {
	    (systfiles[pieces[2]]).push_back(temp);
	  }
	  else
	  {
	    vector<TFile*> tempVec;
	    tempVec.push_back(temp);
	    systfiles.insert(pair<string,vector<TFile*> >(pieces[2],tempVec));
	  }
	}

      }
    }//closes while-loop over drive file's goodness
  }//closes if-statment over drive file's openness

  /////////////////////////////////////////////////////////////////////////
  //Create latex table to store results
  ofstream outfile;
  outfile.open(outfilename.c_str());
  outfile << fixed << setprecision(3);
  outfile << "\\begin{tabular}{|c|c|c|}" << endl;
  outfile << "\\hline" << endl;

  /////////////////////////////////////////////////////////////////////////
  //Create root file to store results
  TFile* outrfile = new TFile(outrootname.c_str(),"RECREATE");

  /////////////////////////////////////////////////////////////////////////
  //If available, grab the data result and store it
  double rawMass = 173.2;
  double rawUncert = 0.2;
  if (datafile != 0)
  {
    TCanvas* dresult = (TCanvas*)datafile->Get("canMass0");
    //    cout << dresult << endl;
    TPaveText* dtext = (TPaveText*)dresult->GetPrimitive("massResult");
    TText* dsnippet = dtext->GetLine(0);
    string dinfo = dsnippet->GetTitle();

    //    outfile << "Uncalibrated Result & " << dinfo << " \\\\" << endl;
    rawMass = makeDouble(dinfo.substr(0,6));
    rawUncert = makeDouble(dinfo.substr(10,4));
    outfile << "\\multicolumn{3}{|c|}{Uncalibrated Result & " << rawMass << " $\\pm " << rawUncert << "}\\\\" << endl;
  }
  else
  {
    outfile << "\\multicolumn{3}{|c|}{Using dummy values: " << rawMass << " " << rawUncert << "}\\\\" << endl;
  }
  /////////////////////////////////////////////////////////////////////////
  //If applicable, create calibration curve, calibrate data result,
  //and estimate systematic from curve.
  if ((int)calibfiles.size() > 1)
  {
    //Loop over calib files and store masses and uncerts in arrays
    int cSize = (int)calibfiles.size();
    double xvals [cSize];
    double xerrs [cSize];
    double yvals [cSize];
    double yerrs [cSize];
    double uvals [cSize];
    double uerrs [cSize];
    map<double,double> nomPts;
    for (int i = 0; i != cSize; ++i)
    {
      //First get the true mass value
      TTree* tree = (TTree*)calibfiles[i]->Get("psexpInfoTree");
      double trueMass;
      tree->SetBranchStatus("*",0);
      tree->SetBranchStatus("trueMass",1);
      tree->SetBranchAddress("trueMass",&trueMass);
      tree->GetEntry(0);

      //Now get the measured mass with resampling uncertainty
      TCanvas* cresult = (TCanvas*)calibfiles[i]->Get(("Mass_"+makeString(trueMass,1)).c_str());
      TPaveStats* ctext = (TPaveStats*)cresult->GetPrimitive("massMeas");
      TText* csnippet = ctext->GetLine(0);
      string cinfo = csnippet->GetTitle();
      vector<string> cpieces = splitString(cinfo);

      xvals[i] = trueMass - OFFSET;
      xerrs[i] = 0.0;
      yvals[i] = makeDouble(cpieces[0]) - OFFSET;
      yerrs[i] = makeDouble(cpieces[2]);
      nomPts.insert(pair<double,double>(trueMass,makeDouble(cpieces[0])));

      //And get the pull widths in case they are needed
      TCanvas* cresult2 = (TCanvas*)calibfiles[i]->Get(("MassPull_"+makeString(trueMass,1)).c_str());
      ctext = (TPaveStats*)cresult2->GetPrimitive("masspullMeas");
      csnippet = ctext->GetLine(1);
      cinfo = csnippet->GetTitle();

      cpieces = splitString(cinfo);
      uvals[i] = makeDouble(cpieces[2]);
      uerrs[i] = makeDouble(cpieces[4]);
    }

    //Make a TGraphErrors with the information, fit a straight line to it, and make
    //a pretty plot
    outrfile->cd();
    string clabel = "calibcurve_nodata";
    TCanvas* cpad = new TCanvas(clabel.c_str(),clabel.c_str(),14,33,600,600);
    cpad->SetTitle(clabel.c_str());
    cpad->SetFillColor(0);
    cpad->SetLeftMargin(0.15);
    cpad->SetBottomMargin(0.14);

    TGraphErrors* cgraph = new TGraphErrors(cSize,xvals,yvals,xerrs,yerrs);
    cgraph->SetTitle("");
    cgraph->SetMarkerColor(kBlack);
    cgraph->SetMarkerStyle(20);
    cgraph->GetXaxis()->SetTitle(("m_{t}^{gen}-"+makeString(OFFSET,1)+" [GeV]").c_str());
    cgraph->GetXaxis()->SetTitleSize(0.05);
    cgraph->GetXaxis()->SetTitleOffset(1.0);
    cgraph->GetYaxis()->SetTitle(("m_{t}^{MEM}-"+makeString(OFFSET,1)+" [GeV]").c_str());
    cgraph->GetYaxis()->SetTitleSize(0.05);
    cgraph->GetYaxis()->SetTitleOffset(1.2);
    cgraph->Draw("ap");

    TF1* func = new TF1("func","[0]+[1]*x",-15.0,15.0);
    func->SetParameter(0,0.0);
    func->SetParameter(1,1.0);
    func->SetLineColor(kBlue);
    func->SetLineWidth(2);

    // TF1* deff = new TF1("deff","[0]+[1]*x",-15.0,15.0);
    // deff->SetParameter(0,0.0);
    // deff->SetParameter(1,1.0);
    // deff->SetLineColor(kBlack);
    // deff->SetLineWidth(2);

    TFitResultPtr cfit = cgraph->Fit(func,"RMSV");
    TPaveText* ctext = new TPaveText(0.3,0.65,0.55,0.8,"NDC");
    ctext->SetFillColor(0);
    ctext->SetLineColor(kBlue);
    ctext->SetLineWidth(1);
    ctext->SetTextSize(0.03);
    ctext->AddText(("f(m_{t}^{MEM}-"+makeString(OFFSET,1)+")=a (m_{t}^{gen}-"+makeString(OFFSET,1)+")+b").c_str());
    double aval = cfit->Value(1);
    double aerr = cfit->ParError(1);
    string bit = "a = "+makeString(aval,3) + "#pm"+makeString(aerr,3);
    ctext->AddText(bit.c_str());
    double bval = cfit->Value(0);
    double berr = cfit->ParError(0);
    bit = "b = "+makeString(bval,3) + "#pm"+makeString(berr,3);
    ctext->AddText(bit.c_str());
    bit = "corr. = " + makeString(cfit->Correlation(0,1), 3);
    ctext->AddText(bit.c_str());
    ctext->Draw();

    cpad->cd();
    cgraph->Draw("ap");
    ctext->Draw("same");
    //    deff->Draw("same");

    outrfile->cd();
    cpad->Write();

    //And made a "calibration curve" for the uncertainties
    string ulabel = "calibcurve_uncert";
    TCanvas* upad = new TCanvas(ulabel.c_str(),ulabel.c_str(),14,33,600,600);
    upad->SetTitle(ulabel.c_str());
    upad->SetFillColor(0);
    upad->SetLeftMargin(0.15);
    upad->SetBottomMargin(0.14);

    TGraphErrors* ugraph = new TGraphErrors(cSize,xvals,uvals,xerrs,uerrs);
    ugraph->SetTitle("");
    ugraph->SetMarkerColor(kBlack);
    ugraph->SetMarkerStyle(20);
    ugraph->GetXaxis()->SetTitle(("m_{t}^{gen}-"+makeString(OFFSET,1)+" [GeV]").c_str());
    ugraph->GetXaxis()->SetTitleSize(0.05);
    ugraph->GetXaxis()->SetTitleOffset(1.0);
    ugraph->GetYaxis()->SetTitle("Pull Width");
    ugraph->GetYaxis()->SetTitleSize(0.05);
    ugraph->GetYaxis()->SetTitleOffset(1.2);
    ugraph->Draw("ap");

    TF1* ufunc = new TF1("ufunc","[0]",-15.0,15.0);
    ufunc->SetParameter(0,0.0);
    ufunc->SetParameter(1,1.0);
    ufunc->SetLineColor(kBlue);
    ufunc->SetLineWidth(2);

    TFitResultPtr ufit = ugraph->Fit(ufunc,"RMS");
    TPaveText* utext = new TPaveText(0.2,0.65,0.4,0.8,"NDC");
    utext->SetFillColor(0);
    utext->SetLineColor(kBlue);
    utext->SetLineWidth(1);
    utext->SetTextSize(0.03);
    double uval = ufit->Value(0);
    double uerr = ufit->ParError(0);
    bit = "a = "+makeString(uval,3) + "#pm"+makeString(uerr,3);
    utext->AddText(bit.c_str());
    utext->Draw();

    // TLine* uline = new TLine(-1.8,1.0,1.8,1.0);
    // uline->SetLineColor(kBlack);
    // uline->SetLineWidth(2);

    upad->cd();
    ugraph->Draw("ap");
    utext->Draw("same");
    //    uline->Draw("same");

    outrfile->cd();
    upad->Write();

    //If data point has also been supplied,
    //Calculate the calibrated mass from the fit and make a pretty plot
    double calMass = ((rawMass) - (1-aval)*OFFSET - bval)/(aval);
    ///////////////////////////////////////////////////////////////////////
    //If requested, calibrate the uncertainty from the pull widths, too.
    double calUncert = rawUncert;
    if (calcUncert)
    {
      calUncert = rawUncert*uval;
      outfile << "\\hline" << endl;
      outfile << "\\multicolumn{3}{|c|}{Calibrated Result: " << calMass << " $\\pm$ " << calUncert << "}\\\\" << endl;
    }
    else
    {
      outfile << "\\hline" << endl;
      outfile << "\\multicolumn{3}{|c|}{Calibrated Result: " << calMass << " $\\pm$ " << rawUncert << "}\\\\" << endl;
    }

    clabel = "calibcurve_data";
    TCanvas* cpad2 = new TCanvas(clabel.c_str(),clabel.c_str(),4,33,600,600);
    cpad2->SetTitle(clabel.c_str());
    cpad2->SetFillColor(0);
    cpad2->SetLeftMargin(0.15);
    cpad2->SetBottomMargin(0.14);
    cgraph->Draw("ap");
    func->Draw("same");
    ctext->Draw("same");
    
    TPaveText* ctext2 = new TPaveText(0.6,0.3,0.8,0.4,"NDC");
    ctext2->SetFillColor(0);
    ctext2->SetTextColor(kGreen+3);
    ctext2->SetLineWidth(1);
    ctext2->SetTextSize(0.03);
    string text2 = "m_{t}^{cal} = " + makeString(calMass,2) + "#pm" + makeString(calUncert,4);
    ctext2->AddText(text2.c_str());
    ctext2->Draw("same");

    TPaveText* ctext3 = new TPaveText(0.6,0.4,0.8,0.5,"NDC");
    ctext3->SetFillColor(0);
    ctext3->SetTextColor(kRed+1);
    ctext3->SetLineWidth(1);
    ctext3->SetTextSize(0.03);
    string text3 = "m_{t}^{MEM} = " + makeString(rawMass,2) + "#pm" + makeString(rawUncert,4);
    ctext3->AddText(text3.c_str());
    ctext3->Draw("same");
    
    TLine* rawline = new TLine(cgraph->GetXaxis()->GetBinLowEdge(1),rawMass-OFFSET,calMass-OFFSET,rawMass-OFFSET);
    rawline->SetLineColor(kRed+1);
    rawline->SetLineWidth(2);
    rawline->Draw("same");
    TLine* calline = new TLine(calMass-OFFSET,cgraph->GetYaxis()->GetBinLowEdge(1),calMass-OFFSET,rawMass-OFFSET);
    calline->SetLineColor(kGreen+3);
    calline->SetLineWidth(2);
    calline->Draw("same");
    
    outrfile->cd();
    cpad2->Write();
    
    //Fluctuate the fit parameters up and down by one uncertainty and calculate
    //the systematic uncertainty
    double aup = aval+aerr;
    double bup = bval+berr;
    double sys_calup = calMass-(rawMass-(1-aup)*OFFSET-bup)/(aup);
    double adown = aval-aerr;
    double bdown = bval-berr;
    double sys_caldown = calMass-(rawMass-(1-adown)*OFFSET-bdown)/(adown);
    
    // outfile << "\\hline" << endl;
    // outfile << "\\multicolumn{2}{|c|}{Systematic Uncertainties}\\\\" << endl;
    // outfile << "Calibration & $+$" << sys_calup <<" $-$" << sys_caldown << "\\\\" << endl;
    
    upSysts.insert(pair<string,double>("CALIBRATION_UP",sys_calup));
    downSysts.insert(pair<string,double>("CALIBRATION_DOWN",sys_caldown));
    double tSysUp = sys_calup;
    double tSysDown = sys_caldown;
    double tSyst = max(fabs(sys_calup),fabs(sys_caldown));
    double tSyst_sym = (fabs(sys_calup)+fabs(sys_caldown))/2;
    
    ///////////////////////////////////////////////////////////////////////
    //If applicable, loop over systematics, creating curves if needed,
    //and calculate differences to nominal curve.
    map<string,vector<TFile*> >::iterator iter;
    for (iter = systfiles.begin(); iter != systfiles.end(); ++iter)
    {
      string syst_name = iter->first;
      vector<TFile*> sfiles = iter->second;
      
      int sSize = (int)sfiles.size();
      double xvals [sSize];
      double xerrs [sSize];
      double yvals [sSize];
      double yerrs [sSize];
      for (int i = 0; i != sSize; ++i)
      {
	//First get the true mass value
	TTree* tree = (TTree*)sfiles[i]->Get("psexpInfoTree");
	double trueMass;
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("trueMass",1);
	tree->SetBranchAddress("trueMass",&trueMass);
	tree->GetEntry(0);
	
	//Now get the measured mass with resampling uncertainty
	TCanvas* sresult = (TCanvas*)sfiles[i]->Get(("Mass_"+makeString(trueMass,1)).c_str());
	TPaveStats* stext = (TPaveStats*)sresult->GetPrimitive("massMeas");
	TText* ssnippet = stext->GetLine(0);
	string sinfo = ssnippet->GetTitle();
	vector<string> spieces = splitString(sinfo);
	
	xvals[i] = trueMass - OFFSET;
	xerrs[i] = 0.0;
	yvals[i] = makeDouble(spieces[0]) - OFFSET;
	yerrs[i] = makeDouble(spieces[2]);
      }
      double syst = 0.0;
      if (sSize == 1)
      {
	//Calculate systematic from point
	syst = yvals[0]+OFFSET - nomPts[xvals[0]+OFFSET];
      }
      else
      {
	//Fit a straight line to systematic points and take difference between it
	//and nominal curve
	TGraphErrors* sgraph = new TGraphErrors(sSize,xvals,yvals,xerrs,yerrs);
	sgraph->SetTitle("");
	sgraph->SetMarkerColor(kBlack);
	sgraph->SetMarkerStyle(20);
	sgraph->GetXaxis()->SetTitle(("m_{t}^{gen}-"+makeString(OFFSET,1)+" [GeV]").c_str());
	sgraph->GetXaxis()->SetTitleSize(0.05);
	sgraph->GetXaxis()->SetTitleOffset(1.0);
	sgraph->GetYaxis()->SetTitle(("m_{t}^{MEM}-"+makeString(OFFSET,1)+" [GeV]").c_str());
	sgraph->GetYaxis()->SetTitleSize(0.05);
	sgraph->GetYaxis()->SetTitleOffset(1.2);
	
	TF1* sfunc = new TF1("sfunc","[0]+[1]*x",160.0,184.0);
	sfunc->SetParameter(0,0.0);
	sfunc->SetParameter(1,1.0);
	sfunc->SetLineColor(kBlue);
	sfunc->SetLineWidth(2);
	
	TFitResultPtr sfit = sgraph->Fit(func,"RMS");
	double aval = sfit->Value(1);
	double bval = sfit->Value(0);
	
	syst = calMass - (rawMass - (1-aval)*OFFSET - bval)/aval;

	//And save the line for visual comparison
	string slabel = "calibcurve_"+syst_name;
	TCanvas* spad = new TCanvas(slabel.c_str(),slabel.c_str(),4,33,600,600);
	spad->SetTitle(slabel.c_str());
	spad->SetFillColor(0);
	spad->SetLeftMargin(0.15);
	spad->SetBottomMargin(0.14);
	sgraph->Draw("ap");

	TText* stop = new TText(0.15,0.92,syst_name.c_str());
	stop->SetNDC();
	stop->Draw();

	TPaveText* stext = new TPaveText(0.3,0.65,0.55,0.8,"NDC");
	stext->SetFillColor(0);
	stext->SetLineColor(kBlue);
	stext->SetLineWidth(1);
	stext->SetTextSize(0.03);
	stext->AddText(("f(m_{t}^{MEM}-"+makeString(OFFSET,1)+")=a (m_{t}^{gen}-"+makeString(OFFSET,1)+")+b").c_str());
	double aerr = sfit->ParError(1);
	string bit = "a = "+makeString(aval,3) + "#pm"+makeString(aerr,3);
	stext->AddText(bit.c_str());
	double berr = sfit->ParError(0);
	bit = "b = "+makeString(bval,3) + "#pm"+makeString(berr,3);
	stext->AddText(bit.c_str());
	bit = "corr. = " + makeString(sfit->Correlation(0,1), 3);
	stext->AddText(bit.c_str());
	stext->Draw("same");

	spad->Write();
      }
      if (syst_name.find("DOWN") != string::npos)
      {
	downSysts.insert(pair<string,double>(syst_name,syst));
      }
      else
      {
	upSysts.insert(pair<string,double>(syst_name,syst));
      }

      if (syst_name.find("UP") != string::npos)
      {
	tSysUp = sqrt(tSysUp*tSysUp + syst*syst);
      }
      else if (syst_name.find("DOWN") != string::npos)
      {
	tSysDown = sqrt(tSysDown*tSysDown + syst*syst);
      }
      else
      {
	tSysUp  = sqrt(tSysUp*tSysUp + syst*syst);
	tSysDown = sqrt(tSysDown*tSysDown + syst*syst);
      }      
    }
    
    //Now write out all the systematics
    map<string,double>::iterator it;
    double jes_pupt_up = 0.0;
    double jes_pupt_down = 0.0;
    for (it = upSysts.begin(); it != upSysts.end(); ++it)
    {
      string syst_name = it->first;
      double syst = it->second;

      if (syst_name.find("UP") != string::npos)
      {
	string syst_stub = syst_name.substr(0,syst_name.find("_UP"));
	if (syst_stub.compare("COLOR") ==  0)
	{
	  //Special case: take difference between tunes, not between tunes and nominal
	  map<string,double>::iterator it2;
	  for (it2 = downSysts.begin(); it2 != downSysts.end(); ++it2)
	  {
	    if ((it2->first).find(syst_stub) != string::npos)
	    {
	      double syst_d = it2->second;
	      double final_syst = syst - syst_d;
	      outfile<<convertLabel(syst_stub)<<" & \\multicolumn{2}{c|}{" << final_syst<<" } \\\\"<<endl;
	      
	      tSyst = sqrt(tSyst*tSyst + pow(max(fabs(syst),fabs(syst_d)),2));
	      tSyst_sym = sqrt(tSyst_sym*tSyst_sym + pow((fabs(syst)+fabs(syst_d))/2,2));
	      tSysUp  = sqrt(tSysUp*tSysUp - syst*syst + final_syst*final_syst);
	      tSysDown = sqrt(tSysDown*tSysDown - syst_d*syst_d + final_syst*final_syst);	    
	    }
	  }
	}//closes color case
	else if (syst_stub.compare("UNDER_EVT") == 0)
	{
	  //Special case: compare variations to COLOR_UP, not NOMINAL
	  //Find COLOR_UP
	  double syst_mid = 0.0;
	  map<string,double>::iterator it3;
	  for (it3 = upSysts.begin(); it3 != upSysts.end(); ++it3)
	  {
	    string temp = it3->first;
	    if (temp.compare("COLOR_UP") == 0)
	    {
	      syst_mid = it3->second;
	    }
	  }
	  //Find UNDER_EVT_DOWN
	  double syst_d = 0.0;
	  map<string,double>::iterator it2;
	  for (it2 = downSysts.begin(); it2 != downSysts.end(); ++it2)
	  {
	    if ((it2->first).find(syst_stub) != string::npos)
	    {
	      syst_d = it2->second;
	    }
	  }
	  outfile<<convertLabel(syst_stub)<<" & " << syst-syst_mid<<" & " << syst_d-syst_mid << "\\\\"<<endl;
	  tSyst = sqrt(tSyst*tSyst + pow(max(fabs(syst-syst_mid),fabs(syst_d-syst_mid)),2));
	  tSyst_sym = sqrt(tSyst_sym*tSyst_sym + pow((fabs(syst-syst_mid)+fabs(syst_d-syst_mid))/2,2));
	  tSysUp  = sqrt(tSysUp*tSysUp - syst*syst + (syst-syst_mid)*(syst-syst_mid));
	  tSysDown = sqrt(tSysDown*tSysDown - syst_d*syst_d + (syst_d-syst_mid)*(syst_d-syst_mid));

	}//closes underlying event case
	else if (syst_stub.find("PUPT") != string::npos)
	{
	  //Special case: JES pile-up pT terms
	  jes_pupt_up = sqrt(pow(jes_pupt_up,2)+pow(syst,2));
	  map<string,double>::iterator it2;
	  for (it2 = downSysts.begin(); it2 != downSysts.end(); ++it2)
	  {
	    if ((it2->first).find(syst_stub) != string::npos)
	    {
	      double syst_d = it2->second;

	      outfile << convertLabel(syst_stub) <<" & " << syst << " & " << syst_d << "\\\\" << endl;

	      jes_pupt_down = sqrt(pow(jes_pupt_down,2)+pow(syst_d,2));
	      tSyst = sqrt(tSyst*tSyst + pow(max(fabs(syst),fabs(syst_d)),2));
	      tSyst_sym = sqrt(tSyst_sym*tSyst_sym + pow((fabs(syst)+fabs(syst_d))/2,2));
	    }
	  }
	}//closes JES pile-up pT case
	else
	{
	  //	  outfile << convertLabel(syst_stub) << " & " << syst << " & ";

	  map<string,double>::iterator it2;
	  bool found = false;
	  for (it2 = downSysts.begin(); it2 != downSysts.end(); ++it2)
	  {
	    if ((it2->first).find(syst_stub) != string::npos)
	    {
	      found = true;
	      double syst_d = it2->second;
	      outfile << convertLabel(syst_stub) << " & " << syst << " & " <<syst_d << "\\\\"<<endl;
	      
	      tSyst = sqrt(tSyst*tSyst + pow(max(fabs(syst),fabs(syst_d)),2));
	      tSyst_sym = sqrt(tSyst_sym*tSyst_sym + pow((fabs(syst)+fabs(syst_d))/2,2));
	    
	    }
	  }
	  if (!found)
	  {
	    outfile << convertLabel(syst_stub) << " & \\multicolumn{2}{c|}{" << syst << " } \\\\"<<endl;
	    tSyst = sqrt(tSyst*tSyst + syst*syst);
	    tSyst_sym = sqrt(tSyst_sym*tSyst_sym + syst*syst);
	  }
	}
	//	outfile << "\\\\" << endl;
      }
      else
      {
	outfile << convertLabel(syst_name) << " & \\multicolumn{2}{c|}{" << syst << " }\\\\" << endl;
	tSyst = sqrt(tSyst*tSyst + syst*syst);
	tSyst_sym = sqrt(tSyst_sym*tSyst_sym + syst*syst);
      }
    }

    //Now print out the JES groups
    outfile << "Pile-Up $p_{T}$~Group & "<<jes_pupt_up<<" & " << jes_pupt_down<<"\\\\"<<endl;

    outfile << "\\hline" << endl;
    outfile <<"\\multicolumn{3}{|c|}{Assymetric: $m_{t}^{meas} =$"<< calMass << "$\\pm$" << calUncert<< "(stat.)$^{+"<< tSysUp << "}_{-"<<tSysDown << "}$ }\\\\" << endl;
    outfile<<"\\multicolumn{3}{|c|}{Conservative: $m_{t}^{meas} =$"<<calMass<<"$\\pm$"<<calUncert<<"(stat.)$\\pm$"<<tSyst<<"(syst.)}\\\\"<<endl;
    outfile<<"\\multicolumn{3}{|c|}{Symmetrized: $m_{t}^{meas} =$"<<calMass<<"$\\pm$"<<calUncert<<"(stat.)$\\pm$"<<tSyst_sym<<"(syst.)}\\\\"<<endl;
  }//closes if-statement over calibration

  outfile << "\\hline" << endl;
  outfile << "\\end{tabular}" << endl;
  outfile.close();
  outrfile->Write();
  outrfile->Close();
  return 0;
}

int main(int argc, char* argv[])
{
  if (argc == 2)
  {
    int result = calcResults((string)argv[1]);
    return result;
  }
  else
  {
    cout << "Please provide a file with calibration inputs." << endl;
    return -1;
  }
}
