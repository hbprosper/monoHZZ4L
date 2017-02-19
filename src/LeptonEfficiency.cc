//-----------------------------------------------------------------------------
// File:        JECUncertainty.cc
// Description: see header
// Created:     18-Feb-2017 Harrison B. Prosper & Sezen Sekmen     
//-----------------------------------------------------------------------------
#include <iostream>
#include <cassert>
#include "LeptonEfficiency.h"
//-----------------------------------------------------------------------------
using namespace std;
//-----------------------------------------------------------------------------
LeptonEfficiency::LeptonEfficiency()
  : rfile(0), table_(0)
	  {}

LeptonEfficiency::LeptonEfficiency(string rootFileName, string histName,
				   bool useAbsEta_)
  : rfile(new TFile(rootFileName.c_str())),
    table_(0),
    useAbsEta(useAbsEta_)
{
  if ( ! rfile->IsOpen() )
    {
      cout << "LeptonEfficiency ** file " << rootFileName 
	   << " not opened " << endl
	   << "check file name"
	   << endl;
      exit(0);
    }
  table_ = dynamic_cast<TH2F*>(rfile->Get(histName.c_str()));
  assert(table_);

  TAxis* xaxis = table_->GetXaxis();
  ptbins = xaxis->GetNbins();
  ptmin  = xaxis->GetBinLowEdge(1);
  ptmax  = xaxis->GetBinUpEdge(ptbins);
  ptstep = (ptmax-ptmin)/ptbins;

  TAxis* yaxis = table_->GetYaxis();
  etabins= yaxis->GetNbins();
  etamin = yaxis->GetBinLowEdge(1);
  etamax = yaxis->GetBinUpEdge(etabins);
  etastep= (etamax-etamin)/etabins;
}
    
LeptonEfficiency::~LeptonEfficiency()
{ 
  if ( rfile ) rfile->Close();
}

double LeptonEfficiency::operator()(double pt, double eta)
{
  if ( ! table_ ) return 1;
  
  if ( useAbsEta ) eta = abs(eta);

  // eta bins should not be a problem
  int ix = table_->GetXaxis()->FindBin(eta);
  if ( ix < 1 ) ix = 1;
  if ( ix > etabins ) ix = etabins;
  
  int iy = table_->GetYaxis()->FindBin(pt);
  if ( iy > ptbins ) iy = ptbins;

  //cout << "\teff bins: pt = " << pt << ", eta = " << eta <<  ", ix = " << ix
  //    << "\tiy = " << iy << endl;
  
  // for pTs lower than lower pT bound linearly extrapolate efficiency
  if ( iy > 0 )
    {
      return table_->GetBinContent(ix, iy);
    }
  else
    {
      double e1 = table_->GetBinContent(ix, 1);      
      // double e2 = table_->GetBinContent(ix, 2);
      // double p1 = ptmin;
      // double p2 = ptmin+ptstep;
      // double m  = (e2-e1)/(p2-p1);
      // double c  = (e1-m*p1);
      // //double eff= m*pt + c;
      // //return eff > 0 ? eff : 0;
      return e1*0.95;
    }
}
