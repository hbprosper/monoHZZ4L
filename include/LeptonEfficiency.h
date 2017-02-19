#ifndef LEPTONEFFICIENCY_H
#define LEPTONEFFICIENCY_H
//-----------------------------------------------------------------------------
// File:        LeptonEfficiency.h
// Description: Return lepton efficiency pT and eta of lepton.
//              Efficiencies are read from a 2-D histogram of pt vs. eta,
//              that is, pt on the y-axis and eta on the x-axis
// Usage:
//   string rooFileName("mytable.root");
//   leff = LeptonEfficiency(rootFilename)
//           :  :
//   eff  = leff(pt, eta)
//
// Created:     18-Feb-2017 Harrison B. Prosper & Sezen Sekmen
//-----------------------------------------------------------------------------
#include <string>
#include "TFile.h"
#include "TH2F.h"
//-----------------------------------------------------------------------------
class LeptonEfficiency
{
 public:
  LeptonEfficiency();

  LeptonEfficiency(std::string rootFileName, std::string histName,
		   bool useAbsEta_=false);
  ~LeptonEfficiency();
  double operator()(double pt, double eta);

  int    ptBins()  { return ptbins; }    
  double ptMin()   { return ptmin; }
  double ptMax()   { return ptmax; }

  int    etaBins() { return etabins; }    
  double etaMin()  { return etamin; }
  double etaMax()  { return etamax; }

  TH2F*  table()   { return table_; }

 private:
  TFile* rfile;
  TH2F*  table_;

  bool   useAbsEta;
  int    ptbins;
  double ptmin;
  double ptmax;
  double ptstep;

  int    etabins;
  double etamin;
  double etamax;
  double etastep;
};

#endif

