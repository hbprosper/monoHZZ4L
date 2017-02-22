// ---------------------------------------------------------------------------
// mono H->ZZ->4L analysis
// Created: Les Houches 2015 Nic & HBP
// Updated: 31-Oct-2015 HBP - clean up, add ntuple output
//          15-Oct-2016 HBP - clean up further
//          29-Oct-2016 HBP - clean up more
//          17-Feb-2016 HBP - simplify
//
// Official cut flow (4e)
// Sample       qqbar -> ZZ  gg -> ZZ  Z + X     Higgs
// ---------------------------------------------------
// Initial      3970         51.5      1.79e+14  52.57
// lepton cuts   116          7.02     3740       4.53
// Z1             93          5.94     2760       4.28
// Z2             21          2.91        1.4     1.45
// ---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <map>

#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TLatex.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "LHParticle.h"
#include "HardScatter.h"

#include "monoHZZ4L.h"
#include "nic.h"
#include "Shrub.h"
#include "LeptonEfficiency.h"

using namespace std;

// ---------------------------------------------------------------------------
namespace {
  const double ZMASS=91.19;
  const int k2E2MU = 1;
  const int k4E    = 2;
  const int k4MU   = 3;

  const int ELECTRON = 11;
  const int MUON     = 13;
  const int TAU      = 15;
  const int ZBOSON   = 23;
  const int HBOSON   = 25;
  int  DEBUG = 0;
};

void copyLeptons(TClonesArray* electrons,
		 TClonesArray* muons,
		 std::vector<LHParticle>& lepton)
{
  if ( electrons )
    {
      if ( DEBUG > 0 )
	cout << "==> number of electrons " << electrons->GetEntriesFast()
	     << endl;
      
      for(int i=0; i < electrons->GetEntriesFast(); i++)
	{
	  Electron* p = static_cast<Electron*>(electrons->At(i));
	  int PID = -ELECTRON * p->Charge;
	  LHParticle q(PID, p->PT, p->Eta, p->Phi);
	  lepton.push_back(q);
	  
	  if ( DEBUG > 2 )
	    cout << i << "\t" << lepton.back() << endl;	  
	}
    }
  
 if ( muons )
   {
      if ( DEBUG > 0 )
	cout << "==> number of muons " << muons->GetEntriesFast()
	     << endl;
          
     for(int i=0; i < muons->GetEntriesFast(); i++)
       {
	 Muon* p = static_cast<Muon*>(muons->At(i));
	 int PID = -MUON * p->Charge;
	 LHParticle q(PID, p->PT, p->Eta, p->Phi);
	 lepton.push_back(q);
	 
	 if ( DEBUG > 2 )
	   cout << i << "\t" << lepton.back() << endl;

       }
   }
 
  // sort in descending order of pT
  sort(lepton.begin(), lepton.end());
}


double copyGenLeptons(std::vector<LHParticle*>& genlepton,
		      std::vector<LHParticle>& lepton,
		      LeptonEfficiency& muonEff,
		      LeptonEfficiency& elecEff)
{
  double weight = 1.0;

  for(size_t c = 0; c < genlepton.size(); c++)
    {
      int ID = abs(genlepton[c]->PID);
      if( (ID == ELECTRON) || (ID == MUON) )
	{
	  lepton.push_back(*genlepton[c]);
	  double w=1;
	  if ( ID == ELECTRON )
	    w = elecEff(lepton.back().Pt(), lepton.back().Eta());
	  else if ( ID == MUON )
	    w = muonEff(lepton.back().Pt(), lepton.back().Eta());
	  weight *= w;
	}
    }
  
  // sort in descending order of pT
  sort(lepton.begin(), lepton.end());
  if ( DEBUG > 0 )
    {
      cout << "LEPTONS" << endl;
      for(size_t c=0; c < lepton.size(); c++)
	cout << lepton[c] << endl;
    }
  return weight;
}

void purgeParticles(std::vector<LHParticle>& particles)
{
  // remove particles flagged as bad
  int c = 0;
  for(size_t i = 0; i < particles.size(); i++)
    {
      if ( particles[i].Skip ) continue;
      particles[c] = particles[i];
      c++;
    }
  particles.resize(c);
}

void filterLeptons(std::vector<LHParticle>& lepton)
{
  // ----------------------------------------------------------
  // Exclude leptons that lie outside the (pT, eta) phase space.
  // ----------------------------------------------------------
  
  for(size_t i = 0; i < lepton.size(); i++)
    {
      LHParticle& p = lepton[i]; // get a reference not a copy
      p.Skip = true;
      
      // cuts differ for electrons and muons
      double PtCut;
      double EtaCut;
      if ( abs(p.PID) == MUON )
	{
	  PtCut  = 5.0;
	  EtaCut = 2.4;
	}
      else
	{
	  PtCut  = 7.0;
	  EtaCut = 2.5;
	}
      
      if ( !(p.Pt() > PtCut) ) continue;
      if ( !(abs(p.Eta()) < EtaCut) ) continue;
      
      p.Skip = false;
    }

  // purge out-of-acceptance leptons from list of leptons
  purgeParticles(lepton);
}


void
listParticles(TClonesArray* particles)
{
  char record[256];
  sprintf(record,"%4s %8s %-10s %8s %8s %8s %4s %4s %4s %4s %4s",
	  "",   "PID", "name",
	  "PT", "Eta", "Phi",
	  "M1", "M2",
	  "D1", "D2", "Stat");
  cout << record << endl;
  for(int i = 0; i < particles->GetEntriesFast(); i++)
    {
      GenParticle* p = static_cast<GenParticle*>(particles->At(i));
      sprintf(record,"%4d %8d %-10s %8.2f %8.3f %8.3f %4d %4d %4d %4d %4d",
	      i, p->PID, LHParticle::name(p->PID).c_str(),
	      p->PT, p->Eta, p->Phi,
	      p->M1, p->M2,
	      p->D1, p->D2, p->Status);
      cout << record << endl;
    }
}


LHParticle* 
getBosons(vector<LHParticle>&  particles,
	  vector<LHParticle*>& genZ,
	  vector<LHParticle*>& genL)
{
  // get Higgs and Z bosons
  LHParticle* H=0;
  vector<LHParticle*> tmpZ;
  vector<pair<float, int> > mass;
  for(size_t c = 0; c < particles.size(); c++)
    {
      if      ( particles[c].PID == HBOSON )
	{
	  H = &particles[c];
	}
      else if ( particles[c].PID == ZBOSON )
	{
	  tmpZ.push_back( &particles[c] );
	  mass.push_back(pair<float, int>(1.0/tmpZ.back()->M(), tmpZ.size()-1));
	}
    }
  // order in descending order of mass so that Z1 is first, then Z2
  sort(mass.begin(), mass.end());
  
  // get leptons from Z bosons
  for(size_t c=0; c < tmpZ.size(); c++)
    {
      pair<float, int> m = mass[c];
      LHParticle* Z = tmpZ[m.second];
      if ( DEBUG > 0 )
	cout << *Z << endl;
      
      int d1 = Z->Daughters[0];
      int d2 = Z->Daughters[1];
      int ID = abs(particles[d1].PID);
      if( (ID == ELECTRON) || (ID == MUON) || (ID == TAU))
	{
	  genZ.push_back(Z);
	  if ( particles[d1].Pt() > particles[d2].Pt() )
	    {
	      genL.push_back(&particles[d1]);
	      genL.push_back(&particles[d2]);
	    }
	  else
	    {
	      genL.push_back(&particles[d2]);
	      genL.push_back(&particles[d1]);
	    }
	  }
      }
  assert(genZ.size()>1);
  assert(genL.size()>3);
  return H;
}


// match objects
void 
matchObjects(vector<LHParticle>&  r,
	     vector<LHParticle*>& p,
	     nic::Match& match,
	     float dRcut=0.05)
{
  for(size_t i = 0; i < r.size(); i++) r[i].ID  = -1;
  for(size_t i = 0; i < p.size(); i++) p[i]->ID = -1;

  for(size_t i = 0; i < r.size(); i++)
    for(size_t j = 0; j < p.size(); j++)
      match.add(i, r[i].Eta(),  r[i].Phi(),
		j, p[j]->Eta(), p[j]->Phi());
  match.run();

  for(size_t i = 0; i < match.order.size(); i++)
    {
      if ( !(match.order[i].first < dRcut) ) continue;
      int ii = match.order[i].second.first;
      int jj = match.order[i].second.second;
      r[ii].ID  = jj;
      p[jj]->ID = ii;
    }
}

void 
findDileptons(vector<LHParticle>& lepton, vector<LHParticle>& dilepton)
{
  dilepton.clear();
  
  // create opposite-sign, same-flavor dileptons
  for(size_t k = 0; k < lepton.size(); k++)
    {
      LHParticle& pk = lepton[k]; // get a reference (i.e., alias), not a copy
      
      for(size_t j = k+1; j < lepton.size(); j++)
	{
	  LHParticle& pj = lepton[j]; // get a reference, not a copy

	  // require same flavor opposite sign
	  if ( !( (pk.PID + pj.PID)==0) ) continue;

	  // require 12 < m(ll) < 120 GeV
	  LHParticle ll  = pk + pj;
	  if ( !(ll.M() > 12) )  continue;
	  if ( !(ll.M() < 120) ) continue;
	  
	  // we have same flavor opposite sign dileptons
	  // set constituents of dilepton
	  ll.PID  = abs(pk.PID);
	  ll.Skip = false;
	  
	  // set leptons so that higher pT
	  // lepton comes first
	  if (pk.Pt() > pj.Pt())
	    {
	      ll.Daughters.push_back(k);
	      ll.Daughters.push_back(j);
	    }
	  else
	    {
	      ll.Daughters.push_back(j);
	      ll.Daughters.push_back(k);
	    }	    
	  dilepton.push_back(ll);
	}
    }
  return;
}

void
purgeDileptons(vector<LHParticle>& dilepton, LHParticle& Z)
{
  if ( dilepton.size() < 1 ) return;

  // purge dileptons that share a lepton with Z candidate
  int j1 = Z.Daughters[0]; // higher pT lepton
  int j2 = Z.Daughters[1]; // lower pT lepton
  for(size_t c=0; c < dilepton.size(); c++)
    {
      int k1 = dilepton[c].Daughters[0];
      int k2 = dilepton[c].Daughters[1];
      if      ( k1 == j1 )
	dilepton[c].Skip = true;
      else if ( k1 == j2 )
	dilepton[c].Skip = true;
      else if ( k2 == j1 )
	dilepton[c].Skip = true;
      else if ( k2 == j2 )
	dilepton[c].Skip = true;      
    }
  purgeParticles(dilepton);
}


int
getFinalState(LHParticle& Z1, LHParticle& Z2)
{
  int Z1PID = abs(Z1.PID);
  int Z2PID = abs(Z2.PID);
  if ( Z1PID != Z2PID )
    return k2E2MU;
  else if ( Z1PID == ELECTRON )
    return k4E;
  else
    return k4MU;
}

int
findZ1(vector<LHParticle>& dilepton)
{
  // Find dilepton of correct flavor closest to the Z pole mass
  // and with mass > 40 GeV
  // Note: findDileptons orders the dileptons in descending pT
  
  int which =-1;
    
  // need at least one dilepton
  if ( dilepton.size() < 1 ) return which;

  double smallest=1e6;
  for(size_t c=0; c < dilepton.size(); c++)
    {
      // find dilepton with the smallest |mass - Zmass|
      double dmass = abs(dilepton[c].M()-ZMASS);
      if ( dmass < smallest )
	{
	  which = c;
	  smallest = dmass;
	}
    }
  if ( which < 0 ) return which;

  if ( dilepton[which].M() < 40 ) return -1;
  
  dilepton[which].Name = "Z1";
  return which;  
}


int 
findZ2(vector<LHParticle>& dilepton)
{
  if ( dilepton.size() < 1 ) return -1;
  dilepton[0].Name = "Z2";
  return 0;
}


bool ghostFree(vector<LHParticle>& lepton)
{
  for(size_t i = 0; i < lepton.size(); i++)
    {
      LHParticle& pi = lepton[i];
      for(size_t j = i+1; j < lepton.size(); j++)
	{
	  LHParticle& pj = lepton[j];

	  double dR = nic::deltaR(pi.Eta(), pi.Phi(),
				  pj.Eta(), pj.Phi());
	  // ghost removal
	  if ( dR < 0.02 ) return false;
	}
    }
  return true;
}
// ---------------------------------------------------------------------------
// inputFile    input file (with .root extension) or a filelist
// numberEvents obvious, no?!
// luminosity   integrated luminosity
// xsection     cross section duh!
// ---------------------------------------------------------------------------
void monoHZZ4L::analysis(string inputFile,
			 int    numberEvents,
			 double luminosity,
			 double xsection,
			 bool useRECO)
{
  cout << endl << "\t== monoHZZ4L::analysis ==" << endl;
  if ( useRECO )
    cout << "\t== use RECO level leptons ==" << endl;
  else
    cout << "\t== use GEN level leptons ==" << endl;

  LeptonEfficiency* muonEff=0;
  LeptonEfficiency* elecEff=0;
  
  if ( ! useRECO )
    {
      muonEff = new LeptonEfficiency("EfficienciesAndSF_BCDEF.root",
				     "TightISO_TightID_pt_eta/"
				     "efficienciesMC/"
				     "abseta_pt_MC",
				     true); // use |eta|
				       
      elecEff = new LeptonEfficiency("egammaEffi.txt_EGM2D.root",
				     "EGamma_EffMC2D");
    }
  
  // -----------------------------------------
  // Objects
  //  muons  - PT > 5 GeV
  //           |Eta| < 2.4
  //
  //  muons  - PT > 7 GeV
  //           |Eta| < 2.7
  //           relisol < 0.4
  //
  //  Z1     - min({mZ - Zmass})
  //           mZ1 > 40 GeV
  //
  // -----------------------------------------
  // try to load libDelphes
  try
    {
      gSystem->Load("libDelphes");
    }
  catch (...)
    {
      nic::ciao("can't load libDelphes");
    }

  string namen = nic::nameonly(inputFile);
  namen = nic::replace(namen, "pp13TeV_", "ntuple_");
  
  // make sure plots and histos directories exist
  nic::shell("mkdir -p plots; mkdir -p histos");
  
  // create empty root file for analysis results
  char filename[256];
  sprintf(filename, "histos/histos_%s.root", namen.c_str());
  TFile* theFile = new TFile(filename, "RECREATE");
  if ( ! theFile->IsOpen() ) nic::ciao(string("can't create ") + filename);

  // create a chain of input files
  cout << endl << "\t=> input files:" << endl;
  TChain* chain = new TChain("Delphes");
  
  // inputFile could be either a root file or a filelist
  std::vector<std::string> inputFiles;
  if ( inputFile.find(".root") != std::string::npos )
    inputFiles.push_back(inputFile);
  else
    {
      // assume this is a filelist
      ifstream inp(inputFile.c_str());
      if ( !inp.good() ) nic::ciao("can't open " + inputFile);
      string infilename;
      while (getline(inp, infilename))
	{
	  if ( infilename.substr(0,1) == "#" ) continue;
	  if ( infilename.substr(0,1) == " " ) continue;
	  if ( infilename.substr(0,1) == "\n" ) continue;
	  inputFiles.push_back(infilename);
	}
    }
  for (size_t i=0; i < inputFiles.size(); i++){
    chain->Add(inputFiles[i].c_str());
    std::cout << "\t=> added input file " << inputFiles[i] << std::endl;
  }

  // -----------------------------------------
  // create A Delphes tree reader
  // -----------------------------------------
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);


  // -----------------------------------------
  // initialize branch pointers
  // -----------------------------------------
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet      = treeReader->UseBranch("Jet");     
  TClonesArray *branchMET      = treeReader->UseBranch("MissingET");  
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon     = treeReader->UseBranch("Muon");

  // -----------------------------------------
  // create empty histograms
  // -----------------------------------------
  nic::setStyle();
  
  vector<TH1F*> h;
  
  cout << endl << "\t=> initialize histograms" << endl;
  TH1F* h_nEvent = new TH1F("h_nEvent", "Cut flow", 9, 0, 9);
  h.push_back(h_nEvent);
  h_nEvent->Sumw2();
  h_nEvent->GetXaxis()->SetBinLabel(1,"no cuts");
  h_nEvent->GetXaxis()->SetBinLabel(2,"lepton >= 4");
  h_nEvent->GetXaxis()->SetBinLabel(3,"dilepton >= 1");
  h_nEvent->GetXaxis()->SetBinLabel(4,"Z1");
  h_nEvent->GetXaxis()->SetBinLabel(5,"dilepton >= 2");
  h_nEvent->GetXaxis()->SetBinLabel(6,"Z2");
  h_nEvent->GetXaxis()->SetBinLabel(7,"pTl1 > 20");
  h_nEvent->GetXaxis()->SetBinLabel(8,"pTl2 > 10");
  h_nEvent->GetXaxis()->SetBinLabel(9,"m4l > 70");
  
  TH1F* h_nleptons = new TH1F("nleptons", "", 10, 0, 10);
  h.push_back(h_nleptons);
  h_nleptons->GetXaxis()->SetTitle("#font[12]{n}_{leptons}");
  
  TH1F* h_njets = new TH1F("njets", "", 10, 0., 10.);
  h.push_back(h_njets);
  h_njets->GetXaxis()->SetTitle("#font[12]{n}_{jets}");

  int maxlep=4;
  TH1F* h_Eta[maxlep];
  TH1F* h_PT[maxlep];
  TH1F* h_genEta[maxlep];
  TH1F* h_genPT[maxlep];
  for(int c=0; c < maxlep; c++)
    {
      char name[20];
      
      sprintf(name, "Eta%d", c+1);
      h_Eta[c] = new TH1F(name, "", 100, -5., 5.);
      h.push_back(h_Eta[c]);
      h_Eta[c]->GetXaxis()->SetTitle("#font[12]{#eta}_{reco}");

      sprintf(name, "PT%d", c+1);
      h_PT[c] = new TH1F(name, "", 100, 0., 200.);
      h.push_back(h_PT[c]);
      h_PT[c]->GetXaxis()->SetTitle("#font[12]{p}_{T,reco} (GeV)");

      sprintf(name, "genEta%d", c+1);
      h_genEta[c] = new TH1F(name, "", 100, -5., 5.);
      h.push_back(h_genEta[c]);
      h_genEta[c]->GetXaxis()->SetTitle("#font[12]{#eta}_{gen}");

      sprintf(name, "genPT%d", c+1);
      h_genPT[c]  = new TH1F(name,  "", 100, 0., 200.);
      h.push_back(h_genPT[c]);
      h_genPT[c]->GetXaxis()->SetTitle("#font[12]{p}_{T,gen} (GeV)");
    }
  TH1F* h_Z1mass = new TH1F("Z1mass", "", 200,  0,  200.);
  h.push_back(h_Z1mass);
  h_Z1mass->GetXaxis()->SetTitle("#font[12]{m}_{Z1,reco} (GeV)");
  
  TH1F* h_Z2mass = new TH1F("Z2mass", "", 200,  0., 200.);
  h.push_back(h_Z2mass);
  h_Z2mass->GetXaxis()->SetTitle("#font[12]{m}_{Z2,reco} (GeV)");
  
  TH1F* h_genZ1mass = new TH1F("genZ1mass", "", 200, 0., 200.);
  h.push_back(h_genZ1mass);
  h_genZ1mass->GetXaxis()->SetTitle("#font[12]{m}_{Z1,gen} (GeV)");
    
  TH1F* h_genZ2mass = new TH1F("genZ2mass", "", 200, 0., 200.);
  h.push_back(h_Z2mass);
  h_genZ2mass->GetXaxis()->SetTitle("#font[12]{m}_{Z2,gen} (GeV)");
  
  TH1F* h_Hmass = new TH1F("Hmass", "", 200, 0., 200.);
  h.push_back(h_Hmass);
  h_Hmass->GetXaxis()->SetTitle("#font[12]{m}_{H,reco} (GeV)");
  
  TH1F* h_HmassLarge = new TH1F("HmassLarge", "", 200, 0., 400.);
  h.push_back(h_HmassLarge);
  h_HmassLarge->GetXaxis()->SetTitle("#font[12]{m}_{H,reco} (GeV)");
  
  TH1F* h_genHmass = new TH1F("genHmass", "", 200, 0., 200.);
  h.push_back(h_genHmass);
  h_genHmass->GetXaxis()->SetTitle("#font[12]{m}_{H,gen} (GeV)");
  
  TH1F* h_MET  = new TH1F("MET", "", 200, 0., 200.);
  h.push_back(h_MET);
  h_MET->GetXaxis()->SetTitle("missing #font[12]{E}_{T} (GeV)");
  
  TH1F* h_dZ1 = new TH1F("dZ1", "", 80,  0., 20.);
  h.push_back(h_dZ1);
  h_dZ1->GetXaxis()->SetTitle("#font[12]{m}_{Z1(gen)-Z1(reco)} (GeV)");
  
  TH1F* h_dZ2 = new TH1F("dZ2", "",  80,  0., 20.);
  h.push_back(h_dZ2);
  h_dZ2->GetXaxis()->SetTitle("#font[12]{m}_{Z2(gen)-Z2(reco)} (GeV)");

  TH1F* h_dRlepton = new TH1F("dRlepton", "", 100,  0., 0.5);
  h.push_back(h_dRlepton);
  h_dRlepton->GetXaxis()->SetTitle("#Delta#font[12]{R}"
				   "(#font[12]{l}, #font[12]{l_{gen}}");  

  TH1F* h_maxmatch = new TH1F("maxmatch", "", 10,  0., 10);
  h.push_back(h_maxmatch);
  h_maxmatch->GetXaxis()->SetTitle("max[match-index]");

  TH1F* h_nummatch = new TH1F("nummatch", "", 10,  0., 10);
  h.push_back(h_nummatch);
  h_nummatch->GetXaxis()->SetTitle("match multiplicity");

  TH1F* h_costhetastar = new TH1F("costhetastar", "", 100, -1, 1);
  h.push_back(h_costhetastar);
  h_costhetastar->GetXaxis()->SetTitle("cos(#Theta_{*})");

  TH1F* h_costheta1 = new TH1F("costheta1", "", 100, -1, 1);
  h.push_back(h_costheta1);
  h_costheta1->GetXaxis()->SetTitle("cos(#theta_{1})");

  TH1F* h_costheta2 = new TH1F("costheta2", "", 100, -1, 1);
  h.push_back(h_costheta2);
  h_costheta2->GetXaxis()->SetTitle("cos(#theta_{2})");

  TH1F* h_cosPhi = new TH1F("cosPhi", "", 100, -1, 1);
  h.push_back(h_cosPhi);
  h_cosPhi->GetXaxis()->SetTitle("cos(#Phi)");

  TH1F* h_cosPhi1 = new TH1F("cosPhi1", "", 100, -1, 1);
  h.push_back(h_cosPhi1);
  h_cosPhi1->GetXaxis()->SetTitle("cos(#Phi_{1})");

  for(size_t c=0; c < h.size(); c++)
    {
      h[c]->SetLineWidth(2);
      h[c]->SetLineColor(kBlack);
      h[c]->SetFillStyle(3001);
      h[c]->SetFillColor(kBlue);
    }

  // ========================================================================
  // EVENT LOOP
  // ========================================================================  

  cout << endl << "\t=> begin event loop" << endl;
    
  long int numberOfEntries = treeReader->GetEntries();

  //numberOfEntries = 100;
  
  if ( numberEvents > 0 ) numberOfEntries = numberEvents;
  std::cout << "\t=> analyzing " << numberOfEntries <<" events" << std::endl;

  // -----------------------------------------
  // cross sections @ 13 TeV
  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CrossSections#
  // Higgs_cross_sections_and_decay_b
  //
  // ggF = 43.920 pb +/- 3.5  pb (mH = 125.0)
  // VBF =  3.748 pb +/- 0.13 pb (mH = 125.0)
  // WH  =  1.380 pb
  // ZH  =  0.870 pb
  // ttH =  0.509 pb
  // bbH =  0.512 pb
  //
  // BR(H->ZZ->e+e-e+e-)   = 0.327e-4
  // BR(H->ZZ->e+e-mu+mu-) = 0.593e-4
  // BR(H->ZZ->4L)         = 1.250e-4 (L=e,mu)
  // BR(H->ZZ->4L)         = 2.760e-4 (L=e,mu,tau)
  //
  // BR(Z->l+l-)           = 3.3658e-2
  // -----------------------------------------
  double eventCount  = xsection * luminosity;
  double eventWeightOrig = eventCount / numberOfEntries;
  double eventWeight = eventWeightOrig;
  
  char preamble[512];
  sprintf(preamble,
	  "========================================\n"
	  "Luminosity:       %10.2f events/fb\n"
	  "Cross section*BR: %10.2f fb\n"
	  "Event count:      %10.2f events\n"
	  "Event weight:     %10.2e events\n"
	  "========================================\n",
	  luminosity, xsection, eventCount, eventWeight);
  cout << endl << preamble << endl;
  
  // -----------------------------------------
  // Set up output tree
  // -----------------------------------------
  string outfile = namen + string(".root");
  Shrub evtTree(outfile, "Analysis", namen);


  int passed = 0;
  double totalPassed = 0.0;
  double totalPassedUnc = 0.0;
  
  HardScatter hardscatter(branchParticle);

  // -----------------------------------------
  // START EVENT LOOP
  // -----------------------------------------
  for(Int_t entry = 0; entry < numberOfEntries; entry++)
    {
      // remember to reset event weight to its original value
      
      eventWeight = eventWeightOrig;
      
      bool printMe = entry % 10000 == 0;
      
      if ( printMe )
	cout << "\t=> processing event "
	     << entry
	     << "\t=> selected " << passed
	     << "\t=> " << totalPassed
	     << endl;
      
      if ( DEBUG > 0 )
	cout << endl
	     << "==> event: " << entry
	     << " =========================================================="
	     << endl;
      
      //Load branches for event
      treeReader->ReadEntry(entry);

      // -----------------------------------------------------
      // first bin of h_nEvent contains the total event weight
      // -----------------------------------------------------

      h_nEvent->Fill(0.1, eventWeight);
    
      // -----------------------------------------------------
      // get hard scatter event
      // the Z bosons are ordered in descending order of mass
      // and the leptons in descending order of pT
      // -----------------------------------------------------
      vector<LHParticle> gparticles;
      hardscatter.get(gparticles, 0);
      
      vector<LHParticle*> genZ;
      vector<LHParticle*> genL;
      LHParticle* genH  = getBosons(gparticles, genZ, genL);

      
      LHParticle* genZ1 = 0;
      LHParticle* genL1 = 0;
      LHParticle* genL2 = 0;
      if ( genZ.size() > 0 )
	{
	  genZ1 = genZ[0];	
	  if ( genL.size() > 0 ) genL1 = genL[0];
	  if ( genL.size() > 1 ) genL2 = genL[1];
	}

      LHParticle* genZ2 = 0;
      LHParticle* genL3 = 0;
      LHParticle* genL4 = 0;
      if ( genZ.size() > 1 )
	{
	  genZ2 = genZ[1];
	  if ( genL.size() > 2 ) genL3 = genL[2];
	  if ( genL.size() > 3 ) genL4 = genL[3];
	}
      
      // Histogram gen-level objects
      if ( genZ1 )
	h_genZ1mass->Fill(genZ1->M(), eventWeight);
      
      if ( genZ2 )
	h_genZ2mass->Fill(genZ2->M(), eventWeight);
   
      if ( genH )
	h_genHmass->Fill(genH->M(), eventWeight);

      if ( genL1 )
	{
	  h_genPT[0]->Fill(genL1->Pt(), eventWeight);
	  h_genEta[0]->Fill(genL1->Eta(), eventWeight);
	}
      
      if ( genL2 )
	{
	  h_genPT[1]->Fill(genL2->Pt(), eventWeight);
	  h_genEta[1]->Fill(genL2->Eta(), eventWeight);
	}
      
      if ( genL3 )
	{
	  h_genPT[2]->Fill(genL3->Pt(), eventWeight);
	  h_genEta[2]->Fill(genL3->Eta(), eventWeight);
	}
      
      if ( genL4 )
	{
	  h_genPT[3]->Fill(genL4->Pt(), eventWeight);
	  h_genEta[3]->Fill(genL4->Eta(), eventWeight);
	}
      
      // -----------------------------------------------------
      // copy data to standard Les Houches objects and
      // apply lepton filter
      // -----------------------------------------------------
      LHParticle::s_UID = 0; // reset unique ID (UID) 
      vector<LHParticle> lepton;
      if ( useRECO )
	{
	  copyLeptons(branchElectron, branchMuon, lepton);
	}
      else
	{
	  if ( DEBUG > 0 )
	    cout << "==> get GEN leptons " << endl;
	  
	  double effWeight = copyGenLeptons(genL, lepton, *muonEff, *elecEff);
	  if ( DEBUG > 0 )
	    cout << "==> lepton efficiency weight: " << effWeight << endl;
	  
	  eventWeight *= effWeight;
	}

      // apply lepton filter
      filterLeptons(lepton);
      
      int nleptons = (int)lepton.size();    
      h_nleptons->Fill(nleptons);
      
      if ( DEBUG > 0 )
	cout << "==> number of selected leptons: "
	     << lepton.size() << endl;

      // ----------------------------------------------------------
      // CUT 1: number of leptons >= 4
      // ----------------------------------------------------------
      if ( !(lepton.size() > 3) )  continue;
      h_nEvent->Fill(1.1, eventWeight);


      // Match gen leptons to reco leptons
      nic::Match match;
      if ( genL.size() > 0 )
	matchObjects(lepton, genL, match);
      
      // histogram dR
      for(size_t c=0; c < match.order.size(); c++)
	{
	  float dR = match.order[c].first;
	  h_dRlepton->Fill(dR);
	}
      
      // ID < 0 if no match
      int maxii  =-1;
      int nmatch = 0;
      for(size_t c=0; c < lepton.size(); c++)
	{
	  if ( lepton[c].ID < 0 ) continue;
	  maxii = c;
	  nmatch++;
	}
      h_maxmatch->Fill(maxii);
      h_nummatch->Fill(nmatch);
      
      // ----------------------------------------------------------
      // CUT 2: number of (opposite-sign, same-flavor dileptons
      //        with 12 < mass < 120 GeV) > 0
      // ----------------------------------------------------------
      vector<LHParticle> dilepton;
      findDileptons(lepton, dilepton);
      if ( !(dilepton.size() > 0) )  continue;
      h_nEvent->Fill(2.1, eventWeight);
      
      if ( DEBUG > 0 )
	{
	  cout << "==> Dileptons"
	       << endl;
	  for (size_t c=0; c < dilepton.size(); c++)
	    cout << dilepton[c] << endl;
	}
      
      // ----------------------------------------------------------
      // CUT 3: find Z1 with mass > 40 GeV
      // ----------------------------------------------------------
      // find dilepton closest to Z pole mass
      int which = findZ1(dilepton);
      if ( which < 0 ) continue;    
      h_nEvent->Fill(3.1, eventWeight);
      
      // cache Z1 and associated leptons
      LHParticle Z1 = dilepton[which];
      LHParticle L1 = lepton[Z1.Daughters[0]];
      LHParticle L2 = lepton[Z1.Daughters[1]];
      assert((L1.PID+L2.PID)==0);
      
      if ( DEBUG > 0 )
	{
	  cout << "==> Z1"
	       << endl << Z1
	       << endl;
	  cout << "\tLeptons"
	       << endl << L1
	       << endl << L2
	       << endl;
	}
       
      // ----------------------------------------------------------
      // CUT 4: find other dileptons
      // ----------------------------------------------------------
      // purge dileptons that share a lepton with Z1
      purgeDileptons(dilepton, Z1);    
      if ( !(dilepton.size() > 0) )  continue;
      h_nEvent->Fill(4.1, eventWeight);
      
      if ( DEBUG > 0 )
	{
	  cout << "==> Purged Dileptons"
	       << endl;
	  for (size_t c=0; c < dilepton.size(); c++)
	    cout << dilepton[c] << endl;
	}
      
      // ----------------------------------------------------------
      // CUT 5: find Z2
      // ----------------------------------------------------------
      which = findZ2(dilepton);
      if ( which < 0 ) continue;
      h_nEvent->Fill(5.1, eventWeight);

      LHParticle Z2 = dilepton[which];
      LHParticle L3 = lepton[Z2.Daughters[0]];
      LHParticle L4 = lepton[Z2.Daughters[1]];
      assert((L3.PID+L4.PID)==0);

      if ( DEBUG > 0 )
	{
	  cout << "==> Z2"
	       << endl << Z2
	       << endl;
	  cout << "\tLeptons"
	       << endl << L3
	       << endl << L4
	       << endl;
	}

      // ----------------------------------------------------------
      // CUT 6: pTl1 > 20 GeV
      // ----------------------------------------------------------
      vector<LHParticle> lep4;
      lep4.push_back(L1);
      lep4.push_back(L2);
      lep4.push_back(L3);
      lep4.push_back(L4);
      sort(lep4.begin(), lep4.end());
      
      if ( !(lep4[0].Pt() > 20) ) continue;
      h_nEvent->Fill(6.1, eventWeight);

      // ----------------------------------------------------------
      // CUT 7: pTl2 > 10 GeV
      // ----------------------------------------------------------
      if ( !(lep4[1].Pt() > 10) ) continue;
      h_nEvent->Fill(7.1, eventWeight);      
      
      // ----------------------------------------------------------
      // CUT 8: m4l > 70 GeV
      // ----------------------------------------------------------
      // compute 4-lepton 4-vector
      LHParticle H = Z1 + Z2;
      if ( !(H.M() > 70) ) continue;
      h_nEvent->Fill(8.1, eventWeight);      			 


            
      // number of events that pass selection criteria
      passed++;
      totalPassed += eventWeight;
      totalPassedUnc += eventWeight*eventWeight;
      
      // ----------------------------------------------------------
      // END OF EVENT SELECTION
      // ----------------------------------------------------------
      
      h_Z1mass->Fill(Z1.M(), eventWeight);
      h_PT[0]->Fill(L1.Pt(), eventWeight);
      h_PT[1]->Fill(L2.Pt(), eventWeight);
      h_PT[2]->Fill(L3.Pt(), eventWeight);

      h_PT[3]->Fill(L4.Pt(), eventWeight);	  
      h_Z2mass->Fill(Z2.M(), eventWeight);
      h_Hmass->Fill(H.M(), eventWeight);
      h_HmassLarge->Fill(H.M(), eventWeight);
      
      // Missing transverse energy
      MissingET* met = static_cast<MissingET*>(branchMET->At(0));    
      h_MET->Fill(met->MET, eventWeight);
      
      // store jet info
      float HT = 0;
      vector<LHParticle> jet;
      for(int i = 0; i < branchJet->GetEntriesFast(); i++)
	{
	  Jet* j = static_cast<Jet*>(branchJet->At(i));    
	  if ( ! (j->PT > 30) ) continue;
	  if ( ! (abs(j->Eta) < 4.7) ) continue;
	  jet.push_back(LHParticle(81,    // ID for jet
				   j->PT,
				   j->Eta,
				   j->Phi,
				   j->Mass));
	  HT += j->PT;
	}
      sort(jet.begin(), jet.end());
      
      h_njets->Fill(jet.size(), eventWeight);
      
      if ( genZ1 )
	{
	  // compute (Z1 - genZ1).M()
	  // this will be zero for a perfect match
	  // between reco Z1 and gen Z1
	  TLorentzVector dZ1 = Z1 - *genZ1;
	  h_dZ1->Fill(dZ1.M(), eventWeight);
	}
      
      if ( genZ2 )
	{
	  // compute (Z2 - genZ2).M()
	  TLorentzVector dZ2 = Z2 - *genZ2;
	  h_dZ2->Fill(dZ2.M(), eventWeight);
	}
      
      // ----------------------------------------------------------
      // write out event quantities
      // ----------------------------------------------------------
      evtTree.Clear();
      
      evtTree.weight = eventWeight;
      evtTree.fstate = getFinalState(Z1, Z2);
      evtTree.njets  = (int)jet.size();
      evtTree.nleps  = (int)lepton.size();
      
      if ( genZ1 )
	{
	  evtTree.genZ1pt   = genZ1->Pt();
	  evtTree.genZ1eta  = genZ1->Eta();
	  evtTree.genZ1phi  = genZ1->Phi();
	  evtTree.genZ1mass = genZ1->M();
	}
      if ( genZ2 )
	{
	  evtTree.genZ2pt   = genZ2->Pt();
	  evtTree.genZ2eta  = genZ2->Eta();
	  evtTree.genZ2phi  = genZ2->Phi();
	  evtTree.genZ2mass = genZ2->M();
	}
      if ( genL1 )
	{
	  evtTree.genl1match= genL1->ID;
	  evtTree.genl1PID  = genL1->PID;
	  evtTree.genl1pt   = genL1->Pt();
	  evtTree.genl1eta  = genL1->Eta();
	  evtTree.genl1phi  = genL1->Phi();
	}
      if ( genL2 )
	{
	  evtTree.genl2match= genL2->ID;
	  evtTree.genl2PID  = genL2->PID;    
	  evtTree.genl2pt   = genL2->Pt();
	  evtTree.genl2eta  = genL2->Eta();
	  evtTree.genl2phi  = genL2->Phi();
	}
      if ( genL3 )
	{
	  evtTree.genl3match= genL3->ID;
	  evtTree.genl3PID  = genL3->PID;    
	  evtTree.genl3pt   = genL3->Pt();
	  evtTree.genl3eta  = genL3->Eta();
	  evtTree.genl3phi  = genL3->Phi();
	}
      if ( genL4 )
	{
	  evtTree.genl4match= genL4->ID;
	  evtTree.genl4PID  = genL4->PID;    
	  evtTree.genl4pt   = genL4->Pt();
	  evtTree.genl4eta  = genL4->Eta();
	  evtTree.genl4phi  = genL4->Phi();    
	}
      
      evtTree.Z1pt   = Z1.Pt();
      evtTree.Z1eta  = Z1.Eta();
      evtTree.Z1phi  = Z1.Phi();
      evtTree.Z1mass = Z1.M();
      
      evtTree.Z2pt   = Z2.Pt();
      evtTree.Z2eta  = Z2.Eta();
      evtTree.Z2phi  = Z2.Phi();
      evtTree.Z2mass = Z2.M();
      
      evtTree.Hpt    = H.Pt();
      evtTree.Heta   = H.Eta();
      evtTree.Hphi   = H.Phi();
      evtTree.Hmass  = H.M();    
      
      evtTree.l1match= L1.ID;
      evtTree.l1PID  = L1.PID;
      evtTree.l1pt   = L1.Pt();
      evtTree.l1eta  = L1.Eta();
      evtTree.l1phi  = L1.Phi();
      
      evtTree.l2match= L2.ID;
      evtTree.l2PID  = L2.PID;    
      evtTree.l2pt   = L2.Pt();
      evtTree.l2eta  = L2.Eta();
      evtTree.l2phi  = L2.Phi();
      
      evtTree.l3match= L3.ID;
      evtTree.l3PID  = L3.PID;    
      evtTree.l3pt   = L3.Pt();
      evtTree.l3eta  = L3.Eta();
      evtTree.l3phi  = L3.Phi();

      evtTree.l4match= L4.ID;
      evtTree.l4PID  = L4.PID;    
      evtTree.l4pt   = L4.Pt();
      evtTree.l4eta  = L4.Eta();
      evtTree.l4phi  = L4.Phi();
  
      if ( jet.size() > 0 )
	{
	  evtTree.j1pt   = jet[0].Pt();
	  evtTree.j1eta  = jet[0].Eta();
	  evtTree.j1phi  = jet[0].Phi();
	  evtTree.j1mass = jet[0].M();
	}
      
      if ( jet.size() > 1 )
	{
	  evtTree.j2pt   = jet[1].Pt();
	  evtTree.j2eta  = jet[1].Eta();
	  evtTree.j2phi  = jet[1].Phi();
	  evtTree.j2mass = jet[1].M();
	  
	  LHParticle dijet = jet[0] + jet[1];
	  evtTree.massjj   = dijet.M();
	  evtTree.deltaetajj = abs(jet[0].Eta() - jet[1].Eta());
	}
      
      evtTree.HT  = HT;
      evtTree.met = met->MET;
      
      // add MELA variables
      float costhetastar=-2;
      float costheta1=-2;
      float costheta2=-2;
      float cosPhi=-2;
      float cosPhi1=-2;

      nic::computeMELAangles(L1, L1.PID,
			     L2, L2.PID,
			     L3, L3.PID,
			     L4, L4.PID,
			     costhetastar,
			     costheta1,
			     costheta2,
			     cosPhi,
			     cosPhi1);
	  
      evtTree.costhetastar = costhetastar;
      evtTree.costheta1    = costheta1;
      evtTree.costheta2    = costheta2;
      evtTree.cosPhi       = cosPhi;
      evtTree.cosPhi1      = cosPhi1;
      
      h_costhetastar->Fill(costhetastar, eventWeight);
      h_costheta1->Fill(costheta1, eventWeight);
      h_costheta2->Fill(costheta2, eventWeight);
      h_cosPhi->Fill(cosPhi, eventWeight);
      h_cosPhi1->Fill(cosPhi1, eventWeight);
  
      evtTree.Fill();
      
    } // END OF EVENT LOOP
  
  evtTree.Close();
  
  // -----------------------------------------------------------
  // Finish up
  // -----------------------------------------------------------
  gStyle->SetOptStat("eimr");

  TCanvas* cPT = new TCanvas("PT", "PT", 10, 10, 800, 800);
  cPT->Divide(2, 2);

  TCanvas* cgenPT = new TCanvas("genPT", "genPT", 100, 100, 800, 800);
  cgenPT->Divide(2, 2);       
  
  for(int c=0; c < 4; c++)
    {
      cPT->cd(c+1);
      h_PT[c]->Draw("hist");
      
      cgenPT->cd(c+1);
      h_genPT[c]->Draw("hist");
    }
  cPT->Update();
  cgenPT->Update();
  
  theFile->cd();
  cPT->Write();
  cgenPT->Write();

  TCanvas* cmass = new TCanvas("masses","masses", 200, 200, 800, 800);
  cmass->Divide(2, 2);
  
  cmass->cd(1); 
  h_Z1mass->Draw("hist"); 
  
  cmass->cd(2); 
  h_Z2mass->Draw("hist"); 
  
  cmass->cd(3); 
  h_Hmass->Draw("hist");
  
  cmass->cd(4); 
  h_HmassLarge->Draw("hist");
  
  cmass->Update();
  
  theFile->cd();
  cmass->Write();

  // -----------------------------------------------   
  // get total weighted event count
  // -----------------------------------------------   
  double summedWeight    = h_nEvent->GetBinContent(1);
  double summedWeightUnc = h_nEvent->GetBinError(1);
  
  char summary[512];
  totalPassedUnc = sqrt(totalPassedUnc);
  sprintf(summary,
	   "======================================================\n"
	   "Event count:      %10.3f +/- %-10.3f events\n"
	   "Selected count:   %10.3f +/- %-10.3f events\n"
	   "======================================================\n",
	   summedWeight, summedWeightUnc,
	   totalPassed, totalPassedUnc);
   
   cout << summary << endl;   
   // -----------------------------------------------
   // cut flow
   // -----------------------------------------------
   gStyle->SetOptStat(0);

   cout << endl << "cut flow" << endl;
   for(int i=0; i < h_nEvent->GetNbinsX(); i++)
     {
       double c = h_nEvent->GetBinContent(i+1);
       char rec[80];
       sprintf(rec, "\t%-20s %10.3f",
	       h_nEvent->GetXaxis()->GetBinLabel(i+1), c);
       cout << rec << endl;
     }
   
   h_nEvent->Scale(1.0/summedWeight);
   h_nEvent->SetMaximum(1.e-3);
   h_nEvent->SetMaximum(1.1);
    
   TCanvas* cuts = new TCanvas("cutflow", "cut flow", 500, 10, 500, 500);
   cuts->SetGridy();
   cuts->cd();
   h_nEvent->Draw("hist");
   cuts->Update();
   sprintf(filename, "plots/%s_cutflow.png", namen.c_str());
   cuts->Print(filename, "png");
   
   theFile->cd();
   theFile->Write();
   cuts->Write();
   
   gSystem->Sleep(5000);


   if ( muonEff ) delete muonEff;
   if ( elecEff ) delete elecEff;
}
