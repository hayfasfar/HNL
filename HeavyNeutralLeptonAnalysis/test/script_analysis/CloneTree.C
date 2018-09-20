//**********************************************************************************************************************************
// Remove some branches + selects the events + add variables -- for muonic channel
//***************************************** To Compile******************************************************************************
// g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
//**********************************************************************************************************************************

#ifndef __CINT__
#include "RooGlobalFunc.h"
//------------------------------------------------                                                                                                           
#endif
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"
#include <vector>
#include <string>
#include <iostream>
#include "RooRandom.h"
#include "RooMinuit.h"
#include "TRandom3.h"
#include <time.h>
#include <TROOT.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <iostream>
#include <TMath.h>
#include "TH1D.h"
#include "TH2.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooArgusBG.h"
#include "TString.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMappedCategory.h"
#include "RooCmdArg.h"
#include "RooChebychev.h"
#include "RooUnblindUniform.h"
#include "RooUnblindPrecision.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooSimWSTool.h"
#include "RooWorkspace.h"
#include <TLatex.h>
#include "RooFit.h"
#include "RooConstVar.h"
#include "RooSimPdfBuilder.h"
#include "RooStringVar.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TLorentzVector.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;




int main(){
  

  //Get old file, old tree and set top branch address
  //TFile *oldfile = new TFile("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu/Background_Analysis.root");
  TFile *oldfile = new TFile("/afs/cern.ch/user/j/jpriscia/work/Background_Analysis.root");
  TTree *oldtree = (TTree*)oldfile->Get("HeavyNeutralLepton/tree_");


  //////selection cuts

  Float_t isoCut = 0.15;


  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  //TChain *oldtree = new TChain("SigmaPMuMuTuple/DecayTree");
  /*oldtree->Add(sigmapmumu2016DownNEW);
  oldtree->Add(sigmapmumu2016UpNEW);
  */

  Long64_t nentries = oldtree->GetEntries();
  cout  << nentries<<endl;

  // These are the variables I cut on 
  Bool_t passIsoMu24All;
  oldtree->SetBranchAddress("passIsoMu24All",&passIsoMu24All);

  vector<Float_t>  *mu_isTightMuon = 0;
  vector<Float_t>  *mu_isLoose = 0;
  vector<Float_t>  *mu_pt = 0;
  vector<Float_t>  *mu_eta =0;
  vector<Float_t>  *mu_phi=0; 
  vector<Float_t>  *mu_charge=0; 
  vector<Float_t>  *mu_en=0;
  vector<Float_t>  *mu_et=0;
  vector<Float_t>  *deltaBeta=0;

  oldtree->SetBranchAddress("mu_isTightMuon",&mu_isTightMuon);
  oldtree->SetBranchAddress("mu_isLoose",&mu_isLoose);
  oldtree->SetBranchAddress("mu_pt",&mu_pt);
  oldtree->SetBranchAddress("mu_eta",&mu_eta);
  oldtree->SetBranchAddress("mu_phi",&mu_phi);
  oldtree->SetBranchAddress("mu_charge",&mu_charge);
  oldtree->SetBranchAddress("mu_en",&mu_en);
  oldtree->SetBranchAddress("mu_et",&mu_et);
  oldtree->SetBranchAddress("mu_recoDeltaBeta",&deltaBeta);


  //Throw away some branches -- thos are all the  one I don't want to keep in the ntuples
  //IMPORTANT: we cannot throw away the branches with the variable we are using in the loop!
  oldtree->SetBranchStatus("ele_*",    0);
  oldtree->SetBranchStatus("mu_ST*", 0);

			
  //Create a new file + a clone of old tree in new file 
  TFile *newfile = new TFile("/eos/cms/store/group/phys_exotica/HNL/Background/crab_Analysis_WZToLLLNu/skimmedSample.root","recreate");
  TTree *newtree = oldtree->CloneTree(0);
  cout<<"cloning done"<<endl;

  
  Float_t mu_promptPt,mu_promptEta, mu_promptPhi, mu_promptCharge, mu_promptEt, mu_promptE;
  TBranch* branch_mu_promptPt = newtree->Branch("mu_promptPt",&mu_promptPt,"mu_promptPt/F");
  TBranch* branch_mu_promptEta = newtree->Branch("mu_promptEta",&mu_promptEta,"mu_promptEta/F");
  TBranch* branch_mu_promptPhi = newtree->Branch("mu_promptPhi",&mu_promptPhi,"mu_promptPhi/F");
  TBranch* branch_mu_promptCharge = newtree->Branch("mu_promptCharge",&mu_promptCharge,"mu_promptCharge/F");
  TBranch* branch_mu_promptE = newtree->Branch("mu_promptE",&mu_promptE,"mu_promptE/F");
  TBranch* branch_mu_promptEt = newtree->Branch("mu_promptEt",&mu_promptEt,"mu_promptEt/F");
  
  
  for (int i=0;i<oldtree->GetEntries(); i++) {
    if (i%10000==0) cout<<i<<endl;
    oldtree->GetEntry(i);

    if (passIsoMu24All==0) {continue;}  // cut on the trigger!
    
    Float_t   minPt = -1000;
    Int_t  myIndex = -1;

    for(size_t t=0; t<mu_isTightMuon->size(); t++){

      if (mu_isTightMuon->at(t)==0. || deltaBeta->at(t)>isoCut) continue;    // I took the tight isolation..can be changed, we have to check the efficiency
      if (mu_pt->at(t)>minPt){
	minPt=mu_pt->at(t);
	myIndex=t;	  
      }
    }	

    //here I save only the info of the prompt muon -- for the sv muon we can clone the branch as it is..
    if (myIndex == -1) continue;
    mu_promptPt = mu_pt->at(myIndex);
    mu_promptEta = mu_eta->at(myIndex);
    mu_promptPhi = mu_phi->at(myIndex);
    mu_promptCharge = mu_charge->at(myIndex);
    mu_promptE = mu_en->at(myIndex);
    mu_promptEt = mu_et->at(myIndex);
    

    newtree->Fill();

  }


  newtree->Print();
  newtree->AutoSave();

  delete newfile;

  return 0;
}
