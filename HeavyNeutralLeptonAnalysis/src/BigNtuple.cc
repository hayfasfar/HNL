#include "FWCore/Framework/interface/Event.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "TTree.h"

void BigNtuple::set_evtinfo(TTree* tree) {
	tree->Branch("run" , &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt" , &evt_, "evt/i");
}

void BigNtuple::fill_evtinfo(const edm::EventID& id) {
	lumi_ = id.run();
  run_  = id.luminosityBlock();
  evt_  = id.event();
}
