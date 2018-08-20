#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "TTree.h"

void BigNtuple::set_evtinfo(TTree* tree) {
  tree->Branch("run" , &run_, "run/I");
  tree->Branch("lumi", &lumi_, "lumi/I");
  tree->Branch("evt" , &evt_, "evt/L");
  //tree->Branch("NbGoodMuons",&NbGoodMuons);
}

void BigNtuple::fill_evtinfo(const edm::EventID& id) {
	lumi_ = id.run();
        run_  = id.luminosityBlock();
        evt_  = id.event();
}
