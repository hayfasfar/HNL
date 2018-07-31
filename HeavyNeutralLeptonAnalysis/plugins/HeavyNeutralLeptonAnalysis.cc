// -*- C++ -*-
//
// Package:    HNL/HeavyNeutralLeptonAnalysis
// Class:      HeavyNeutralLeptonAnalysis
// 
/**\class HeavyNeutralLeptonAnalysis HeavyNeutralLeptonAnalysis.cc HNL/HeavyNeutralLeptonAnalysis/plugins/HeavyNeutralLeptonAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  root
//         Created:  Fri, 23 Mar 2018 12:09:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for the plotting
#include "HNL/HeavyNeutralLeptonAnalysis/interface/SmartSelectionMonitor.h"

// root includes
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include <Math/VectorUtil.h>



//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"

using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HeavyNeutralLeptonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HeavyNeutralLeptonAnalysis(const edm::ParameterSet&);
      ~HeavyNeutralLeptonAnalysis();
      
      bool PrimaryVertex( const reco::VertexCollection &vtx);
      //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);




   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void initialize(const edm::Event&); 
      virtual void endJob() override;



 edm::Service<TFileService> fs;
 
  TTree * tree_;
 BigNtuple ntuple_;


//--------------Variables------------------------


  // ----------  Files  ---------------------------


//--------------template-------------------------

//std::vector<bool>  GoodPV;
//
//std::vector<int> NbGoodMuons;





      // ----------member data ---------------------------
      bool   isMC;

      edm::EDGetTokenT         < reco::VertexCollection > vtxMiniAODToken_;
      edm::EDGetTokenT                         < double > rhoToken_;
      edm::EDGetTokenT            < pat::MuonCollection > muonsMiniAODToken_;
      edm::EDGetTokenT        < pat::ElectronCollection > electronsMiniAODToken_;
      edm::EDGetTokenT           < EcalRecHitCollection > recHitEBToken_;
      edm::EDGetTokenT           < EcalRecHitCollection > recHitEEToken_;
      edm::EDGetTokenT             < pat::TauCollection > tausMiniAODToken_;
      edm::EDGetTokenT < pat::PackedCandidateCollection > packedCandidateToken_;
      edm::EDGetTokenT             < pat::JetCollection > jetsMiniAODToken_;
      edm::EDGetTokenT             < pat::METCollection > pfMETAODToken_;
      edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;
      edm::EDGetTokenT            < edm::TriggerResults > metFilterResultsToken_;
      edm::EDGetTokenT    < reco::GenParticleCollection > genParticleToken_;
      edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
      edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
      edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;

    protected:
      edm::Handle         < reco::VertexCollection > vtxHandle;
      edm::Handle                         < double > rhoHandle;
      edm::Handle            < pat::MuonCollection > muonsHandle;
      edm::Handle        < pat::ElectronCollection > electronsHandle;
      edm::Handle           < EcalRecHitCollection > recHitCollectionEBHandle;
      edm::Handle           < EcalRecHitCollection > recHitCollectionEEHandle;
      edm::Handle             < pat::TauCollection > tausHandle;
      edm::Handle < pat::PackedCandidateCollection > pfCandidatesHandle;
      edm::Handle             < pat::JetCollection > jetsHandle;
      edm::Handle             < pat::METCollection > metsHandle;
      edm::Handle            < edm::TriggerResults > triggerResultsHandle;
      edm::Handle            < edm::TriggerResults > metFilterResultsHandle;
      /* Only for MC */
      edm::Handle    < reco::GenParticleCollection > genHandle;
      edm::Handle            < GenEventInfoProduct > genEventInfoHandle;
      edm::Handle < std::vector<PileupSummaryInfo> > puInfoH;
      edm::Handle                < LHEEventProduct > lheEPHandle;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HeavyNeutralLeptonAnalysis::HeavyNeutralLeptonAnalysis(const edm::ParameterSet& iConfig):
  ntuple_(),
  isMC(iConfig.getParameter<bool>("isMC")),
  vtxMiniAODToken_(mayConsume<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
  rhoToken_(mayConsume<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonsMiniAODToken_(mayConsume<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronsMiniAODToken_(mayConsume<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  recHitEBToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEBSrc"))),
  recHitEEToken_(mayConsume< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("recHitCollectionEESrc"))),
  tausMiniAODToken_(mayConsume< pat::TauCollection >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
  packedCandidateToken_(mayConsume< pat::PackedCandidateCollection >(iConfig.getParameter<edm::InputTag>("packCandSrc"))),
  jetsMiniAODToken_(mayConsume< pat::JetCollection >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  pfMETAODToken_(mayConsume<pat::METCollection>(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultSrc"))),
  metFilterResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterResultSrc"))),
  genParticleToken_(mayConsume<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  genEventInfoToken_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts")))
{
   //now do what ever initialization is needed
  usesResource("TFileService");
  
  tree_ = fs->make<TTree>("tree_", "tree_");
  ntuple_.set_evtinfo(tree_);
  
}


HeavyNeutralLeptonAnalysis::~HeavyNeutralLeptonAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

void HeavyNeutralLeptonAnalysis::initialize(const edm::Event& iEvent){
  iEvent.getByToken(vtxMiniAODToken_, vtxHandle);
  iEvent.getByToken(rhoToken_, rhoHandle);
  iEvent.getByToken(muonsMiniAODToken_, muonsHandle);
  iEvent.getByToken(electronsMiniAODToken_, electronsHandle);
  iEvent.getByToken(recHitEBToken_, recHitCollectionEBHandle);
  iEvent.getByToken(recHitEEToken_, recHitCollectionEEHandle);
  iEvent.getByToken(tausMiniAODToken_, tausHandle);
  iEvent.getByToken(packedCandidateToken_, pfCandidatesHandle);
  iEvent.getByToken(jetsMiniAODToken_, jetsHandle);
  iEvent.getByToken(pfMETAODToken_, metsHandle);
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle);
  iEvent.getByToken(metFilterResultsToken_, metFilterResultsHandle);

  if (isMC){
    iEvent.getByToken(genParticleToken_, genHandle);
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    iEvent.getByToken(PUInfoToken_, puInfoH);
    iEvent.getByToken(lheEventProductToken_, lheEPHandle);
  }
}



//
// member functions

bool HeavyNeutralLeptonAnalysis::PrimaryVertex( const reco::VertexCollection &vtx)
{
  int nbGoodPv = 0;
  bool result = false;

     for(reco::VertexCollection::const_iterator PV = vtx.begin(); PV!=vtx.end();++PV) 
    {
       if(!PV->isFake()) {
         if(PV->ndof() > 4 && fabs(PV->position().z()) <= 24 && fabs(PV->position().rho()) <= 2 ) nbGoodPv++;
       }
    }
 
  if( nbGoodPv>= 1) result = true;
  return result;
}







// ------------ method called for each event  ------------
void
HeavyNeutralLeptonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){



   using namespace edm;
   
   //initialize(iEvent);
   
   ntuple_.fill_evtinfo(iEvent.id());
    
    

   
 //  GoodPV.clear();
   reco::VertexCollection vtx;
   if(vtxHandle.isValid()){ vtx = *vtxHandle;
   bool GoodPv = PrimaryVertex(vtx);
  // GoodPV.push_back(GoodPv);
   }



   /*
   double rho = 0;
   if(rhoHandle.isValid()){ rho = *rhoHandle;}
   */


  // NbGoodMuons.clear();
   int NbMuons = 0;
   pat::MuonCollection muons;
   if(muonsHandle.isValid()){ muons = *muonsHandle;

   for (const pat::Muon mu : muons) {
     if( mu.pt() < 0.0 ) continue;
      if (!( fabs(mu.eta()) < 2.4 && mu.pt() > 5. )) continue;
       ++NbMuons;
  }
    //  NbGoodMuons.push_back(NbMuons);
      cout<<"nb of muons"<< NbMuons<<endl;
}



   EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}

   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}

   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}

   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

   
   pat::ElectronCollection electrons;
   if(electronsHandle.isValid()){ electrons = *electronsHandle;}

   pat::JetCollection jets;
   if(jetsHandle.isValid()){ jets = *jetsHandle;}

   pat::METCollection mets;
   if(metsHandle.isValid()){ mets = *metsHandle;}
   pat::MET met = mets[0];

}


// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNeutralLeptonAnalysis::beginJob()
{



}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNeutralLeptonAnalysis::endJob() 
{

//outputFile_->cd();
//tree_->Write();
//outputFile_->Close(); 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
HeavyNeutralLeptonAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}*/
//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNeutralLeptonAnalysis);
