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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
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
 TH1F * hist_ ; 

//--------------Variables------------------------
    reco::Vertex pv;
    int nbmuons = -1;
//-------------- Templates________________________

std::vector<bool>  GoodPV;
std::vector<int> NbGoodMuons;
std::vector<bool> tightmuons;
std::vector<bool> pfmuons;
std::vector<bool> loosemuons;
std::vector<float> muonsp;
std::vector<float> muonspt;
std::vector<float> muonspx;
std::vector<float> muonspy;
std::vector<float> muonspz;
std::vector<float> muonseta;
std::vector<float> muonstheta;
std::vector<float> muonsphi;
std::vector<float> muonsvx;
std::vector<float> muonsvy;
std::vector<float> muonsvz;
std::vector<float> muonsdxy;

std::vector<float> muonsPtError;
std::vector<float> muonsEtaError;
std::vector<float> muonsThetaError;
std::vector<float> muonsPhiError;
std::vector<string>HLTnames;
std::vector<bool>  HLTresult;
std::vector<int> pu_info; 





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
      edm::EDGetTokenT                <pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT         < std::vector<pat::TriggerObjectStandAlone>> triggetObjects_; 
      const std::vector<std::string> bDiscriminators_;

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
      edm::Handle                < pat::PackedTriggerPrescales> PrescaleHandle;
      edm::Handle              <std::vector< pat::TriggerObjectStandAlone>> triggerObjectsHandle;
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
  lheEventProductToken_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
  triggerPrescales_(mayConsume<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  triggetObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  bDiscriminators_(iConfig.getParameter<std::vector<std::string> >("bDiscriminators")) {
   //now do what ever initialization is needed
   usesResource("TFileService");
  
   tree_ = fs->make<TTree>("tree", "tree");
 
  
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
  iEvent.getByToken(triggerPrescales_, PrescaleHandle);
  iEvent.getByToken(triggetObjects_, triggerObjectsHandle);

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
   
   initialize(iEvent);
   
   ntuple_.fill_evtinfo(iEvent.id());
   reco::VertexCollection vtx;
   if(vtxHandle.isValid()){ vtx = *vtxHandle;
   pv = vtx.front();
   bool GoodPv = PrimaryVertex(vtx);
   GoodPV.push_back(GoodPv);
   }




   
//************************************* ElectronCollection Study***************************************

   pat::ElectronCollection electrons;
   if(electronsHandle.isValid()){ electrons = *electronsHandle;}

//************************************* MuonCollection Study***************************************



   pat::MuonCollection muons;
   if(muonsHandle.isValid()){ muons = *muonsHandle;
 
   nbmuons=0;  
   for (const pat::Muon mu : muons) {
     if( mu.pt() < 0.0 ) continue;
      if (!( fabs(mu.eta()) < 2.4)) continue;
      ++nbmuons;

      tightmuons.push_back(mu.isTightMuon(pv)); 
      pfmuons.push_back(mu.isPFMuon());   
      loosemuons.push_back(mu.isLooseMuon());    
      muonsp.push_back(mu.muonBestTrack()->p());
      muonspt.push_back(mu.muonBestTrack()->pt());
      muonspx.push_back(mu.muonBestTrack()->px());
      muonspy.push_back(mu.muonBestTrack()->py());
      muonspz.push_back(mu.muonBestTrack()->pz());
      muonseta.push_back(mu.muonBestTrack()->eta());
      muonsphi.push_back(mu.muonBestTrack()->phi());
      muonstheta.push_back(mu.muonBestTrack()->theta());
      muonsvx.push_back(mu.muonBestTrack()->vx());
      muonsvy.push_back(mu.muonBestTrack()->vy());
      muonsvz.push_back(mu.muonBestTrack()->vz());
      muonsdxy.push_back(mu.muonBestTrack()->dxy(pv.position()));

      muonsPtError.push_back(mu.muonBestTrack()->ptError());
      muonsEtaError.push_back(mu.muonBestTrack()->etaError());
      muonsThetaError.push_back(mu.muonBestTrack()->thetaError());
      muonsPhiError.push_back(mu.muonBestTrack()->phiError());
  }
  NbGoodMuons.push_back(nbmuons);
}

//************************************* TauCollection Study***************************************

   pat::TauCollection taus;
   if(tausHandle.isValid()){ taus = *tausHandle;}

   EcalRecHitCollection recHitCollectionEB;
   if(recHitCollectionEBHandle.isValid()){ recHitCollectionEB = *recHitCollectionEBHandle;}

   EcalRecHitCollection recHitCollectionEE;
   if(recHitCollectionEEHandle.isValid()){ recHitCollectionEE = *recHitCollectionEEHandle;}


   pat::PackedCandidateCollection pfCandidates;
   if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }


//*************************************JetCollection Study***************************************

   pat::JetCollection jets;
   if(jetsHandle.isValid()){ jets = *jetsHandle;
    for(const pat::Jet jet : jets ){

        cout<<"Jets pt"<<jet.pt()<<endl;
    for( const std::string &bDiscr : bDiscriminators_ )
    cout<<"B discriminator: "<<jet.bDiscriminator(bDiscr)<<endl;  
}
    



}

//************************************* METCollection Study***************************************

   pat::METCollection mets;
   if(metsHandle.isValid()){ mets = *metsHandle;}
   pat::MET met = mets[0];



//************************************* HLT Paths and Objects Study***************************************


  edm::TriggerResults  HLTresults;
  if(triggerResultsHandle.isValid()){HLTresults = *triggerResultsHandle; 

   const edm::TriggerNames &names =  iEvent.triggerNames(HLTresults);
           //cout<< "\n == TRIGGER PATHS= "<<endl;

     for(unsigned int i = 0, n = HLTresults.size(); i<n ; ++i ) {
           //cout<< "Trigger" << names.triggerName(i) <<", prescale  = " <<PrescaleHandle->getPrescaleForIndex(i) << ":" << (HLTresults.accept(i) ? "PASS" : "fail (or not run)")<<endl;

           //add the name of the trigger and check the trigger result->accept(i) variable.
           //cout<<"trigger name : "<<names.triggerName(i)<<endl;  
          // cout<<"Trigger result is : "<< HLTresults.accept(i)<<endl; 
          // if(HLTresults.accept(i)==1) cout<<"trigger name : "<<names.triggerName(i)<<endl;
        HLTnames.push_back(names.triggerName(i));
        HLTresult.push_back(HLTresults.accept(i));
    }




 std::vector<pat::TriggerObjectStandAlone>  HLTObject;

    if(triggerObjectsHandle.isValid()){ 

        HLTObject = *triggerObjectsHandle;
 
    for(unsigned i =0 , n= HLTObject.size() ; i< n ; ++i){
 
      pat::TriggerObjectStandAlone src = HLTObject.at(i);

   }


    for (pat::TriggerObjectStandAlone obj : HLTObject) {
       obj.unpackPathNames(names);
      // std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

                // Print trigger object collection and type
             
        //         std::cout << "\t   Collection: " << obj.collection() << std::endl;
        //         std::cout << "\t   Type IDs:   ";
       //              for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds( )[h] ;
       //          std::cout << std::endl;

         }

  
} 

}


//************************************* Pile Up Collection info***************************************

std::vector<PileupSummaryInfo> puInfo ; 

if(puInfoH.isValid()){

for(vector<PileupSummaryInfo> ::const_iterator pu = puInfoH->begin(); pu != puInfoH->end() ; ++pu ) {
//cout<< "Number of interactions is : "<< pu->getPU_NumInteractions()<<endl;
pu_info.push_back(pu->getPU_NumInteractions());

}


}

   /*
   double rho = 0;
   if(rhoHandle.isValid()){ rho = *rhoHandle;}
   */

     tree_->Fill();      
}


// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNeutralLeptonAnalysis::beginJob()
{

ntuple_.set_evtinfo(tree_);
tree_->Branch("NbGoodMuons",&NbGoodMuons);
tree_->Branch("GoodPV",&GoodPV);
tree_->Branch("tightmuons",&tightmuons);
tree_->Branch("pfmuons",&pfmuons);
tree_->Branch("loosemuons",&loosemuons);
tree_->Branch("muonsp",&muonsp);
tree_->Branch("muonspt",&muonspt);
tree_->Branch("muonspx",&muonspx);
tree_->Branch("muonspy",&muonspy);
tree_->Branch("muonspz",&muonspz);
tree_->Branch("muonseta",&muonseta);
tree_->Branch("muonsphi",&muonsphi);
tree_->Branch("muonstheta",&muonstheta);
tree_->Branch("muonsvx",&muonsvx);
tree_->Branch("muonsvy",&muonsvy);
tree_->Branch("muonsvz",&muonsvz);
tree_->Branch("muonsdxy",&muonsdxy);

//tree_->Branch("muonsPresol",&muonsPresol);
tree_->Branch("muonsPtError",&muonsPtError);
tree_->Branch("muonsThetaError",&muonsThetaError);
tree_->Branch("muonsEtaError",&muonsEtaError);
tree_->Branch("muonsPhiError",&muonsPhiError);

tree_->Branch("muons",&muonsPhiError);
tree_->Branch("muonsPhiError",&muonsPhiError);
tree_->Branch("muonsPhiError",&muonsPhiError);

tree_->Branch("HLTnames",&HLTnames);
tree_->Branch("HLTresult",&HLTresult);
tree_->Branch("pu_info", &pu_info);

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
