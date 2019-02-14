#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for the plotting
 #include "HNL/HeavyNeutralLeptonAnalysis/interface/SmartSelectionMonitor.h"
//
// // root includes
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

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
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
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"

#include "HNL/DisplacedSVAssociator/interface/VertexAssociation.h"

#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"



class HeavyNeutralLeptonAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HeavyNeutralLeptonAnalysis(const edm::ParameterSet&);
  ~HeavyNeutralLeptonAnalysis();

  reco::VertexCollection getMatchedVertex_Muon(const pat::Muon & mu, const reco::VertexCollection& vertexCollection);
  reco::VertexCollection PrimaryVertex( const reco::VertexCollection &vtx);
    bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate* particle);
    double MatchGenMuon(const edm::Event&,  reco::TrackRef BestTrack , int pdgId);
//  double MatchGenElectron(const edm::Event& iEvent, const pat::Electron& ele_, int pdgId);
  double MatchGenVertex(const edm::Event& iEvent, reco::Vertex vertex, int pdgId);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void initialize(const edm::Event&);
  virtual void endJob() override;

  edm::Service<TFileService> fs;
  TTree * tree_;
  BigNtuple ntuple_;


  //--------------Variables------------------------
  int    debug;
  float pvCompatibilityScore = .05;
  float indexMu1 = -10 ;
  float indexMu2 = -10 ; 
  float npT = -1;
  float npIT = -1;

  edm::LumiReWeighting Lumiweights_;
  edm::LumiReWeighting LumiweightsUp_;
  edm::LumiReWeighting LumiweightsDown_;
  bool   isMC;
  bool   isElectonStudy;
  bool   isMCSignal;
  bool   isData ;

  math::XYZPoint gen_PvPosition ;
  TLorentzVector  gen_pv ; 
  SmartSelectionMonitor mon;

  //--------------template-------------------------
   std::vector<float> dRStudy ;
  // ----------member data ---------------------------

  edm::EDGetTokenT         < reco::VertexCollection > vtxMiniAODToken_;
  edm::EDGetTokenT                         < double > rhoToken_;
  edm::EDGetTokenT            < pat::MuonCollection > muonsMiniAODToken_;
  edm::EDGetTokenT        < pat::ElectronCollection > electronsMiniAODToken_;
  edm::EDGetTokenT           < EcalRecHitCollection > recHitEBToken_;
  edm::EDGetTokenT           < EcalRecHitCollection > recHitEEToken_;
  edm::EDGetTokenT             < pat::TauCollection > tausMiniAODToken_;
  edm::EDGetTokenT < pat::PackedCandidateCollection > packedCandidateToken_;
  edm::EDGetTokenT <pat::PackedCandidateCollection >  lostTracks_  ; 
  edm::EDGetTokenT            <pat::JetCollection>    jetsMiniAODToken_;
  edm::EDGetTokenT             < pat::METCollection > pfMETAODToken_;
  edm::EDGetTokenT            < edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT            < edm::TriggerResults > metFilterResultsToken_;
  edm::EDGetTokenT    < reco::GenParticleCollection > genParticleToken_;
  edm::EDGetTokenT    < pat::PackedGenParticleCollection > genPackedParticleToken_;
  edm::EDGetTokenT            < GenEventInfoProduct > genEventInfoToken_;
  edm::EDGetTokenT < std::vector<PileupSummaryInfo> > PUInfoToken_;
  edm::EDGetTokenT                < LHEEventProduct > lheEventProductToken_;
  edm::EDGetTokenT         < reco::VertexCollection > inclusiveSecondaryVertices_;
  edm::EDGetTokenT          < std::vector<reco::VertexCompositePtrCandidate>> slimmedSecondaryVertices_ ;

  const std::vector        <std::string> bDiscriminators_;

  /*edm::EDGetTokenT         <edm::ValueMap<float>>     eleMvaToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleVetoToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleLooseToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleMediumToken_;
  edm::EDGetTokenT         <edm::ValueMap<bool>>      eleTightToken_;*/
  edm::EDGetTokenT         <reco::GenJetCollection>  GenJetsMiniAODToken_ ;


protected:
  edm::Handle         < reco::VertexCollection > vtxHandle;
  edm::Handle                         < double > rhoHandle;
  edm::Handle            < pat::MuonCollection > muonsHandle;
  edm::Handle       <std::vector<pat::Electron>> electronsHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEBHandle;
  edm::Handle           < EcalRecHitCollection > recHitCollectionEEHandle;
  edm::Handle             < pat::TauCollection > tausHandle;
  edm::Handle < pat::PackedCandidateCollection > pfCandidatesHandle;
  edm::Handle < pat::PackedCandidateCollection > lostTracksHandle;
  edm::Handle             < pat::JetCollection > jetsHandle;
  edm::Handle             < pat::METCollection > metsHandle;
  edm::Handle            < edm::TriggerResults > triggerResultsHandle;
  edm::Handle            < edm::TriggerResults > metFilterResultsHandle;
  edm::Handle         < reco::VertexCollection > secondaryVertexHandle;
  edm::Handle         <std::vector<reco::VertexCompositePtrCandidate>> slimmedSecondaryVerticesHandle;

  /*edm::Handle         <edm::ValueMap<float>> electronsMva;
  edm::Handle         <edm::ValueMap<bool>> electronsVeto;
  edm::Handle         <edm::ValueMap<bool>> electronsLoose;
  edm::Handle         <edm::ValueMap<bool>> electronsMedium;
  edm::Handle         <edm::ValueMap<bool>> electronsTight;*/


  /* Only for MC */
  edm::Handle    < reco::GenParticleCollection > genHandle;
  edm::Handle    < pat::PackedGenParticleCollection > genPackedHandle;

  edm::Handle            < GenEventInfoProduct > genEventInfoHandle;
  edm::Handle < std::vector<PileupSummaryInfo> > puInfoH;
  edm::Handle                < LHEEventProduct > lheEPHandle;
  edm::Handle            <reco::GenJetCollection>  GenJetHandle;
};


