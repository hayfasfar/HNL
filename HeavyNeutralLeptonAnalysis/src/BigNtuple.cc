#include "FWCore/Framework/interface/Event.h"
#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
#include "TTree.h"

void BigNtuple::set_evtInfo(TTree* tree) {
	tree->Branch("run" , &run_, "run/i");
	tree->Branch("lumi", &lumi_, "lumi/i");
	tree->Branch("evt" , &evt_, "evt/i");
}

void BigNtuple::fill_evtInfo(const edm::EventID& id) {
	lumi_ = id.run();
	run_  = id.luminosityBlock();
	evt_  = id.event();
}


void BigNtuple::set_pvInfo(TTree* tree){
       tree->Branch("pvX" , &pvX_, "pvX/i");
       tree->Branch("pvY" , &pvY_, "pvY/i");
       tree->Branch("pvZ" , &pvZ_, "pvZ/i");
       tree->Branch("pvXErr" , &pvXErr_, "pvXErr/i");
       tree->Branch("pvYErr" , &pvYErr_, "pvYErr/i");
       tree->Branch("pvZErr" , &pvZErr_, "pvZErr/i");
       tree->Branch("pvMass" , &pvMass_, "pvMass/i");
       tree->Branch("pvLxy" , &pvLxy_, "pvLxy/i");
       tree->Branch("pvLxyz" , &pvLxyz_, "pvLxyz/i");
       tree->Branch("pvLxySig" , &pvLxySig_, "pvLxySig/i");
       tree->Branch("pvLxyzSig" , &pvLxyzSig_, "pvLxyzSig/i");
       tree->Branch("pvChi2" , &pvChi2_, "pvChi2/i");
       tree->Branch("pvNTrack" , &pvNTrack_, "pvNTrack/i");
       tree->Branch("pvSumPtSq" , &pvSumPtSq_, "pvSumPtSq/i");
       tree->Branch("numberPV" , &numberPV_, "numberPV/i");
       
}

void BigNtuple::fill_pvInfo(const std::vector<reco::VertexCollection>& pvs){
       numberPV_ = pvs->size();
       pv = pvs -> front();

       float x  = pv->x(), y = pv->y(), z = pv->z();
       float xE = pv->xError(), yE = pv->yError(), zE = pv->zError();
       
       pvX_ = x;
       pvY_ = y;
       pvZ_ = z;
       pvXErr_ = xE;
       pvYErr_ = yE;
       pvZErr_ = zE;
       pvMass_ = 0;
       pvLxy_ = std::sqrt( x * x + y * y );
       pvLxyz_ = std::sqrt( x * x + y * y + z * z);
       pvLxySig_ = std::sqrt( x * x + y * y ) / std::sqrt(xE * xE + yE * yE);
       pvLxyzSig_ = std::sqrt( x * x + y * y + z * z )/ std::sqrt(xE * xE + yE * yE + zE * zE);
       pvChi2_ = pv->chi2();

       reco::Vertex::trackRef_iterator vtxIter = pv->tracks_begin();
       float  SumPtSq =  0;
       int NTrack = 0;
       for(; vtxIter != pv->tracks_end(); ++vtxIter) {
         NTrack++;
         SumPtSq += (*vtxIter)->pt() * (*vtxIter)->pt();
       }
       pvNTrack_.push_back(NTrack);
       pvSumPtSq_.push_back(SumPtSq);       
}



void BigNtuple::set_svInfo(TTree* tree){

    tree->Branch("sv_TrackSize" , &sv_TrackSize_, "sv_TrackSize/i");
    tree->Branch("sv_LXYSig" , &sv_LXYSig_, "sv_LXYSig/i");
    tree->Branch("sv_LXYZSig" , &sv_LXYZSig_, "sv_LXYZSig/i");
    tree->Branch("sv_LXY" , &sv_LXY_, "sv_LXY/i");
    tree->Branch("sv_LXYZ" , &sv_LXYZ_, "sv_LXYZ/i");
    tree->Branch("sv_mass" , &sv_mass_, "sv_mass/i");
    tree->Branch("sv_charge" , &sv_charge_, "sv_charge/i");
    tree->Branch("sv_eta" , &sv_eta_, "sv_eta/i");
    tree->Branch("sv_phi" , &sv_phi_, "sv_phi/i");
    tree->Branch("sv_pt" , &sv_pt_, "sv_pt/i");
    tree->Branch("sv_p" , &sv_p_, "sv_p/i");
    tree->Branch("sv_Beta" , &sv_Beta_, "sv_Beta/i");
    tree->Branch("sv_Gamma" , &sv_Gamma_, "sv_Gamma/i");
    tree->Branch("sv_CTau0" , &sv_CTau0_, "sv_CTau0/i");
    tree->Branch("sv_NDof" , &sv_NDof_, "sv_NDof/i");
    tree->Branch("sv_Chi2" , &sv_Chi2_, "sv_Chi2/i");
    tree->Branch("sv_Angle3D" , &sv_Angle3D_, "sv_Angle3D/i");
    tree->Branch("sv_Angle2D" , &sv_Angle2D_, "sv_Angle2D/i");
    tree->Branch("sv_bestScore" , &sv_bestScore_, "sv_bestScore/i");
    tree->Branch("sv2_Match" , &sv2_Match_, "sv2_Match/i");
    tree->Branch("sv_tracks_charge" , &sv_tracks_charge_, "sv_tracks_charge/i");
    tree->Branch("sv_tracks_eta" , &sv_tracks_eta_, "sv_tracks_eta/i");
    tree->Branch("sv_tracks_phi" , &sv_tracks_phi_, "sv_tracks_phi/i");
    tree->Branch("sv_tracks_pt" , &sv_tracks_pt_, "sv_tracks_pt/i");
    tree->Branch("sv_tracks_dxySig" , &sv_tracks_dxySig_, "sv_tracks_dxySig/i");
    tree->Branch("sv_tracks_dxy" , &sv_tracks_dxy_, "sv_tracks_dxy/i");
    tree->Branch("sv_tracks_dxyz" , &sv_tracks_dxyz_, "sv_tracks_dxyz/i");    
    tree->Branch("sv_tracks_Sumcharge" , &sv_tracks_Sumcharge_, "sv_tracks_Sumcharge/i");
    tree->Branch("sv_tracks_Sumpt" , &sv_tracks_Sumpt_, "sv_tracks_Sumpt/i");

}


void BigNtuple::fill_svInfo(const reco::Vertex bestVertex, const float bestVertexScore){


  float  svChi2 = bestVertex.chi2();
  float  svNDof = bestVertex.ndof();

  //flight distance from the firstPV                                                                                                                                                        
  float x  = bestVertex.x(), y = bestVertex.y(), z = bestVertex.z();
  float dx = x - selPV.x() , dy = y - selPV.y(), dz = z - selPV.z();

  // set the compatibility score                                                                                                                                                            
  float  selIVFIsPVScore = std::sqrt((dx/x)*(dx/x) + (dy/y)*(dy/y) + (dz/z)*(dz/z));
  selIVFIsPV = selIVFIsPVScore < pvCompatibilityScore; // default 5% consitency check                                                                                                       

  //build the total error                                                                                                                                                                   
  float svxE = bestVertex.xError(), svyE = bestVertex.yError(), svzE = bestVertex.zError();
  float pvxE = selPV.xError(), pvyE = selPV.yError(), pvzE = selPV.zError();
  float xE   = std::sqrt(svxE * svxE + pvxE * pvxE), yE = std::sqrt(svyE * svyE + pvyE * pvyE), zE = std::sqrt(svzE * svzE + pvzE * pvzE);

  // mother beta, gamma, ctau                                                                                                                                                               
  float   beta_mom  = bestVertex.p4().P() / bestVertex.p4().energy();
  float   gamma_mom = bestVertex.p4().energy() / bestVertex.p4().mass();

  TVector3 pvVector3D(selPV.x(), selPV.y(), selPV.z());
  TVector3 pvVector2D(selPV.x(), selPV.y(), 0);
  TVector3 svVector3D(x, y, z);
  TVector3 svVector2D(x, y, 0);

  // line pointing form the primary vertex through the sceondary vertex                                                                                                                     
  TVector3 svMom3D( bestVertex.p4().x(), bestVertex.p4().y(), bestVertex.p4().z());
  TVector3 svMom2D( bestVertex.p4().x(), bestVertex.p4().y(), 0);

  // you want the negative when the momentum and sv are in the same                                                                                                                         
  // direction relative to the PV                                                                                                                                                           
  // this makes sure the angle is not pi when the vertex is fit behind                                                                                                                      
  // the primary vertex                                                                                                                                                                     
  float sign2D =  (svMom2D * (svVector2D - pvVector2D)) > 0 ? -1: 1;
  float sign3D =  (svMom3D * (svVector3D - pvVector3D)) > 0 ? -1: 1;

  TVector3 pvToVertex3D( sign3D * dx, sign3D * dy, sign3D * dz);
  TVector3 pvToVertex2D( sign2D * dx, sign2D * dy, 0);

  float  svAngle3D = pvToVertex3D.Angle(svMom3D);
  float  svAngle2D = pvToVertex2D.Angle(svMom2D);

  sv_TrackSize_ = selIVFIsPV ? 0 : bestVertex.nTracks();
  sv_LXY_ = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy );
  sv_LXYZ_ = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz );
  sv_LXYSig_ = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE);
  sv_LXYZSig_ = selIVFIsPV ? 0 :  std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE);
  sv_mass_ = selIVFIsPV ? 0 : bestVertex.p4().mass();
  sv_eta_ = selIVFIsPV ? 0 : bestVertex.p4().eta();
  sv_phi_ = selIVFIsPV ? 0 : bestVertex.p4().phi();
  sv_pt_ = selIVFIsPV ? 0 : bestVertex.p4().pt();
  sv_p_ = selIVFIsPV ? 0 : bestVertex.p4().P();
  sv_Beta_ = selIVFIsPV ? 0 : beta_mom;
  sv_Gamma_ = selIVFIsPV ? 0 : gamma_mom;
  sv_CTau0_ = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz) / (beta_mom * gamma_mom);
  sv_NDof_ = selIVFIsPV ? 0 : svNDof;
  sv_Chi2_ = selIVFIsPV ? 0 : svChi2;
  sv_Angle3D_ = selIVFIsPV ? -3 : svAngle3D;
  sv_Angle2D_ = selIVFIsPV ? -3 : svAngle2D;

  sv_bestScore_.push_back(bestVertexScore);

  int ch = 0;
  float pt = 0;

  if(!selIVFIsPV && bestVertexScore > 0 ){ //what is this selIVFIsPV??
    
    reco::Vertex::trackRef_iterator tt = bestVertex.tracks_begin();
    for(; tt != bestVertex.tracks_end(); ++tt) {
      
      sv_tracks_charge_.push_back((*tt)->charge());
      sv_tracks_eta_.push_back((*tt)->eta());
      sv_tracks_phi_.push_back((*tt)->phi());
      sv_tracks_pt_.push_back((*tt)->pt());
      sv_tracks_dxySig_.push_back(fabs((*tt)->dxy(selPV.position()))/fabs((*tt)->dxyError()));
      sv_tracks_dxy_.push_back((*tt)->dxy(selPV.position()));
      
      ROOT::Math::SVector<double, 3> lxyz1((*tt)->vx()-selPV.position().x(), (*tt)->vy()-selPV.position().y(), (*tt)->vz()-selPV.position().z());
      float dxyz = (float)ROOT::Math::Mag(lxyz1); // magntude of the vector                                                                                                                 
      sv_tracks_dxyz_.push_back(dxyz);
      ch+=(*tt)->charge();
      pt+=(*tt)->pt();
    }
    sv_tracks_Sumcharge_.push_back(ch);
    sv_tracks_Sumpt_.push_back(pt);
  }



	 
}


 void BigNtuple::set_muInfo(TTree* tree){
   tree->Branch("mu_nbMuon" , &mu_nbMuon_, "mu_nbMuon/i");
   tree->Branch("mu_en" , &mu_en_, "mu_en/i");
   tree->Branch("mu_pt" , &mu_pt_, "mu_pt/i");
   tree->Branch("mu_eta" , &mu_eta_, "mu_eta/i");
   tree->Branch("mu_phi" , &mu_phi_, "mu_phi/i");
   tree->Branch("mu_et" , &mu_et_, "mu_et/i");
   tree->Branch("mu_charge" , &mu_charge_, "mu_charge/i");
   tree->Branch("mu_trackiso" , &mu_trackiso_, "mu_trackiso/i");
   tree->Branch("mu_pfSumChargedHadronPt" , &mu_pfSumChargedHadronPt_, "mu_pfSumChargedHadronPt/i");
   tree->Branch("mu_pfSumNeutralHadronEt" , &mu_pfSumNeutralHadronEt_, "mu_pfSumNeutralHadronEt/i");
   tree->Branch("mu_PFSumPhotonEt" , &mu_PFSumPhotonEt_, "mu_PFSumPhotonEt/i");
   tree->Branch("mu_pfSumPUPt" , &mu_pfSumPUPt_, "mu_pfSumPUPt/i");
   tree->Branch("mu_numberOfValidMuonHits" , &mu_numberOfValidMuonHits_, "mu_numberOfValidMuonHits/i");
   tree->Branch("mu_emIso" , &mu_emIso_, "mu_emIso/i");
   tree->Branch("mu_hadIso" , &mu_hadIso_, "mu_hadIso/i");
   tree->Branch("mu_normalizedChi2" , &mu_normalizedChi2_, "mu_normalizedChi2/i");
   tree->Branch("mu_numberOfMatchedStations" , &mu_numberOfMatchedStations_, "mu_numberOfMatchedStations/i");
   tree->Branch("mu_numberOfValidPixelHits" , &mu_numberOfValidPixelHits_, "mu_numberOfValidPixelHits/i");
   tree->Branch("mu_numberOftrackerLayersWithMeasurement" , &mu_numberOftrackerLayersWithMeasurement_, "mu_numberOftrackerLayersWithMeasurement/i");
   tree->Branch("mu_numberOfpixelLayersWithMeasurement" , &mu_numberOfpixelLayersWithMeasurement_, "mu_numberOfpixelLayersWithMeasurement/i");
   tree->Branch("mu_TrackQuality" , &mu_TrackQuality_, "mu_TrackQuality/i");
   tree->Branch("mu_InnerTrackQuality" , &mu_InnerTrackQuality_, "mu_InnerTrackQuality/i");
   tree->Branch("mu_pxTunePMuonBestTrack" , &mu_pxTunePMuonBestTrack_, "mu_pxTunePMuonBestTrack/i");
   tree->Branch("mu_pyTunePMuonBestTrack" , &mu_pyTunePMuonBestTrack_, "mu_pyTunePMuonBestTrack/i");
   tree->Branch("mu_pzTunePMuonBestTrack" , &mu_pzTunePMuonBestTrack_, "mu_pzTunePMuonBestTrack/i");
   tree->Branch("mu_pTunePMuonBestTrack" , &mu_pTunePMuonBestTrack_, "mu_pTunePMuonBestTrack/i");
   tree->Branch("mu_etaTunePMuonBestTrack" , &mu_etaTunePMuonBestTrack_, "mu_etaTunePMuonBestTrack/i");
   tree->Branch("mu_LXYZ" , &mu_LXYZ_, "mu_LXYZ/i");
   tree->Branch("mu_LXY" , &mu_LXY_, "mu_LXY/i");
   tree->Branch("mu_ptTunePMuonBestTrack" , &mu_ptTunePMuonBestTrack_, "mu_ptTunePMuonBestTrack/i");
   tree->Branch("mu_phiTunePMuonBestTrack" , &mu_phiTunePMuonBestTrack_, "mu_phiTunePMuonBestTrack/i");
   tree->Branch("mu_thetaTunePMuonBestTrack" , &mu_thetaTunePMuonBestTrack_, "mu_thetaTunePMuonBestTrack/i");
   tree->Branch("mu_chargeTunePMuonBestTrack" , &mu_chargeTunePMuonBestTrack_, "mu_chargeTunePMuonBestTrack/i");
   tree->Branch("mu_dPToverPTTunePMuonBestTrack" , &mu_dPToverPTTunePMuonBestTrack_, "mu_dPToverPTTunePMuonBestTrack/i");
   tree->Branch("mu_absdxyTunePMuonBestTrack" , &mu_absdxyTunePMuonBestTrack_, "mu_absdxyTunePMuonBestTrack/i");
   tree->Branch("mu_absdxyErrorTunePMuonBestTrack" , &mu_absdxyErrorTunePMuonBestTrack_, "mu_absdxyErrorTunePMuonBestTrack/i");
   tree->Branch("mu_absdxySigTunePMuonBestTrack" , &mu_absdxySigTunePMuonBestTrack_, "mu_absdxySigTunePMuonBestTrack/i");
   tree->Branch("mu_absdzTunePMuonBestTrack" , &mu_absdzTunePMuonBestTrack_, "mu_absdzTunePMuonBestTrack/i");
   tree->Branch("mu_absdzErrorTunePMuonBestTrack" , &mu_absdzErrorTunePMuonBestTrack_, "mu_absdzErrorTunePMuonBestTrack/i");
   tree->Branch("mu_absdzSigTunePMuonBestTrack" , &mu_absdzSigTunePMuonBestTrack_, "mu_absdzSigTunePMuonBestTrack/i");
   tree->Branch("mu_recoDeltaBeta" , &mu_recoDeltaBeta_, "mu_recoDeltaBeta/i");
   tree->Branch("mu_recoiso" , &mu_recoiso_, "mu_recoiso/i");
   tree->Branch("mu_isGlobalMuon" , &mu_isGlobalMuon_, "mu_isGlobalMuon/i");
   tree->Branch("mu_isStandAloneMuon" , &mu_isStandAloneMuon_, "mu_isStandAloneMuon/i");
   tree->Branch("mu_isPF" , &mu_isPF_, "mu_isPF/i");
   tree->Branch("mu_isRPCMuon" , &mu_isRPCMuon_, "mu_isRPCMuon/i");
   tree->Branch("mu_isTrackerMuon" , &mu_isTrackerMuon_, "mu_isTrackerMuon/i");
   tree->Branch("mu_isGoodMuon" , &mu_isGoodMuon_, "mu_isGoodMuon/i");
   tree->Branch("mu_isSoftMuon" , &mu_isSoftMuon_, "mu_isSoftMuon/i");
   tree->Branch("mu_isLoose" , &mu_isLoose_, "mu_isLoose/i");
   tree->Branch("mu_isTightMuon" , &mu_isTightMuon_, "/i");
   tree->Branch("mu_STAnHits" , &mu_STAnHits_, "mu_STAnHits/i");
   tree->Branch("mu_STAnLost" , &mu_STAnLost_, "mu_STAnLost/i");
   tree->Branch("mu_STAnStationsWithAnyHits" , &mu_STAnStationsWithAnyHits_, "mu_STAnStationsWithAnyHits/i");
   tree->Branch("mu_STAnCscChambersWithAnyHits" , &mu_STAnCscChambersWithAnyHits_, "mu_STAnCscChambersWithAnyHits/i");
   tree->Branch("mu_STAnDtChambersWithAnyHits" , &mu_STAnDtChambersWithAnyHits_, "mu_STAnDtChambersWithAnyHits/i");
   tree->Branch("mu_STAnRpcChambersWithAnyHits" , &mu_STAnRpcChambersWithAnyHits_, "mu_STAnRpcChambersWithAnyHits/i");
   tree->Branch("mu_STAinnermostStationWithAnyHits" , &mu_STAinnermostStationWithAnyHits_, "mu_STAinnermostStationWithAnyHits/i");
   tree->Branch("mu_STAoutermostStationWithAnyHits" , &mu_STAoutermostStationWithAnyHits_, "mu_STAoutermostStationWithAnyHits/i");
   tree->Branch("mu_STAnStationsWithValidHits" , &mu_STAnStationsWithValidHits_, "mu_STAnStationsWithValidHits/i");
   tree->Branch("mu_STAnCscChambersWithValidHits" , &mu_STAnCscChambersWithValidHits_, "mu_STAnCscChambersWithValidHits/i");
   tree->Branch("mu_STAnDtChambersWithValidHit" , &mu_STAnDtChambersWithValidHit_, "mu_STAnDtChambersWithValidHit/i");
   tree->Branch("mu_STAnRpcChambersWithValidHits" , &mu_STAnRpcChambersWithValidHits_, "mu_STAnRpcChambersWithValidHits/i");
   tree->Branch("mu_STAnValidMuonHits" , &mu_STAnValidMuonHits_, "mu_STAnValidMuonHits/i");
   tree->Branch("mu_STAnValidCscHits" , &mu_STAnValidCscHits_, "mu_STAnValidCscHits/i");
   tree->Branch("mu_STAnValidDtHits" , &mu_STAnValidDtHits_, "mu_STAnValidDtHits/i");
   tree->Branch("mu_STAnValidRpcHits" , &mu_STAnValidRpcHits_, "mu_STAnValidRpcHits/i");
   tree->Branch("mu_STAinnermostStationWithValidHits" , &mu_STAinnermostStationWithValidHits_, "mu_STAinnermostStationWithValidHits/i");
   tree->Branch("mu_STAoutermostStationWithValidHits" , &mu_STAoutermostStationWithValidHits_, "mu_STAoutermostStationWithValidHits/i");
   tree->Branch("mu_STATofDirection" , &mu_STATofDirection_, "mu_STATofDirection/i");
   tree->Branch("mu_STATofNDof" , &mu_STATofNDof_, "mu_STATofNDof/i");
   tree->Branch("mu_STATofTimeAtIpInOut" , &mu_STATofTimeAtIpInOut_, "mu_STATofTimeAtIpInOut/i");
   tree->Branch("mu_STATofTimeAtIpInOutErr" , &mu_STATofTimeAtIpInOutErr_, "mu_STATofTimeAtIpInOutErr/i");
   tree->Branch("mu_STATofTimeAtIpOutIn" , &mu_STATofTimeAtIpOutIn_, "mu_STATofTimeAtIpOutIn/i");
   tree->Branch("mu_STATofTimeAtIpOutInErr" , &mu_STATofTimeAtIpOutInErr_, "mu_STATofTimeAtIpOutInErr/i");
   tree->Branch("mu_SecondGenMatch" , &mu_SecondGenMatch_, "mu_SecondGenMatch/i");
   tree->Branch("mu_FirstGenMatch" , &mu_FirstGenMatch_, "mu_FirstGenMatch/i");
   tree->Branch("" , &_, "/i");


 }



 void BigNtuple::fill_muInfo(TTree* tree){


 }


void BigNtuple::set_jetInfo(TTree* tree){

  tree->Branch("jet_nb" , &jet_nb_, "jet_nb/i");
  tree->Branch("jet_charge" , &jet_charge_, "jet_charge/i");
  tree->Branch("jet_et" , &jet_et_, "jet_et/i");
  tree->Branch("jet_pt" , &jet_pt_, "jet_pt/i");
  tree->Branch("jet_eta" , &jet_eta_, "jet_eta/i");
  tree->Branch("jet_phi" , &jet_phi_, "jet_phi/i");
  tree->Branch("jet_theta" , &jet_theta_, "jet_theta/i");
  tree->Branch("jet_en" , &jet_en_, "jet_en/i");
  tree->Branch("jet_chargedEmEnergy" , &jet_chargedEmEnergy_, "jet_chargedEmEnergy/i");
  tree->Branch("jet_neutralEmEnergyFraction" , &jet_neutralEmEnergyFraction_, "jet_neutralEmEnergyFraction/i");
  tree->Branch("jet_chargedHadronEnergy" , &jet_chargedHadronEnergy_, "jet_chargedHadronEnergy/i");
  tree->Branch("jet_neutralHadronEnergyFraction" , &jet_neutralHadronEnergyFraction_, "jet_neutralHadronEnergyFraction/i");
  tree->Branch("jet_chargedMuEnergy" , &jet_chargedMuEnergy_, "jet_chargedMuEnergy/i");
  tree->Branch("jet_chargedMuEnergyFraction" , &jet_chargedMuEnergyFraction_, "jet_chargedMuEnergyFraction/i");
  tree->Branch("jet_numberOfDaughters" , &jet_numberOfDaughters_, "jet_numberOfDaughters/i");
  tree->Branch("jet_muonEnergy" , &jet_muonEnergy_, "jet_muonEnergy/i");
  tree->Branch("jet_muonEnergyFraction" , &jet_muonEnergyFraction_, "jet_muonEnergyFraction/i");
  tree->Branch("jet_muonMultiplicity" , &jet_muonMultiplicity_, "jet_muonMultiplicity/i");
  tree->Branch("jet_neutralEmEnergy" , &jet_neutralEmEnergy_, "jet_neutralEmEnergy/i");
  tree->Branch("jet_neutralHadronEnergy" , &jet_neutralHadronEnergy_, "jet_neutralHadronEnergy/i");
  tree->Branch("jet_neutralHadronMultiplicity" , &jet_neutralHadronMultiplicity_, "jet_neutralHadronMultiplicity/i");
  tree->Branch("jet_neutralMultiplicity" , &jet_neutralMultiplicity_, "jet_neutralMultiplicity/i");
  tree->Branch("" , &_, "/i");
}



void BigNtuple::fill_jetInfo(const pat::JetCollection jets){
  int jetnumber=0;
  
  for (const pat::Jet JET : jets) {
    if( JET.pt() < 0.0 ) continue;
    if (!( fabs(JET.eta()) < 3 && JET.pt() > 5. )) continue;
    jetnumber++;

    jet_charge_.push_back(JET.charge());
    jet_et_.push_back(JET.et());
    jet_pt_.push_back(JET.pt());
    jet_eta_.push_back(JET.eta());
    jet_phi_.push_back(JET.phi());
    jet_theta_.push_back(JET.theta());
    jet_en_.push_back(JET.energy());
    jet_chargedEmEnergy_.push_back(JET.chargedEmEnergy());
    jet_neutralEmEnergyFraction_.push_back(JET.neutralEmEnergyFraction());
    jet_chargedHadronEnergy_.push_back(JET.chargedHadronEnergy());
    jet_neutralHadronEnergyFraction_.push_back(JET.neutralHadronEnergyFraction());
    jet_chargedMuEnergy_.push_back(JET.chargedMuEnergy());
    jet_chargedMuEnergyFraction_.push_back(JET.chargedMuEnergyFraction());
    jet_chargedMultiplicity_.push_back(JET.chargedMultiplicity());
    jet_numberOfDaughters_.push_back(JET.numberOfDaughters());
    jet_muonEnergy_.push_back(JET.muonEnergy());
    jet_muonEnergyFraction_.push_back(JET.muonEnergyFraction());
    jet_muonMultiplicity_.push_back(JET.muonMultiplicity());
    jet_neutralEmEnergy_.push_back(JET.neutralEmEnergy());
    jet_neutralHadronEnergy_.push_back(JET.neutralHadronEnergy());
    jet_neutralHadronMultiplicity_.push_back(JET.neutralHadronMultiplicity());
    jet_neutralMultiplicity_.push_back(JET.neutralMultiplicity());

  }

  jet_nb_ = jetnumber;

}




void BigNtuple::set_metInfo(TTree* tree){

  tree->Branch("pfMet_et" , &pfMet_et_, "pfMet_et/i");
  tree->Branch("pfMet_pt" , &pfMet_pt_, "pfMet_pt/i");
  tree->Branch("pfMet_phi" , &pfMet_phi_, "pfMet_phi/i");
  tree->Branch("pfMet_en" , &pfMet_en_, "pfMet_en/i");
  tree->Branch("pfMet_px" , &pfMet_px_, "pfMet_px/i");
  tree->Branch("pfMet_py" , &pfMet_py_, "pfMet_py/i");
  tree->Branch("pfMet_pz" , &pfMet_pz_, "pfMet_pz/i");
  tree->Branch("pfMet_sumEt" , &pfMet_sumEt_, "pfMet_sumEt/i");
  tree->Branch("caloMet_pt" , &caloMet_pt_, "caloMet_pt/i");
  tree->Branch("caloMet_phi" , &caloMet_phi_, "caloMet_phi/i");   
  
}



void BigNtuple::fill_metInfo(const pat::MET met){

  pfMet_et_ = met.et();
  pfMet_pt_ = met.pt();
  pfMet_phi_ = met.phi();
  pfMet_en_ = met.energy();
  pfMet_px_ = met.px();
  pfMet_py_ = met.py();
  pfMet_pz_ = met.pz();
  pfMet_sumEt_ = met.sumEt();

  caloMet_pt_ = met.caloMETPt();
  caloMet_phi_ = met.caloMETPhi();


}


