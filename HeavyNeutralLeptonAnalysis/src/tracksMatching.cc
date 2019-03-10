#include "HNL/HeavyNeutralLeptonAnalysis/interface/BigNtuple.h"
using namespace std ;
void  BigNtuple::set_sv_genInfo(TTree* tree){
//tree->Branch("partPosition",&partPosition_);
//tree->Branch("mothPosition",&mothPosition_);
tree->Branch("resolutionAllTracks", & resolutionTrack_);
tree->Branch("resolutionMatchedTracks", & reso_Matched_); 
tree->Branch("gen_chargedPrtMult",&chargedPrtMult_);
//tree->Branch("gen_PartMatchingDr", &trackDr) ; 
tree->Branch("gen_mu1Pt", & gen_mu1Pt_);
tree->Branch("gen_mu1Eta", & gen_mu1Eta_);
tree->Branch("gen_mu1Phi", & gen_mu1Phi_);
tree->Branch("gen_mu1dxy", & gen_mu1dxy_);
tree->Branch("gen_mu2Pt", & gen_mu2Pt_);
tree->Branch("gen_mu2Eta", & gen_mu2Eta_);
tree->Branch("gen_mu2Phi", & gen_mu2Phi_);
tree->Branch("gen_mu2dxy", & gen_mu2dxy_);
tree->Branch("gen_mu2dxyzPv", & gen_mu2dxyzPv_);
tree->Branch("gen_mu2dxyPv", & gen_mu2dxyPv_) ;
tree->Branch("gen_partId23",&gen_partId23_);
tree->Branch("gen_partPt",&gen_partPt_);
tree->Branch("gen_partPz", &gen_partPz_);
tree->Branch("gen_partEta", &gen_partEta_);
tree->Branch("gen_partPhi", &gen_partPhi_);
tree->Branch("gen_partVx" , &gen_partVx_) ; 
tree->Branch("gen_partVy", &gen_partVy_);
tree->Branch("gen_partVz", &gen_partVz_);
tree->Branch("gen_partDz", &gen_partDz_) ;
tree->Branch("gen_partdxy0", &gen_partdxy0_) ;
tree->Branch("gen_partdxyz0",  &gen_partdxyz0_) ;
tree->Branch("gen_partDxyPv", &gen_partDxy_);
tree->Branch("gen_partDxyzPv", &gen_partDxyz_);
tree->Branch("gen_partVertex", &gen_partVertex_) ; 
tree->Branch("genMuMult", &genMuMult_) ;
tree->Branch("gen_partMu2dR", &gen_partMu2dR_);
tree->Branch("gen_partMu1dR", &gen_partMu1dR_);
tree->Branch("gen_dRMu1Mu2", &dr_Mu1Mu2Gen_);
tree->Branch("gen_dEtaMu1Mu2", &gen_dEtaMu1Mu2_);
tree->Branch("gen_dPhiMu1Mu2", &gen_dPhiMu1Mu2_); 
tree->Branch("gen_dPzMu1Mu2", &gen_dPzMu1Mu2_);
tree->Branch("reco_nbOfTracks", &reco_nbOfTracks_);
tree->Branch("reco_TrackPt",  &reco_TrackPt_) ;
tree->Branch("deltaRMatchedTr", &deltaRMatchedTr_);
tree->Branch("reco_TrackEta",  &reco_TrackEta_);
tree->Branch("reco_TrackPhi",  &reco_TrackPhi_);
tree->Branch("reco_TrackVx", &reco_TrackVx_); 
tree->Branch("reco_TrackVy", &reco_TrackVy_);
tree->Branch("reco_TrackVz", &reco_TrackVz_);
tree->Branch("reco_TrackDxyPV", &reco_TrackDxy_);
tree->Branch("reco_TrackDzPv",&reco_TrackDzPv_);
tree->Branch("reco_TrackNbofhits",&reco_TrackNbofhits_);
tree->Branch("reco_TrackNbofPixelhits" , &reco_TrackNbofPixelhits_) ;
tree->Branch("reco_TrackDz", &reco_TrackDz_) ;
tree->Branch("reco_TrackHighPurity", &reco_TrackHighPurity_);
tree->Branch("dR_trackGen" , &dR_trackGen_ ) ;
tree->Branch("losttracks", &losttrack_);
tree->Branch("recoMuMult_", &recoMuMult_);
tree->Branch("reco_mu1Pt", &reco_mu1Pt_);
tree->Branch("reco_mu1Eta", &reco_mu1Eta_);
tree->Branch("reco_mu1dxy0", &reco_mu1dxy0_);
tree->Branch("reco_mu1Dxy", &reco_mu1Dxy_) ;
tree->Branch("reco_mu1vertex", &reco_mu1vertex_); 
tree->Branch("reco_mu2Pt", &reco_mu2Pt_);
tree->Branch("reco_mu2Eta", &reco_mu2Eta_);
tree->Branch("reco_mu2_vx", &reco_mu2_vx_);
tree->Branch("reco_mu2_vy",&reco_mu2_vy_);
tree->Branch("reco_mu2_vz", &reco_mu2_vz_);
tree->Branch("reco_mu2dxy0", &reco_mu2dxy0_);
tree->Branch("reco_mu2Dxy", &reco_mu2Dxy_) ;
tree->Branch("reco_mu2vertex", &reco_mu2vertex_); 
tree->Branch("itMatchesMu2", &itMatchesMu2_);
tree->Branch("itMatchesTrack", &itMatchesTrack_);
tree->Branch("gen_pv_position", &gen_pv_position_);
tree->Branch("gen_pv_x", &gen_pv_x_);
tree->Branch("gen_pv_y", &gen_pv_x_);
tree->Branch("gen_pv_z", &gen_pv_x_);


}

void BigNtuple::muonMCmatching(  edm::Handle<pat::PackedGenParticleCollection> genPackedHandle , edm::Handle<pat::MuonCollection> muonsHandle, float indexMu1 , float indexMu2 , math::XYZPoint gen_PvPosition  ){
nbRecoMu = muonsHandle->size() ;

for (size_t m = 0 ; m < muonsHandle->size() ; m ++ ){

     if((*muonsHandle)[m].bestTrack() != nullptr && fabs((*genPackedHandle)[indexMu1].eta()) < 2.4 && fabs((*genPackedHandle)[indexMu2].eta()) <2.4 ) {
        if((*genPackedHandle)[indexMu1].charge() == (*muonsHandle)[m].bestTrack()->charge() ){
        mu1_Matching = deltaR( (*genPackedHandle)[indexMu1].eta(),(*genPackedHandle)[indexMu1].phi(), (*muonsHandle)[m].bestTrack()->eta(), (*muonsHandle)[m].bestTrack()->phi());
       
        if(mu1_Matching < 0.1 && mu1_Matching < drminMu1 ){

            indexRecoMu1 = m ;
            drminMu1 = mu1_Matching ;
       	}
     } 
        if((*genPackedHandle)[indexMu1].charge() == (*muonsHandle)[m].bestTrack()->charge() ){
        mu2_Matching = deltaR( (*genPackedHandle)[indexMu2].eta(),(*genPackedHandle)[indexMu2].phi(), (*muonsHandle)[m].bestTrack()->eta(), (*muonsHandle)[m].bestTrack()->phi());
        if(mu2_Matching < 0.1 && mu2_Matching < drminMu2 ){
            indexRecoMu2 = m ;
            drminMu2 = mu2_Matching ; 
       	}
     
        }
     }  
        
}

recoMuMult_.push_back(nbRecoMu);

if(indexRecoMu1 != -10 ){

reco_mu1Pt_.push_back( (*muonsHandle)[indexRecoMu1].bestTrack()->pt() ) ;
reco_mu1Eta_.push_back((*muonsHandle)[indexRecoMu1].bestTrack()->eta()) ;
reco_mu1Dxy_.push_back((*muonsHandle)[indexRecoMu1].bestTrack()->dxy(gen_PvPosition)) ; 
float mu_x = (*muonsHandle)[indexRecoMu1].bestTrack()->vx()  ; 
float mu_y = (*muonsHandle)[indexRecoMu1].bestTrack()->vy() ;
float mu_z = (*muonsHandle)[indexRecoMu1].bestTrack()->vz() ; 
reco_mu1dxy0_.push_back(sqrt(mu_x *mu_x + mu_y*mu_y)) ; 
reco_mu1vertex_.push_back(sqrt(mu_x *mu_x + mu_y*mu_y + mu_z*mu_z)) ;


}
if(indexRecoMu2 == -10) itMatchesMu2_.push_back(0);
if(indexRecoMu2 != -10 ){
itMatchesMu2_.push_back(1);
reco_mu2Pt_.push_back((*muonsHandle)[indexRecoMu2].bestTrack()->pt());
reco_mu2Eta_.push_back((*muonsHandle)[indexRecoMu2].bestTrack()->eta()) ;
reco_mu2Dxy_.push_back((*muonsHandle)[indexRecoMu2].bestTrack()->dxy(gen_PvPosition));



float mu2_x = (*muonsHandle)[indexRecoMu2].bestTrack()->vx()  ; 
float mu2_y = (*muonsHandle)[indexRecoMu2].bestTrack()->vy() ;
float mu2_z = (*muonsHandle)[indexRecoMu2].bestTrack()->vz() ;
reco_mu2_vx_.push_back(mu2_x) ;
reco_mu2_vy_.push_back(mu2_y) ;
reco_mu2_vz_.push_back(mu2_z) ;   
reco_mu2dxy0_.push_back(sqrt(mu2_x*mu2_x + mu2_y*mu2_y )) ; 
reco_mu2vertex_.push_back(sqrt(mu2_x*mu2_x + mu2_y*mu2_y + mu2_z*mu2_z)) ; 

}

}

void  BigNtuple::fill_sv_genInfo(edm::Handle<reco::GenParticleCollection>  genHandle , edm::Handle<pat::PackedGenParticleCollection>  genPackedHandle , edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle , edm::Handle<pat::PackedCandidateCollection> lostTracks , float * indexMu1 , float* indexMu2 , math::XYZPoint * gen_PvPosition ){
std::vector<float> matched ;

 cout<<"first point " <<endl ;
 for (size_t j = 0 ; j <  genPackedHandle->size() ; j++ ){
   if(fabs((*genPackedHandle)[j].pdgId()) == 13 ) nbGenMu +=1 ;  
   if(fabs((*genPackedHandle)[j].pdgId()) == 13 && (*genPackedHandle)[j].mother(0) != nullptr)
   {
    // nbGenMu +=1 ;
     bool  motherFound = false ; 
     const reco::Candidate * mother  = (*genPackedHandle)[j].mother(0) ;
     while(motherFound == false) {
       
       if(mother == nullptr) break ;

       if((fabs(mother->pdgId())== 9900012 || fabs(mother->pdgId()) == 24) ) { 
         motherFound = true ;
       }    
       if (motherFound == false ) mother = mother->mother(0) ; 
     }

    if(mother != nullptr ){    
    if( fabs(mother->pdgId()) == 24){
     *indexMu1 = j ;
    }

    if(mother->pdgId() == 9900012 ){
     *indexMu2 = j ; 
    }
    }   
    
   if(*indexMu1 != -10 && *indexMu2 != -10 ) {
       gen_mu1Pt_.push_back((*genPackedHandle)[*indexMu1].pt()); 
       gen_mu1Eta_.push_back((*genPackedHandle)[*indexMu1].eta()); 
       gen_mu1Phi_.push_back((*genPackedHandle)[*indexMu1].phi()); 
       gen_mu2Pt_.push_back((*genPackedHandle)[*indexMu2].pt()) ;
       gen_mu2Eta_.push_back((*genPackedHandle)[*indexMu2].eta()); 
       gen_mu2Phi_.push_back((*genPackedHandle)[*indexMu2].phi()); 

       float dr_genMu1Mu2 = deltaR((*genPackedHandle)[*indexMu1].eta(),(*genPackedHandle)[*indexMu1].phi(), (*genPackedHandle)[*indexMu2].eta(),(*genPackedHandle)[*indexMu2].phi());


       gen_dPhiMu1Mu2_.push_back(fabs((*genPackedHandle)[*indexMu1].phi() - (*genPackedHandle)[*indexMu2].phi()));
       gen_dEtaMu1Mu2_.push_back(fabs((*genPackedHandle)[*indexMu1].eta() -(*genPackedHandle)[*indexMu2].eta())); 
       gen_dPzMu1Mu2_.push_back(fabs((*genPackedHandle)[*indexMu1].pz() - (*genPackedHandle)[*indexMu2].pz()));
       dr_Mu1Mu2Gen_.push_back(dr_genMu1Mu2); 
       break ; 
   }
  }
 }
    genMuMult_.push_back(nbGenMu);

   // Looking for the Gen final state particles from Signal hadronisation and match them to reco tracks.
  float prt_x =0., prt_y=0. , prt_z =0. ; 
  float deltaDxy = 0. , deltaDxyz = 0. ;

  if(*indexMu2 != -10 && *indexMu1 != -10) { // index mu1 and index mu2 should be always different than -10, but this is only for safety.


   for (size_t j = 0 ; j < genPackedHandle->size() ; j++ ){


      if(j == *indexMu2 || j == *indexMu1)  continue ;       

     if(fabs((*genPackedHandle)[j].pdgId() ) > 100 ) {
 //      if((*genPackedHandle)[j].charge()!= 0){

     for (size_t i = 0 ; i< genHandle->size(); i++){

        if(fabs((*genHandle)[i].pdgId()) == 24 && (*genHandle)[i].isLastCopy() && ! pvFound){
                  *gen_PvPosition =  (*genHandle)[i].vertex() ;
		   gen_pv_x = (*genHandle)[i].vx() ; gen_pv_y = (*genHandle)[i].vy() ; gen_pv_z = (*genHandle)[i].vz() ;
                   pvFound = true ;     
                   gen_pv_x_.push_back(gen_pv_x);  
                   gen_pv_y_.push_back(gen_pv_y);  
                   gen_pv_z_.push_back(gen_pv_z);  
        }
        if((*genHandle)[i].mother(0) != nullptr && (*genHandle)[i].mother(0)->pdgId() == 9900012 && fabs((*genHandle)[i].status() ) == 23 && fabs((*genHandle)[i].pdgId()) < 5){
          const reco::Candidate * HNL =  (*genHandle)[i].mother(0) ; 
          const reco::Candidate * motherInPrunedCollection = (*genPackedHandle)[j].mother(0);
          if(motherInPrunedCollection != nullptr && isAncestor( HNL , motherInPrunedCollection)) {
                      if((*genPackedHandle)[j].charge() != 0 ){
                         chParticles += 1;
                         gen_partId23_.push_back((*genPackedHandle)[j].pdgId());
			 gen_partPt_.push_back((*genPackedHandle)[j].pt());
                         gen_partPz_.push_back((*genPackedHandle)[j].pz()) ;
                         gen_partEta_.push_back((*genPackedHandle)[j].eta());
			 gen_partPhi_.push_back((*genPackedHandle)[j].phi());
			 gen_partDz_.push_back((*genPackedHandle)[j].dz(*gen_PvPosition));
                         prt_x = (*genHandle)[i].vx() ;
                         prt_y = (*genHandle)[i].vy() ;
                         prt_z = (*genHandle)[i].vz() ;
			 gen_partVx_.push_back(prt_x) ;  
			 gen_partVy_.push_back(prt_y) ;  
			 gen_partVz_.push_back(prt_z) ;
 			 gen_partdxy0_.push_back(sqrt(prt_x*prt_x + prt_y*prt_y)); 
 			 gen_partdxyz0_.push_back(sqrt(prt_x*prt_x + prt_y*prt_y + prt_z*prt_z));
                         deltaDxy = sqrt(pow(prt_x - gen_pv_x,2) + pow(prt_y - gen_pv_y,2) ) ;
                         deltaDxyz = sqrt(pow(prt_x - gen_pv_x,2) + pow(prt_y - gen_pv_y,2) + pow(prt_z - gen_pv_z,2));
 			 gen_partDxy_.push_back(deltaDxy);
 			 gen_partDxyz_.push_back(deltaDxyz);
                       
                         float gendRmu1part = deltaR((*genPackedHandle)[j].eta(), (*genPackedHandle)[j].phi(), (*genPackedHandle)[*indexMu1].eta(), (*genPackedHandle)[*indexMu1].phi());
                         float gendRmu2part = deltaR((*genPackedHandle)[j].eta(), (*genPackedHandle)[j].phi(), (*genPackedHandle)[*indexMu2].eta(), (*genPackedHandle)[*indexMu2].phi()); 

                            gen_partMu2dR_.push_back(gendRmu2part);
                            gen_partMu1dR_.push_back(gendRmu1part);
// here we do match the tracks :
//   		            
			    int indexRecoTrack = -10 ; 
			    double drmintrack = 10. ; 
			    int losttrack = -10 ; 
			    double resolution = 10000;
 
 			    for (size_t k = 0 ; k < pfCandidatesHandle->size() ; k++ )
                            {
		//this to avoid double counting.
 				      for(size_t m = 0 ; m < matched.size() ; m ++ ) { 
                                         if(matched.size() == 0)  continue ; 
					 if(matched[m] == k && losttrack_[m] == 0 ) continue ; 
                                      }
				//if there is best track  i take it, if not i just access directly  eta and phi of the pf Candidates.
				//
                                   if((*pfCandidatesHandle)[k].bestTrack() == nullptr){ //cout << "yes it is nullptr"<<endl;

                                       float track_Matching = deltaR((*pfCandidatesHandle)[k].eta(),(*pfCandidatesHandle)[k].phi(),(*genPackedHandle)[j].eta(),(*genPackedHandle)[j].phi());


                                       if(track_Matching < drmintrack && track_Matching < 0.15 && (*pfCandidatesHandle)[k].charge() == (*genPackedHandle)[j].charge() ){ 

  			                 	drmintrack = track_Matching ;
 			                	indexRecoTrack = k ;
  		 				dR_trackGen_.push_back(drmintrack ) ;
						losttrack = 0 ;
 						 
                                       }
                                   }
                   
                                   if((*pfCandidatesHandle)[k].bestTrack()!= nullptr ){
                                     float track_Matching = deltaR((*pfCandidatesHandle)[k].bestTrack()->eta(),(*pfCandidatesHandle)[k].bestTrack()->phi(),(*genPackedHandle)[j].eta(),(*genPackedHandle)[j].phi());
 
                                   if(track_Matching < drmintrack && track_Matching < 0.15 && (*pfCandidatesHandle)[k].bestTrack()->charge() == (*genPackedHandle)[j].charge() ){					   
  			                 	drmintrack = track_Matching ;
						cout<< "dr min is "<<drmintrack <<endl ; 
 			                	indexRecoTrack = k ;
     						dR_trackGen_.push_back(drmintrack ) ; 
						losttrack = 0 ;	
                                      }
 			//		cout<<"point 7" <<endl ;     
                                 }    
                          }
			if(drmintrack > 0.02){
		       	cout<<"it enter the lost loop "<<endl ;

                          for(size_t l = 0 ; l < lostTracks->size() ; l++){
		//this to avoid double counting.
                           for(size_t m = 0 ; m < matched.size() ; m ++ ) {
                                         if(matched.size() == 0)  continue ;
  					 if(matched[m] == l &&  losttrack_[m] == 1 ) continue ; 

                           }
				if((*lostTracks)[l].bestTrack() != nullptr ){
                                float track_Matching = deltaR((*lostTracks)[l].bestTrack()->eta(),(*lostTracks)[l].bestTrack()->phi(),(*genPackedHandle)[j].eta(),(*genPackedHandle)[j].phi());
					  if(track_Matching < drmintrack &&  track_Matching < 0.15 && (*lostTracks)[l].charge() == (*genPackedHandle)[j].charge()){

                                                drmintrack = track_Matching ;

                                                indexRecoTrack = l ;
                                                dR_trackGen_.push_back(drmintrack ) ;
                                                losttrack = 1 ;
                                }

			        }
  
			        if((*lostTracks)[l].bestTrack() == nullptr ){

                               float track_Matching = deltaR((*lostTracks)[l].eta(),(*lostTracks)[l].phi(),(*genPackedHandle)[j].eta(),(*genPackedHandle)[j].phi());
				 if(track_Matching < drmintrack && track_Matching < 0.15 && (*lostTracks)[l].charge() == (*genPackedHandle)[j].charge()){
					indexRecoTrack = l ;
					dR_trackGen_.push_back(drmintrack ) ;
					losttrack = 1 ;
				 }
   			          	
				}

                              //Do the  same but with  == nullptr  !!!! 
			   }
			 }

//   cout<< "losttrack "<< losttrack << "resolution "<< resolution <<"ndexRecoTrack "<<indexRecoTrack << endl ; 
//   cout<<"next gen particle : "<<endl ;
// here you add things about the reconstructed tracks from pfCandidates. 

     if(indexRecoTrack == -10) itMatchesTrack_.push_back(false);
     if(losttrack == 0 ){
     if(indexRecoTrack != -10 && (*pfCandidatesHandle)[indexRecoTrack].bestTrack() ==nullptr ){
     itMatchesTrack_.push_back(true);
     matched.push_back(indexRecoTrack);
     deltaRMatchedTr_.push_back(drmintrack); 
     losttrack_.push_back(losttrack); 
     reco_TrackHighPurity_.push_back((*pfCandidatesHandle)[indexRecoTrack].trackHighPurity() ) ; 
     reco_TrackEta_.push_back((*pfCandidatesHandle)[indexRecoTrack].eta());
     reco_TrackPhi_.push_back((*pfCandidatesHandle)[indexRecoTrack].phi());
     reco_TrackPt_.push_back((*pfCandidatesHandle)[indexRecoTrack].pt());
     reco_TrackVx_.push_back((*pfCandidatesHandle)[indexRecoTrack].vx());
     reco_TrackVy_.push_back((*pfCandidatesHandle)[indexRecoTrack].vy());
     reco_TrackVz_.push_back((*pfCandidatesHandle)[indexRecoTrack].vz());
     reco_TrackDxy_.push_back((*pfCandidatesHandle)[indexRecoTrack].dxy(*gen_PvPosition));
     reco_TrackDxy0_.push_back((*pfCandidatesHandle)[indexRecoTrack].dxy());
     reco_TrackNbofhits_.push_back((*pfCandidatesHandle)[indexRecoTrack].numberOfHits());
     reco_TrackNbofPixelhits_.push_back((*pfCandidatesHandle)[indexRecoTrack].numberOfPixelHits () ) ; 
  
     reco_TrackDz_.push_back((*pfCandidatesHandle)[indexRecoTrack].dz());
     reco_TrackDzPv_.push_back((*pfCandidatesHandle)[indexRecoTrack].dz(*gen_PvPosition));
     
     nbOftracks += 1 ;
    }   
    if(indexRecoTrack != -10 && (*pfCandidatesHandle)[indexRecoTrack].bestTrack() !=nullptr){
     matched.push_back(indexRecoTrack);
     losttrack_.push_back(losttrack); 
     deltaRMatchedTr_.push_back(drmintrack); 
     itMatchesTrack_.push_back(true);
     reco_TrackHighPurity_.push_back((*pfCandidatesHandle)[indexRecoTrack].trackHighPurity() ) ; 
     reco_TrackEta_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->eta());
     reco_TrackPhi_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->phi());
     reco_TrackPt_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->pt());
     reco_TrackVx_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->vx());
     reco_TrackVy_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->vy());
     reco_TrackVz_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->vz());
     reco_TrackDxy_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->dxy(*gen_PvPosition));
     reco_TrackDxy0_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->dxy());
     reco_TrackDz_.push_back((*pfCandidatesHandle)[indexRecoTrack].bestTrack()->dz());
     reco_TrackDzPv_.push_back((*pfCandidatesHandle)[indexRecoTrack].dz(*gen_PvPosition));
     reco_TrackNbofhits_.push_back((*pfCandidatesHandle)[indexRecoTrack].numberOfHits());
     reco_TrackNbofPixelhits_.push_back((*pfCandidatesHandle)[indexRecoTrack].numberOfPixelHits () ) ; 
     nbOftracks += 1 ; 
    }              
    }
   if(losttrack == 1 )	{
   if(indexRecoTrack != -10 && (*lostTracks)[indexRecoTrack].bestTrack() ==nullptr){
     itMatchesTrack_.push_back(true);
  
     matched.push_back(indexRecoTrack);
     deltaRMatchedTr_.push_back(drmintrack); 
    

     losttrack_.push_back(losttrack); 
     reco_TrackHighPurity_.push_back((*lostTracks)[indexRecoTrack].trackHighPurity() ) ; 

     reco_TrackEta_.push_back((*lostTracks)[indexRecoTrack].eta());
     reco_TrackPhi_.push_back((*lostTracks)[indexRecoTrack].phi());
     reco_TrackPt_.push_back((*lostTracks)[indexRecoTrack].pt());
     reco_TrackVx_.push_back((*lostTracks)[indexRecoTrack].vx());
     reco_TrackVy_.push_back((*lostTracks)[indexRecoTrack].vy());
     reco_TrackVz_.push_back((*lostTracks)[indexRecoTrack].vz());
     reco_TrackDxy_.push_back((*lostTracks)[indexRecoTrack].dxy(*gen_PvPosition));
     reco_TrackDxy0_.push_back((*lostTracks)[indexRecoTrack].dxy());
     reco_TrackNbofhits_.push_back((*lostTracks)[indexRecoTrack].numberOfHits());
     reco_TrackNbofPixelhits_.push_back((*pfCandidatesHandle)[indexRecoTrack].numberOfPixelHits () ) ; 
     reco_TrackDz_.push_back((*lostTracks)[indexRecoTrack].dz());
     reco_TrackDzPv_.push_back((*lostTracks)[indexRecoTrack].dz(*gen_PvPosition));
     
     nbOftracks += 1 ;
    }   
    if(indexRecoTrack != -10 && (*lostTracks)[indexRecoTrack].bestTrack() !=nullptr ){
     matched.push_back(indexRecoTrack);
     losttrack_.push_back(losttrack); 
     reco_TrackHighPurity_.push_back((*lostTracks)[indexRecoTrack].trackHighPurity() ) ; 
     deltaRMatchedTr_.push_back(drmintrack); 
     itMatchesTrack_.push_back(true);
     reco_TrackEta_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->eta());
     reco_TrackPhi_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->phi());
     reco_TrackPt_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->pt());
     reco_TrackVx_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->vx());
     reco_TrackVy_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->vy());
     reco_TrackVz_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->vz());
     reco_TrackDxy_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->dxy(*gen_PvPosition));
     reco_TrackDxy0_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->dxy());
     reco_TrackNbofhits_.push_back((*lostTracks)[indexRecoTrack].numberOfHits());
     reco_TrackNbofPixelhits_.push_back((*lostTracks)[indexRecoTrack].numberOfPixelHits () ) ; 
     reco_TrackDz_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->dz());
     reco_TrackDzPv_.push_back((*lostTracks)[indexRecoTrack].bestTrack()->dz(*gen_PvPosition));
     nbOftracks += 1 ; 
    } 

    }

          }
        }
break;
     } 
   }
  }
 }
//Here you put push_back  reco_Mu1/2 variables and signal tracks variables. 
if(*indexMu2 != -10 && *indexMu2 != -10)
chargedPrtMult_.push_back(chParticles);                  


//if(nbOftracks!=0)
reco_nbOfTracks_.push_back(nbOftracks);
cout<<"last point " <<endl ; 
 }
 }


bool BigNtuple::isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle)
{ 
  if(ancestor == particle ) return true;
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
     if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  return false;
}

