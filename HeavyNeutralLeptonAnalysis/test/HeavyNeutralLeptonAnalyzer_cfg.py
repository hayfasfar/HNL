from os import path as path
import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#import PhysicsTools.PythonAnalysis.LumiList as LumiList

#LumiList.LumiList().getVLuminosityBlockRange()
import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = ''

process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_Zpt-200toInf_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/04B4B947-A1E4-E611-BF2C-ECF4BBE16468.root'
#'file:/afs/cern.ch/work/h/hrejebsf/public/HLTReco/CMSSW_8_0_21/src/miniAOD.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
))
process.TFileService = cms.Service("TFileService", fileName = cms.string("HeavyNeutralLepton.root"))

#################################################
## Update PAT jets
#################################################

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

## b-tag discriminators
bTagDiscriminators = [
    'pfTrackCountingHighEffBJetTags',
    'pfTrackCountingHighPurBJetTags',
    'pfJetProbabilityBJetTags',
    'pfJetBProbabilityBJetTags',
    'pfSimpleSecondaryVertexHighEffBJetTags',
    'pfSimpleSecondaryVertexHighPurBJetTags',
    'pfCombinedSecondaryVertexV2BJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfCombinedMVAV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = bTagDiscriminators
)


process.HeavyNeutralLepton = cms.EDAnalyzer('HeavyNeutralLeptonAnalysis',
                                            isMC = cms.bool(True),
                                            vtxSrc                = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            rho                   = cms.InputTag("fixedGridRhoFastjetAll"),
                                            muonSrc               = cms.InputTag("slimmedMuons"),
                                            electronSrc           = cms.InputTag("slimmedElectrons"),
                                            recHitCollectionEBSrc = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                            recHitCollectionEESrc = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                            tauSrc                = cms.InputTag("slimmedTaus"),
                                            packCandSrc           = cms.InputTag("packedPFCandidates"),
                                            jetSrc                = cms.InputTag("slimmedJets"),
                                            pfMETSrc              = cms.InputTag("slimmedMETs"),
                                            triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
                                            metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
					    prescales             = cms.InputTag("patTrigger","","PAT"),
					    objects               = cms.InputTag("slimmedPatTrigger","","PAT"),
                                            genParticleSrc        = cms.InputTag("prunedGenParticles"),
                                            genEventInfoProduct   = cms.InputTag("generator"),
                                            PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
                                            lheEventProducts      = cms.InputTag("externalLHEProducer"),
				            bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
        				      'pfTrackCountingHighEffBJetTags',
 				              'pfTrackCountingHighPurBJetTags',
				              'pfJetProbabilityBJetTags',
				              'pfJetBProbabilityBJetTags',
				              'pfSimpleSecondaryVertexHighEffBJetTags',
				              'pfSimpleSecondaryVertexHighPurBJetTags',
   				              'pfCombinedSecondaryVertexV2BJetTags',
			                      'pfCombinedInclusiveSecondaryVertexV2BJetTags',
                                              'pfCombinedMVAV2BJetTags'
					    )				 
)

#process.HeavyNeutralLepton.outputFileName  = cms.untracked.string('HeavyNeutralLepton.root')
#process.HeavyNeutralLepton.genTreeName     = cms.untracked.string('Tree')


process.p = cms.Path(process.HeavyNeutralLepton)
