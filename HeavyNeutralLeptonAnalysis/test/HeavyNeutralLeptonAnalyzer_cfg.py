from os import path as path
import FWCore.ParameterSet.Config as cms

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

import PhysicsTools.PythonAnalysis.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = ''

process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_105.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_106.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_108.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_109.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_110.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_111.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_113.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_115.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_116.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_117.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_118.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_119.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_120.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_121.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_124.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_129.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_131.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_135.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_136.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_137.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_138.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_139.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_140.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_141.root',
                                                               'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-2_V-0.00836660026534_mu_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_142.root') )

process.TFileService = cms.Service("TFileService", fileName = "HeavyNeutralLepton.root")

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
                                            genParticleSrc        = cms.InputTag("prunedGenParticles"),
                                            genEventInfoProduct   = cms.InputTag("generator"),
                                            PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
                                            lheEventProducts      = cms.InputTag("externalLHEProducer"),
                                            )

process.p = cms.Path(process.HeavyNeutralLepton)
