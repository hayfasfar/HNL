import sys, os
#from os import path as path
import FWCore.ParameterSet.Config as cms

#function to return JSON file
def getJSON(is2017):
    if is2017: return "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
    else: return "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"


debugLevel    = -1 

isMC_         = False
isMCSignal_   = False
isData_       = True
is2017        = False 

algorithm     = "AK4PFchs"

GT_MC   = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
GT_DATA = '80X_dataRun2_2016SeptRepro_v7'
GT      =  GT_MC if isMC_ else GT_DATA

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1
if isData_:
 import FWCore.PythonUtilities.LumiList as LumiList
 LumiList.LumiList(filename = "../data/JSON/" + getJSON(is2017)).getVLuminosityBlockRange()

 #from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/50000/EE9CC3E0-0DED-E711-BCAC-00E081CB560C.root'
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-2_V-0.00316227766017_mu_massiveAndCKM_LO/heavyNeutrino_40.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_mu_massiveAndCKM_LO/heavyNeutrino_96.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-1_V-0.00836660026534_e_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_96.root'
))


outputFile      = 'data.root'     # trilep    --> skim three leptons (basic pt/eta criteria)
                                 # dilep     --> skim two leptons
                                 # singlelep --> skim one lepton
                                 # ttg       --> skim two leptons + one photon
                                 # fakerate  --> not implemented

def getVal(arg):
    return arg.split('=')[-1]

# Loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i]
    if "outputFile"  in sys.argv[i]: outputFile = getVal(sys.argv[i])
    elif "inputFile" in sys.argv[i]: inputFile  = getVal(sys.argv[i])
    elif "events" in sys.argv[i]: nEvents = int(getVal(sys.argv[i]))



process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))
process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')

from HNL.HeavyNeutralLeptonAnalysis.ele_Sequence_cff import addElectronSequence

addElectronSequence(process)

process.get = cms.EDAnalyzer("EventSetupRecordDataGetter",
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        data = cms.vstring('JetCorrectorParametersCollection/AK4PFchs')
    )
    ),
    verbose = cms.untracked.bool(True)
)


process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute',
              'L2L3Residual'],
    payload = 'AK4PFchs' ) 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )


#process.ak4PFJetsCHSL1FastL2L3 = process.ak4PFCHSJetsL1.clone(correctors = ['ak4PFL1FastL2L3Corrector'])
#process.correctJets           = cms.Sequence(process.ak4PFL1FastL2L3CorrectorChain * process.ak4PFJetsCHSL1FastL2L3)

process.HeavyNeutralLepton = cms.EDAnalyzer('HeavyNeutralLeptonAnalysis',
                                            debugLevel            = cms.int32(debugLevel),
                                            isMC                  = cms.bool(isMC_),
                                            isMCSignal            = cms.bool(isMCSignal_),
                                            vtxSrc                = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            rho                   = cms.InputTag("fixedGridRhoFastjetAll"),
                                            muonSrc               = cms.InputTag("slimmedMuons"),
                                            electronSrc           = cms.InputTag("slimmedElectrons"),
                                            recHitCollectionEBSrc = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                            recHitCollectionEESrc = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                            tauSrc                = cms.InputTag("slimmedTaus"),
                                            packCandSrc           = cms.InputTag("packedPFCandidates"),
                                            jetSrc                = cms.InputTag("slimmedJetsJEC"),
                                            pfMETSrc              = cms.InputTag("slimmedMETs"),
                                            triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
                                            metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
                                            genParticleSrc        = cms.InputTag("prunedGenParticles"),
                                            genEventInfoProduct   = cms.InputTag("generator"),
                                            PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
                                            lheEventProducts      = cms.InputTag("externalLHEProducer"),
                                            SecondaryVertices     = cms.InputTag("displacedInclusiveSecondaryVertices"),
                                            bDiscriminators       = cms.vstring("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                            electronsMva          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                            electronsVeto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                            electronsLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                            electronsMedium= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                            electronsTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
                                            )

process.p = cms.Path(
    process.displacedInclusiveVertexing 
    *process.ele_Sequence
    *process.jetCorrFactors
    *process.slimmedJetsJEC
    *process.HeavyNeutralLepton
    )
