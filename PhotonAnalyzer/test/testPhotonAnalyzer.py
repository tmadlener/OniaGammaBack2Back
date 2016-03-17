import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonAnalyzerTest")

process.load('FWCore.MessageService.MessageLogger_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:myOutputFile.root'
                                # 'file:skimming.root'
                            )
)

process.photonAnalysis = cms.EDAnalyzer('PhotonAnalyzerPAT',
                        outputFileName = cms.string('photonAnalysis.root'),
                        pFlowPhotons = cms.InputTag('OGB2BPhotonProducer:PFlowPhotons'),
                        photons = cms.InputTag('OGB2BPhotonProducer:photons'),
                        conversions = cms.InputTag('OGB2BPhotonProducer:convertedPhotons'),
                        dimuons = cms.InputTag('onia2MuMuPAT'),
                        muons = cms.InputTag('REPLACEME'),
                        primaryVertices = cms.InputTag('offlinePrimaryVertices'),
                        TriggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
                        pi0WideWindow = cms.vdouble(0.110, 0.160),
                        pi0NarrowWindow = cms.vdouble(0.130, 0.140),
                        oniaMassMax = cms.double(4.0),
                        oniaMassMin = cms.double(2.2),
                        storeOnlyBest = cms.bool(False), # store only the best dimuon candidate in each event
)

process.p = cms.Path(process.photonAnalysis)
