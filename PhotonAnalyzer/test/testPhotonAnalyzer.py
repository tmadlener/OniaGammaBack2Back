import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonAnalyzerTest")

process.load('FWCore.MessageService.MessageLogger_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:myOutputFile.root'
                            )
)

process.photonAnalysis = cms.EDAnalyzer('PhotonAnalyzerPAT',
                        outputFileName = cms.string('photonAnalysis.root'),
                        pFlowPhotons = cms.InputTag('OGB2BPhotonProducer:PFlowPhotons'),
                        photons = cms.InputTag('OGB2BPhotonProducer:photons'),
                        conversions = cms.InputTag('OGB2BPhotonProducer:convertedPhotons'),
                        pi0WideWindow = cms.vdouble(0.110, 0.160),
                        pi0NarrowWindow = cms.vdouble(0.130, 0.140),
)

process.p = cms.Path(process.photonAnalysis)
