import FWCore.ParameterSet.Config as cms

process = cms.Process("ConversionPhotonsProdTest")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/scratch/tmadlener/data/Run2015D/Charmonium/AOD/PromptReco-v3/000/256/673/00000/04BC2F52-F55E-E511-BB01-02163E0143A2.root'
        # 'root://xrootd-cms.infn.it//store/data/Run2015D/Charmonium/AOD/PromptReco-v3/000/256/673/00000/04BC2F52-F55E-E511-BB01-02163E0143A2.root'
    )
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# make patCandidates, select and clean them
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')
process.patMuons.embedTrack  = True

process.selectedPatMuons.cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    )


#make patTracks
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.7',                    # preselection cut on candidates
                    selection    = 'pt > 0.7',                    # selection on PAT Layer 1 objects
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None                           # replicate MC match as the one used for Muons
)
process.patTrackCands.embedTrack = True

###############################################
# HERE STARTS MY PART                         #
###############################################
process.OGB2BPhotonProducer = cms.EDProducer('ConversionPhotonProducer',
                conversions = cms.InputTag('allConversions'),
                allPhotons = cms.InputTag('photons'),
                pfcandidates = cms.InputTag('particleFlow'),
                beamspot = cms.InputTag('offlineBeamSpot'),
                primaryVtxTag = cms.InputTag('offlinePrimaryVertices'),
                convSelection = cms.string('conversionVertex.position.rho > 1.5'), # hardcoded at the moment
                convAlgorithm = cms.string('undefined'),
                convQuality = cms.vstring(['highPurity', 'generalTracksOnly']),
                tkVtxCompSigma = cms.double(5.0),
                pfCandSelection = cms.string('pt > 4'),
                photonSelection = cms.string(''),
                vertexChi2ProbCut = cms.double(0.0005),
                trackChi2Cut = cms.double(10),
                trackMinNDOF = cms.double(3.),
                minDistanceOfApproachMaxCut = cms.double(1.00),
                minDistanceOfApproachMinCut = cms.double(-0.25),
                pi0NarrowWindow = cms.vdouble(0.130, 0.140),
                pi0WideWindow = cms.vdouble(0.110, 0.160),
                convFlags = cms.string('10101'), # NOTE: shortened here! internally will be expanded to the right length
                photonFlags = cms.string('0'), # Only one flag at the moment, which is tied to the cut selection
                pFlowFlags = cms.string('0'), # Only one flag at the moment, which is tied to the cut selection
)

process.p = cms.Path(process.OGB2BPhotonProducer)
################################
# AND HERE IT ENDS (BASICALLY) #
################################

process.load('HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi')
process.onia2MuMuPAT.muons = cms.InputTag('cleanPatMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')

process.onia2MuMuPATCounter = cms.EDFilter('CandViewCountFilter',
      src = cms.InputTag('OGB2BPhotonProducer', 'convertedPhotons'),
      minNumber = cms.uint32(0),
      filter = cms.bool(True)
   )

# reduce MC genParticles a la miniAOD
process.load('PhysicsTools.PatAlgos.slimming.genParticles_cff')
process.packedGenParticles.inputVertices = cms.InputTag('offlinePrimaryVertices')

# make photon candidate conversions for P-wave studies
process.load('HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi')

# copied from CompactSkim_cff.py
SlimmedEventContent = [
    'keep recoVertexs_offlinePrimaryVertices_*_*',
    'keep *_inclusiveSecondaryvertices_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep *_gtDigis_*_RECO',
    'keep *_cleanPatTrackCands_*_*',
    'keep *_onia2MuMuPAT_*_*',
    'keep *_generalV0Candidates_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep *_OGB2BPhotonProducer_*_*',
    'keep *_PhotonCandidates_*_*',
]

# copied from CompactSkim_cff.py
from PhysicsTools.PatAlgos.tools.coreTools import runOnData
runOnData(process, outputModules = [] )

process.FilterOutput = cms.Path(process.onia2MuMuPATCounter)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root'),
    outputCommands = cms.untracked.vstring('drop *', *SlimmedEventContent),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('FilterOutput')),
    # SelectEvents = cms.untracked.PSet(),
)

process.e = cms.EndPath(process.out)
