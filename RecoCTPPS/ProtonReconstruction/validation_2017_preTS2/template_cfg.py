import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("CTPPSTestProtonReconstruction", eras.ctpps_2016)

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
  statistics = cms.untracked.vstring(),
  destinations = cms.untracked.vstring('cout'),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string('WARNING')
  )
)

# raw data source
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
  #lumisToProcess = cms.untracked.VLuminosityBlockRange("$run:1-$run:max")
)
$input

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

# provide LHCInfo
process.load("RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsLHCInfoESSourceJSON_cfi")

# proton reconstruction
from RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsProtonReconstructionOF_cfi import *
process.ctppsProtonReconstructionOFDB = ctppsProtonReconstructionOFDB
process.ctppsProtonReconstructionOFDB.alignmentFiles = cms.vstring("RecoCTPPS/ProtonReconstruction/data/alignment/2017_preTS2/collect_alignments_$alignment.out")

# reconstruction plotter
process.ctppsProtonReconstructionPlotter = cms.EDAnalyzer("CTPPSProtonReconstructionPlotter",
    tagTracks = cms.InputTag("ctppsLocalTrackLiteProducer"),
    tagRecoProtons = cms.InputTag("ctppsProtonReconstructionOFDB"),

    rpId_45_F = cms.uint32(23),
    rpId_45_N = cms.uint32(3),
    rpId_56_N = cms.uint32(103),
    rpId_56_F = cms.uint32(123),

    outputFile = cms.string("$output"),
    maxNonEmptyEvents = cms.untracked.int32(100000)
)

process.p = cms.Path(
    process.ctppsProtonReconstructionOFDB
    * process.ctppsProtonReconstructionPlotter
)
