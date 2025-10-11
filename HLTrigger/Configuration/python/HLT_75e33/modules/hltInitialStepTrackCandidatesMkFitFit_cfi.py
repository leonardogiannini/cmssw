import FWCore.ParameterSet.Config as cms

# InitialStepTrackCandidatesMkFit options

hltInitialStepTrackCandidatesMkFitFit = cms.EDProducer("MkFitFitProducer",
   config = cms.ESInputTag("","hltInitialStepTrackCandidatesMkFitConfig"),
   eventOfHits = cms.InputTag("hltMkFitEventOfHits"),
   limitConcurrency = cms.untracked.bool(True),
   mightGet = cms.optional.untracked.vstring,
   mkFitSilent = cms.untracked.bool(True),
   tracks = cms.InputTag("hltInitialStepTrackCandidatesMkFit"),
   mkFitPixelHits = cms.InputTag("hltMkFitSiPixelHits"),
   pixelCPE = cms.string('PixelCPEGeneric'),
   pixelHits = cms.InputTag("hltMkFitSiPixelHits"),
   stripHits = cms.InputTag("hltMkFitSiPhase2Hits")
)
