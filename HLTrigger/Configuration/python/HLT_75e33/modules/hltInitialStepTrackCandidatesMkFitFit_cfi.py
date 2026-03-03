import FWCore.ParameterSet.Config as cms

hltInitialStepTrackCandidatesMkFitFit = cms.EDProducer("MkFitFitProducer",
    eventOfHits = cms.InputTag("hltMkFitEventOfHits"),
    config = cms.ESInputTag("","hltInitialStepTrackCandidatesMkFitConfig"),
    pixelCPE = cms.string('PixelCPEGeneric'),
    mkFitPixelHits = cms.InputTag("hltMkFitSiPixelHits"),
    tracks = cms.InputTag("hltInitialStepTrackCandidatesMkFit"),
    mkFitSilent = cms.untracked.bool(True),
    limitConcurrency = cms.untracked.bool(False),
    candCutSel = cms.bool(True),
    candMinPtCut = cms.double(0.8),
    candMinPtCutOuter = cms.double(0),
    candEtaRegion = cms.int32(0),
    candMinNHitsCut = cms.int32(3),
    mightGet = cms.optional.untracked.vstring
)

_hltInitialStepTrackCandidatesMkFitFitLSTSeeds = _hltInitialStepTrackCandidatesMkFit.clone(
    candMinPtCut = 0.9,
    candMinPtCutOuter = 0.8,
    candEtaRegion = 1.4
)

from Configuration.ProcessModifiers.singleIterPatatrack_cff import singleIterPatatrack
from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
from Configuration.ProcessModifiers.seedingLST_cff import seedingLST
from Configuration.ProcessModifiers.hltTrackingMkFitInitialStep_cff import hltTrackingMkFitInitialStep
(singleIterPatatrack & seedingLST & trackingLST & hltTrackingMkFitInitialStep).toReplaceWith(hltInitialStepTrackCandidatesMkFitFit, _hltInitialStepTrackCandidatesMkFitFitLSTSeeds)
