import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from Configuration.StandardSequences.Eras import eras

from RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsProtonReconstructionOF_cfi import *



ctppsProtonTable = cms.EDProducer("SimpleProtonTrackFlatTableProducer",
    src = cms.InputTag("ctppsProtonReconstructionOFDB"),
    cut = cms.string(""), # filtered already above
    name = cms.string("Proton"),
    doc  = cms.string("bon"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        
        xi = Var("xi",float,doc="xi or dp/p",precision=10),  
        xiError = Var("xiErr",float,doc="error on xi or dp/p",precision=10),
        
        vx = Var("vertex().x()",float,doc="vx",precision=10),
        vy = Var("vertex().y()",float,doc="vy",precision=10),
        vz = Var("vertex().z()",float,doc="vz",precision=10),
        
        #x,y,z
        dirx = Var("direction().x()",float,doc="dir x",precision=10),
        diry = Var("direction().y()",float,doc="dir y",precision=10),
        dirz = Var("direction().z()",float,doc="dir z",precision=10),
        
        #const double th_x = p.direction().x() / p.direction().mag();
        #const double th_y = p.direction().y() / p.direction().mag();
        dirMag = Var("direction().mag()",float,doc="dir mag",precision=10),
        thx = Var("direction().x() / direction().mag()",float,doc="th x",precision=10),
        thy = Var("direction().y() / direction().mag()",float,doc="th y",precision=10),
        
        #th,phi,eta ??????
        #sono di tipo  Geom::Phi<T>, Geom::Theta<T> e non numerico
        #dirTheta = Var("direction().theta().value()",float,doc="dir theta",precision=10),
        #dirPhi = Var("direction().phi().value(); ",float,doc="dir phi",precision=10),
        dirEta = Var("direction().eta()",float,doc="dir eta",precision=10),
        
        rmSingleRP = Var("singleRP", bool,doc="is rmSingleRP"),  
        rmMultiRP = Var("multiRP", bool,doc="is rmMultiRP"),
        sector45 = Var("isSector45", bool,doc="is sector45"),
        sector56 = Var("isSector56", bool,doc="is sector56"),
        
        pixelFlag = Var("pixelFlag", bool, doc="is pixel"),  
        armFlag = Var("armFlag", bool, doc="arm where"),
        stationFlag = Var("statFlag", bool, doc="station"),
        rpFlag = Var("rpFlag", bool, doc="rp is 3"),
        rpArmStatId = Var("rpArmStatId", int, doc="all info as in plotter"),        
        
        fitNDF = Var("fitNDF",int ,doc="fitNDF"),  
        fitChiSq = Var("fitChiSq",float,doc="fitChiSquare",precision=10),
        
        #halfCrossingAngleSector45=0., halfCrossingAngleSector56=0.; 
        xangle45 = Var("halfCrossingAngleSector45",float,doc="halfCrossingAngleSector45",precision=10),         
        xangle56 = Var("halfCrossingAngleSector56",float,doc="halfCrossingAngleSector56",precision=10),

        valid = Var("valid",bool,doc="valid"),  
        
    ),

)


nmuTable = cms.EDProducer("GlobalVariablesTableProducer",
    variables = cms.PSet(
        #MuonSlimmedMuonsHT = ExtVar( cms.InputTag("slimmedMuons"), "candidatescalarsum", doc = "number of of all the slimmed muons" ),
        SlimmedMuonsN = ExtVar( cms.InputTag("slimmedMuons"), "candidatesize", doc = "number of of all the slimmed muons" ),

    )
)


ctppsProtonTables = cms.Sequence(ctppsProtonReconstructionOFDB+ctppsProtonTable+nmuTable )
