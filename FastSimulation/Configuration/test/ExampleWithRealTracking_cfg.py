import FWCore.ParameterSet.Config as cms

#########################################################################################################################
#
# Example to show how to run the real tracking instead of the emulated one after having created the tracker hits
#
#########################################################################################################################

process = cms.Process("PROD")

# Include DQMStore, needed by the famosSimHits
process.DQMStore = cms.Service( "DQMStore")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
    
# For valgrind studies
# process.ProfilerService = cms.Service("ProfilerService",
#    lastEvent = cms.untracked.int32(13),
#    firstEvent = cms.untracked.int32(3),
#    paths = cms.untracked.vstring('p1')
#)

# Include the RandomNumberGeneratorService definition
process.load("IOMC.RandomEngine.IOMC_cff")

# Generate H -> ZZ -> l+l- l'+l'- (l,l'=e or mu), with mH=200GeV/c2
process.load("Configuration.Generator.H200ZZ4L_cfi")

# Common inputs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('FastSimulation.Configuration.Geometries_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['mc']

# Parametrized magnetic field (new mapping, 4.0 and 3.8T)
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

# If you want to turn on/off pile-up
#process.load('FastSimulation.PileUpProducer.PileUpSimulator_2012_Startup_inTimeOnly_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')

# Detector simulation with FastSim
process.load("FastSimulation.EventProducer.FamosSimHits_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.simulationSequence = cms.Sequence(
    process.offlineBeamSpot+
    process.famosMixing+
    process.famosSimHits+
    process.mix)


# Needed to run the tracker digitizers
process.load('Configuration.StandardSequences.Digi_cff')
process.load('SimGeneral.MixingModule.pixelDigitizer_cfi')
process.pixelDigitizer.hitsProducer =  'famosSimHitsTrackerHits'
process.pixelDigitizer.makeDigiSimLinks = False # if you set to True, you need some more replacements
process.load('SimGeneral.MixingModule.stripDigitizer_cfi')
process.stripDigitizer.hitsProducer =  'famosSimHitsTrackerHits'
process.stripDigitizer.ROUList = ['famosSimHitsTrackerHits']
process.load('SimTracker.SiStripDigitizer.SiStripDigiSimLink_cfi')
process.simSiStripDigiSimLink.ROUList = ['famosSimHitsTrackerHits']

# Extend the MixingModule parameterset that tells which digitizers must be executed
from SimGeneral.MixingModule.pixelDigitizer_cfi import *
from SimGeneral.MixingModule.stripDigitizer_cfi import *

from SimGeneral.MixingModule.aliases_cfi import simEcalUnsuppressedDigis, simHcalUnsuppressedDigis, simSiPixelDigis, simSiStripDigis

from SimGeneral.MixingModule.ecalDigitizer_cfi import *
from SimCalorimetry.EcalSimProducers.ecalDigiParameters_cff import *
simEcalUnsuppressedDigis.hitsProducer = cms.string('famosSimHits')
ecal_digi_parameters.hitsProducer = cms.string('famosSimHits')
ecalDigitizer.hitsProducer = cms.string('famosSimHits')

import SimCalorimetry.HcalSimProducers.hcalUnsuppressedDigis_cfi
hcalSimBlockFastSim = SimCalorimetry.HcalSimProducers.hcalUnsuppressedDigis_cfi.hcalSimBlock.clone()
hcalSimBlockFastSim.hitsProducer = cms.string('famosSimHits')
hcalDigitizer = cms.PSet(
        hcalSimBlockFastSim,
        accumulatorType = cms.string("HcalDigiProducer"),
        makeDigiSimLinks = cms.untracked.bool(False))

process.mix.digitizers = cms.PSet(pixel = cms.PSet(pixelDigitizer),
                                  strip = cms.PSet(stripDigitizer),
                                  ecal = cms.PSet(ecalDigitizer),
                                  hcal = cms.PSet(hcalDigitizer))

# Needed to run the tracker reconstruction
process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
process.siPixelClusters.src = cms.InputTag("mix")
process.siStripZeroSuppression.RawDigiProducersList = cms.VInputTag( cms.InputTag('mix','VirginRaw'),
                                                                     cms.InputTag('mix','ProcessedRaw'),
                                                                     cms.InputTag('mix','ScopeMode'))
process.siStripZeroSuppression.DigisToMergeZS = cms.InputTag('siStripDigis','ZeroSuppressed')
process.siStripZeroSuppression.DigisToMergeVR = cms.InputTag('siStripVRDigis','VirginRaw')
##### THE FOLLOWING MUST BE FIXED SOMEHOW:
process.siStripClusters.DigiProducersList = DigiProducersList = cms.VInputTag(
    cms.InputTag('mix','ZeroSuppressed'),
    cms.InputTag('mix','VirginRaw'),
    cms.InputTag('mix','ProcessedRaw'),
    cms.InputTag('mix','ScopeMode'))

process.load('RecoTracker.Configuration.RecoTracker_cff')

# Temporary, for debugging
process.dumpContent = cms.EDAnalyzer('EventContentAnalyzer')

# Produce Tracks and Clusters
process.source = cms.Source("EmptySource")
process.gensim_step = cms.Path(process.generator*process.simulationSequence) # choose any sequence that you like in FamosSequences_cff
process.reconstruction_step = cms.Path(process.dumpContent+process.trackerlocalreco)

# To write out events
process.o1 = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("OutputFileWithRealTracks.root"),
    outputCommands = cms.untracked.vstring("keep *",
                                           "drop *_mix_*_*")
    )

process.outpath = cms.EndPath(process.o1)

process.schedule = cms.Schedule(process.gensim_step,process.reconstruction_step,process.outpath)

# Keep the logging output to a nice level #

#process.Timing =  cms.Service("Timing")
#process.load("FWCore/MessageService/MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt","cout")
#process.MessageLogger.categories.append("FamosManager")
#process.MessageLogger.cout = cms.untracked.PSet(threshold=cms.untracked.string("INFO"),
#                                                default=cms.untracked.PSet(limit=cms.untracked.int32(0)),
#                                                FamosManager=cms.untracked.PSet(limit=cms.untracked.int32(100000)))

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )
