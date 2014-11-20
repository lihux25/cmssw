# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:upgradePLS3 -n 10 --eventcontent FEVTDEBUGHLT,DQM -s RAW2DIGI,L1Reco,RECO,VALIDATION,DQM --datatier GEN-SIM-RECO,DQMIO --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023HGCalMuon --geometry Extended2023HGCalV4Muon,Extended2023HGCalV4MuonReco --magField 38T_PostLS1 --io RecoFull_Extended2023HGCalV4.io --python RecoFull_Extended2023HGCalV4.py --no_exec --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/067472B3-BF5F-E411-8A52-0025905B8576.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/08E4D96C-C35F-E411-B5A3-002618943900.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/0C2FC8FD-B85F-E411-B9F2-0025905B858A.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/426B9D2E-BE5F-E411-91A5-0025905A6136.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/8428E767-C35F-E411-898F-002354EF3BDD.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/CE8D2C41-BE5F-E411-90EA-00261894395B.root',
       '/store/relval/CMSSW_6_2_0_SLHC20/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPG2023SHNoTaper-v1/00000/D050D768-B85F-E411-BCE0-0025905A6138.root'
#	'/store/relval/CMSSW_6_2_0_SLHC19/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/1856D652-0154-E411-AD05-0025905A7786.root',
#	'/store/relval/CMSSW_6_2_0_SLHC19/RelValSinglePiE50HCAL/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/BA4F854E-0154-E411-AB15-003048FF86CA.root'
#        '/store/relval/CMSSW_6_2_0_SLHC19/RelValQCD_Pt_80_120_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/36FB7BB4-1354-E411-B2CF-0025905A60DE.root',
#        '/store/relval/CMSSW_6_2_0_SLHC19/RelValQCD_Pt_80_120_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/684182D4-1A54-E411-8FD6-0025905964B2.root',
#        '/store/relval/CMSSW_6_2_0_SLHC19/RelValQCD_Pt_80_120_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/7E4E1700-0F54-E411-907E-0026189438EB.root',
#        '/store/relval/CMSSW_6_2_0_SLHC19/RelValQCD_Pt_80_120_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/94EC63D5-1A54-E411-B3A4-0025905964A2.root',
#        '/store/relval/CMSSW_6_2_0_SLHC19/RelValQCD_Pt_80_120_14TeV/GEN-SIM-DIGI-RAW/DES23_62_V1_UPGHGCalV4-v1/00000/C878E064-1254-E411-9E0F-0025905A6132.root'

#       '/store/relval/CMSSW_6_2_0_SLHC19/RelValSinglePiE50HCAL/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV4-v1/00000/42AC872B-0654-E411-876B-0025905A6080.root',
#       '/store/relval/CMSSW_6_2_0_SLHC19/RelValSinglePiE50HCAL/GEN-SIM-RECO/DES23_62_V1_UPGHGCalV4-v1/00000/EC2EEA2D-0654-E411-BA43-0025905A6104.root' 
                                                                                                                                                          )  
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:step3.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.DQMEventContent.outputCommands,
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DQMIO')
    )
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.prevalidation_step = cms.Path(process.prevalidation)
process.dqmoffline_step = cms.Path(process.DQMOffline)
process.validation_step = cms.EndPath(process.validation)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.prevalidation_step,process.validation_step,process.dqmoffline_step,process.DQMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023SHCal

#call to customisation function cust_2023SHCal imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023SHCal(process)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

file = open('allDump_step3_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
