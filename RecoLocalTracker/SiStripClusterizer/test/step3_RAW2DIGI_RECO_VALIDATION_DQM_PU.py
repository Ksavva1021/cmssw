# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2021_realistic -n 5000 --era Run3 --eventcontent RECOSIM,DQM --runUnscheduled -s RAW2DIGI,RECO:reconstruction_trackingOnly,VALIDATION:@trackingOnlyValidation,DQM:@trackingOnlyDQM --datatier GEN-SIM-RECO,DQMIO --geometry DB:Extended --pileup_input das:/RelValMinBias_13/CMSSW_10_6_0_pre3-105X_postLS2_realistic_v6-v1/GEN-SIM --pileup AVE_35_BX_25ns --filein file:step2.root --fileout step3.root --nThreads 20 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('reReco', Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/eos/user/a/abrusamo/TTbar_13TeV_generation_CMSSW_12_1_0_pre2/step2_'+str(i)+'.root' for i in range(1, 6)),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:5000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:/eos/user/a/abrusamo/siStripClusters.root'),
    outputCommands = cms.untracked.vstring(
      'keep *',
      'keep *_siStripClusters*_*_*'
    ),
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(#'file:/eos/user/a/abrusamo/originalTTbar35PU_filter_pixelLessOnly_inDQM.root'
                                    'file:approxClusters.root' 
    ),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '120X_upgrade2018_realistic_v1', '')


#this module gets a collection of standard strip clusters (list of amplitudes) as input
#it computes the barycenter using the strip amplitudes as weights, 
#the avg charge (tot charge/width).
#the output is a collection of SiStripApproximatedClusters; each cluster is described through its
#barycenter, avg charge and width (number of strips) and each cluster is of type SiStripApproximateCluster

process.siStripClusters2ApproxClusters = cms.EDProducer("SiStripClusters2ApproxClusters",
	inputClusters = cms.InputTag("siStripClusters")
)

#this module allows to have a different resolution on the barycenter and on
#the charge (see RecoLocalTracker/SiStripClusterizer/pluging/SiStripApprox2Approx
#for the different resolutions)
#the input of this module is a collection of SiStripApproximateClusters and the output is 
#still a collection of SiStripApproximateClusters
process.siStripApprox2ApproxClusters = cms.EDProducer("SiStripApprox2ApproxClusters",
	inputApproxClusters = cms.InputTag("siStripClusters2ApproxClusters"),
    approxVersion= cms.string("ORIGINAL")
)

#process.siStripApproximatedClustersDump = cms.EDAnalyzer("SiStripApproximatedClustersDump",
    #approximatedClustersTag = cms.InputTag("SiStripClusters2ApproxClusters")
#    approximatedClustersTag = cms.InputTag("siStripApprox2ApproxClusters")
#)

#now the collection of SiStripApproximatedClusters is converted back to
#a collection of SiStripClusters; we use the SiStripCluster(SiStripApproximateCluster cluster)
#constructor

process.load('RecoLocalTracker.SiStripClusterizer.SiStripApprox2Clusters_cfi')
process.siStripConvertedClusters = process.SiStripApprox2Clusters.clone(
  inputApproxClusters = 'siStripClusters2ApproxClusters'
)

#remove materialDumpAnalyzer
process.DQMOfflineTracking = cms.Sequence(process.TrackingDQMSourceTier0+process.DQMOfflineVertex)

process.SiStripDQMSource = cms.Sequence(
  process.APVPhases+
  process.consecutiveHEs+
  process.siStripFEDCheck+
  process.siStripFEDMonitor+
  process.SiStripMonitorDigi+
  process.SiStripMonitorClusterBPTX+
  process.SiStripMonitorTrackCommon+
  process.refittedForPixelDQM+
  process.MonitorTrackResiduals+
  process.dqmInfoSiStrip
)

#Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.siStripClusters_step = cms.Path(process.siStripClusters)
process.siStripClusters2ApproxClusters_step = cms.Path(process.siStripClusters2ApproxClusters)
process.approxToApproxClusters_step = cms.Path(process.siStripApprox2ApproxClusters)
process.reconstruction_step = cms.Path(process.reconstruction_trackingOnly)
process.striptrackerlocalrecoTask = cms.Task(process.siStripConvertedClusters,process.siStripMatchedRecHits)
process.striptrackerlocalreco_step = cms.Path(process.striptrackerlocalrecoTask)
process.prevalidation_step = cms.Path(process.globalPrevalidationTrackingOnly) #tracking only validation; no infos from other subdetectors
process.validation_step = cms.EndPath(process.globalValidationTrackingOnly)     
process.dqmoffline_step = cms.EndPath(process.SiStripDQMSource+
  process.DQMOfflineTracking)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)



# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,
  process.L1Reco_step, 
  process.siStripClusters_step, #clusterizer 
  process.siStripClusters2ApproxClusters_step, #convert SiStripCluster to SiStripApproximateCluster
  process.approxToApproxClusters_step, #change resolution on barycenter and charge (optional; see SiStripApprox2ApproxCluster) 
  process.striptrackerlocalreco_step, #convert SiStripApproximateCluster to SiStripCluster
  process.reconstruction_step, #standard reconstruction
  process.prevalidation_step,
  process.validation_step,
  process.dqmoffline_step,process.dqmofflineOnPAT_step,
  #process.RECOSIMoutput_step,
  process.DQMoutput_step
)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(30)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.
process.MeasurementTrackerEventPreSplitting.stripClusterProducer = "siStripConvertedClusters"
process.MeasurementTrackerEvent.stripClusterProducer = "siStripConvertedClusters"

from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
MassReplaceInputTag(process, new='siStripConvertedClusters', old='siStripClusters') #use siStripConvertedCluster instead of siStripCluster everywhere in the reconstruction steps


process.siStripMatchedRecHits.ClusterProducer = "siStripConvertedClusters" 

process.siStripClusters2ApproxClusters.inputClusters = "siStripClusters" #MassReplaceInputTag changes it, but we need SiStripCluster in this step

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
