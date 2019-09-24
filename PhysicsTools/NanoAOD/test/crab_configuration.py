from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '1NANOREG16UL'
config.General.workArea = '../../../../../crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myNanoProdMc2016_NANO.py'

config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 80000
config.Data.totalUnits = 80000
config.Data.outLFNDirBase = '/store/user/%s/Nano16legacyReg/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'NANO16legacy_RegInputs'

config.Site.storageSite = 'T2_IT_Legnaro'
