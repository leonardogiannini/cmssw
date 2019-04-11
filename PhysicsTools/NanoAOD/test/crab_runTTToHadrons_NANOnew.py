from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'NANOREG_2'
config.General.workArea = '../../../../../crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myNanoProdMc_NANO.py'

config.JobType.allowUndistributedCMSSW = True

#config.JobType.outputFiles = ['lzma.root']

#/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD
#/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD
#/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD
#/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD
#/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD

#config.Data.lumiMask = '/scratch/lgiannini/ctpps/slc6/CMSSW_10_2_0/src/PhysicsTools/NanoAOD/test/combined_RPIN_CMS.json'

#ALTERNATIVE
#TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
#TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8
#TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8
#TTToHadronic_TuneCP5_13TeV-powheg-pythia8
config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM'
#config.Data.inputDataset = '/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 80000
#config.Data.totalUnits = 8000
config.Data.outLFNDirBase = '/store/user/%s/NanoREGInputs/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'NANO_withRegInputs'

config.Site.storageSite = 'T2_IT_Legnaro'
