# HLTTutorial

## Setup
Setup the release, import HLTrigger package, compile
```
cmsrel CMSSW_10_6_1_patch3 
cd CMSSW_10_6_1_patch3/src
cmsenv
git cms-addpkg HLTrigger/Configuration
scram b -j 8
```
Creat a proxy (to access remote files) 

```
voms-proxy-init --voms cms --valid 168:00 
```

Clone this repository, and compile (N.B. some unused variables in the code introduce warnings)

```
git clone https://github.com/lathomas/HLTTutorial.git
USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b
```
Now, fetch the needed trigger configuration associated to the B-tagged Jet + soft muon path: 

```
hltGetConfiguration --cff --offline /dev/CMSSW_10_6_0/GRun --paths HLTriggerFirstPath,HLTriggerFinalPath --unprescale > HLT_TutoEffcySession_cff.py
hltGetConfiguration --cff /users/lathomas/HLTTuto2019/HLTTuto/V5   \
--globaltag auto:run2_hlt_GRun   \
--unprescale >> HLT_TutoEffcySession_cff.py
 ```
 You then need to modify HLT_TutoEffcySession_cff.py : search for these two lines, they appear twice in the config file. Comment out their second appearance (not the first one). 
 ```
#/users/HLTTutoOct2017/ExamplePath/Mu3Jet200Btag/V1 (CMSSW_9_2_10)   (already commented)                                                                                                                                       
#import FWCore.ParameterSet.Config as cms                                                                                                                                                                   
#fragment = cms.ProcessFragment( "HLT" )  
 ```
 
 We will also play a bit with the usual single electron path, we therefore need a config for that one too. 
```
hltGetConfiguration --cff /dev/CMSSW_10_6_0/GRun   \
--globaltag auto:run2_hlt_GRun   \
--path HLTriggerFirstPath,HLTriggerFinalPath,HLT_Ele35_WPTight_Gsf_v9 \
--unprescale > HLT_TutoEle35WPTight_cff.py
```
  
Copy the files to the python subdirectory
```
cp HLT_TutoEffcySession_cff.py ../python/.
cp HLT_TutoEle35WPTight_cff.py  ../python/.
```

 Produce a HLT2_HLT.py config file to be run with cmsRun HLT2_HLT.py. We will need to rerun HLT using RAW data and and then access offline information from MINIAOD. In order to do that, we need to specify the MINIAOD file as the main input file and the list of all its parent RAW files as secondary input files: 
```
cmsDriver.py HLT2 --step=HLT:TutoEffcySession --era=Run2_2018 --data --conditions auto:run2_hlt_GRun --filein root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/6307D2AA-47CD-A849-9F7F-DD2D6D183E5E.root --secondfilein file:root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/RAW/v1/000/321/164/00000/BAF0A515-439E-E811-838D-02163E019F14.root,root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/RAW/v1/000/321/164/00000/407FB415-439E-E811-9977-FA163E2FE09D.root,root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/RAW/v1/000/321/164/00000/00833922-439E-E811-9D43-02163E01A04C.root --processName=HLT2 -n 100 --no_exec 
```


So far the HLT2_HLT.py file contains a configuration to rerun HLT on top of RAW. 
We need to add our EDAnalyzer as a second module in the process. After: 
```
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
```
Add: 
```
process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD')
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( "out.root" )
                                   )
process.demo_step = cms.EndPath(process.demo)
```
Replace the line:
```
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])
```
by:
```
process.schedule.extend([process.endjob_step, process.demo_step])
```
N.B.: Keeping RECOSIMoutput_step in the process would create an updated (big!) RAW file 
with the information related to the trigger menu you reran. 
This file is typically quite big so we will drop it here. As an exercise, you may wish to produce 
it and take a look at its content, for example by doing: 
```
edmDumpEventContent --regex HLT HLT2_HLT.root
```
Test the code: 
```
cmsRun HLT2_HLT.py 
```

## Exercices

All the code (and exercises) related to this tutorial is to be found in:
https://github.com/lathomas/HLTTutorial/blob/master/TriggerAnalyzerRAWMiniAOD/plugins/TriggerAnalyzerRAWMiniAOD.cc

You will need to modify this file and compile (in CMSSW_9_2_12/src). The exercises are inserted as comments in this EDAnalyzer. After each modification of the analyzer, do: 
```
cd $CMSSW_BASE/src
USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j 6
cd HLTrigger/Configuration/test
cmsRun HLT2_HLT.py
```

There are six exercises: 


 - Exercise 1: learn how to access the trigger decision (boolean) associated to any path (in RAW, AOD, MINIAOD).
 - Exercise 2: learn how to access the trigger objects stored in MINIAOD.
 - Exercise 3: learn how to access the trigger objects stored in RAW.
 - Exercise 4: perform an efficiency measurement and define a proper denominator. NB: switch "maxEvents" to 10000 instead of 100 for exercise 4. (this may take some time...)
 - Exercise 5: Check some efficiency plots similar to those obtained in ex 4 but with higher stats, study some interesting behaviours.
 - Exercise 6: Perform a measurement using Tag and Probe with dielectron events and access the HLT variables used in the various filters associated to a single electron path. One can process 1000 events in this exercise. 

