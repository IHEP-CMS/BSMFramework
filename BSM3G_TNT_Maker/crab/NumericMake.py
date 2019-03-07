import os
import re
import sys
import commands
import subprocess

#########################
############### Information provided by user
#################################

# Path you run this script
workpath = "/publicfs/cms/user/libh/Test/Rootplizer/analyzer"
#Specify needed variables
#  Type   Reading  Saving  Example
#  CaseA   Yes      Yes     Jet_pt: you want to save it and need also to read it (using SetBranchStatus+CopyEntries would require more processing time) 
#  CaseB   Yes      No      Jet_muonEnergyFraction: you want to read it to define jet ID, but may not want to save it  
#  CaseC   No       Yes     Jet_LooseID: you can not read it from ntuplas if there is not, but you want to save in rootplas 

# Variable Case
#Case = "CaseA"
#Case = "CaseB"
Case = "CaseC"
# Variable Definition

VariableType  = "double"
VariableNames = [
####### Data ############
'Fall17V1_SEleBlockB',
'Fall17V1_SEleBlockC',
'Fall17V1_SEleBlockD',
'Fall17V1_SEleBlockE',
'Fall17V1_SEleBlockF',
'Fall17V1_SMuBlockB',
'Fall17V1_SMuBlockC',
'Fall17V1_SMuBlockD',
'Fall17V1_SMuBlockE',
'Fall17V1_SMuBlockF',
###### MC ###########
'Fall17V1_TTHnobb',
'Fall17V1_TTWToLNu',
'Fall17V1_TTW_PSwgt_ToLNu',
'Fall17V1_TTZToLLNuNu_M10',
'Fall17V1_TTZToLL_M1to10',
'Fall17V1_TTWW',
'Fall17V1_DY_M10to50',
'Fall17V1_DY_M50',
'Fall17V1_DY_ext_M50',
'Fall17V1_WJets',
'Fall17V1_WWTo2L2Nu',
'Fall17V1_WZTo3LNu',
'Fall17V1_ZZTo4L',
'Fall17V1_ZZ_ext_To4L',
'Fall17V1_TT_PSwgt_To2L2Nu',
'Fall17V1_TTTo2L2Nu',
'Fall17V1_TT_PSwgt_ToSemiLep',
'Fall17V1_TTToSemiLep',
'Fall17V1_TT_PSwgt_ToHadron',
'Fall17V1_TTToHadron',
'Fall17V1_ST_tW_top',
'Fall17V1_ST_tW_antitop',
'Fall17V1_STt_top',
'Fall17V1_STt_antitop',
'Fall17V1_STs',
'Fall17V1_TTGJets',
'Fall17V1_tZq',
'Fall17V1_WW_DoubleScatter',
'Fall17V1_WW_DS_To2L2Nu',
'Fall17V1_WWW',
'Fall17V1_WWZ',
'Fall17V1_WZZ',
'Fall17V1_ZZZ',
'Fall17V1_TTTT_Tune',
]

#VariableType  = "int"
#VariableNames = [
# CaseA Variables
#"Muon_soft","Muon_loose","Muon_medium","Muon_tight","Muon_isHighPt",
#"Muon_POGisGood","Muon_pdgId","Muon_pf","Muon_isGlobal","Muon_isTrackerMuon",
#"Muon_tunePBestTrackType","Muon_gen_pdgId","Muon_gen_isPromptFinalState","Muon_gen_isDirectPromptTauDecayProductFinalState"

# CaseB Variables

# CaseC Variables
#]

# Create the variable file
Vectorname = VariableType+Case+"_Numeric.cc"
vector     = file (Vectorname,"w")

#ReadTreeptr & WriteTreeptr
RTreeptr = "readingtree"
WTreeptr = "newtree"

#Name of Current Entry
ParEntry = "tentry"
#Name of index in Push_back
ParSel = "gen_en"

###################
### Script itself
#################
if Case == "CaseA":
 print >> vector, "//This is Case A"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector,VariableType+" r"+Variable+"; TBranch* b_r"+Variable+" =0;" 
 
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector,VariableType+" "+Variable+";" 
 
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Initialize Number"
 for Variable in VariableNames:
     print >> vector, " "+Variable+"= -999;"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, " b_r"+Variable+"->GetEntry("+ParEntry+");"
 
 print >> vector, "//Write_variables"
 for Variable in VariableNames:
     print >> vector, "  "+Variable+" = r"+Variable+";"

elif Case == "CaseB":
 print >> vector, "//This is Case B"
 print >> vector, "//Head file declaration"
 print >> vector, "//variables to be read"
 for Variable in VariableNames:
     print >> vector,VariableType+" r"+Variable+"; TBranch* b_r"+Variable+" =0;" 
 
 print >> vector, "//source file definition"
 print >> vector, "//read setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+RTreeptr+'->SetBranchAddress("'+Variable+'",&r'+Variable+",&b_r"+Variable+");"
 
 print >> vector, "//GetEntry"
 for Variable in VariableNames:
     print >> vector, " b_r"+Variable+"->GetEntry("+ParEntry+");"

 print >> vector, "//This is for MVA"
 for Variable in VariableNames:
     print >> vector, ' factory->AddVariable("'+Variable+'",'+"'D');"
 
 print >> vector, "//This is for dasclient"
 for Variable in VariableNames:
     print >> vector, './das_client.py --query="file dataset='+ Variable+' | sum(file.nevents)"'
 
elif Case == "CaseC":
 print >> vector, "//This is Case C"
 print >> vector, "//Head file declaration"
 
 print >> vector, "//variables to be written"
 for Variable in VariableNames:
     print >> vector,VariableType+" "+Variable+";" 
 
 
 print >> vector, "//source file definition"
 print >> vector, "//Write setbranchaddress"
 for Variable in VariableNames:
     print >> vector, " "+WTreeptr+'->Branch("'+Variable+'",&'+Variable+");"
 
 print >> vector, "//Initialize Number"
 for Variable in VariableNames:
     print >> vector, " "+Variable+"= -999;"
 
 print >> vector, "//Write_variables"
 for Variable in VariableNames:
     print >> vector, "  "+Variable+" = r"+Variable+";"

 print >> vector, "//this is for crab"
 for Variable in VariableNames:
     print >> vector, "crab resubmit -d "+Variable+"/crab_"+Variable
else :
 print >> vector, "//!!!!!!Please specify a correct Case!!!!!!//"
