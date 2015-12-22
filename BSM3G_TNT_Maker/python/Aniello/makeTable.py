import ROOT, sys, math

tf = ROOT.TFile("/afs/cern.ch/work/a/aspiezia/ttH/CMSSW_7_4_12/src/BSMFramework/BSM3G_TNT_Maker/python/OutTree_tth.root")
tt = tf.Get("TNT/BOOM")
data = False

print "run,lumi,event,\
is_SL,is_DL,\
lep1_pt,lep1_eta,lep1_phi,lep1_iso,lep1_pdgId,\
lep2_pt,lep2_eta,lep2_phi,lep2_iso,lep2_pdgId,\
mll,mll_passed,\
jet1_pt,jet2_pt,jet3_pt,jet4_pt,\
jet1_CSVv2,jet2_CSVv2,jet3_CSVv2,jet4_CSVv2,\
jet1_JesSF,jet2_JesSF,jet3_JesSF,jet4_JesSF,\
jet1_JerSF,jet2_JerSF,jet3_JerSF,jet4_JerSF,\
MET_pt,MET_phi,\
met_passed,\
n_jets,n_btags,\
bWeight,\
ttHFCategory,\
final_discriminant1,final_discriminant2,\
n_fatjets,\
pt_fatjet_1,pt_fatjet_2,\
pt_nonW_1,pt_nonW_2,\
pt_W1_1,pt_W1_2,\
pt_W2_1,pt_W2_2,\
m_top_1,m_top_2,\
higgstag_fatjet_1,higgstag_fatjet_2"


for ev in range(tt.GetEntries()):
    tt.GetEntry(ev)

    ##INITIALIZE LEPTON VARIABLES
    muons_pt = tt.Muon_pt
    electrons_pt = tt.patElectron_pt
    SelLepton_pt = []
    SelLepton_id = []
    SelLepton_eta = []
    SelLepton_phi = []
    SelLepton_iso = []
    SelLepton_energy = []

    ##MUON SELECTION
    for i in range(len(muons_pt)):
        if (muons_pt[i]>15 and math.fabs(tt.Muon_eta[i])<2.4 and tt.Muon_tight[i]==1 and tt.Muon_relIsoDeltaBetaR04[i]<0.25):
            SelLepton_pt.append(tt.Muon_pt[i])
            SelLepton_id.append(tt.Muon_pdgId[i])
            SelLepton_eta.append(tt.Muon_eta[i])
            SelLepton_phi.append(tt.Muon_phi[i])
            SelLepton_iso.append(tt.Muon_relIsoDeltaBetaR04[i])
            SelLepton_energy.append(tt.Muon_energy[i])

    #ELECTRON SELECTION
    for i in range(len(electrons_pt)):
        if (tt.patElectron_pt[i]>15 and math.fabs(tt.patElectron_eta[i])<2.4 and tt.patElectron_inCrack[i]==0 and tt.patElectron_isPassMvatrig[i]==1 and tt.patElectron_relIsoRhoEA[i]<0.15 and ((math.fabs(tt.patElectron_SCeta[i])<1.4442 and tt.patElectron_full5x5_sigmaIetaIeta[i]<0.012 and tt.patElectron_hOverE[i]<0.09 and (tt.patElectron_ecalPFClusterIso[i]/tt.patElectron_pt[i])<0.37 and (tt.patElectron_hcalPFClusterIso[i]/tt.patElectron_pt[i])<0.25 and (tt.patElectron_isolPtTracks[i]/tt.patElectron_pt[i])<0.18 and math.fabs(tt.patElectron_dEtaIn[i])<0.0095 and math.fabs(tt.patElectron_dPhiIn[i])<0.065) or (math.fabs(tt.patElectron_SCeta[i])>1.5660 and tt.patElectron_full5x5_sigmaIetaIeta[i]<0.033 and tt.patElectron_hOverE[i]<0.09 and (tt.patElectron_ecalPFClusterIso[i]/tt.patElectron_pt[i])<0.45 and (tt.patElectron_hcalPFClusterIso[i]/tt.patElectron_pt[i])<0.28 and (tt.patElectron_isolPtTracks[i]/tt.patElectron_pt[i])<0.18))):
            SelLepton_pt.append(tt.patElectron_pt[i])
            SelLepton_id.append(tt.patElectron_pdgId[i])
            SelLepton_eta.append(tt.patElectron_eta[i])
            SelLepton_phi.append(tt.patElectron_phi[i])
            SelLepton_iso.append(tt.patElectron_relIsoRhoEA[i])
            SelLepton_energy.append(tt.patElectron_energy[i])

    #INITIALIZE JET VARIABLES
    jet0_pt = 0
    jet0_eta = 0
    jet0_phi = 0
    jet0_csv = 0  
    jet0_JesSF = 0    
    jet0_JerSF = 0    
    jet1_pt = 0
    jet1_eta = 0
    jet1_phi = 0
    jet1_csv = 0
    jet1_JesSF = 0    
    jet1_JerSF = 0    
    jet2_pt = 0
    jet2_eta = 0
    jet2_phi = 0
    jet2_csv = 0
    jet2_JesSF = 0    
    jet2_JerSF = 0    
    jet3_pt = 0
    jet3_eta = 0
    jet3_phi = 0
    jet3_csv = 0
    jet3_JesSF = 0    
    jet3_JerSF = 0    
    SelJet_pt = []
    SelJet_eta = []
    SelJet_phi = []
    SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags = []   
    SelJet_JesSF = []   
    SelJet_JerSF = []
    SelTightJet_pt = []
    SelTightJet_eta = []
    SelTightJet_phi = []
    SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags = []  
    SelTightJet_JesSF = []   
    SelTightJet_JerSF = []
    jets_pt = tt.Jet_pt

    ##JET SELECTION FOR DILEPTON EVENTS
    for i in range(len(jets_pt)):
        #jet_pt=(tt.Jet_Uncorr_pt[i]*tt.Jet_JesSF[i]*tt.Jet_JerSF[i])
        jet_pt=tt.Jet_pt[i]*tt.Jet_JerSF[i]
        if (jet_pt>20 and math.fabs(tt.Jet_eta[i])<2.4 and tt.Jet_neutralHadEnergyFraction[i]<0.99 
            and tt.Jet_chargedEmEnergyFraction[i]<0.99 and tt.Jet_neutralEmEnergyFraction[i]<0.99 and tt.Jet_numberOfConstituents[i]>1 
            and tt.Jet_chargedHadronEnergyFraction[i]>0.0 and tt.Jet_chargedMultiplicity[i]>0.0):
            deltaRJetLepBoolean = False
            for j in range(len(SelLepton_pt)):
                deltaEta = SelLepton_eta[j]-tt.Jet_eta[i];
                deltaPhi = math.fabs(SelLepton_phi[j]-tt.Jet_phi[i]);
                if(deltaPhi > math.pi): 
                    deltaPhi = 2*math.pi - deltaPhi;
                if(math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)<0.4): 
                    deltaRJetLepBoolean = True
            if(deltaRJetLepBoolean==False):
                SelJet_pt.append(jet_pt)
                SelJet_eta.append(tt.Jet_eta[i])
                SelJet_phi.append(tt.Jet_phi[i])
                SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.append(tt.Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i])
                SelJet_JesSF.append(tt.Jet_JesSF[i])
                SelJet_JerSF.append(tt.Jet_JerSF[i])
        #if(jet_pt>15 and tt.EVENT_event==904458):
        #    print ""
        #    print "%i tt.Jet_pt[i] %f tt.Jet_JerSF[i] %f jet_pt %f" %(tt.EVENT_event,tt.Jet_pt[i],tt.Jet_JerSF[i],jet_pt)
        #    print "%i tt.Jet_eta[i] %f math.fabs(tt.Jet_eta[i]) %f" %(tt.EVENT_event,tt.Jet_eta[i],math.fabs(tt.Jet_eta[i]))
        #    print "%i tt.Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i] %f" %(tt.EVENT_event,(tt.Jet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i]))
        #    print "%i jet_pt>20 %i" %(tt.EVENT_event,jet_pt>20)
        #    print "%i math.fabs(tt.Jet_eta[i])<2.4 %i" %(tt.EVENT_event,math.fabs(tt.Jet_eta[i])<2.4)
        #    print "%i tt.Jet_neutralHadEnergyFraction[i]<0.99 %i" %(tt.EVENT_event,tt.Jet_neutralHadEnergyFraction[i]<0.99)
        #    print "%i tt.Jet_chargedEmEnergyFraction[i]<0.99 %i" %(tt.EVENT_event,tt.Jet_chargedEmEnergyFraction[i]<0.99)
        #    print "%i tt.Jet_neutralEmEnergyFraction[i]<0.99 %i" %(tt.EVENT_event,tt.Jet_neutralEmEnergyFraction[i]<0.99)
        #    print "%i tt.Jet_numberOfConstituents[i]>1 %i" %(tt.EVENT_event,tt.Jet_numberOfConstituents[i]>1)
        #    print "%i tt.Jet_chargedHadronEnergyFraction[i]>0.0 %i" %(tt.EVENT_event,tt.Jet_chargedHadronEnergyFraction[i]>0.0)
        #    print "%i tt.Jet_chargedMultiplicity[i]>0.0 %i" %(tt.EVENT_event,tt.Jet_chargedMultiplicity[i]>0.0)
        #    print "%i deltaRJetLepBoolean==0 %i" %(tt.EVENT_event,deltaRJetLepBoolean==0)

    ##JET SELECTION FOR SINGLE LEPTON EVENTS
    for i in range(len(SelJet_pt)):
        if(SelJet_pt[i]>30):
            SelTightJet_pt.append(SelJet_pt[i])
            SelTightJet_eta.append(SelJet_eta[i])
            SelTightJet_phi.append(SelJet_phi[i])
            SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags.append(SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i])
            SelTightJet_JesSF.append(SelJet_JesSF[i])
            SelTightJet_JerSF.append(SelJet_JerSF[i])

    #BTAG FOR SINGLE AND DILEPTON EVENTS
    nBCSVM_SL = 0
    for i in range(len(SelTightJet_pt)):
        if (SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i]>0.89):
            nBCSVM_SL=nBCSVM_SL+1
    nBCSVM_DL = 0
    for i in range(len(SelJet_pt)):
        if (SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[i]>0.89):
            nBCSVM_DL=nBCSVM_DL+1

    ##MET CORRECTION
    deltaJetPx = 0
    deltaJetPy = 0
    for i in range(len(jets_pt)):
        if(tt.Jet_Uncorr_pt[i]>10 and ((tt.Jet_isPFJet[i]==0 and tt.Jet_emEnergyFraction[i] < 0.9) or (tt.Jet_isPFJet[i]==1 and tt.Jet_neutralEmEnergyFraction[i]+tt.Jet_chargedEmEnergyFraction[i]<0.9))):
            deltaJetPx = deltaJetPx + (tt.Jet_px[i]*tt.Jet_JerSF[i] - tt.Jet_px[i])
            deltaJetPy = deltaJetPy + (tt.Jet_py[i]*tt.Jet_JerSF[i] - tt.Jet_py[i])
    correctedMetPx = tt.Met_type1PF_px - deltaJetPx;
    correctedMetPy = tt.Met_type1PF_py - deltaJetPy;
    METp4 = ROOT.TLorentzVector(0, 0, 0, 0)
    METp4.SetPxPyPzE(correctedMetPx,correctedMetPy,0,math.sqrt(correctedMetPx*correctedMetPx+correctedMetPy*correctedMetPy))
    #MET = METp4.Pt()
    #MET_phi = METp4.Phi()
    MET = tt.Met_type1PF_pt
    MET_phi = tt.Met_type1PF_phi

    ##DEFINE LEPTON VARIABLES
    lep0_pt = 0
    lep0_id = 0
    lep0_eta = 0
    lep0_phi = 0
    lep0_iso = 0
    lep1_pt = 0
    lep1_id = 0
    lep1_eta = 0
    lep1_phi = 0
    lep1_iso = 0
    isSL = 0
    isDL = 0
    mll = -99
    mll_pass = 1
    met_pass = 1

    ##DILEPTON SELECTION
    if (len(SelLepton_pt) == 2):
        dilepSelection = True
        isSL = 0
        isDL = 1
        lepCouplePt=-99
        for i in range(len(SelLepton_pt)):
            for j in range(i+1,len(SelLepton_pt)):
                lep0 = ROOT.TLorentzVector(0, 0, 0, 0)
                lep0.SetPtEtaPhiE(SelLepton_pt[i], SelLepton_eta[i], SelLepton_phi[i], SelLepton_energy[i])
                lep1 = ROOT.TLorentzVector(0, 0, 0, 0)
                lep1.SetPtEtaPhiE(SelLepton_pt[j], SelLepton_eta[j], SelLepton_phi[j], SelLepton_energy[j])
                mll = (lep0+lep1).M()
                if(SelLepton_pt[i]<20 and SelLepton_pt[j]<20):
                    dilepSelection=False
                if(SelLepton_id[i]*SelLepton_id[j]>0):
                    dilepSelection=False
                if((lep0+lep1).M()<20): 
                    mll_pass=0
                #    dilepSelection=False
                if(math.fabs(SelLepton_id[i])==math.fabs(SelLepton_id[j]) and (lep0+lep1).M()>76 and (lep0+lep1).M()<106): 
                    mll_pass=0
                #    dilepSelection=False
                if(math.fabs(SelLepton_id[i])==math.fabs(SelLepton_id[j]) and MET<40): 
                    met_pass=0
                #    dilepSelection=False
                if(dilepSelection==True and lepCouplePt<SelLepton_pt[i]+SelLepton_pt[j]):
                    lep0_pt = SelLepton_pt[i]
                    lep0_id = SelLepton_id[i]
                    lep0_eta = SelLepton_eta[i]
                    lep0_phi = SelLepton_phi[i]
                    lep0_iso = SelLepton_iso[i]
                    lep0_energy = SelLepton_energy[i]
                    lep1_pt = SelLepton_pt[j]
                    lep1_id = SelLepton_id[j]
                    lep1_eta = SelLepton_eta[j]
                    lep1_phi = SelLepton_phi[j]
                    lep1_iso = SelLepton_iso[j]
                    lep1_energy = SelLepton_energy[i]
                    lepCouplePt = SelLepton_pt[i]+SelLepton_pt[j]
    
    ##SELECT SINGLE-LEPTON EVENT
    NumberOfJets=0
    SingleLeptonEvent=False
    ##MC
    if(data==False and len(SelLepton_pt)==1 and ((math.fabs(SelLepton_id[0])==11 and SelLepton_pt[0]>30 and tt.HLT_Ele27_WP85_Gsf==1) or (math.fabs(SelLepton_id[0])==13 and SelLepton_pt[0]>25 and SelLepton_iso[0]<0.15 and tt.HLT_IsoMu17_eta2p1==1)) and math.fabs(SelLepton_eta[0])<2.1 and len(SelTightJet_pt)>=4 and nBCSVM_SL>=2):
        isSL = 1
        isDL = 0
        SingleLeptonEvent=True
        nBCSVM=nBCSVM_SL
        NumberOfJets=len(SelTightJet_pt)
        lep0_pt = SelLepton_pt[0]
        lep0_id = SelLepton_id[0]
        lep0_eta = SelLepton_eta[0]
        lep0_phi = SelLepton_phi[0]
        lep0_iso = SelLepton_iso[0]
        lep0_energy = SelLepton_energy[0]
        jet0_pt = SelTightJet_pt[0]
        jet0_eta = SelTightJet_eta[0] 
        jet0_phi = SelTightJet_phi[0] 
        jet0_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[0]
        jet0_JesSF = SelTightJet_JesSF[0]
        jet0_JerSF = SelTightJet_JerSF[0]
        jet1_pt = SelTightJet_pt[1]
        jet1_eta = SelTightJet_eta[1] 
        jet1_phi = SelTightJet_phi[1] 
        jet1_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[1]
        jet1_JesSF = SelTightJet_JesSF[1]
        jet1_JerSF = SelTightJet_JerSF[1]
        jet2_pt = SelTightJet_pt[2]
        jet2_eta = SelTightJet_eta[2] 
        jet2_phi = SelTightJet_phi[2] 
        jet2_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[2]
        jet2_JesSF = SelTightJet_JesSF[2]
        jet2_JerSF = SelTightJet_JerSF[2]
        jet3_pt = SelTightJet_pt[3]
        jet3_eta = SelTightJet_eta[3] 
        jet3_phi = SelTightJet_phi[3] 
        jet3_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[3] 
        jet3_JesSF = SelTightJet_JesSF[3]
        jet3_JerSF = SelTightJet_JerSF[3]
    ##DATA
    if(data==True and len(SelLepton_pt)==1 and ((math.fabs(SelLepton_id[0])==11 and SelLepton_pt[0]>30 and tt.HLT_Ele27_eta2p1_WPLoose_Gsf==1) or (math.fabs(SelLepton_id[0])==13 and SelLepton_pt[0]>25 and SelLepton_iso[0]<0.15 and tt.HLT_IsoMu18==1)) and math.fabs(SelLepton_eta[0])<2.1 and len(SelTightJet_pt)>=4 and nBCSVM_SL>=2):
        isSL = 1
        isDL = 0
        SingleLeptonEvent=True
        nBCSVM=nBCSVM_SL
        NumberOfJets=len(SelTightJet_pt)
        lep0_pt = SelLepton_pt[0]
        lep0_id = SelLepton_id[0]
        lep0_eta = SelLepton_eta[0]
        lep0_phi = SelLepton_phi[0]
        lep0_iso = SelLepton_iso[0]
        lep0_energy = SelLepton_energy[0]
        jet0_pt = SelTightJet_pt[0]
        jet0_eta = SelTightJet_eta[0] 
        jet0_phi = SelTightJet_phi[0] 
        jet0_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[0]
        jet0_JesSF = SelTightJet_JesSF[0]
        jet0_JerSF = SelTightJet_JerSF[0]
        jet1_pt = SelTightJet_pt[1]
        jet1_eta = SelTightJet_eta[1] 
        jet1_phi = SelTightJet_phi[1] 
        jet1_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[1]
        jet1_JesSF = SelTightJet_JesSF[1]
        jet1_JerSF = SelTightJet_JerSF[1]
        jet2_pt = SelTightJet_pt[2]
        jet2_eta = SelTightJet_eta[2] 
        jet2_phi = SelTightJet_phi[2] 
        jet2_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[2]
        jet2_JesSF = SelTightJet_JesSF[2]
        jet2_JerSF = SelTightJet_JerSF[2]
        jet3_pt = SelTightJet_pt[3]
        jet3_eta = SelTightJet_eta[3] 
        jet3_phi = SelTightJet_phi[3] 
        jet3_csv = SelTightJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[3] 
        jet3_JesSF = SelTightJet_JesSF[3]
        jet3_JerSF = SelTightJet_JerSF[3]
    
    ##SELECT DOUBLE-LEPTON EVENT
    DoubleLeptonEvent=False
    if((tt.HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ==1 or tt.HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL==1 or tt.HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL==1 or tt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ==1 or tt.HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ==1) and len(SelLepton_pt)==2 and len(SelJet_pt)>=2 and SelJet_pt[0]>30 and SelJet_pt[1]>30 and nBCSVM_DL>=1 and dilepSelection==True):
        DoubleLeptonEvent=True
        nBCSVM=nBCSVM_DL
        NumberOfJets=len(SelJet_pt)
        jet0_pt = SelJet_pt[0]
        jet0_eta = SelJet_eta[0] 
        jet0_phi = SelJet_phi[0] 
        jet0_JesSF = SelJet_JesSF[0]
        jet0_JerSF = SelJet_JerSF[0]
        jet0_csv = SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[0]
        jet1_pt = SelJet_pt[1]
        jet1_eta = SelJet_eta[1] 
        jet1_phi = SelJet_phi[1] 
        jet1_csv = SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[1]
        jet1_JesSF = SelJet_JesSF[1]
        jet1_JerSF = SelJet_JerSF[1]
        if(len(SelJet_pt)>2):
            jet2_pt = SelJet_pt[2]
            jet2_eta = SelJet_eta[2] 
            jet2_phi = SelJet_phi[2] 
            jet2_csv = SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[2]
            jet2_JesSF = SelJet_JesSF[2]
            jet2_JerSF = SelJet_JerSF[2]
        if(len(SelJet_pt)>3):
            jet3_pt = SelJet_pt[3]
            jet3_eta = SelJet_eta[3] 
            jet3_phi = SelJet_phi[3] 
            jet3_csv = SelJet_pfCombinedInclusiveSecondaryVertexV2BJetTags[3] 
            jet3_JesSF = SelJet_JesSF[3]
            jet3_JerSF = SelJet_JerSF[3]
        
    #if(tt.EVENT_event==18923434 or tt.EVENT_event==18946508 or tt.EVENT_event==19029883 or tt.EVENT_event==19063781 or tt.EVENT_event==19178549 or tt.EVENT_event==19394723):
    #    print ""
    #    print "SINGLE LEPTON"
    #    print "%i tt.triggerSL==1                     %i" %(tt.EVENT_event, tt.triggerSL==1)
    #    print "%i len(SelLepton_pt)==1                %i" %(tt.EVENT_event, len(SelLepton_pt)==1)
    #    if(len(SelLepton_pt)>0):
    #        print "%i pt cut                              %i" %(tt.EVENT_event, ((math.fabs(SelLepton_id[0])==11 and SelLepton_pt[0]>30) or (math.fabs(SelLepton_id[0])==13 and SelLepton_pt[0]>25)))
    #        print "%i math.fabs(SelLepton_eta[0])<2.1     %i" %(tt.EVENT_event, math.fabs(SelLepton_eta[0])<2.1)
    #    print "%i len(SelTightJet_pt)>=4              %i" %(tt.EVENT_event, len(SelTightJet_pt)>=4)
    #    print "%i nBCSVM_SL>=2                        %i" %(tt.EVENT_event, nBCSVM_SL>=2)
    #    print "DOUBLE LEPTON"
    #    print "%i tt.triggerDL==1                     %i" %(tt.EVENT_event, tt.triggerDL==1)
    #    print "%i len(SelLepton_pt)==2                %i" %(tt.EVENT_event, len(SelLepton_pt)==2)
    #    print "%i len(SelJet_pt)>=2                   %i" %(tt.EVENT_event, len(SelJet_pt)>=2)
    #    print "%i nBCSVM_DL>=1                        %i" %(tt.EVENT_event, nBCSVM_DL>=1)
    #    print "%i dilepSelection                      %i" %(tt.EVENT_event, dilepSelection)
    #    print ""
    #    print "%i len(SelLepton_pt)                   %i" %(tt.EVENT_event, len(SelLepton_pt))
    #    print ""

    ##ORDERING LEPTON PT
    if(DoubleLeptonEvent):
        lep_PROV_pt  = lep0_pt 
        lep_PROV_eta = lep0_eta
        lep_PROV_phi = lep0_phi
        lep_PROV_iso = lep0_iso
        lep_PROV_id  = lep0_id 
        if(lep0_pt<lep1_pt):
            lep0_pt  = lep1_pt 
            lep0_eta = lep1_eta
            lep0_phi = lep1_phi
            lep0_iso = lep1_iso
            lep0_id  = lep1_id 
            lep1_pt  = lep_PROV_pt 
            lep1_eta = lep_PROV_eta
            lep1_phi = lep_PROV_phi
            lep1_iso = lep_PROV_iso
            lep1_id  = lep_PROV_id 
       
    ttHFCategory = -99
    if(data==False):
        ttHFCategory =(tt.ttHFCategory)

    #if(SingleLeptonEvent or DoubleLeptonEvent):
    if(SingleLeptonEvent):
    #if(DoubleLeptonEvent):
        arr = [int(tt.EVENT_run), int(tt.EVENT_lumiBlock), int(tt.EVENT_event),
               int(isSL), int(isDL),
               lep0_pt, lep0_eta, lep0_phi, lep0_iso, lep0_id,
               lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_id,
               mll, mll_pass,
               jet0_pt, jet1_pt, jet2_pt, jet3_pt,
               jet0_csv, jet1_csv, jet2_csv, jet3_csv,
               jet0_JesSF, jet1_JesSF, jet2_JesSF, jet3_JesSF,
               jet0_JerSF, jet1_JerSF, jet2_JerSF, jet3_JerSF,
               MET, MET_phi,
               met_pass,
               int(NumberOfJets), int(nBCSVM),
               float(tt.bWeight),
               int(ttHFCategory), int(-99), int(-99), int(-99),
               int(-99), int(-99), int(-99), int(-99), int(-99), int(-99),
               int(-99), int(-99), int(-99), int(-99), int(-99), int(-99)
               ]
        
        s = ""
        for i in range(len(arr)):
            if not isinstance(arr[i], int):
                s += str(round(arr[i], 4)) + ","
            else:
                s += str(arr[i]) + ","
        print s[:-1]
