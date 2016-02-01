#!/bin/bash
#Specify needed variables
vartype=int #double, int, float
obj=patElectron_  #rMet_type1PF_pt #Tau_, Jet_, Muon_
ob= #met
varList=(isHEEPId)
#_Declare variables for reading
echo " "
pos=0
for count in ${varList[@]};
do
 echo "  vector<$vartype>* r$obj${varList[$pos]}; r$obj${varList[$pos]} = 0; TBranch* b_r$obj${varList[$pos]} = 0; readingtree->SetBranchAddress(\"$obj${varList[$pos]}\",&r$obj${varList[$pos]},&b_r$obj${varList[$pos]});"
 let pos=pos+1
done

#Declare variables for writing
echo " "
pos=0
for count in ${varList[@]};
do
 echo "  vector<$vartype>* $obj${varList[$pos]} = new std::vector<$vartype>; newtree->Branch(\"$obj${varList[$pos]}\",&$obj${varList[$pos]});"  
 let pos=pos+1
done

#Read branches
echo " "
pos=0
for count in ${varList[@]};
do
 echo "   b_r$obj${varList[$pos]}->GetEntry(en);"
 let pos=pos+1
done

#Clear new branches
echo " "
pos=0
for count in ${varList[@]};
do
 echo "   $obj${varList[$pos]}->clear();"  
 let pos=pos+1
done

#Save new branches content
echo " "
pos=0
for count in ${varList[@]};
do
 echo "    $obj${varList[$pos]}->push_back(r$obj${varList[$pos]}->at("$ob"_en));"
 let pos=pos+1
done
#Muon
#pt eta phi energy pTErrOVpT_it charge soft loose medium tight isHighPt POGisGood pdgId pf isGlobal isTrackerMuon trackIso ecalIso hcalIso caloIso isoSum pfEcalEnergy chi2 chi2LocalPosition matchedStat validHits validHitsInner TLayers ndof validFraction pixelLayersWithMeasurement qualityhighPurity trkKink segmentCompatibility dz_pv dxy_pv miniIsoRel miniIsoCh miniIsoNeu miniIsoPUsub  jetdr jetpt jetptratio jetcsv ptrel IP3Dsig_it pvass etarel ptOVen mujet_pfJetProbabilityBJetTag mujet_pfCombinedMVABJetTags mujetmass mujetWmass mujetTopmass mujetWTopmass IP3D_val IP3D_err IP3D_sig IP2D_val IP2D_err IP2D_sig sIP3D_val sIP3D_err sIP3D_sig sIP2D_val sIP2D_err sIP2D_sig IP1D_val IP1D_err IP1D_sig sIP1D_val sIP1D_err sIP1D_sig lepjetMaxIP3D_val lepjetMaxIP3D_sig lepjetMaxsIP3D_val lepjetMaxsIP3D_sig lepjetMaxIP2D_val lepjetMaxIP2D_sig lepjetMaxsIP2D_val lepjetMaxsIP2D_sig lepjetMaxIP1D_val lepjetMaxIP1D_sig lepjetMaxsIP1D_val lepjetMaxsIP1D_sig lepjetAvIP3D_val lepjetAvIP3D_sig lepjetAvsIP3D_val lepjetAvsIP3D_sig lepjetAvIP2D_val lepjetAvIP2D_sig lepjetAvsIP2D_val lepjetAvsIP2D_sig lepjetAvIP1D_val lepjetAvIP1D_sig lepjetAvsIP1D_val lepjetAvsIP1D_sig lepjetchtrks lepjetpvchtrks lepjetnonpvchtrks lepjetndaus lepjetpvchi2 lepjetnumno2trk gen_pt gen_eta gen_phi gen_en gen_pdgId gen_isPromptFinalState gen_isDirectPromptTauDecayProductFinalState
#Electron
#pt eta phi energy Et SCeta inCrack charge passVetoId passLooseId passMediumId passTightId passHEEPId pdgId isEcalDriven isoChargedHadrons isoNeutralHadrons isoPhotons isoPU relIsoDeltaBeta relIsoRhoEA isolPtTracks dEtaIn dPhiIn full5x5_sigmaIetaIeta full5x5_e2x5Max full5x5_e5x5 full5x5_e1x5 hOverE ooEmooP passConversionVeto_ expectedMissingInnerHits gsfTrack_ndof gsfTrack_normChi2 gsfTrack_dz_pv gsfTrack_dxy_pv d0 miniIsoRel miniIsoCh miniIsoNeu miniIsoPUsub jetdr jetpt jetptratio jetcsv ptrel IP3Dsig eleMVASpring15NonTrig25ns eleMVASpring15NonTrig25ns_VL gen_pt gen_eta gen_phi gen_en gen_pdgId gen_isPromptFinalState gen_isDirectPromptTauDecayProductFinalState
#Tau
#pt eta phi energy charge againstMuonLoose3 againstMuonTight3 againstElectronVLooseMVA5 againstElectronLooseMVA5 againstElectronMediumMVA5 againstElectronTightMVA5 againstElectronVTightMVA5 againstElectronMVA5raw byLooseCombinedIsolationDeltaBetaCorr3Hits leadChargedCandDz_pv leadChargedCandDxy_pv)
#Jet
#pt eta phi energy mass px py pz Uncorr_pt pfCombinedInclusiveSecondaryVertexV2BJetTags pfCombinedMVABJetTags pfJetProbabilityBJetTags pileupId isPFJet isCaloJet neutralHadEnergyFraction neutralEmEnergyFraction chargedHadronEnergyFraction chargedEmEnergyFraction muonEnergyFraction electronEnergy photonEnergy emEnergyFraction numberOfConstituents chargedMultiplicity vtxMass vtxNtracks vtx3DVal vtx3DSig JesSF JesSFup JesSFdown JerSF JerSFup JerSFdown partonFlavour hadronFlavour
