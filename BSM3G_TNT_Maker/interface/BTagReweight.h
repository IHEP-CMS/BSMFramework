#ifndef __BTAGREWEIGHT_HE_H_
#define __BTAGREWEIGHT_HE_H_
/////
//   Include files and namespaces
/////
#include <memory>
#include <map>
#include <algorithm>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <TTree.h>
#include "TLorentzVector.h"
#include "TMath.h"
#include <Math/VectorUtil.h>
#include "baseTree.h"
using namespace std;
using namespace edm;
/////
//   Class declaration
/////
class BTagReweight : public baseTree{
 public:
  BTagReweight(std::string name, TTree* tree, bool debug, const edm::ParameterSet& cfg, edm::ConsumesCollector && ic);
  ~BTagReweight();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void JECInitialization();
  void GetJER(pat::Jet jet, float JesSF, float rhoJER, bool AK4PFchs, float &JERScaleFactor, float &JERScaleFactorUP, float &JERScaleFactorDOWN);
  void GetLeptonsForDeltaRWithJets(vector<TLorentzVector> &LeptonsForDeltaRWithJets, const edm::Event& iEvent);
  double get_effarea(double eta);
  bool isGoodVertex(const reco::Vertex& vtx);
 private:
  BTagReweight(){};
  /////
  //   Config variables
  /////
  bool _is_data;
  double bWeight;
  double bWeightLFup;
  double bWeightLFdown;
  double bWeightHFup;
  double bWeightHFdown;
  double bWeightHFStats1up;
  double bWeightHFStats1down;
  double bWeightLFStats1up;
  double bWeightLFStats1down;
  double bWeightHFStats2up;
  double bWeightHFStats2down;
  double bWeightLFStats2up;
  double bWeightLFStats2down;
  double bWeightCErr1up;
  double bWeightCErr1down;
  double bWeightCErr2up;
  double bWeightCErr2down;
  double bWeightJESup;
  double bWeightJESdown;
  edm::EDGetTokenT<reco::VertexCollection> vtx_h_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electron_pat_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muon_h_;
  edm::EDGetTokenT<pat::JetCollection> jets_;
  edm::EDGetTokenT<double> rhopogHandle_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVATrigIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>  > eleMVAnonTrigIdMap_;
  edm::EDGetTokenT<double> rhoJERHandle_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC1_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC2_;
  edm::FileInPath jecPayloadNamesAK4PFchsMC3_;
  edm::FileInPath jecPayloadNamesAK4PFchsMCUnc_;
  std::string jerAK4PFchs_;
  std::string jerAK4PFchsSF_;
  boost::shared_ptr<FactorizedJetCorrector>   jecAK4PFchsMC_;
  boost::shared_ptr<JetCorrectionUncertainty> jecAK4PFchsMCUnc_;
  int    _vtx_ndof_min;
  int    _vtx_rho_max;
  double _vtx_position_z_max;
  /////
  //   IHEP methods/variables
  /////
  edm::FileInPath BTAGReweightfile1_;
  edm::FileInPath BTAGReweightfile2_;
  void fillCSVhistos(TFile *fileHF, TFile *fileLF);
  double get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
  // CSV reweighting (2 x 2D array, 1 x 3D array)
  TH1D* h_csv_wgt_hf[9][6];
  TH1D* c_csv_wgt_hf[9][6];
  TH1D* h_csv_wgt_lf[9][4][3];
};
#endif
