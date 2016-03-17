// -*- C++ -*-
//
// Package:    LowEnergyPhotons
// Class:      PhotonAnalyzer
//
/**\class ConversionPhotonProducer OniaGammaBack2Back/LowEnergyPhotons/interface/PhotonAnalyzer.h

 Description: Class for Analyzing Photons

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Madlener
//         Created:  Thu Mar 7  09:52:14 CET 2016
// $Id$
//
//

#ifndef __OniaGamma_PhotonAnalyzer_h_
#define __OniaGamma_PhotonAnalyzer_h_

// system include files
#include <memory>

#include "OniaGammaBack2Back/PhotonAnalyzer/interface/PhotonAnalysisRootVariables.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// #include "DataFormats/VertexReco/interface/VertexFwd.h"

// #include "DataFormats/Common/interface/TriggerResults.h"

#include <TFile.h>
#include <TTree.h>

#include <string>
#include <vector>

//
// class declaration
//

class PhotonAnalyzerPAT : public edm::EDAnalyzer {
public:
  explicit PhotonAnalyzerPAT(const edm::ParameterSet&);
  ~PhotonAnalyzerPAT();

  static void fillDescription(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------input parameters-----------------------
  std::string m_rootFileName; /**< File name of the output root file. */
  edm::EDGetTokenT<pat::PFParticleCollection> m_pfCollTok; /**< Token to retrieve the pat::PFParticle photons. */
  edm::EDGetTokenT<pat::PhotonCollection> m_photonCollTok; /**< Token to retrieve the pat::Photon collection. */
  edm::EDGetTokenT<pat::CompositeCandidateCollection> m_convCollTok; /**< Token to retrieve the conversions. */
  edm::EDGetTokenT<pat::CompositeCandidateCollection> m_diMuonTok; /**< Token to retrieve the dimuon candidates. */
  edm::EDGetTokenT<pat::MuonCollection> m_muonTok; /**< Token to retrieve the pat::Muons. */
  edm::EDGetTokenT<reco::VertexCollection> m_primVertTok; /**< Token to retrieve the primary vertices. */
  edm::EDGetTokenT<edm::TriggerResults> m_triggerResultsTok; /**< Token to retrieve the trigger results. */

  double m_oniaMassMax; /**< Max. mass a dimuon candidate is allowed to have for storage. */
  double m_oniaMassMin; /**< Min. mass a dimuon candidate is allowed to have for storage. */
  bool m_storeOnlyBest; /**< Store only the 'best' dimuon candidate. */

  // ----------handles to used data ------------------------------
  edm::Handle<pat::PFParticleCollection> m_pfCandColl; /**< Handle to the PFParticle photons. */
  edm::Handle<pat::PhotonCollection> m_photonColl; /**< Handle to the pat Photons. */
  edm::Handle<pat::CompositeCandidateCollection> m_convColl; /**< Handle to the conversions. */
  edm::Handle<pat::CompositeCandidateCollection> m_diMuonColl; /**< Handle to the dimuons. */
  edm::Handle<pat::MuonCollection> m_muonColl; /**< Handle to the muon collection. */
  edm::Handle<reco::VertexCollection> m_primVertices; /**< Handle to the primary vertices. */
  edm::Handle<edm::TriggerResults> m_triggerResults; /**< Handle to the trigger results. */

  // ----------member data ---------------------------
  TFile* m_outputFile; /**< The root output file. */
  TTree* m_tree; /**< The tree in the root file. */

  PARootVariables m_rootVariables; /**< All Variables that are written to the root file. */

  // ----------counter variables -------------------------
  unsigned m_convPi0{}; /**< conversions with pi0 veto counter. */

  // ----------member functions -------------------------
  /** pack the trigger results into 1 bitset */
  unsigned getTriggerBits(const edm::Event& triggerResults);

  /** Get the information available from the dimuon candidates and put them to the rootVariables.
   * Returns the number of stored dimuon candidates.
   */
  size_t getDiMuonInfo(const pat::CompositeCandidateCollection& dimuonColl);

  /** Get the information available from muon candidates and put them to the rootVariables.
   * Returns the number of stored muons.
   * NOTE: This is not implemented yet, as I see no point in doing it right now. Returns 0 always!
   */
  size_t getSingleMuonInfo(const pat::MuonCollection& muonColl);

  /** Get the information avialable from the photon candidates. */
  template<typename PATType>
  size_t getPhotonInfo(const std::vector<PATType>& photonColl);

  /** Get the source category from the passed Type for storing. */
  template<typename PATType>
  short getPhotonCat(const std::vector<PATType>& coll) const;
};

//
// constants, enums and typedefs
//
/** names of the triggers that are checked and put into the triggerBits with v1 and v2. */
const std::vector<std::string> CheckTriggerNames = {"HLT_Dimuon16_Jpsi_v",
                                                    "HLT_Dimuon13_PsiPrime_v",
                                                    "HLT_Dimuon13_Upsilon_v",
                                                    "HLT_Dimuon10_Jpsi_Barrel_v",
                                                    "HLT_Dimuon8_PsiPrime_Barrel_v",
                                                    "HLT_Dimuon8_Upsilon_Barrel_v",
                                                    "HLT_Dimuon20_Jpsi_v",
                                                    "HLT_Dimuon0_Phi_Barrel_v",
                                                    "HLT_HIL1DoubleMu0_v",
                                                    "HLT_HIL2DoubleMu0_v",
                                                    "HLT_HIL2Mu3_v",
                                                    "HLT_HIL3Mu3_v"};


//
// static data member definitions
//

#endif
