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
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
  std::vector<double> m_pi0WideWindow; /**< wide mass window for pi0 veto. */
  std::vector<double> m_pi0NarrowWindow; /**< narrow mass window for pi0 veto. */

  // ----------handles to used data ------------------------------
  edm::Handle<pat::PFParticleCollection> m_pfCandColl; /**< Handle to the PFParticle photons. */
  edm::Handle<pat::PhotonCollection> m_photonColl; /**< Handle to the pat Photons. */
  edm::Handle<pat::CompositeCandidateCollection> m_convColl; /**< Handle to the conversions. */

  // ----------member data ---------------------------
  TFile* m_outputFile; /**< The root output file. */
  TTree* m_tree; /**< The tree in the root file. */

  PARootVariables m_rootVariables; /**< All Variables that are written to the root file. */
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
