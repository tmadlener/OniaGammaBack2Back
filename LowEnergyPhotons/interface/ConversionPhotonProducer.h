// -*- C++ -*-
//
// Package:    LowEnergyPhotons
// Class:      LowEnergyPhotons
//
/**\class ConversionPhotonProducer OniaGammaBack2Back/LowEnergyPhotons/interface/LowEnergyPhotons.h

 Description: Class for Producing Conversion Photons

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Madlener
//         Created:  Thu Feb  4 09:52:14 CET 2016
// $Id$
//
//

#ifndef __OniaGamma_ConversionPhotonProducer_h_
#define __OniaGamma_ConversionPhotonProducer_h_

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// actually needed objects
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// stl
#include <vector>
#include <string>

/**
 * ConversionPhotonProducer prducing pat::CompositeCandidate.
 */
class ConversionPhotonProducer : public edm::EDProducer {
public:
  explicit ConversionPhotonProducer(const edm::ParameterSet&);
  ~ConversionPhotonProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::ConversionCollection> m_convCollTok; /**< token to get the ConversionCollection */
  edm::EDGetTokenT<reco::PhotonCollection> m_photonCollTok; /**< token to get the photon collection */
  edm::EDGetTokenT<edm::View<reco::PFCandidate> > m_pfCandViewTok; /**< token to get the pf candidate view */
  edm::EDGetTokenT<reco::BeamSpot> m_beamSpotTok; /**< token to get the beamspot */
  edm::EDGetTokenT<reco::VertexCollection> m_vertexCollTok; /**< token to get the vertex collection */
  StringCutObjectSelector<reco::Conversion> m_conversionCutSel; /**< string cut selector for conversions */
  double m_TkVtxCompSigma; /**< standard deviations in distance in z a track can be displaced to a vertex for compatiblity */

  reco::Conversion::ConversionAlgorithm m_convAlgo; /**< desired conversion algorithm */
  std::vector<reco::Conversion::ConversionQuality> m_convQualities; /**< desired conversion qualities */

  // ---------- data handles -------------------------
  edm::Handle<reco::ConversionCollection> m_convColl; /**< handle to the ConversionCollection */
  edm::Handle<reco::PhotonCollection> m_photonColl; /**< handle to the PhotonCollection */
  edm::Handle<edm::View<reco::PFCandidate> > m_pfCandView; /**< handle to the PFCandidates (view) */
  edm::Handle<reco::BeamSpot> m_beamSpotHand; /**< handle to the beamspot */
  edm::Handle<reco::VertexCollection> m_vertexColl; /**< handle to the vertex collection */

  // ---------- counter variables ------------------------
  unsigned m_patConvCtr{}; /**< counter for created conversions */
  unsigned m_patPfPartCtr{}; /**< counter for created particle flow particles */

  // --------- private member functions --------------
  /** get all ParticleFlow candidates with particleId 'gamma' from the passed PFCandidateCollection */
  const pat::PFParticleCollection getPFPhotons(const edm::Handle<edm::View<reco::PFCandidate> >& pfCands);

  /** collect all reco::Photons and create pat::Photons from them */
  const pat::PhotonCollection getPhotons(const edm::Handle<reco::PhotonCollection>& photons);

  /** collect all conversions and put them into a CompositeCandidateCollection and add some additional info */
  const pat::CompositeCandidateCollection getConversions(const edm::Handle<reco::ConversionCollection>& convs);

  /** create a pat::CompositeCandidate from a reco::Conversion */
  pat::CompositeCandidate makePhotonCandidate(const reco::Conversion& conversion);

  /** convert from one LorentzVectorType to the other */
  inline reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v) {
    return reco::Candidate::LorentzVector(v.x(), v.y(), v.z(), v.t());
  }

  /** add some flag stuff to the pat::CompositeCandidate */
  void annotate(pat::CompositeCandidate& patConv, const reco::Conversion& recoConv) const;

  /** check if the conversion candidate has all desired quality flags */
  bool checkConversionQuality(const reco::Conversion& conv) const
  {
    for(const auto& qual : m_convQualities) {
      if( !conv.quality(qual)) return false;
    }
    return true;
  }

  /** get the bits that have to be set to true for the conversion for storing the bitset */
  const std::vector<unsigned short> getConversionFlagBits(const reco::Conversion& conv) const;

  /** check the compatibility of inner hits (i.e. hits closer to the beamspot than the vertex).
   * TODO: CHECK IF THIS IS WHAT IS ACTUALLY DONE HERE
   */
  bool checkCompatibleInnerHits(const reco::Conversion& conv) const;

  /** TODO: documentation */
  bool foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) const;

  /** check if the tracks of the conversion are compatible with one of the primary vertices. */
  bool checkTkVtxCompatibility(const reco::Conversion& conv) const;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

#endif
