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
#include <bitset>

/**
 * ConversionPhotonProducer prducing pat::CompositeCandidate.
 */
class ConversionPhotonProducer : public edm::EDProducer {
public:
  explicit ConversionPhotonProducer(const edm::ParameterSet&);
  ~ConversionPhotonProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  /** define the number of bits globally for this module */
  const static size_t nBits = 8; // at the moment need only 8
  /** typedef for consistent usage inside the module */
  using bitsetT = std::bitset<nBits>;

  /** small helper struct to group together any object with flags that have to be set in the corresponding pat object. */
  template<typename T>
  struct AnnotatedT {
    /** default constructor setting all flags to 0. Deliberately declared non-explicit! */
    AnnotatedT(T* obj, bitsetT f = bitsetT()) : object(obj), flags(f) {}
    T* object; /**< Object that shall be annotated with the flags. */
    bitsetT flags; /**< the bitset holding the flags. */
  };

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
  StringCutObjectSelector<reco::PFCandidate> m_pfCandCutSel; /**< string cut selector for PFCandidates */
  StringCutObjectSelector<reco::Photon> m_photonCutSel; /**< string cut selector for reco::Photons */
  double m_TkVtxCompSigma; /**< standard deviations in distance in z a track can be displaced to a vertex for compatiblity */
  double m_vertexChi2ProbCut; /**< chi2 probability cut that has to be fullfilled by the conversion vertex */
  double m_trackChi2Cut; /**< chi2 value cut for each of the tracks of the conversion */
  double m_minDistanceOfApproachMinCut; /**< min cut for the distance of minimum approach for the conversion */
  double m_minDistanceOfApproachMaxCut; /**< max cut for the distance of minimum approach for the conversion */
  double m_trackMinNDOF; /**< minimum numbers of degree of freedom for each of the tracks of the conversion */
  std::vector<double> m_pi0NarrowWindow; /**< narrow invariant mass window for the pi0 veto */
  std::vector<double> m_pi0WideWindow; /**< wide invariant mass window for the pi0 veto */

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
  unsigned m_photonCtr{}; /**< counter for created pat::Photons */

  // --------- private member functions --------------
  /** get all ParticleFlow candidates with particleId 'gamma' from the passed PFCandidateCollection*/
  const pat::PFParticleCollection getPFPhotons(const edm::Handle<edm::View<reco::PFCandidate> >& pfCands);

  /** collect all reco::Photons and create pat::Photons from them */
  const pat::PhotonCollection getPhotons(const edm::Handle<reco::PhotonCollection>& photons);

  /** collect all conversions and put them into a CompositeCandidateCollection and add some additional info */
  const pat::CompositeCandidateCollection getConversions(const edm::Handle<reco::ConversionCollection>& convs);

  /** create a pat::CompositeCandidate from a reco::Conversion */
  pat::CompositeCandidate makePhotonCandidate(const AnnotatedT<const reco::Conversion>& conversion);

  /** convert from one LorentzVectorType to the other */
  inline reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v) {
    return reco::Candidate::LorentzVector(v.x(), v.y(), v.z(), v.t());
  }

  /** add some flags to the pat Particle using data from the reco Particle (which can already have some other flags set)*/
  template<typename PatType, typename RecoType>
  inline void annotate(PatType& patPart, const AnnotatedT<RecoType>& recoPart) const
  {
     patPart.addUserInt("flags", recoPart.flags.to_ulong());
  }

  /** check if the conversion candidate has all desired quality flags */
  bool checkConversionQuality(const reco::Conversion& conv) const
  {
    for(const auto& qual : m_convQualities) {
      if( !conv.quality(qual)) {return false;}
    }
    return true;
  }

  /** get the bits that have to be set to true for the conversion for storing the bitset */
  const std::vector<unsigned short> getFlagBits(const reco::Conversion& conv) const;

  /** get the bits that have to be set to true for the photon for storing the bitset */
  const std::vector<unsigned short> getFlagBits(const reco::Photon& photon) const;

  /** get the bits that have to be set to true for the pfCandidate for storing the bitset */
  const std::vector<unsigned short> getFlagBits(const reco::PFCandidate& pfCand) const;

  /** check the compatibility of inner hits (i.e. hits closer to the beamspot than the vertex).
   * TODO: CHECK IF THIS IS WHAT IS ACTUALLY DONE HERE
   */
  bool checkCompatibleInnerHits(const reco::Conversion& conv) const;

  /** TODO: documentation */
  bool foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) const;

  /** check if the tracks of the conversion are compatible with one of the primary vertices. */
  bool checkTkVtxCompatibility(const reco::Conversion& conv) const;

  /** check if the conversion is follows a high purity definition
   * parameters used: trackMinNDOF, trackChi2Cut, vertexChi2ProbCut, minDistanceOfApproach[Min|Max]Cut
   */
  bool checkHighPuritySubset(const reco::Conversion& conv, const reco::VertexCollection& vtxs) const;

  /** check if two conversions share tracks.
   * If two conversions share a track only the one with the best chi2 will not get the corresponding flag set
   */
  void checkTrackSharing(std::vector<AnnotatedT<const reco::Conversion> >& conversions);

  /** Remove all objects from the collection that have any of the flags set in flags. */
  template<typename RecoType>
  void removeFlagged(std::vector<AnnotatedT<const RecoType> >& coll, const bitsetT& flags);

  /** Store the flags into an unsigned and store it in the patCandidate. */
  template<typename patType>
  void setFlags(patType& patCand, const bitsetT& bits);
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

#endif
