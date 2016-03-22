#include "OniaGammaBack2Back/PhotonsProducer/interface/OGB2BPhotonProducer.h"
#include "OniaGammaBack2Back/PhotonsProducer/interface/PiZeroChecker.h"

#include "CommonTools/Utils/interface/StringToEnumValue.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// root stuff
#include <TMath.h>

// stl stuff
#include <vector>
#include <typeinfo>
#include <sstream>
#include <map>

// helper stuff for dev
#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"
#include "/afs/hephy.at/user/t/tmadlener/snippets/type_deduction_helper.h"

// ============================== constructor / destructor ==============================
OGB2BPhotonProducer::OGB2BPhotonProducer(const edm::ParameterSet& iConfig) :
  m_convCollTok( consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions")) ),
  m_photonCollTok( consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("allPhotons")) ),
  m_pfCandViewTok( consumes<edm::View<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("pfcandidates"))  ),
  m_beamSpotTok( consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")) ),
  m_vertexCollTok( consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVtxTag")) ),
  m_conversionCutSel(iConfig.getParameter<std::string>("convSelection")),
  m_pfCandCutSel(iConfig.getParameter<std::string>("pfCandSelection")),
  m_photonCutSel(iConfig.getParameter<std::string>("photonSelection")),
  m_TkVtxCompSigma(iConfig.getParameter<double>("tkVtxCompSigma")),
  m_vertexChi2ProbCut(iConfig.getParameter<double>("vertexChi2ProbCut")),
  m_trackChi2Cut(iConfig.getParameter<double>("trackChi2Cut")),
  m_minDistanceOfApproachMinCut(iConfig.getParameter<double>("minDistanceOfApproachMinCut")),
  m_minDistanceOfApproachMaxCut(iConfig.getParameter<double>("minDistanceOfApproachMaxCut")),
  m_trackMinNDOF(iConfig.getParameter<double>("trackMinNDOF")),
  m_pi0NarrowWindow(iConfig.getParameter<std::vector<double> >("pi0NarrowWindow")),
  m_pi0WideWindow(iConfig.getParameter<std::vector<double> >("pi0WideWindow"))
{
  std::string algo = iConfig.getParameter<std::string>("convAlgorithm");
  // convert the returned int into the enum-value in the constructor to avoid having to do it later
  m_convAlgo = (reco::Conversion::ConversionAlgorithm)StringToEnumValue<reco::Conversion::ConversionAlgorithm>(algo);

  // store the desired qualities (input by string) into the internally used vector
  std::vector<std::string> qualities = iConfig.getParameter<std::vector<std::string> >("convQuality");
  for(const auto& qual : qualities) {
    m_convQualities.push_back( (reco::Conversion::ConversionQuality) StringToEnumValue<reco::Conversion::ConversionQuality>(qual) );
  }


  std::string cutStr = iConfig.getParameter<std::string>("convSelection");
  m_convFlags = createFlags(iConfig.getParameter<std::string>("convFlags"));
  if (!cutStr.empty()) m_convFlags.set(0); // set flag if a cut selection is specified

  cutStr = iConfig.getParameter<std::string>("photonSelection");
  m_photonFlags = createFlags(iConfig.getParameter<std::string>("photonFlags"));
  if (!cutStr.empty()) m_photonFlags.set(0); // set flag if a cut selection is specified

  cutStr = iConfig.getParameter<std::string>("pfCandSelection");
  m_pfPhotonFlags = createFlags(iConfig.getParameter<std::string>("pFlowFlags"));
  if (!cutStr.empty()) m_pfPhotonFlags.set(0); // set flag if a cut selection is specified

  produces<pat::CompositeCandidateCollection>("convertedPhotons");
  produces<pat::PhotonCollection>("photons");
  produces<pat::PFParticleCollection>("PFlowPhotons");
}

OGB2BPhotonProducer::~OGB2BPhotonProducer()
{
  // TODO
}

// ============================== produce ==============================
void OGB2BPhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get the data to the handles (COULDDO: put this into a separate function?)
  iEvent.getByToken(m_convCollTok, m_convColl);
  iEvent.getByToken(m_photonCollTok, m_photonColl);
  iEvent.getByToken(m_pfCandViewTok, m_pfCandView);
  iEvent.getByToken(m_beamSpotTok, m_beamSpotHand);
  iEvent.getByToken(m_vertexCollTok, m_vertexColl);

  // only provide the photons to the pi0checker not all PFCandidates
  std::vector<reco::PFCandidate> recoPfPhotons;
  for(const auto& part : *m_pfCandView) {
    if (part.particleId() == reco::PFCandidate::ParticleType::gamma) recoPfPhotons.push_back(part);
  }

  PiZeroCheckFunction piZeroWindows(m_pi0WideWindow, m_pi0NarrowWindow);
  PiZeroChecker pi0Checker(piZeroWindows);
  pi0Checker.addCollections(*m_photonColl, recoPfPhotons, *m_convColl);
  pi0Checker.check();

  // process the conversions
  std::auto_ptr<pat::CompositeCandidateCollection> patConvOutColl(new pat::CompositeCandidateCollection);
  const auto convColl = getConversions(m_convColl, pi0Checker.getCheckIndices<decltype(*m_convColl)>());
  for(const auto& patConv : convColl) {
    m_patConvCtr++;
    patConvOutColl->push_back(patConv);
  }

  // process the PFCandidates
  std::auto_ptr<pat::PFParticleCollection> pfPartOutColl(new pat::PFParticleCollection);
  pat::PFParticleCollection pfPhotons = getPFPhotons(m_pfCandView, pi0Checker.getCheckIndices<decltype(recoPfPhotons)>());
  for(const auto& pfPhoton : pfPhotons) {
      m_patPfPartCtr++;
      pfPartOutColl->push_back(pfPhoton);
  }

  // process the photons
  std::auto_ptr<pat::PhotonCollection> patPhotonOutColl(new pat::PhotonCollection);
  auto photonColl = getPhotons(m_photonColl, pi0Checker.getCheckIndices<decltype(*m_photonColl)>());
  for(const auto& photon : photonColl) {
    m_photonCtr++;
    patPhotonOutColl->push_back(photon);
  }

  iEvent.put(patPhotonOutColl, "photons");
  iEvent.put(pfPartOutColl, "PFlowPhotons");
  iEvent.put(patConvOutColl, "convertedPhotons");
}

// ============================== GET PFPHOTONS ==============================
const pat::PFParticleCollection
OGB2BPhotonProducer::getPFPhotons(const edm::Handle<edm::View<reco::PFCandidate> >& pfCands,
                                       const std::vector<DoubleMassWindowRT>& piZeroResults)
{
  // NOTE: need the edm::View to have a RefToBase for constructing the pat::PFParticle, this makes this
  // function more messy than it should be.
  std::vector<AnnotatedT<const reco::PFCandidate> > annotatedCands;
  pat::PFParticleCollection photons;
  for(size_t iCand = 0; iCand < pfCands->size(); ++iCand) {
    const auto& cand = (*pfCands)[iCand];
    if(cand.particleId() != reco::PFCandidate::ParticleType::gamma) continue;

    m_pfCandCtr++;
    AnnotatedT<const reco::PFCandidate> annCand(&cand, getFlags(cand));
    annotatedCands.push_back(annCand);

    pat::PFParticle photon( pfCands->refAt(iCand) ); // construct from RefToBase to PFCandidate
    photons.push_back(photon);
  }

  setPiZeroFlags(annotatedCands, piZeroResults);

  for (size_t iCand = 0; iCand < photons.size(); ++iCand) {
    const auto& annCand = annotatedCands[iCand];
    auto& photon = photons[iCand];
    if (!checkFlags(annCand, m_pfPhotonFlags)) {
      annotate(photon, annCand);
    }
  }

  // now loop again over all patCandidates and remove all that have no flags
  auto patIt = photons.begin();
  while(patIt != photons.end()) {
    if(patIt->hasUserInt("flags")) ++patIt;
    else photons.erase(patIt);
  }

  return photons;
}

// ============================= GET PHOTONS ==============================
const pat::PhotonCollection
OGB2BPhotonProducer::getPhotons(const edm::Handle<reco::PhotonCollection>& photons,
                                     const std::vector<DoubleMassWindowRT>& piZeroResults)
{
  std::vector<AnnotatedT<const reco::Photon> > allPhotons;
  for(const auto& phot : *photons) {
    m_photonCandCtr++;
    allPhotons.push_back(AnnotatedT<const reco::Photon>(&phot, getFlags(phot)));
  }

  setPiZeroFlags(allPhotons, piZeroResults);
  removeFlagged(allPhotons, m_photonFlags);

  pat::PhotonCollection outPhotons;
  for (const auto & phot : allPhotons) {
    pat::Photon patPhoton(*(phot.object));
    annotate(patPhoton, phot);
    outPhotons.push_back(patPhoton);
  }

  return outPhotons;
}

// ============================== GET CONVERSIONS ==============================
const pat::CompositeCandidateCollection
OGB2BPhotonProducer::getConversions(const edm::Handle<reco::ConversionCollection>& conversions,
                                         const std::vector<DoubleMassWindowRT>& piZeroResults)
{
  pat::CompositeCandidateCollection collection;
  std::vector<AnnotatedT<const reco::Conversion> > convCands;

  // collect the flags for all conversions
  for(const auto& conv : *conversions) {
    m_recoConvCtr++;
    convCands.push_back( AnnotatedT<const reco::Conversion>(&conv, getFlags(conv)) );
  }

  setPiZeroFlags(convCands, piZeroResults);
  removeFlagged(convCands, m_convFlags);

  checkTrackSharing(convCands);

  for(const auto& conv : convCands) {
    if (conv.flags.test(7)) continue; // do not store overlapping conversion
    auto patCand = makePhotonCandidate(conv);
    annotate(patCand, conv);
    collection.push_back(patCand);
  }

  return collection;
}

// ============================== MAKE PHOTON CANDIDATE ==============================
pat::CompositeCandidate
OGB2BPhotonProducer::makePhotonCandidate(const AnnotatedT<const reco::Conversion>& conversion)
{
  pat::CompositeCandidate candidate;
  candidate.setP4( convertVector(conversion.object->refittedPair4Momentum()) );
  candidate.setVertex( conversion.object->conversionVertex().position() );

  const auto& convTracks = conversion.object->tracks();
  candidate.addUserData<reco::Track>("track0", *convTracks[0]);
  candidate.addUserData<reco::Track>("track1", *convTracks[1]);

  return candidate;
}

template<typename RecoType>
OGB2BPhotonProducer::bitsetT OGB2BPhotonProducer::getFlags(const RecoType& recoPart)
{
  bitsetT flags;
  for ( const auto bit : getFlagBits(recoPart)) flags.set(bit);
  return flags;
}

// ============================== GET CONVERSION FLAG BITS ==============================
const std::vector<unsigned short> OGB2BPhotonProducer::getFlagBits(const reco::Conversion& conv) const
{
  std::vector<unsigned short> flagBits;
  if(!m_conversionCutSel(conv)) flagBits.push_back(0);
  if(m_convAlgo != reco::Conversion::ConversionAlgorithm::undefined && m_convAlgo != conv.algo() ) flagBits.push_back(1);
  if(!checkConversionQuality(conv)) flagBits.push_back(2);

  // if the conversion does not have two tracks attatched to it, there is no point in doing the following checks!
  if(conv.tracks().size() == 2) {
    if(!checkTkVtxCompatibility(conv)) flagBits.push_back(3);
    if(!checkCompatibleInnerHits(conv)) flagBits.push_back(4);
  } else {
    flagBits.push_back(3); flagBits.push_back(4);
  }
  if(!checkHighPuritySubset(conv, *m_vertexColl.product())) flagBits.push_back(5);

  return flagBits;
}

const std::vector<unsigned short> OGB2BPhotonProducer::getFlagBits(const reco::Photon& photon) const
{
  std::vector<unsigned short> flagBits;
  if(!m_photonCutSel(photon)) flagBits.push_back(0);

  return flagBits;
}

const std::vector<unsigned short> OGB2BPhotonProducer::getFlagBits(const reco::PFCandidate& pfCand) const
// const std::vector<unsigned short> OGB2BPhotonProducer::getFlagBits(const pat::PFParticle& pfCand) const
{
  std::vector<unsigned short> flagBits;
  if(!m_pfCandCutSel(pfCand)) flagBits.push_back(0);

  return flagBits;
}

// ============================== CHECK TRACK VERTEX COMPATIBILITY ==============================
bool OGB2BPhotonProducer::checkTkVtxCompatibility(const reco::Conversion& conv) const
{
  std::array<std::vector<std::pair<double, unsigned short> >,2> vtxIdcs; // NOTE: the size has been checked prior to this!
  for(size_t iTk = 0; iTk < 2; ++iTk ) {
    const auto& track = conv.tracks()[iTk];
    for(unsigned short iVtx = 0; iVtx < m_vertexColl->size(); ++iVtx) {
      const auto& vertex = (*m_vertexColl)[iVtx];
      double dz = fabs(track->dz(vertex.position()));
      double dzErr = track->dzError();
      dzErr = sqrt(dzErr * dzErr + vertex.covariance(2,2));
      if(dz / dzErr > m_TkVtxCompSigma) continue;
      vtxIdcs[iTk].push_back(std::make_pair(dz, iVtx));
    }
    if(vtxIdcs[iTk].empty()) return false;

    // define a lambda function for ordering by the .first of the pair (aside from ignoring .second this is what the standard does!)
    // COULDDO: remove this since this is covered by the standard.
    auto ltLam = [] (const std::pair<double, unsigned short>& a, const std::pair<double, unsigned short>& b)
      { return a.first < b.first; };
    std::stable_sort(vtxIdcs[iTk].begin(), vtxIdcs[iTk].end(), ltLam);
  }

  if( vtxIdcs[0][0].second == vtxIdcs[1][0].second ||
      vtxIdcs[0][1].second == vtxIdcs[1][0].second ||
      vtxIdcs[0][0].second == vtxIdcs[1][1].second ) {
    return true;
  }

  return false;
}

// ============================== CHECK COMPATIBLE INNER HITS ==============================
bool OGB2BPhotonProducer::checkCompatibleInnerHits(const reco::Conversion& conv) const
{
  const reco::HitPattern& hitPatA = conv.tracks()[0]->hitPattern();
  const reco::HitPattern& hitPatB = conv.tracks()[1]->hitPattern();

  return foundCompatibleInnerHits(hitPatA, hitPatB) && foundCompatibleInnerHits(hitPatB, hitPatA);
}

// ============================== FOUND COMPATIBLE INNER HITS ==============================
bool OGB2BPhotonProducer::foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) const
{
  reco::HitPattern::HitCategory tkHits = reco::HitPattern::HitCategory::TRACK_HITS;

  size_t count{};
  uint32_t oldSubStr{};
  for(int iHit = 0; iHit < hitPatA.numberOfHits(tkHits) && count < 2; ++iHit) {
    uint32_t hitA = hitPatA.getHitPattern(tkHits, iHit);
    if(!hitPatA.validHitFilter(hitA) || !hitPatA.trackerHitFilter(hitA)) continue;

    if(hitPatA.getSubStructure(hitA) == oldSubStr && hitPatA.getLayer(hitA) == oldSubStr) continue;

    if(hitPatB.getTrackerMonoStereo(tkHits, hitPatA.getSubStructure(hitA), hitPatA.getLayer(hitA))) return true;

    oldSubStr = hitPatA.getSubStructure(hitA);
    count++;
  }

  return false;
}

// ============================== CHECK HIGH PURITY SUBSET ==============================
bool OGB2BPhotonProducer::checkHighPuritySubset(const reco::Conversion& conv, const reco::VertexCollection& vtxColl) const
{
  if(ChiSquaredProbability(conv.conversionVertex().chi2(), conv.conversionVertex().ndof()) < m_vertexChi2ProbCut) {
    return false;
  }

  size_t clVtxIdx = 0; // index of closest vertex to the conversion
  for(size_t iVtx = 0; iVtx < vtxColl.size(); ++iVtx) {
    if (conv.zOfPrimaryVertexFromTracks(vtxColl[iVtx].position()) <
        conv.zOfPrimaryVertexFromTracks(vtxColl[clVtxIdx].position()) ) {
      clVtxIdx = iVtx;
    }
  }

  for(const auto& trackRef : conv.tracks()) {
    // check the impact parameter w.r.t. the closest vertex found prior
    if (-trackRef->dxy(vtxColl[clVtxIdx].position()) * trackRef->charge() / trackRef->dxyError() < 0) return false;
    // chi2 of single tracks
    if (trackRef->normalizedChi2() > m_trackChi2Cut) return false;
    // dof freedom of single tracks
    if (trackRef->ndof() < m_trackMinNDOF) return false;
  }

  double minApp = conv.distOfMinimumApproach();
  if (minApp < m_minDistanceOfApproachMinCut || minApp > m_minDistanceOfApproachMaxCut) return false;

  return true;
}

// ============================== CHECK TRACK SHARING ==============================
void OGB2BPhotonProducer::checkTrackSharing(std::vector<AnnotatedT<const reco::Conversion> >& convs)
{
  if(convs.size() < 2) return; // no purpose in checking if there is only one conversion

  auto convLessChi2 = [](const AnnotatedT<const reco::Conversion>& c1, const AnnotatedT<const reco::Conversion>& c2)
    {return TMath::Prob(c1.object->conversionVertex().chi2(), c1.object->conversionVertex().ndof()) >
     TMath::Prob(c2.object->conversionVertex().chi2(), c2.object->conversionVertex().ndof());};

  std::sort(convs.begin(), convs.end(), convLessChi2); // sort by chi2 values of track vertices

  auto commonTrack = [](const AnnotatedT<const reco::Conversion>&c1, const AnnotatedT<const reco::Conversion>& c2)
    {
      for(const auto& tk1 : c1.object->tracks()) {
        for(const auto& tk2 : c2.object->tracks()) {
          if(tk1 == tk2) { return true; }
        }
      }
     return false;
    };

  for(size_t iCa = 0; iCa < convs.size() - 1; ++iCa) {
    for(size_t iCb = iCa + 1; iCb < convs.size(); ++iCb) {
      if(commonTrack(convs[iCa], convs[iCb])) {
        if (!convs[iCb].flags.test(7)) { // if it is not already marked as check if the better ranke is marked already
          if (convs[iCa].flags.test(7)) {
            convs[iCb].flags.set(8);
          } else {
            convs[iCb].flags.set(7);
          }
        }
      }
    }
  }
}

// ============================== REMOVE FLAGGED ==============================
template<typename RecoType>
void OGB2BPhotonProducer::removeFlagged(std::vector<AnnotatedT<const RecoType> >& coll, const bitsetT& flags)
{
  if (coll.empty()) return;
  auto collIt = coll.cbegin();
  while (collIt != coll.cend()) {
    if ( checkFlags(*collIt, flags) ) { collIt = coll.erase(collIt); }
    else { ++collIt; }
  }
}

// ============================== CREATE FLAGS ==============================
OGB2BPhotonProducer::bitsetT OGB2BPhotonProducer::createFlags(const std::string& flagStr)
{
  for (const char c : flagStr) {
    if (c != '1' && c != '0') {
      std::cerr << "Passed string for flags contains characters other than 0 and 1! Unsetting all flags!" << std::endl;
      return bitsetT{};
    }
  }

  if (flagStr.length() > nBits) {
    std::cerr << "Passed too many bits for the internally used bitset. Using only as many flags as can be fit into bitset!" << std::endl;
  }

  return bitsetT(flagStr);
}

// ============================== SET PI0 FLAGS ==============================
template<typename RecoType>
void OGB2BPhotonProducer::setPiZeroFlags(std::vector<AnnotatedT<const RecoType> >& cands,
                                              const std::vector<DoubleMassWindowRT>& pi0Results)
{
  for (const auto& res : pi0Results) {
    if (res.m_massWindow1) cands[res.m_index].flags.set(8);
    if (res.m_massWindow2) cands[res.m_index].flags.set(9);
  }
}

// ============================== begin / end Job ==============================
void OGB2BPhotonProducer::beginJob()
{
}

void OGB2BPhotonProducer::endJob()
{
  std::cout << "number of conversions (stored / presented):  " << m_patConvCtr << " / " << m_recoConvCtr << std::endl;
  std::cout << "number of PFlow photons (stored /presented): " << m_patPfPartCtr << " / " << m_pfCandCtr << std::endl;
  std::cout << "number of Photons: (stored / presented):     " << m_photonCtr << " / " << m_photonCandCtr << std::endl;
}

// ============================== fill descriptions ==============================
void OGB2BPhotonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(OGB2BPhotonProducer);
