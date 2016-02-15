#include "OniaGammaBack2Back/LowEnergyPhotons/interface/ConversionPhotonProducer.h"

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
#include <bitset>
#include <map>

// helper stuff for dev
#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"

// ============================== constructor / destructor ==============================
ConversionPhotonProducer::ConversionPhotonProducer(const edm::ParameterSet& iConfig) :
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
  m_trackMinNDOF(iConfig.getParameter<double>("trackMinNDOF"))
{
  std::string algo = iConfig.getParameter<std::string>("convAlgorithm");
  // convert the returned int into the enum-value in the constructor to avoid having to do it later
  m_convAlgo = (reco::Conversion::ConversionAlgorithm)StringToEnumValue<reco::Conversion::ConversionAlgorithm>(algo);

  // store the desired qualities (input by string) into the internally used vector
  std::vector<std::string> qualities = iConfig.getParameter<std::vector<std::string> >("convQuality");
  for(const auto& qual : qualities) {
    m_convQualities.push_back( (reco::Conversion::ConversionQuality) StringToEnumValue<reco::Conversion::ConversionQuality>(qual) );
  }

  produces<pat::CompositeCandidateCollection>("convertedPhotons");
  produces<pat::PhotonCollection>("photons");
  produces<pat::PFParticleCollection>("PFlowPhotons");
}

ConversionPhotonProducer::~ConversionPhotonProducer()
{
  // TODO
}

// ============================== produce ==============================
void ConversionPhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // get the data to the handles (COULDDO: put this into a separate function?)
  iEvent.getByToken(m_convCollTok, m_convColl);
  iEvent.getByToken(m_photonCollTok, m_photonColl);
  iEvent.getByToken(m_pfCandViewTok, m_pfCandView);
  iEvent.getByToken(m_beamSpotTok, m_beamSpotHand);
  iEvent.getByToken(m_vertexCollTok, m_vertexColl);

  // process the conversions
  std::auto_ptr<pat::CompositeCandidateCollection> patConvOutColl(new pat::CompositeCandidateCollection);
  auto convColl = getConversions(m_convColl); // need this to be able to modify the elements of the vector
  for(const auto& patConv : convColl) {
    m_patConvCtr++;
    patConvOutColl->push_back(patConv);
  }
  iEvent.put(patConvOutColl, "convertedPhotons");

  // process the PFCandidates
  std::auto_ptr<pat::PFParticleCollection> pfPartOutColl(new pat::PFParticleCollection);
  const pat::PFParticleCollection pfPhotons = getPFPhotons(m_pfCandView);
  for(const auto& pfPhoton : pfPhotons) {
    m_patPfPartCtr++;
    pfPartOutColl->push_back(pfPhoton);
  }
  iEvent.put(pfPartOutColl, "PFlowPhotons");

  // process the photons
  std::auto_ptr<pat::PhotonCollection> patPhotonOutColl(new pat::PhotonCollection);
  auto photonColl = getPhotons(m_photonColl);
  for(const auto& photon : photonColl) {
    m_photonCtr++;
    patPhotonOutColl->push_back(photon);
  }
  iEvent.put(patPhotonOutColl, "photons");

}

// ============================== GET PFPHOTONS ==============================
const pat::PFParticleCollection
ConversionPhotonProducer::getPFPhotons(const edm::Handle<edm::View<reco::PFCandidate> >& pfCands)
{
  pat::PFParticleCollection photons;
  for(size_t iCand = 0; iCand < pfCands->size(); ++iCand) {
    const auto& cand = (*pfCands)[iCand];
    if(cand.particleId() != reco::PFCandidate::ParticleType::gamma) continue;
    pat::PFParticle photon( pfCands->refAt(iCand) ); // construct from RefToBase to PFCandidate
    annotate(photon, cand);
    photons.push_back(photon);
  }
  return photons;
}

// ============================= GET PHOTONS ==============================
const pat::PhotonCollection
ConversionPhotonProducer::getPhotons(const edm::Handle<reco::PhotonCollection>& photons)
{
  pat::PhotonCollection outPhotons;
  for(const auto& phot : *photons) {
    pat::Photon patPhoton(phot);
    annotate(patPhoton, phot);
    outPhotons.push_back(patPhoton);
  }
  return outPhotons;
}

// ============================== GET CONVERSIONS ==============================
const pat::CompositeCandidateCollection
ConversionPhotonProducer::getConversions(const edm::Handle<reco::ConversionCollection>& conversions)
{
  pat::CompositeCandidateCollection collection;

  for(const auto& conv : *conversions) {
    pat::CompositeCandidate cand = makePhotonCandidate(conv);
    annotate(cand, conv);
    collection.push_back(cand);
  }

  return collection;
}

// ============================== MAKE PHOTON CANDIDATE ==============================
pat::CompositeCandidate
ConversionPhotonProducer::makePhotonCandidate(const reco::Conversion& conversion)
{
  pat::CompositeCandidate candidate;
  candidate.setP4( convertVector(conversion.refittedPair4Momentum()) );
  candidate.setVertex( conversion.conversionVertex().position() );

  const auto& convTracks = conversion.tracks();
  candidate.addUserData<reco::Track>("track0", *convTracks[0]);
  candidate.addUserData<reco::Track>("track1", *convTracks[1]);

  return candidate;
}

// ============================== ANNOTATE PAT CANDIDATE  ==============================
template<typename PatType, typename RecoType>
void ConversionPhotonProducer::annotate(PatType& patPart, const RecoType& recoPart) const
{
  const unsigned short nFlags = 10;
  std::bitset<nFlags> flags;
  for (const auto& bit : getFlagBits(recoPart)) flags.set(bit);
  patPart.addUserInt("flags", flags.to_ulong());
}

// ============================== GET CONVERSION FLAG BITS ==============================
const std::vector<unsigned short> ConversionPhotonProducer::getFlagBits(const reco::Conversion& conv) const
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
  // TODO pi0 stuff

  return flagBits;
}

const std::vector<unsigned short> ConversionPhotonProducer::getFlagBits(const reco::Photon& photon) const
{
  std::vector<unsigned short> flagBits;
  if(!m_photonCutSel(photon)) flagBits.push_back(0);

  return flagBits;
}

const std::vector<unsigned short> ConversionPhotonProducer::getFlagBits(const reco::PFCandidate& pfCand) const
{
  std::vector<unsigned short> flagBits;
  if(!m_pfCandCutSel(pfCand)) flagBits.push_back(0);

  return flagBits;
}

// ============================== CHECK TRACK VERTEX COMPATIBILITY ==============================
bool ConversionPhotonProducer::checkTkVtxCompatibility(const reco::Conversion& conv) const
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
bool ConversionPhotonProducer::checkCompatibleInnerHits(const reco::Conversion& conv) const
{
  const reco::HitPattern& hitPatA = conv.tracks()[0]->hitPattern();
  const reco::HitPattern& hitPatB = conv.tracks()[1]->hitPattern();

  return foundCompatibleInnerHits(hitPatA, hitPatB) && foundCompatibleInnerHits(hitPatB, hitPatA);
}

// ============================== FOUND COMPATIBLE INNER HITS ==============================
bool ConversionPhotonProducer::foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) const
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
bool ConversionPhotonProducer::checkHighPuritySubset(const reco::Conversion& conv, const reco::VertexCollection& vtxColl) const
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

// ============================== begin / end Job ==============================
void ConversionPhotonProducer::beginJob()
{
}

void ConversionPhotonProducer::endJob()
{
  std::cout << "number of conversions: " << m_patConvCtr << std::endl;
  std::cout << "number of PFlow photons: " << m_patPfPartCtr << std::endl;
  std::cout << "number of Photons: " << m_photonCtr << std::endl;
}

// ============================== fill descriptions ==============================
void ConversionPhotonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ConversionPhotonProducer);
