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
  m_conversionCutSel( iConfig.getParameter<std::string>("convSelection") ),
  m_TkVtxCompSigma( iConfig.getParameter<double>("tkVtxCompSigma"))
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

  // collect some counting variables
  m_pfPhotonsCtr.push_back(pfPhotons.size());
  m_convCtr.push_back(m_convColl->size());
  m_photonCtr.push_back(m_photonColl->size());
  m_pfCandCtr.push_back(m_pfCandView->size());
  double ratio = (double) m_pfPhotonsCtr.back() / (double) m_pfCandCtr.back();
  m_pfPhotonRatio.push_back(ratio);

}

// ============================== GET PFPHOTONS ==============================
const pat::PFParticleCollection
ConversionPhotonProducer::getPFPhotons(const edm::Handle<edm::View<reco::PFCandidate> >& pfCands)
{
  pat::PFParticleCollection photons;
  for(size_t iCand = 0; iCand < pfCands->size(); ++iCand) {
    pat::PFParticle photon( pfCands->refAt(iCand) ); // construct from RefToBase to PFCandidate
    photons.push_back(photon);
  }
  return photons;
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

// ============================== ANNOTATE CONVERSION ==============================
void ConversionPhotonProducer::annotate(pat::CompositeCandidate& patConv, const reco::Conversion& recConv) const
{
  patConv.addUserInt("photon-number", m_patConvCtr); // probably use-less but still for checking its nice

  const unsigned short nFlags = 10; // defined at prominent place for dev (10 should be enough right now)
  std::bitset<nFlags> flags; // create a bitset with all flags set to false
  // set all the flags for this conversion
  for(const auto& bit : getConversionFlagBits(recConv) ) flags.set(bit);
  patConv.addUserInt("flags", flags.to_ulong());

  // dev output
  std::cout << "annotate_flags " << flags << " " << flags.to_ulong() << std::endl;
}

// ============================== GET CONVERSION FLAG BITS ==============================
const std::vector<unsigned short> ConversionPhotonProducer::getConversionFlagBits(const reco::Conversion& conv) const
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
  // TODO pi0 stuff

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
// ============================== PRINT CONVERSION INFO ==============================
std::string ConversionPhotonProducer::printConversionInfo(const reco::Conversion& conv)
{
  std::stringstream str{};
  str << "pair momentum: " << conv.refittedPairMomentum() << std::endl;
  str << "invariant mass: " << conv.pairInvariantMass() << std::endl;
  str << "energy: " << conv.refittedPair4Momentum().E() << std::endl;

  if ( ConversionTools::isGoodConversion(conv, m_beamSpotHand->position()) ) {
    str << "good conversion\n";
  } else {
    unsigned failbits{}; // checking what failed in the check by default settings

    const reco::Vertex& vtx = conv.conversionVertex();

    str << "valid Vertex: " << vtx.isValid() << "\n";
    str << "fit probability: " << TMath::Prob( vtx.chi2(), vtx.ndof() ) << "\n";

    if(!vtx.isValid()) failbits += 1;
    if( TMath::Prob( vtx.chi2(), vtx.ndof()) < 1e-6 ) failbits += 2;

    math::XYZVector mom(conv.refittedPairMomentum());
    double dbsx = vtx.x() - m_beamSpotHand->position().x();
    double dbsy = vtx.y() - m_beamSpotHand->position().y();
    double lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
    str << "lxy: " << lxy << "\n";
    if( lxy < 2.0 ) failbits += 4;

    const std::vector<uint8_t>& nHitsBeforeVtx = conv.nHitsBeforeVtx();
    std::vector<unsigned int> nHitsVec;
    for(const auto& hits : nHitsBeforeVtx) { nHitsVec.push_back((unsigned int) hits);}
    str << "hits before vertex: " << nHitsVec << "\n";

    for(const auto& hits : nHitsBeforeVtx ) {
      if(hits > 1 ) {
        failbits += 8;
        break;
      }
    }
    str << "failbits " << failbits << "\n";
  }

  return str.str();
}

// ============================== begin / end Job ==============================
void ConversionPhotonProducer::beginJob()
{
}

void ConversionPhotonProducer::endJob()
{
  std::cout << "number of conversions: " << m_patConvCtr << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "average number of PFCandidates assigned with a photon label " << mean(m_pfPhotonsCtr) << std::endl;
  std::cout << "average number of PFCandidates " << mean(m_pfCandCtr) << std::endl;
  std::cout << "average number of conversions " << mean(m_convCtr) << std::endl;
  std::cout << "average number of photons " << mean(m_photonCtr) << std::endl;
  std::cout << "mean ratio of PFPhotons in PFCandidates " << mean(m_pfPhotonRatio) << std::endl;
  std::cout << "average number of overlaps btwn PFCandidates and Photons " << mean(m_pfPhotonOverlapCtr) << std::endl;
  std::cout << "average number of overlaps btwn PFCandidates and Conversions " << mean(m_pfConversionOverlapCtr) << std::endl;
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
