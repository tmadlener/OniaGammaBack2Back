#include "OniaGammaBack2Back/PhotonAnalyzer/interface/PhotonAnalyzerPAT.h"

// #include "OniaGammaBack2Back/LowEnergyPhotons/interface/PiZeroChecker.h"
// #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include <TVector3.h>

#include <limits>
#include <bitset>
#include <algorithm>
#include <typeinfo>

PhotonAnalyzerPAT::PhotonAnalyzerPAT(const edm::ParameterSet& iConfig) :
  m_rootFileName(iConfig.getParameter<std::string>("outputFileName")),
  m_pfCollTok( consumes<pat::PFParticleCollection>(iConfig.getParameter<edm::InputTag>("pFlowPhotons")) ),
  m_photonCollTok( consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons")) ),
  m_convCollTok( consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("conversions")) ),
  m_diMuonTok( consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dimuons")) ),
  m_muonTok( consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")) ),
  m_primVertTok( consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices")) ),
  m_triggerResultsTok( consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults")) ),
  m_oniaMassMax( iConfig.getParameter<double>("oniaMassMax")),
  m_oniaMassMin( iConfig.getParameter<double>("oniaMassMin")),
  m_storeOnlyBest( iConfig.getParameter<bool>("storeOnlyBest")),
  m_outputFile(nullptr),
  m_tree(nullptr)
{

}

PhotonAnalyzerPAT::~PhotonAnalyzerPAT()
{
  // TODO
}

void PhotonAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  m_rootVariables.clear();

  // static unsigned nEv;
  // std::cout << "== " << ++nEv << " ==" << std::endl; // DEBUGGING only
  iEvent.getByToken(m_pfCollTok, m_pfCandColl);
  iEvent.getByToken(m_photonCollTok, m_photonColl);
  iEvent.getByToken(m_convCollTok, m_convColl);
  iEvent.getByToken(m_diMuonTok, m_diMuonColl);
  iEvent.getByToken(m_muonTok, m_muonColl);
  iEvent.getByToken(m_primVertTok, m_primVertices);
  iEvent.getByToken(m_triggerResultsTok, m_triggerResults);

  m_rootVariables.run = iEvent.id().run();
  m_rootVariables.event = iEvent.id().event();
  m_rootVariables.lumiblock = iEvent.id().luminosityBlock();

  m_rootVariables.triggerBits = getTriggerBits(iEvent);
  if (m_primVertices.isValid()) m_rootVariables.nPrimVertices = m_primVertices->size();

  size_t nDiMuons{};
  if (m_diMuonColl.isValid()) {
    nDiMuons = getDiMuonInfo(*m_diMuonColl);
  }
  if (!nDiMuons && m_muonColl.isValid()) {
    getSingleMuonInfo(*m_muonColl); // NOTE: NO EFFECT! (see .h file)
  }

  getPhotonInfo(*m_photonColl);
  getPhotonInfo(*m_pfCandColl);
  getPhotonInfo(*m_convColl);

  m_tree->Fill();
}

void PhotonAnalyzerPAT::beginJob()
{
  m_outputFile = new TFile(m_rootFileName.c_str(), "RECREATE");
  m_tree = new TTree("photonAnalysis", "Photon Analysis Content");
  if (!m_rootVariables.createSetBranches(m_tree)) {
    std::cerr << "Could not create the Branches in root file: " << m_rootFileName << std::endl;
  }
}

void PhotonAnalyzerPAT::endJob()
{
  if(m_outputFile && m_tree) {
    m_outputFile->cd();
    m_tree->Write();
    m_outputFile->Close();
  }

}

void PhotonAnalyzerPAT::fillDescription(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

unsigned PhotonAnalyzerPAT::getTriggerBits(const edm::Event& iEvent)
{
  if (!m_triggerResults.isValid()) {
    std::cerr << "PhotonAnalyzerPAT::getTriggerResults(): Handle to TriggerResults is not valid" << std::endl;
    return std::numeric_limits<unsigned>::max(); // set all bits in this case
  }

  constexpr int storableBits = std::numeric_limits<unsigned>::digits;
  std::bitset<storableBits> bits;

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*m_triggerResults);

  for (const char v : std::vector<char>{'1', '2'}) { // different versions of triggers
    for (size_t iName = 0; iName < CheckTriggerNames.size(); ++iName) {
      if (iName >= storableBits) { // probably unnecessary safety check
        std::cerr << "cannot store Trigger: " << CheckTriggerNames[iName] << " because there are not enough bits" << std::endl;
        continue;
      }

      unsigned bit = triggerNames.triggerIndex( CheckTriggerNames[iName] + v );
      if (bit < m_triggerResults->size()) {
        if (m_triggerResults->accept(bit) && !m_triggerResults->error(bit)) {
          bits.set(iName);
        }
      }
    } // end loop trigger names
  } // end loop versions

  return bits.to_ulong();
}

size_t PhotonAnalyzerPAT::getDiMuonInfo(const pat::CompositeCandidateCollection& dimuonColl)
{
  for (const auto& cand : dimuonColl) { // loop over all candidates
    std::cout << "in getDiMuonInfo " << cand.mass() << " " << cand.charge() << std::endl;
    if (cand.mass() > m_oniaMassMin && cand.mass() < m_oniaMassMax && cand.charge() == 0) {
      std::cout << "starting to store a dimuon candidate" << std::endl;
      m_rootVariables.dimuonP4s.push_back(cand.p4());

      auto pMuonP4 = cand.daughter("muon1")->p4();
      auto nMuonP4 = cand.daughter("muon2")->p4();
      if (cand.daughter("muon1")->charge() < 0) {
        std::swap(pMuonP4, nMuonP4);
      }
      m_rootVariables.pMuonP4s.push_back(pMuonP4);
      m_rootVariables.nMuonP4s.push_back(nMuonP4);

      m_rootVariables.massErr.push_back(cand.userFloat("MassErr"));
      m_rootVariables.vProb.push_back(cand.userFloat("vProb"));
      double dca = cand.hasUserFloat("DCA") ? cand.userFloat("DCA") : -1.;
      m_rootVariables.DCA.push_back(dca);

      double ppdlPV = cand.userFloat("ppdlPV");
      m_rootVariables.ppdlPV.push_back(ppdlPV);
      m_rootVariables.ppdlErrPV.push_back(cand.userFloat("ppdlErrPV"));

      double ppdlBS = cand.userFloat("ppdlBS");
      m_rootVariables.ppdlBS.push_back(ppdlBS);
      m_rootVariables.ppdlErrBS.push_back(cand.userFloat("ppdlErrBS"));

      m_rootVariables.cosAlpha.push_back(cand.userFloat("cosAlpha"));
      m_rootVariables.charge.push_back(cand.charge()); // NOTE: useless info since everything else is cut!

      TVector3 pperp(cand.px(), cand.py(), 0);
      m_rootVariables.lxyPV.push_back(ppdlPV * pperp.Perp() / cand.mass());
      m_rootVariables.lxyBS.push_back(ppdlBS * pperp.Perp() / cand.mass());

      const auto* vertexPos = cand.userData<math::XYZPoint>("vertexPosition");
      if(!vertexPos) vertexPos = new math::XYZPoint{};
      m_rootVariables.dimuonVertices.push_back(*vertexPos);

      if (m_storeOnlyBest) break; // stop here if only the best dimuon candidate is desired
    }
  } // end dimuon cand loop

  return m_rootVariables.dimuonP4s.size(); // should be the same for all variables
}

size_t PhotonAnalyzerPAT::getSingleMuonInfo(const pat::MuonCollection& muonColl)
{
  // TODO: Implement this, but at the moment I see no point in doing so!
  return 0;
}

template<typename PATType>
size_t PhotonAnalyzerPAT::getPhotonInfo(const std::vector<PATType>& photonColl)
{
  size_t preSize = m_rootVariables.photonP4s.size();
  short cat = getPhotonCat(photonColl);

  for (const auto& cand : photonColl) {
    m_rootVariables.photonP4s.push_back(cand.p4());
    m_rootVariables.photonCat.push_back(cat);
    m_rootVariables.photonVertices.push_back(cand.vertex());
  }

  return m_rootVariables.photonP4s.size() - preSize;
}

template<typename PATType>
short PhotonAnalyzerPAT::getPhotonCat(const std::vector<PATType>& coll) const
{
  if (typeid(PATType) == typeid(pat::Photon)) return 1;
  if (typeid(PATType) == typeid(pat::PFParticle)) return 2;
  if (typeid(PATType) == typeid(pat::CompositeCandidate)) return 3;

  return -1;
}

// define plug in
DEFINE_FWK_MODULE(PhotonAnalyzerPAT);
