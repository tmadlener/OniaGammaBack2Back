#include "OniaGammaBack2Back/PhotonAnalyzer/interface/PhotonAnalyzerPAT.h"

#include "OniaGammaBack2Back/LowEnergyPhotons/interface/PiZeroChecker.h"
// #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

PhotonAnalyzerPAT::PhotonAnalyzerPAT(const edm::ParameterSet& iConfig) :
  m_rootFileName(iConfig.getParameter<std::string>("outputFileName")),
  m_pfCollTok( consumes<pat::PFParticleCollection>(iConfig.getParameter<edm::InputTag>("pFlowPhotons")) ),
  m_photonCollTok( consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons")) ),
  m_convCollTok( consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("conversions")) ),
  m_pi0WideWindow( iConfig.getParameter<std::vector<double> >("pi0WideWindow")),
  m_pi0NarrowWindow( iConfig.getParameter<std::vector<double> >("pi0NarrowWindow")),
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
  iEvent.getByToken(m_pfCollTok, m_pfCandColl);
  iEvent.getByToken(m_photonCollTok, m_photonColl);
  iEvent.getByToken(m_convCollTok, m_convColl);

  PiZeroChecker pi0Checker{m_pi0WideWindow, m_pi0NarrowWindow};
  pi0Checker.addCollections(*m_pfCandColl, *m_photonColl, *m_convColl);
  pi0Checker.check();

  const auto pi0Reports = pi0Checker.getPiZeroReport<pat::PFParticle>(); // get the pi0 stuff for PFParticles

  // m_rootVariables.nConversions = m_convColl->size();
  // m_rootVariables.nPhotons = m_photonColl->size();
  // m_rootVariables.nPFPhotons = m_pfCandColl->size();

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

// define plug in
DEFINE_FWK_MODULE(PhotonAnalyzerPAT);
