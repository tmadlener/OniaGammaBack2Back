#ifndef __OniaGamma_PhotonAnalysisRootVariables_h_
#define __OniaGamma_PhotonAnalysisRootVariables_h_

#include <vector>

#include <TTree.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

/**
 * struct containing all the variables that should be written to the root file.
 *
 * Provides 'sane' defaults and manages the TBranch initialization (hopefully)
 */
struct PARootVariables {
  bool createSetBranches(TTree* tree, bool mc = false); /**< create the branches and set the addresses to the variables. */

  void clear(); /**< reset everything for the recording of a new event. */

  // ============================== MEMBER DATA ==============================

  std::vector<math::XYZTLorentzVector> photonP4s{}; /**< the photon four momenta. */
  std::vector<math::XYZTLorentzVector> dimuonP4s{}; /**< the dimuon four momemta. */
  std::vector<math::XYZTLorentzVector> pMuonP4s{}; /**< positively charged muon four momenta. */
  std::vector<math::XYZTLorentzVector> nMuonP4s{}; /**< negatively charged muon four momenta. */

  std::vector<math::XYZPoint> photonVertices{}; /**< the photon vertices. */
  std::vector<math::XYZPoint> dimuonVertices{}; /**< the dimuon vertices. */
  std::vector<short> photonCat{}; /**< the photon collection. -1: unknown, 1 - pat::Photon, 2 - pat::PFParticle, 3 - pat::CompositeCandidate */

  unsigned event{}; /**< the event number. */
  unsigned lumiblock{}; /**< the luminosity block. */
  unsigned run{}; /**< the run number. */

  unsigned triggerBits{}; /**< the trigger bits. */
  unsigned nPrimVertices{}; /**< the number of primary vertices. */

  std::vector<double> massErr{}; /**< the error of the mass estimation. */
  std::vector<double> vProb{}; /**< the p-value of the vertex fit. (TOOD: check if this is true) */
  std::vector<double> DCA{}; /**< the distance of closest approach. */
  std::vector<double> ppdlPV{}; /**< TODO: documentation. */
  std::vector<double> ppdlErrPV{}; /**< TODO: documentation. */
  std::vector<double> ppdlBS{}; /**< TODO: documentation. */
  std::vector<double> ppdlErrBS{}; /**< TODO: documentation. */
  std::vector<double> cosAlpha{}; /**< TODO: documentation. */
  std::vector<double> charge{}; /** the dimuon charge. */
  std::vector<double> lxyPV{}; /**< TODO: documentation. */
  std::vector<double> lxyBS{}; /**< TODO: documentation. */

  bool m_useMC{false};
};

bool PARootVariables::createSetBranches(TTree* tree, bool mc)
{
  m_useMC = mc;

  if(!tree->Branch("photonP4s", &photonP4s)) return false;
  if(!tree->Branch("dimuonP4", &dimuonP4s)) return false;
  if(!tree->Branch("nMuonP4s", &nMuonP4s)) return false;
  if(!tree->Branch("pMuonP4s", &pMuonP4s)) return false;

  if(!tree->Branch("dimuonVertices", &dimuonVertices)) return false;
  if(!tree->Branch("photonVertices", &photonVertices)) return false;
  if(!tree->Branch("photonCat", &photonCat)) return false;

  if(!tree->Branch("event", &event)) return false;
  if(!tree->Branch("lumiblock", &lumiblock)) return false;
  if(!tree->Branch("run", &run)) return false;

  if(!tree->Branch("triggerBits", &triggerBits)) return false;
  if(!tree->Branch("nPrimVertics", &nPrimVertices)) return false;

  if(!tree->Branch("massErr", &massErr)) return false;
  if(!tree->Branch("vProb", &vProb)) return false;
  if(!tree->Branch("DCA", &DCA)) return false;
  if(!tree->Branch("ppdlPV", &ppdlPV)) return false;
  if(!tree->Branch("ppdlErrPV", &ppdlErrPV)) return false;
  if(!tree->Branch("ppdlBS", &ppdlBS)) return false;
  if(!tree->Branch("ppdlErrBS", &ppdlErrBS)) return false;
  if(!tree->Branch("cosAlpha", &cosAlpha)) return false;
  if(!tree->Branch("charge", &charge)) return false;
  if(!tree->Branch("lxyPV", &lxyPV)) return false;
  if(!tree->Branch("lxyBS", &lxyBS)) return false;

  return true;
}

void PARootVariables::clear()
{
  photonP4s.clear();
  dimuonP4s.clear();
  pMuonP4s.clear();
  nMuonP4s.clear();

  dimuonVertices.clear();
  photonVertices.clear();
  photonCat.clear();

  event = lumiblock = run = 0; // not entirely the best style
  triggerBits = 0;
  nPrimVertices = 0;

  massErr.clear();
  vProb.clear();
  DCA.clear();
  ppdlPV.clear();
  ppdlErrPV.clear();
  ppdlBS.clear();
  ppdlErrBS.clear();
  cosAlpha.clear();
  charge.clear();
  lxyPV.clear();
  lxyBS.clear();
}

#endif
