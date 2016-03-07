#ifndef __OniaGamma_PhotonAnalysisRootVariables_h_
#define __OniaGamma_PhotonAnalysisRootVariables_h_

#include <TTree.h>
#include <TLorentzVector.h>

/**
 * struct containing all the variables that should be written to the root file.
 *
 * Provides 'sane' defaults and manages the TBranch initialization (hopefully)
 */
struct PARootVariables {
  bool createSetBranches(TTree* tree); /**< create the branches and set the addresses to the variables. */

  TLorentzVector photonP{}; /**< the photon four momentum. */

  // unsigned short sourceCat{}; /**< the source category. 1 - pat::Photon, 2 - pat::PFParticle, 3 - pat::CompositeCandidate, 0 - not set. */
  // unsigned nConversions{}; /**< number of conversions in the event. */
  // unsigned nPhotons{}; /**< number of pat::Photons in the event. */
  // unsigned nPFPhotons{}; /**< number of pat::PFParticle photons in the event. */
};

bool PARootVariables::createSetBranches(TTree* tree)
{
  if(!tree->Branch("photonP", &photonP)) return false;
  // if(!tree->Branch("sourceCategory", &sourceCat)) return false;
  // if(!tree->Branch("nConversions", &nConversions)) return false;
  // if(!tree->Branch("nPhotons", &nPhotons)) return false;
  // if(!tree->Branch("nPFPhotons", &nPFPhotons)) return false;


  return true;
}

#endif
