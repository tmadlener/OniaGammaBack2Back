#ifndef __OniaGamma_PiZeroChecker_h_
#define __OniaGamma_PiZeroChecker_h_

// #include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
// #include "DataFormats/EgammaCandidates/interface/Conversion.h"
// #include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
// #include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"

#include <vector>
#include <unordered_map>
#include <typeinfo>
#include <iterator>
#include <algorithm>

#include "/afs/hephy.at/user/t/tmadlener/snippets/type_deduction_helper.h"

/** return type kind of thing */
struct PiZeroCheckReport {
  PiZeroCheckReport() = default;
  std::pair<size_t, size_t> overlap;
  bool narrow{false};
};

/** Class that takes different collections of photons (Conversions, PFlow, Photons) and checks if any combination
 * of two photons possibly originates from a \f$\pi^{0}\f$-decay.
 *
 * At the moment this is just a declaration....
 * At the moment planned to not be done in the PAT creation but in the analyzing stage
 * Still TODO
 */
class PiZeroChecker {
public:
  /** construction only possible from mass windows. WARNING: working with a lot of assumptions at the moment! */
  PiZeroChecker(const std::vector<double>& wideWindow, const std::vector<double>& narrowWindow);
  PiZeroChecker() = delete; /**< delete default constructor. */
  ~PiZeroChecker() = default; /**< nothing is allocated, hence default destructor will work. */

  /** add a collection that should be checked.
   * NOTE: somewhat misleading title due to variadic templatization
   */
  template<typename PATType>
  void addCollections(const std::vector<PATType>& collection);

  /** add an arbitrary number of collections that should be checked. */
  template<typename PATType, typename... PATTypes>
  void addCollections(const std::vector<PATType>& coll, PATTypes... others)
  {
    addCollections(coll);
    addCollections(others...);
  }

  /** check for possible pi0 */
  void check();

  template<typename PATType>
  std::vector<PiZeroCheckReport> getPiZeroReport() const;

private:
  /** internally used */
  struct CheckObject {
    template<typename PATType>
    CheckObject(const PATType& object) : p4(object.p4()) {}
    reco::Candidate::LorentzVector p4;
    std::vector<size_t> partnerInds;
    std::vector<bool> narrow;
  };

  /** internally used */
  struct TypeRange {
    TypeRange(std::vector<CheckObject>::const_iterator f) : first(f) {}
    std::vector<CheckObject>::const_iterator first; /**< iterator pointing to the first element. */
    std::vector<CheckObject>::const_iterator last; /**< iterator pointing to the element after the last!*/
    // bool operator<(const TypeRange& rhs) { return first < rhs.first; }
  };

  struct MassWindow {
    // MassWindow(double l, double h) : low(l), high(h) {;}
    double low, high;
    inline bool contains(double mass) { return (mass > low && mass < high); }
  };

  // ============================== MEMBER DATA ==============================
  std::vector<CheckObject> m_allCandidates;
  std::unordered_map<size_t, std::vector<TypeRange> > m_typeMap;
  MassWindow m_massWide;
  MassWindow m_massNarrow;
  std::vector<size_t> m_candidateInds; /**< vector containing the lower index of an pi0 candidate pair. */
};

// ============================== IMPLEMENTATION ==============================
PiZeroChecker::PiZeroChecker(const std::vector<double>& wideWindow, const std::vector<double>& narrowWindow)
{
  if (wideWindow.size() < 2 || narrowWindow.size() < 2) {
    std::cerr << "PiZeroChecker::PiZeroChecker(): passed vectors are too small! " << std::endl;
  }
  m_massWide.low = wideWindow[0]; m_massWide.high = wideWindow[1];
  m_massNarrow.low = narrowWindow[0]; m_massNarrow.high = narrowWindow[1];
}

template<typename PATType>
void PiZeroChecker::addCollections(const std::vector<PATType>& collection)
{
  if (collection.empty()) return;

  TypeRange range(m_allCandidates.end()); // this will be the first place where a new type starts
  for (const auto& obj : collection) { m_allCandidates.push_back(CheckObject(obj)); }
  range.last = m_allCandidates.end();

  m_typeMap[typeid(PATType).hash_code()].push_back(range);
  std::cout << m_typeMap.size() << std::endl;
}

void PiZeroChecker::check()
{
  if (m_allCandidates.size() < 2) {
    std::cout << "nothing to do in this event" << std::endl;
    return;
  }

  for (auto firstIt = m_allCandidates.begin(); firstIt != m_allCandidates.end() - 1; ++firstIt) {
    for (auto secIt = firstIt + 1; secIt != m_allCandidates.end(); ++secIt) {
      std::cout << "Checking combination: " << std::distance(m_allCandidates.begin(), firstIt) << " | " << std::distance(m_allCandidates.begin(), secIt) << std::endl;
      std::cout << "4-vectors: " << firstIt->p4 << " | " << secIt->p4 << std::endl;

      double invMass = (firstIt->p4 + secIt->p4).M();
      std::cout << "invMass: " << invMass << std::endl;

      if (m_massWide.contains(invMass)) {
        std::cout << "contained in wide mass window" << std::endl;

        size_t idx1 = std::distance(m_allCandidates.begin(), firstIt);
        size_t idx2 = std::distance(m_allCandidates.begin(), secIt);
        firstIt->partnerInds.push_back(idx2);
        secIt->partnerInds.push_back(idx1);
        bool narrow = m_massNarrow.contains(invMass);
        firstIt->narrow.push_back(narrow);
        secIt->narrow.push_back(narrow);

        m_candidateInds.push_back(idx1);

        if (m_massNarrow.contains(invMass)) {
          std::cout << "contained in narrow mass window" << std::endl;
        }
      }
    }
    std::cout << "-------------------------" << std::endl;
  }
}

template<typename PATType>
std::vector<PiZeroCheckReport> PiZeroChecker::getPiZeroReport() const
{
  std::vector<PiZeroCheckReport> reports;

  const auto typeRangeIt = m_typeMap.find(typeid(PATType).hash_code());
  if(typeRangeIt != m_typeMap.end()) {
    for (const auto& range : typeRangeIt->second) {
      size_t firstIdx = std::distance(m_allCandidates.cbegin(), range.first);
      size_t lastIdx = std::distance(m_allCandidates.cbegin(), range.last);
      for (const size_t& idx : m_candidateInds) {
        if (idx >= firstIdx && idx < lastIdx) {
          // PiZeroCheckReport report;
          // TODO
        }
      }
    }
  }

  return reports;
}

#endif
