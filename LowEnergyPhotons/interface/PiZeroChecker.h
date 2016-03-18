#ifndef __OniaGamma_PiZeroChecker_h_
#define __OniaGamma_PiZeroChecker_h_

// #include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
// #include "DataFormats/EgammaCandidates/interface/Conversion.h"
// #include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
// #include "DataFormats/EgammaCandidates/interface/Photon.h"
// #include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"

#include <vector>
// #include <unordered_map>
#include <map>
#include <typeinfo>
#include <iterator>
#include <algorithm>
#include <bitset>

#include "/afs/hephy.at/user/t/tmadlener/snippets/type_deduction_helper.h"
#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"
#include "/afs/hephy.at/user/t/tmadlener/snippets/map_stuff.h"

/** return type kind of thing */
struct PiZeroCheckReport {
  PiZeroCheckReport() = default;
  // std::pair<size_t, size_t> overlap;
  size_t index;
  bool narrow{false};
  // std::bitset<3> collection{}; /**< 0 - PF, 1 - Conversion, 2 - Photon. */
  bool sameCollection{true}; /**< is partner in same collection? */
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
    bool veto() const { return !partnerInds.empty(); }
    bool veto_narrow() const { return /*any_equal(narrow, [](const bool& b) { return b == true; });*/ any(narrow); }
  };

  /** internally used */
  struct TypeRange {
    TypeRange(std::vector<CheckObject>::const_iterator f) : first(f) {}
    TypeRange() = default; /**< needed for operating with unordered_map. */
    std::vector<CheckObject>::const_iterator first; /**< iterator pointing to the first element. */
    std::vector<CheckObject>::const_iterator last; /**< iterator pointing to the element after the last!*/
    // bool operator<(const TypeRange& rhs) { return first < rhs.first; }
    /** check if the passed index is within this range.
     * NOTE: the casting is done to suppress compiler warnings, which is seldomly a good idea. */
    bool contains(size_t idx, const std::vector<CheckObject>& vec) const {
      return (static_cast<size_t>(std::distance(vec.cbegin(), first)) <= idx &&
              static_cast<size_t>(std::distance(vec.cbegin(), last)) > idx);
    }
  };

  struct MassWindow {
    // MassWindow(double l, double h) : low(l), high(h) {;}
    double low, high;
    inline bool contains(double mass) { return (mass > low && mass < high); }
  };

  // ============================== MEMBER DATA ==============================
  // std::vector<CheckObject> m_allCandidates;
  // std::unordered_map<size_t, TypeRange> m_typeMap;
  std::multimap<size_t, CheckObject> m_allCandidates;
  MassWindow m_massWide;
  MassWindow m_massNarrow;
  mutable std::vector<size_t> m_candidateInds; /**< cache all the candidateInds. */
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

// ============================== ADD COLLECTIONS ==============================
template<typename PATType>
void PiZeroChecker::addCollections(const std::vector<PATType>& collection)
{
  std::cout << "adding collection: " << typeid(PATType).name() << std::endl;
  if (collection.empty()) return;

  size_t typeHash = typeid(PATType).hash_code();
  for (const auto& obj : collection) { m_allCandidates.insert(std::make_pair(typeHash, CheckObject(obj))); }

  // bool firstAddition = m_allCandidates.empty(); // empty vectors iterator gets invalidated after first push_back
  // TypeRange range(m_allCandidates.cend()); // this will be the first place where a new type starts
  // for (const auto& obj : collection) { m_allCandidates.push_back(CheckObject(obj)); }
  // range.last = m_allCandidates.cend();
  // if (firstAddition) range.first = m_allCandidates.cbegin(); // set appropriate iterator for first collection

  // std::cout << "range: " << std::distance(m_allCandidates.cbegin(), range.first) << " " << std::distance(m_allCandidates.cbegin(), range.last) <<  std::endl;
  // std::cout << "distance: " << std::distance(range.first, range.last) << std::endl;
  // m_typeMap[typeid(PATType).hash_code()] = range;
  // // std::cout << m_typeMap.size() << std::endl;
}

// ============================== CHECK ==============================
void PiZeroChecker::check()
{
  if (m_allCandidates.size() < 2) {
    std::cout << "nothing to do in this event" << std::endl;
    return;
  }

  // NOTE: can't use iterator offsets since they are not defined on forward_iterators!
  // Simply loop over all candidates and make sure to not combine identical objects.
  for (auto firstIt = m_allCandidates.begin(); firstIt != m_allCandidates.end(); ++firstIt) {
    auto idx1 = std::distance(m_allCandidates.begin(), firstIt);
    CheckObject& obj1 = firstIt->second;
    for (auto secIt = firstIt; secIt != m_allCandidates.end(); ++secIt) {
      auto idx2 = std::distance(m_allCandidates.begin(), secIt);
      if (idx1 == idx2) { continue; }
      CheckObject& obj2 = secIt->second;
      double invMass = (obj1.p4 + obj2.p4).M();

      if (m_massWide.contains(invMass)) {
        obj1.partnerInds.push_back(idx2);
        obj2.partnerInds.push_back(idx1);

        bool narrow = m_massNarrow.contains(invMass);
        obj1.narrow.push_back(narrow);
        obj2.narrow.push_back(narrow);
      }
    }
  }
}

template<typename PATType>
std::vector<PiZeroCheckReport> PiZeroChecker::getPiZeroReport() const
{
  std::cout << m_allCandidates << std::endl;

  std::vector<PiZeroCheckReport> reports;
  std::cout << "obtaining piZeroReports for " << typeid(PATType).name() << std::endl;

  std::cout << "this event has " << m_allCandidates.size() << " records" << std::endl;

  const auto typeRange = m_allCandidates.equal_range(typeid(PATType).hash_code());
  std::cout << "got " << std::distance(typeRange.first, typeRange.second) << " entries for this type." << std::endl;

  std::multimap<size_t, CheckObject> typeCands(typeRange.first, typeRange.second);
  std::cout << typeCands << std::endl;

  const auto pfRange = m_allCandidates.equal_range(typeid(pat::PFParticle).hash_code());
  std::cout << "checking against " << std::distance(pfRange.first, pfRange.second) << " PFlow" << std::endl;

  for (auto typeIt = typeRange.first; typeIt != typeRange.second; ++typeIt) {
    if (typeIt->second.veto()) std::cout << typeIt->second.partnerInds << " " << typeIt->second.narrow << std::endl;
    if (typeIt->second.veto_narrow()) {
      PiZeroCheckReport report;
      report.index = std::distance(m_allCandidates.begin(), typeIt); // TODO: return the right value here!
      reports.push_back(report);
    }
  }

  // // if not already done cache all candidate indices
  // if (m_candidateInds.empty()) {
  //   for (auto it = m_allCandidates.cbegin(); it != m_allCandidates.cend(); ++it) {
  //     if (it->veto()) m_candidateInds.push_back(std::distance(m_allCandidates.cbegin(), it));
  //   }
  // }

  // std::cout << m_candidateInds << std::endl;

  // const auto typeRangeIt = m_typeMap.find(typeid(PATType).hash_code());
  // if (typeRangeIt != m_typeMap.end()) {
  //   TypeRange typeRange = typeRangeIt->second;
  //   std::cout << "found type in typeMap for this event." << std::endl;
  //   std::cout << "range: " << std::distance(m_allCandidates.cbegin(), typeRange.first) << " - " <<
  //     std::distance(m_allCandidates.cbegin(), typeRange.last) << ", " << m_allCandidates.size() << std::endl;
  //   for (auto it = typeRange.first; it != typeRange.last; ++it) { // iterate over all check objects in range
  //     if (it->veto()) std::cout << std::distance(m_allCandidates.cbegin(), it) << " ";
  //     // std::cout << it->veto() << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  // if(typeRangeIt != m_typeMap.end()) {
  //   std::cout << "found type in type_map" << std::endl;
  //   for (const auto& range : typeRangeIt->second) {
  //     size_t firstIdx = std::distance(m_allCandidates.cbegin(), range.first);
  //     size_t lastIdx = std::distance(m_allCandidates.cbegin(), range.last);
  //     std::cout << firstIdx << " " << lastIdx << std::endl;
  //     std::cout << "range contains " << lastIdx - firstIdx << " candidates" << std::endl;
  //     for (const size_t& idx : m_candidateInds) {
  //       if (range.contains(idx, m_allCandidates)) {
  //         std::cout << "idx " << idx << " in m_candidateInds" << std::endl;
  //         std::cout << "partner inidces: " << m_allCandidates[idx].partnerInds << std::endl;

  //         std::vector<bool> sameType;
  //         for (const size_t& idx2 : m_allCandidates[idx].partnerInds) {
  //           sameType.push_back(range.contains(idx2, m_allCandidates));
  //         }
  //         std::cout << "same types: " << sameType << std::endl;
  //       }
  //       // if (idx >= firstIdx && idx < lastIdx) {
  //       //   std::cout << "found " << idx << " in pi0 Candidates list, compatible with range" << std::endl;
  //       //   const auto& checkObj = m_allCandidates[idx];``
  //       //   PiZeroCheckReport report;
  //       //   report.index = idx;

  //       //   // report.narrow = m_allCandidates[idx].narrow;
  //         // reports.push_back(report);
  //       // }
  //     }
  //   }
  // }

  return reports;
}

#endif
