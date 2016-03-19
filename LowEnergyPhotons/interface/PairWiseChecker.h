#ifndef __OniaGamma_PairWiseChecker_h_
#define __OniaGamma_PairWiseChecker_h_

#include <map>
#include <utility>
#include <typeinfo>
#include <vector>
#include <iterator>

template<typename CheckObject, typename CheckFunction, typename ReturnType>
class PairWiseChecker {
public:
  PairWiseChecker(const CheckFunction& f) : m_checkFunc(f) {;} /**< Constructor taking a function object. */
  PairWiseChecker() = delete; /**< CheckFunction has to be provided at construction. */
  ~PairWiseChecker() = default; /**< default constructor. Everything that is allocated is managed, so this should work. */

  using mapT = std::multimap<size_t, CheckObject>; /**< typedef of internal used map-type*/

  /** add a collection that should be checked. */
  template<typename Collection>
  void addCollections(const Collection& collection);

  /** add more than one collection of different type to be checked. */
  template<typename Collection, typename... MoreCollections>
  void addCollections(const Collection& collection, MoreCollections... collections)
  {
    addCollections(collection);
    addCollections(collections...);
  }

  /** perform the operation defined by CheckFunction on each pair-wise combination.*/
  void check();

  /** get the indices of the objects in the collection that failed the check (as specified by the CheckFunction). */
  template<typename Collection>
  std::vector<ReturnType> getCheckIndices() const;

private:

  mapT m_allObjects; /**< map that stores all information that is needed to perform the checks. */

  CheckFunction m_checkFunc; /**< the check function to be used. */
};

template<typename CheckObject, typename CheckFunction, typename ReturnType>
template<typename Collection>
void PairWiseChecker<CheckObject, CheckFunction, ReturnType>::addCollections(const Collection& collection)
{
  size_t typeHash = typeid(Collection).hash_code();
  for (const auto& obj : collection) { m_allObjects.insert(std::make_pair(typeHash, CheckObject(obj))); }
}

template<typename CheckObject, typename CheckFunction, typename ReturnType>
void PairWiseChecker<CheckObject, CheckFunction, ReturnType>::check()
{
  // NOTE: cannot use iterator offsets on std::multimap since they only define forward_iterators!
  // Solution: loop over all candidates and make sure to not combine identical objects
  for (auto firstIt = m_allObjects.begin(); firstIt != m_allObjects.end(); ++firstIt) {
    for (auto secIt = m_allObjects.begin(); secIt != m_allObjects.end(); ++secIt) {
      if (firstIt == secIt) continue;

      m_checkFunc(firstIt->second, secIt->second); // COULDDO: ipmlement usage of return type for something.
    }
  }
}

template<typename CheckObject, typename CheckFunction, typename ReturnType>
template<typename Collection>
std::vector<ReturnType> PairWiseChecker<CheckObject, CheckFunction, ReturnType>::getCheckIndices() const
{
  std::vector<ReturnType> returnVec{};
  const auto typeRange = m_allObjects.equal_range(typeid(Collection).hash_code());
  if (!std::distance(typeRange.first, typeRange.second)) return returnVec; // return empty if type not present

  mapT typeObjects(typeRange.first, typeRange.second);
  for (auto typeIt = typeRange.first; typeIt != typeRange.second; ++typeIt) {
    if (typeIt->second.veto()) {
      size_t ind = std::distance(typeRange.first, typeIt);
      returnVec.push_back(ReturnType(ind, typeIt->second));
    }
  }

  return returnVec;
}

#endif
