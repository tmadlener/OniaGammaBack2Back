#ifndef __OniaGamma_PairWiseChecker_h_
#define __OniaGamma_PairWiseChecker_h_

#include <map>
#include <utility>
#include <typeinfo>
#include <vector>
#include <iterator>

/** Templated class to perform checks on pair of object.
 *
 * Different collections can be added via addCollections() (takes more than one collection, which need not be of
 * the same type). Calling check() performs the check specified by CheckFunction on each possible combination
 * of two objects. The results are stored in CheckObject for each passed object (together with the information
 * of which type the object is for returning purposes). Calling getCheckIndices() returns a vector of ReturnType
 * for all objects of a given type. In ReturnType the index of the object in the original collection gets stored
 * as well as arbitrary user defined information.
 *
 * Requirements on the template parameters:
 * @param CheckObject - internally used Object that needs to store all information needed to do the checks. Needs
 * to be constructible from all objects of all collections that are passed into addCollections(). Needs to provide
 * bool veto() to indicate whether or not a check failed on an object.
 * This is essentially used to get a homogenous interface to all objects despite them being of different types with
 * probably different interfaces.
 * @param CheckFunction - function object taking two CheckObject (references) as arguments to operator(). The
 * return type of operator() is not used and can thus be declared as void.
 * @param ReturnType - User defined ReturnType that contains the index of the original object in the original
 * collection. Needs to provide a constructor that takes a size_t as first and a CheckObject as second argument.
 */
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
