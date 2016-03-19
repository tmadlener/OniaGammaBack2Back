#ifndef __OniaGamma_PiZeroChecker_h_
#define __OniaGamma_PiZeroChecker_h_

#include "OniaGammaBack2Back/LowEnergyPhotons/interface/PairWiseChecker.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include <vector>

#include "/afs/hephy.at/user/t/tmadlener/snippets/type_deduction_helper.h"
#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"
#include "/afs/hephy.at/user/t/tmadlener/snippets/map_stuff.h"

/** type that can be used to store the information needed for checking two mass-windows at once.
 * NOTE: at the moment a bit taylored to the needs of PiZero Checking, but should be general enough for it to still
 * work in other circumstances.
 * LorentzVectorT must provide: operator+ and double M()
 */
template<typename LorentzVectorT>
struct DoubleMassWindowCO {
  /** constructor from any type that provides a ::p4() method. */
  template<typename AnyType>
  DoubleMassWindowCO(const AnyType& object) : m_p4(object.p4()) {;}

  /** constructor from a reco::Conversion, which does not provide a ::p4() method. */
  DoubleMassWindowCO(const reco::Conversion& conversion) : m_p4(conversion.refittedPair4Momentum()) {;}

  DoubleMassWindowCO() = delete; /**< should not be default constructible. */
  ~DoubleMassWindowCO() = default; /**< should do. TODO: revise this again after object is finished. */

  /** MANDATORY veto function. true if any of the two mass windows vetoes. */
  bool veto() const { return ( this->veto1() || this->veto2() ); }

  bool veto1() const { return any(m_massWindow1); } /**< veto from mass window 1.*/
  bool veto2() const { return any(m_massWindow2); } /**< veto from mass window 2. */

  LorentzVectorT m_p4; /**< the four momentum. For the mass only this is needed! */

  std::vector<bool> m_massWindow1; /**< overlaps with mass window 1. */

  std::vector<bool> m_massWindow2; /**< overlaps with mass window 2. */
};

/** typedef for PiZero Checking.*/
using PiZeroCheckObject = DoubleMassWindowCO<reco::Candidate::LorentzVector>;

/** functor that does the double mass window checking, provided it gets two appropriate check objects. */
template<typename CheckObject>
struct DoubleMassWindowChecker {
  /** constructor using two mass windows. */
  DoubleMassWindowChecker(const std::vector<double>& w1, const std::vector<double> w2) :
    m_window1(w1), m_window2(w2) {;}
  DoubleMassWindowChecker() = delete; /**< no use for default constructor. */

  void operator()(CheckObject& obj1, CheckObject& obj2);

private:
  std::vector<double> m_window1;
  std::vector<double> m_window2;
};

template<typename CheckObject>
void DoubleMassWindowChecker<CheckObject>::operator()(CheckObject& obj1, CheckObject& obj2)
{
  double invMass = (obj1.m_p4 + obj2.m_p4).M();
  if (invMass > m_window1[0]  && invMass < m_window1[1]) {
    obj1.m_massWindow1.push_back(true);
    obj2.m_massWindow1.push_back(true);
  }
  if (invMass > m_window2[0] && invMass < m_window2[1]) {
    obj1.m_massWindow2.push_back(true);
    obj2.m_massWindow2.push_back(true);
  }
}

/** typedef for PiZero Checking. */
using PiZeroCheckFunction = DoubleMassWindowChecker<PiZeroCheckObject>;


/** struct holding the information needed at return from the DoubleMassWindow checking.*/
struct DoubleMassWindowRT {
  /** */
  template<typename LorentzVectorT>
  DoubleMassWindowRT(size_t ind, const DoubleMassWindowCO<LorentzVectorT>& checkObj) :
    m_index(ind),
    m_massWindow1(checkObj.veto1()),
    m_massWindow2(checkObj.veto2()) {;}

  size_t m_index{};
  bool m_massWindow1{false};
  bool m_massWindow2{false};
};


/** Class that takes different collections of photons and checks if any combination fo two photons could possibly
 * originiate from a \f$\pi^{0}\f$-decay. Does so by summing the four momenta of the two photons and checking if
 * the mass is compatible with two spezialized mass windows (which have to be provided at the construction of the
 * CheckFunction).
 */
using PiZeroChecker = PairWiseChecker<PiZeroCheckObject, PiZeroCheckFunction, DoubleMassWindowRT>;

#endif
