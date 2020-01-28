#ifndef DataFormats_PixelClusterCountsPerEvent_h
#define DataFormats_PixelClusterCountsPerEvent_h
/** \class reco::PixelClusterCountsPerEvent
 *  
 * Reconstructed PixelClusterCountsPerEvent object that will contain the moduleID, bxID, and counts per event.
 *
 * \authors: Sam Higginbotham (shigginb@cern.ch), Chris Palmer (capalmer@cern.ch), Attila Radl (attila.radl@cern.ch) 
 * 
 *
 */
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Luminosity/interface/LumiConstants.h"

namespace reco {
  class PixelClusterCountsPerEvent {
  public:
    void increment(int mD, int count) {
      size_t modIndex = std::distance(m_ModID.begin(), std::find(m_ModID.begin(), m_ModID.end(), mD));
      if (modIndex == m_ModID.size()) {
        m_ModID.push_back(mD);
        m_counts.push_back(0);
      }
      m_counts.at(modIndex) += count;
    }

    void setbxID(unsigned int inputbxID) { bxID = inputbxID; }

    std::vector<int> const& readCounts() const { return (m_counts); }

    std::vector<int> const& readModID() const { return (m_ModID); }

  private:
    std::vector<int> m_counts;
    std::vector<int> m_ModID;
    unsigned int bxID;
  };

}  // namespace reco
#endif

