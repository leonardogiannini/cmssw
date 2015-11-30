#ifndef MuonIsolation_TrackSelector_H
#define MuonIsolation_TrackSelector_H

#include "RecoMuon/MuonIsolation/interface/Range.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/View.h"
#include <list>

namespace muonisolation {

  class TrackSelector {
  public:

    typedef muonisolation::Range<float> Range;
    typedef std::list<const reco::Track*> result_type;
    typedef reco::TrackCollection input_type;
    typedef reco::TrackBase::Point BeamPoint;
  
    //!config parameters
    struct Parameters {
      Parameters()
	:  zRange(-1e6,1e6), rRange(-1e6,1e6), dir(0,0), drMax(1e3), beamPoint(0,0,0),
	   nHitsMin(0), chi2NdofMax(1e64), chi2ProbMin(-1.), ptMin(-1) {}
      Parameters(const Range& dz, const double d0Max, const reco::isodeposit::Direction& dirC, double rMax, 
		 const BeamPoint& point  = BeamPoint(0,0,0))
	: zRange(dz), rRange(Range(0,d0Max)), dir(dirC), drMax(rMax), beamPoint(point),
	  nHitsMin(0), chi2NdofMax(1e64), chi2ProbMin(-1.), ptMin(-1) {}
      Range        zRange;  //! range in z
      Range        rRange;  //! range in d0 or dxy (abs value)
      reco::isodeposit::Direction       dir;  //! direction of the selection cone
      double        drMax;  //! cone size
      BeamPoint beamPoint;  //! beam spot position
      unsigned int       nHitsMin;  //! nValidHits >= nHitsMin
      double  chi2NdofMax;  //! max value of normalized chi2
      double  chi2ProbMin;  //! ChiSquaredProbability( chi2, ndf ) > chi2ProbMin
      double        ptMin;  //! tk.pt>ptMin
    };


    TrackSelector(const Parameters& pars) : thePars(pars) { } 

      result_type operator()(const input_type & tracks) const;


  private:
      Parameters thePars;
  }; 

}

#endif 
