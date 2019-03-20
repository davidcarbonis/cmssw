#ifndef __kalmanSeed_H__
#define __kalmanSeed_H__

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/Utility.h"
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/Sector.h"

#include <vector>
#include <utility>

using namespace std;

//=== Full track finding and fitter KF seed.
//=== Gives access to all stubs used for the seeding and remaining stubs and initial helix estimates.
//=== Also calculates & gives access to associated truth particle (Tracking Particle) if any.

namespace TMTT {

  class KalmanSeed {

  public: 

    KalmanSeed (const Settings* settings, const <Stub*>& seedStubs, vector<Stub*>& otherStubs, pair<float, float> helixRphi, pair<float, float> helixRz, unsigned int iPhiSec, unsigned int iEtaReg, unsigned int optoLinkID) :
        settings_(settings),
        seedStubs_(seedStubs),
        otherStubs_(otherStubs),
        helixRphi_(helixRphi),
        helixRz_(helixRz),
        iPhiSec_(iPhiSec),
        iEtaSec_(iEtaSec),
        optoLinkID_(optoLinkID)
    {
        nSeedLayers_ = Utility::countLayers(settings, seedStubs); // count the tracker layers that the seed stub(s) are in
        matchedTP_   = Utility::matchingTP(settings, seedStubs, matchedSeedLayers_, matchedSeedStubs_);
    }

    ~KalmanSeed() {}

    // Get the initial stubs to be used as seeds by the KF and the remaining ones to be used in the track building
    const vector< const Stub*>&        getSeedStubs()              const  {return seedStubs_;}
    const vector< const Stub*>&        getOtherStubs()             const  {return otherStubs_;}

    // Get the total number of seed stubs
    unsigned int                       getNumSeedStubs()           const  {return seedStubs_.size();}
    // Get the total number of stubs to be considered (i.e. minus seeds) during the track building
    unsigned int                       getNumOtherStubs()          const  {return otherStubs_.size();}

    // Return the conventionally agreed track helix parameters relevant in the r-phi and r-z planes
    pair<float, float>                 getHelixRphi()              const  {return helixRphi_;}
    pair<float, float>                 getHelixRz()                const  {return helixRZ_;}

    // User-friendly access to the helix parameters
    float   charge()     const  {return (this->qOverPt() > 0  ?  1  :  -1);} 
    float   invPt()      const  {return fabs(this->qOverPt());}
    float   pt()         const  {return 1./(1.0e-6 + this->invPt());} // includes protection against 1/pt = 0.
    float   d0()         const  {return 0.;} // Hough transform assumes d0 = 0.
    float   phi0()       const  {return helixRphi_.second;}
    float   z0()         const  {return helixRz_.first;}
    float   tanLambda()  const  {return helixRz_.second;}
    float   theta()      const  {return atan2(1., this->tanLambda());} // Use atan2 to ensure 0 < theta < pi.
    float   eta()        const  {return -log(tan(0.5*this->theta()));}
    float   qOverPt()    const  {return helixRphi_.first;}

    //--- Get phi sector and eta region used by track finding code that this track is in.
    unsigned int iPhiSec() const  {return iPhiSec_;}
    unsigned int iEtaReg() const  {return iEtaReg_;}

    //--- Opto-link ID used to send this track from HT to Track Fitter
    unsigned int optoLinkID() const {return optoLinkID_;}

    //--- Get info about the seed's association (if any) to a truth Tracking Particle
    // Get the best matching TP (=nullptr if none).
    const TP*                  getMatchedTP()          const   {return matchedTP_;}
    // Get the matched seed stubs with this TP
    const vector<const Stub*>& getMatchedSeedStubs()       const   {return matchedSeedStubs_;}
    // Get the number of matched seed stubs with this TP
    unsigned int               getNumMatchedSeedStubs()    const   {return matchedSeedStubs_.size();}
    // Get number of tracker layers with matched seed stubs with this TP
    unsigned int               getNumMatchedSeedLayers()   const   {return nMatchedSeedLayers_;}
    // Get purity of seed stubs on track candidate (i.e. fraction matching best TP)
    float                      getPurity()             const   {return getNumMatchedSeedStubs()/float(getNumSeedStubs());}


    private:

    //--- Configuration parameters
    const Settings*                    settings_;
  
    //--- Seed information, including stubs used as seeds
    vector<const Stub*>                seedStubs_;
    unsigned int                       nSeedLayers_;
    vector<const Stub*>                otherStubs_;
    pair<float, float> >               helixRphi_;
    pair<float, float> >               helixRz_;
    unsigned int                       iPhiSec_;
    unsigned int                       iEtaReg_; 
    unsigned int                       optoLinkID_;


    //--- Information about its association (if any) to a truth Tracking Particle.  
    const TP*                          matchedTP_;
    vector<const Stub*>                matchedSeedStubs_;
    unsigned int                       nMatchedSeedLayers_;

  }
}



