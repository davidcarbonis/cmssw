///=== This is the base class for the Full Kalman Combinatorial Filter track finding and fitting algorithm.

///=== Written by: A. D. Morton, adapted from L1KalmanComb class

#include "L1Trigger/TrackFindingTMTT/interface/FullKalmanComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/Utility.h"

#include <TMatrixD.h> 
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
//#include "L1Trigger/TrackFindingTMTT/interface/kalmanSeed.h"
#include "L1Trigger/TrackFindingTMTT/interface/kalmanState.h"
#include "L1Trigger/TrackFindingTMTT/interface/StubCluster.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <algorithm>
#include <functional>
#include <fstream>
#include <iomanip>
#include <TH2F.h>

//#define CKF_DEBUG
// Enable debug printout to pair with that in Histos.cc enabled by recalc_debug.
//#define RECALC_DEBUG

// Enable merging of nearby stubs.
//#define MERGE_STUBS

namespace TMTT {

void FullKalmanComb::init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
		  float etaMinSector, float etaMaxSector, float phiCentreSector) {

  settings_ = settings;

  // Sector number
  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;

    // Rapidity range of sector.
  etaMinSector_ = etaMinSector;
  etaMaxSector_ = etaMaxSector;

  invPtToDphi_   = settings->invPtToDphi();  // B*c/2E11
 
  chosenRofPhi_       = settings->chosenRofPhi();
  phiCentreSector_    = phiCentreSector; // Centre of phiTrk sector.

  // Used to kill excess stubs or tracks that can't be transmitted within time-multiplexed period.
  nReceivedStubs_ = 0;
  busyInputSectorKill_     = settings_->busyInputSectorKill();   // Kill excess stubs going fron GP to HT?
  busyInputSectorNumStubs_ = settings_->busyInputSectorNumStubs(); // Max. num. of stubs that can be sent from GP to HT within TM period
  busySectorKill_          = settings_->busySectorKill();        // Kill excess tracks flowing out of HT?
  busySectorNumStubs_      = settings_->busySectorNumStubs();    // Max. num. of stubs that can be sent out of HT within TM period

  // Full CKF Options
  seedingOption_           = settings_->kalmanSeedingOption();
}

void FullKalmanComb::stubBuffer( const Stub* stub ) {

  //=== Filter stubs
  // Stub Bend Filter?

  // Optionally, only store stubs that can be sent from GP to KF within TM period.
  if ( ( ! busyInputSectorKill_) || (nReceivedStubs_ <  busyInputSectorNumStubs_) ) {
    nReceivedStubs_++;

    unsigned int stubLayer = stub->layerIdReduced();

    //Seeding options:
    //Default option L0 + beamspot
    if ( seedingOption_ == 0 ) {
      if ( stubLayer == 1 ) vSeedStubs_.push_back(stub);
      else vOtherStubs_.push_back(stub);
    }
  }
}

void FullKalmanComb::createSeeds() {

  vector<L1track3D> trackCands3D;

  // Read in stubs from initial buffer and create seeds from the seed stub collection	

  const vector<const Stub*>& seedStubs = vSeedStubs_;


  //=== Create Seeds  
  // SPIT OUT A KalmanSeed or L1track3D to be used by the KF fit

//    KalmanSeed (const Settings* settings, const vector<const Stub*>& seedStubs, pair<float, float> helixRphi,
//                pair<float, float> helixRz, unsigned int iPhiSec, unsigned int iEtaReg, unsigned int optoLinkID) :

  // If seeding with just a single stub from Layer 0
  if (seedingOption_ ==0) {
    
  }

  for ( auto stub : seedStubs ) {

    vector<const Stub*> stubs {stub};
//    stubs.push_back(stub);
    stubs.insert( stubs.end(), vOtherStubs_.begin(), vOtherStubs_.end() );

    const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location

    float qOverPt = stub->qOverPt();
    float phi0 = stub->beta();
    float z0 = 0;
    float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector_))) + 1/tan(2*atan(exp(-etaMaxSector_))));

    const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
    const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda

    unsigned int optoLinkID = this->calcOptoLinkID();    

    L1track3D l1Trk3D(settings_, stubs, cellLocation, helixParamsRphi, helixParamsRz, iPhiSec_, iEtaReg_, optoLinkID, false);

    trackCands3D_.push_back(l1Trk3D);
  }


/*  L1track3D l1Trk3D(settings_, stubs, cellLocation, helixParamsRphi, helixParamsRz, iPhiSec, iEtaReg, optoLinkID, false);

  float deltaPhi = wrapRadian( innerStub->phi() - outerStub->phi() );
  float displacement = sqrt( pow( outerStub->r(), 2 ) + pow( innerStub->r(), 2 )
                        - 2*outerStub->r()*innerStub->r()*cos(deltaPhi));
  float qOverPt = 2*sin(deltaPhi)/displacement;

  if ( endcap stubs && innerStub->() > outerStub->r() ) qOverPt = -qOverPt;

  float phi0 = innerStub->phi() - this->secPhiMin() + asin(0.5*innerStub->r()*rInv);
  float z0 = innerStub->z() - tanLambda * ( 2 * ( asin( 0.5 * rInv * innerStub->r() ) ) ) / rInv;
  float tan_lambda = ( innerStub->z() - outerStub->z() ) * rInv / ( 2 * ( asin( 0.5*rInv*innerStub->r() ) - asin( 0.5*rInv*outerStub->r() ) ) );


  // Store all this info about the track ...
*/
// Something other seeding option?
}

void FullKalmanComb::run() {}


}
