///=== This is the base class for seeding the  Kalman Combinatorial Filter track finding and fitting algorithm.

///=== Written by: A. D. Morton, adapted from L1KalmanComb class

#include "L1Trigger/TrackFindingTMTT/interface/KalmanCombSeeder.h"
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

static double wrapRadian( double t ){
	
  if( t > 0 ){
    while( t > M_PI ) t-= 2*M_PI; 
  }
  else{
    while( t < - M_PI ) t+= 2*M_PI; 
  }
  return t;
}

void KalmanCombSeeder::init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
		  float etaMinSector, float etaMaxSector, float phiCentreSector) {

  settings_ = settings;

  // Sector number
  iPhiSec_ = iPhiSec;
  iEtaReg_ = iEtaReg;

    // Rapidity range of sector.
  etaMinSector_ = etaMinSector;
  etaMaxSector_ = etaMaxSector;

  invPtToInvR_   = settings->invPtToInvR();
  invPtToDphi_   = settings->invPtToDphi();  // B*c/2E11
 
  chosenRofPhi_       = settings->chosenRofPhi();
  phiCentreSector_    = phiCentreSector; // Centre of phiTrk sector.

  // Used to kill excess stubs or tracks that can't be transmitted within time-multiplexed period.
  nReceivedStubs_ = 0;
  busyInputSectorKill_     = settings_->busyInputSectorKill();   // Kill excess stubs going fron GP to HT?
  busyInputSectorNumStubs_ = settings_->busyInputSectorNumStubs(); // Max. num. of stubs that can be sent from GP to HT within TM period
  busySectorKill_          = settings_->busySectorKill();        // Kill excess tracks flowing out of HT?
  busySectorNumStubs_      = settings_->busySectorNumStubs();    // Max. num. of stubs that can be sent out of HT within TM period

  // CKF Seeder Options
  seedingOption_           = settings_->kalmanSeedingOption();
}

void KalmanCombSeeder::stubBuffer( const Stub* stub ) {

  //=== Filter stubs
  // Stub Bend Filter?

  // Optionally, only store stubs that can be sent from GP to KF within TM period.
  if ( ( ! busyInputSectorKill_) || (nReceivedStubs_ <  busyInputSectorNumStubs_) ) {
    nReceivedStubs_++;

    unsigned int stubLayer = stub->layerIdReduced();

    //Seeding options:
    // Search for pairs of stubs, thinking about missing layers
    if ( seedingOption_ == 1 ) {
      if ( stubLayer == 1 ) vSeedStubs_.push_back(stub);
      else if ( stubLayer == 2 ) vSeedStubs_.push_back(stub);
      else vOtherStubs_.push_back(stub);

      if ( stubLayer == 1 ) vLayer1Stubs_.push_back(stub);
      if ( stubLayer == 2 ) vLayer2Stubs_.push_back(stub);
    }
    //Default option L0 + beamspot
    else {
      if ( stubLayer == 1 ) vSeedStubs_.push_back(stub);
      else vOtherStubs_.push_back(stub);
    }
  }
}

void KalmanCombSeeder::createSeeds() {

  vector<L1track3D> trackCands3D;

  //=== Create Seeds  
  // SPIT OUT A KalmanSeed or L1track3D to be used by the KF fit

//    KalmanSeed (const Settings* settings, const vector<const Stub*>& seedStubs, pair<float, float> helixRphi,
//                pair<float, float> helixRz, unsigned int iPhiSec, unsigned int iEtaReg, unsigned int optoLinkID) :

  // Seeding from layers 1+2 (or just 2)
  if ( seedingOption_ == 1 ) {

    const vector<const Stub*>& layer1Stubs = vLayer1Stubs_;
    const vector<const Stub*>& layer2Stubs = vLayer2Stubs_;
    vector<const Stub*>& backupStubs = vLayer2Stubs_;
    const vector<const Stub*>& otherStubs = vOtherStubs_;

    for ( auto innerStub : layer1Stubs ) {
      if ( layer2Stubs.size() == 0 ) continue;
      for ( auto outerStub : layer2Stubs ) {

        // compataibility test 1
        if ( (innerStub->r() - (outerStub->z()-innerStub->z())*innerStub->r()/(outerStub->r()-innerStub->r())) > 30. ) continue; 

        double deltaPhi = wrapRadian( innerStub->phi() - outerStub->phi() );
        double displacement = sqrt( pow( outerStub->r(), 2 ) + pow( innerStub->r(), 2 )
                              - 2*outerStub->r()*innerStub->r()*cos(deltaPhi));
  

        double rInv = 2*sin(deltaPhi)/displacement;
        // compataibility test 2
        if ( rInv > 0.0057 ) continue;

        float phi0 = innerStub->beta();
        double tan_lambda = ( innerStub->z() - outerStub->z() ) * rInv / ( 2 * ( asin( 0.5*rInv*innerStub->r() ) - asin( 0.5*rInv*outerStub->r() ) ) );
        double z0 = innerStub->z() - tan_lambda * ( 2 * ( asin( 0.5 * rInv * innerStub->r() ) ) ) / rInv;
        // compataibility test 3
        if ( z0 > 15. ) continue;

        float qOverPt = rInv*invPtToInvR_;

        const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
        const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda

        unsigned int optoLinkID = this->calcOptoLinkID();    

        vector<const Stub*> stubs {innerStub,outerStub};
        stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );

        const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location

        L1track3D l1Trk3D(settings_, stubs, cellLocation, helixParamsRphi, helixParamsRz, iPhiSec_, iEtaReg_, optoLinkID, false);

        // push back candidate
        trackCands3D_.push_back(l1Trk3D);      

        // erase Layer 2 stub from 2nd pass seed collcetion 
//        vector<const Stub*>::iterator position = std::find(layer2Stubs.begin(), layer2Stubs.end(), outerStub)
//        if ( position != layer2Stubs.end() ) backupStubs.erase(position);
      }
    }

  // Create seeds from layer 2 stubs not used to make pairs with layer 1

  }

  // If seeding with just a single stub from Layer 1
  else {
    // Read in stubs from initial buffer and create seeds from the seed stub collection	
    const vector<const Stub*>& seedStubs = vSeedStubs_;
    const vector<const Stub*>& otherStubs = vOtherStubs_;

    createSingleStubSeeds( seedStubs, otherStubs );

  }

}

void KalmanCombSeeder::createSingleStubSeeds( const vector<const Stub*>& seedStubs, const vector<const Stub*>& otherStubs ) {

  for ( auto stub : seedStubs ) {
    
    vector<const Stub*> stubs {stub};
    //    stubs.push_back(stub);
    stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );
    
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
}

}
