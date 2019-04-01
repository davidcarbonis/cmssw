///=== This is the base class for the Full Kalman Combinatorial Filter track finding and fitting algorithm.
 
#ifndef __FULL_KALMAN_COMB__
#define __FULL_KALMAN_COMB__
 
#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/kalmanState.h"
#include <map>
#include <vector>
#include <fstream>
#include <TString.h>

class TH1F;
class TH2F;

namespace TMTT {

class TP; 
class kalmanState;
//class KalmanSeed;
class StubCluster;
class Stub;

class FullKalmanComb {
 
    public:

	FullKalmanComb() {}
	virtual ~FullKalmanComb() {}


        void init (const Settings* settings, unsigned int iPhiSec, unsigned int iEtaReg, 
                  float etaMinSector, float etaMaxSector, float phiCentreSector);
        void stubBuffer (const Stub* stub);
        void createSeeds();
//        void run();
        virtual const vector<L1track3D>& trackCands3D() const {return trackCands3D_;}


        const Settings* getSettings() const{return settings_;}


    protected:

        // Configuration parameters
        const Settings* settings_;
        // NEW STUFF

        unsigned int iPhiSec_;
        unsigned int iEtaReg_;
        float etaMinSector_;
        float etaMaxSector_;
        float phiCentreSector_;

        unsigned int numSubSecs_;

        float invPtToDphi_;
        float chosenRofPhi_;  

        // Options for killing stubs/tracks that cant be sent within time-multiplexed period.
        bool busyInputSectorKill_; 
        unsigned int busyInputSectorNumStubs_; 
        bool busySectorKill_;     
        unsigned int busySectorNumStubs_; 

        // Full CKF options
        unsigned int seedingOption_;

        // Number of stubs received from GP, irrespective of the number used in the KF (e.g.filtered, truncation)
        unsigned int nReceivedStubs_;

        // input stub data
        vector<const Stub*> vLayer1Stubs_; // input seed stubs from Layer 1	
        vector<const Stub*> vLayer2Stubs_; // input seed stubs from Layer 2	
        vector<const Stub*> vLayer3Stubs_; // input seed stubs from Layer 3	
        vector<const Stub*> vOtherStubs_; // input non-seed stubs	

        // List of all the Kalman Seeds constructed
//        vector<KalmanSeed> kalmanSeeds_;

        // List of all track candidates found and fitted by the KF
        vector<L1track3D> trackCands3D_;

    private:
        // Calculate output opto-link ID from HT, assuming there is no MUX stage.
        virtual unsigned int calcOptoLinkID() const {
        unsigned int numPhiSecPerOct =  settings_->numPhiSectors() / settings_->numPhiOctants();
        return (iEtaReg_ * numPhiSecPerOct + iPhiSec_);
        }
};

}

#endif




