///=== This is the base class for the Full Kalman Combinatorial Filter track finding and fitting algorithm.
 
#ifndef __FULL_KALMAN_COMB__
#define __FULL_KALMAN_COMB__
 
#include <TMatrixD.h>
#include "L1Trigger/TrackFindingTMTT/interface/TrackFitGeneric.h"
#include "L1Trigger/TrackFindingTMTT/interface/Stub.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1track3D.h"
#include "L1Trigger/TrackFindingTMTT/interface/L1fittedTrack.h"
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
class StubCluster;

class FullKalmanComb : public TrackFitGeneric{
 
    public:
	enum OVERLAP_TYPE { TYPE_NORMAL, TYPE_V2, TYPE_NOCLUSTERING, TYPE_TP };
    public:
	FullKalmanComb(const Settings* settings, const uint nPar, const string &kfName="", const uint nMeas=2 );

	virtual ~FullKalmanComb() {this->resetStates();}

        void hitBuffer (const Stub* stub, const vector<bool>& inSubSecs);
        void createSeeds();
        void run ();

	L1fittedTrack fit(const L1track3D& l1track3D);
	void bookHists();

        //=== Get results
        //=== Most of these functions operate on stubs that have been filtered (e.g. have they passed the bend filter)
        //=== If no filters were requested, they are identical to the unfiltered stubs.)

        // Get the filtered stubs to be used in the KF
        const vector<const Stub*>& stubs() const { return vFilteredStubs_; }

       // Return info useful for deciding if the track candidate meets various criteria
       unsigned int numStubs()         const { return vFilteredStubs_.size(); }      // Number of filtered stubs 
       unsigned int numLayers()        const { return numFilteredLayersInCell_; }    // Number of tracker layers with filtered stubs
       unsigned int numLayersSubSec()  const { return numFilteredLayersInCellBestSubSec_; }  // Number of tracker layers with filtered stubs,  requiring all stubs to be in same subsector to be counted. The number returned is the highest layer count found in any of the subsectors in this sector. If subsectors are not used, it is equal to numLayers().

       // Useful for debugging.
       unsigned int numUnfilteredStubs()   const { return vStubs_.size(); }    // Number of unfiltered stubs 

    protected:
	static  std::map<std::string, double> getTrackParams( const FullKalmanComb *p, const kalmanState *state );
	virtual std::map<std::string, double> getTrackParams( const kalmanState *state ) const=0;

        // Get track params with beam-spot constraint & chi2 increase from applying it..
        virtual std::map<std::string, double> getTrackParams_BeamConstr( const kalmanState *state, double& deltaChi2 ) const {
          deltaChi2 = 0.0;
          return (this->getTrackParams(state)); // Returns unconstrained result, unless derived class overrides it.
        }

	double sectorPhi()const
	{
	    return 2.*M_PI * (0.5 + float(iCurrentPhiSec_)) / float(getSettings()->numPhiSectors()) - M_PI; // Centre of sector in phi
	}
        //bool kalmanUpdate( const StubCluster *stubCluster, kalmanState &state, kalmanState &new_state, const TP *tpa );
	virtual const kalmanState *kalmanUpdate( unsigned skipped, unsigned layer, const StubCluster* stubCluster, const kalmanState &state, const TP *);
	void resetStates();
	void deleteStubClusters();
	const kalmanState *mkState( const L1track3D &candidate, unsigned skipped, unsigned layer, unsigned layerId, const kalmanState *last_state, 
				    const std::vector<double> &x, const TMatrixD &pxx, const TMatrixD &K, const TMatrixD &dcov, const StubCluster* stubCluster, double chi2 );

    protected:
	/* Methods */
	std::vector<double> Hx( const TMatrixD &pH, const std::vector<double> &x )const;
	std::vector<double> Fx( const TMatrixD &pF, const std::vector<double> &x )const;
	TMatrixD HxxH( const TMatrixD &pH, const TMatrixD &xx )const;
	double Chi2( const TMatrixD &dcov, const std::vector<double> &delta, bool debug = false )const;
	TMatrixD GetKalmanMatrix( const TMatrixD &h, const TMatrixD &pxcov, const TMatrixD &dcov )const;
	void GetAdjustedState( const TMatrixD &K, const TMatrixD &pxcov, 
			       const std::vector<double> &x, const StubCluster *stubCluster, const std::vector<double>& delta,  
			       std::vector<double> &new_x, TMatrixD &new_xcov )const;


	virtual std::vector<double> seedx(const L1track3D& l1track3D)const=0;
	virtual TMatrixD seedP(const L1track3D& l1track3D)const=0;
	virtual void barrelToEndcap( double r, const StubCluster *stubCluster, std::vector<double> &x, TMatrixD &cov_x )const{}
	virtual std::vector<double> d(const StubCluster* stubCluster )const=0;
	virtual TMatrixD H(const StubCluster* stubCluster)const=0;
	virtual TMatrixD F(const StubCluster* stubCluster=0, const kalmanState *state=0 )const=0;
	virtual TMatrixD PxxModel( const kalmanState *state, const StubCluster *stubCluster )const=0; 
  	virtual TMatrixD PddMeas(const StubCluster* stubCluster, const kalmanState *state )const=0;

        virtual std::vector<double> residual(const StubCluster* stubCluster, const std::vector<double> &x, double candQoverPt )const;
	virtual const kalmanState *updateSeedWithStub( const kalmanState &state, const StubCluster *stubCluster ){ return 0; }
	virtual bool isGoodState( const kalmanState &state )const{ return true; }

	double calcChi2( const kalmanState &state )const;
	void printTP( std::ostream &os, const TP *tp )const;


	virtual double getRofState( unsigned layerId, const vector<double> &xa )const{ return 0;}
	std::vector<const kalmanState *> doKF( const L1track3D &l1track3D, const std::vector<const StubCluster *> &stubClusters, const TP *tpa );

	void fillSeedHists( const kalmanState *state, const TP *tpa );
	void fillCandHists( const kalmanState &state, const TP *tpa=0 );
	void fillStepHists( const TP *tpa, unsigned nItr, const kalmanState *new_state );

	double DeltaRphiForClustering( unsigned layerId, unsigned endcapRing );
	double DeltaRForClustering( unsigned endcapRing );
	bool isOverlap( const Stub* a, const Stub*b, OVERLAP_TYPE type );

	set<unsigned> getKalmanDeadLayers( unsigned layerMap[18][8], bool& remove2PSCut ) const;

        // Function to calculate approximation for tilted barrel modules (aka B) copied from Stub class.
        float getApproxB(float z, float r) const;

        // Is this HLS code?
        virtual bool isHLS() {return false;};

    protected:
	unsigned nPar_;
	unsigned nMeas_;
	std::vector<kalmanState *> state_list_;
	std::vector<StubCluster *> stbcl_list_;

        std::vector<const Stub*> vStubs_;

	std::vector<double> hxaxtmin;
	std::vector<double> hxaxtmax;
	std::vector<double> hxmin;
	std::vector<double> hxmax;
	std::vector<double> hymin;
	std::vector<double> hymax;
	std::vector<double> hdxmin;
	std::vector<double> hdxmax;
	std::vector<double> hresmin;
	std::vector<double> hresmax;
	std::vector<double> hddMeasmin;
	std::vector<double> hddMeasmax;

	TH1F * hTrackEta_;
	TH1F * hUniqueTrackEta_;
	std::map<TString, TH2F*> hBarrelStubMaxDistanceMap;
	std::map<TString, TH2F*> hEndcapStubMaxDistanceMap;
	std::map<TString, TH2F*> hphiErrorRatioMap;
	std::map<TString, TH1F*> hstubCombMap;

	TH1F*           hndupStub_;
	TH1F*           hnmergeStub_;
	std::map<TString, TH1F*> hytMap;
	std::map<TString, TH1F*> hy0Map;
	std::map<TString, TH1F*> hyfMap;
	std::map<TString, TH1F*> hxMap;
	std::map<TString, TH1F*> hxcovMap;
	std::map<TString, TH1F*> hkMap;
	std::map<TString, TH1F*> hresMap;
	std::map<TString, TH1F*> hmcovMap;

      	double hchi2min;
	double hchi2max;

	unsigned maxNfitForDump_;
	bool     dump_;
	unsigned int      iCurrentPhiSec_;
	unsigned int      iCurrentEtaReg_;
	unsigned int      iLastPhiSec_;
	unsigned int      iLastEtaReg_;

        unsigned int      minStubLayersRed_;

        unsigned int      numUpdateCalls_;

       const TP* tpa_;


       // NEW STUFF

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
       vector<const Stub*> vStubs_; // input stubs	
       vector<const Stub*> vFilteredStubs_; // input stubs that have passed the stub filters
       map<const Stub*, vector<bool>> subSectors_; // Indicate which subsectors within the sector this stub is consistent with.

       // List of all track candidates found and fitted by the KF
       vector<L1track3D> trackCands3D_;
};

}

#endif




