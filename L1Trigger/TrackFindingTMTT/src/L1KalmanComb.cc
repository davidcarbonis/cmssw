
///=== This is the base class for the Kalman Combinatorial Filter track fit algorithm.

///=== Written by: S. Summers, K. Uchida, M. Pesaresi
///=== Modifications by A. D. Morton for track finding + fitting

#include "L1Trigger/TrackFindingTMTT/interface/L1KalmanComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/Utility.h"

#include <TMatrixD.h> 
#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
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

unsigned LayerId[16] = { 1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25 };


static double wrapRadian( double t ){
	
  if( t > 0 ){
    while( t > M_PI ) t-= 2*M_PI; 
  }
  else{
    while( t < - M_PI ) t+= 2*M_PI; 
  }
  return t;
}

static bool orderStubsByLayer(const Stub* a, const Stub* b){
  return (a->layerId() < b->layerId());
}

static bool orderStubsByZ(const Stub* a, const Stub* b){
  return (a->z() < b->z());
}

static bool orderStubsByR(const Stub* a, const Stub* b){
  return (a->r() < b->r());
}

void printTPSummary( std::ostream &os, const TP *tp, bool addReturn=true ){
	
  if( tp ){
		
    os << "TP ";
    //  os << "addr=" << tp << " ";
    os << "index=" << tp->index() << " ";
    os << "qOverPt=" << tp->qOverPt() << " ";
    os << "phi0=" << tp->phi0() << " ";
    os << "z0=" << tp->z0() << " ";
    os << "t=" << tp->tanLambda() << " ";
    os << "d0=" << tp->d0();
    if( addReturn ) os << endl;
    else os << " | ";
  }
}

void L1KalmanComb::printTP( std::ostream &os, const TP *tp )const{
        
  std::map<std::string, double> tpParams;
  bool useForAlgEff(false);
  if( tp ){
    useForAlgEff = tp->useForAlgEff();
    tpParams["qOverPt"] = tp->qOverPt();
    tpParams["phi0"] = tp->phi0();
    tpParams["z0"] = tp->z0();
    tpParams["t"] = tp->tanLambda();
    tpParams["d0"] = tp->d0();
  }
  if( tp ){
    os << "  TP index = " << tp->index() << " useForAlgEff = " << useForAlgEff << " ";
    for( auto pair : tpParams ){
      os << pair.first << ":" << pair.second << ", "; 
    }
    os << "  inv2R = " << tp->qOverPt() * getSettings()->invPtToInvR() * 0.5; 
  }
  else{
    os << "  Fake"; 
  }
  os << endl;
}

static void printStubLayers( std::ostream &os, std::vector<const Stub *> &stubs ){

  if( stubs.size() == 0 ) os << "stub layers = []" << endl;
  else{
    os << "stub layers = [ ";
    for( unsigned i=0; i<stubs.size()-1; i++ ) os << stubs[i]->layerId() << ", ";
    os << stubs.back()->layerId() << " ]" << endl;
  }
}

static void printStubCluster( std::ostream &os, const StubCluster * stubCluster, bool addReturn=true ){
  os << "stub: ";
  //   os << "addr=" << stub << " "; 
  os << "layer=" << stubCluster->layerId() << " ";
  os << "ring=" << stubCluster->endcapRing() << " ";
  os << "r=" << stubCluster->r() << " ";
  os << "phi=" << stubCluster->phi() << " ";
  os << "z=" << stubCluster->z() << " ";
  os << "sigmaX=" << stubCluster->sigmaX() << " ";
  os << "sigmaZ=" << stubCluster->sigmaZ() << " ";
  os << "dphi_dr=" << stubCluster->dphi_dr() << " ";
  os << "#stubs= " << stubCluster->nStubs() << " ";
  os << "TPids="; 
  std::set<const TP*> tps = stubCluster->assocTPs();
  for( auto tp : tps ) os << tp->index() << ","; 
  if( addReturn ) os << endl;
  else os << " | ";
}

static void printStubClusters( std::ostream &os, std::vector<const StubCluster *> &stubClusters ){

  for( auto &stubcl : stubClusters ){
    printStubCluster( os, stubcl );
  }
}
static void printStub( std::ostream &os, const Stub * stub, bool addReturn=true ){
  os << "stub ";
  //   os << "addr=" << stub << " "; 
  os << "index=" << stub->index() << " ";
  os << "layerId=" << stub->layerId() << " ";
  os << "endcapRing=" << stub->endcapRing() << " ";
  os << "r=" << stub->r() << " ";
  os << "phi=" << stub->phi() << " ";
  os << "z=" << stub->z() << " ";
  os << "sigmaX=" << stub->sigmaX() << " ";
  os << "sigmaZ=" << stub->sigmaZ() << " ";
  os << "TPids="; 
  std::set<const TP*> tps = stub->assocTPs();
  for( auto tp : tps ) os << tp->index() << ","; 
  if( addReturn ) os << endl;
  else os << " | ";

}

static void printStubs( std::ostream &os, std::vector<const Stub *> &stubs ){

  for( auto &stub : stubs ){
    printStub( os, stub );
  }
}




L1KalmanComb::L1KalmanComb(const Settings* settings, const uint nPar, const string &fitterName, const uint nMeas ) : TrackFitGeneric(settings, fitterName ){
  nPar_ = nPar;
  nMeas_ = nMeas;
  hymin = vector<double>( nPar_, -1 );
  hymax = vector<double>( nPar_,  1 );
  hymin[0] = -0.05;
  hymax[0] = +0.05;
  hymin[1] = -3.2;
  hymax[1] = +3.2;
  hymin[2] = -20;
  hymax[2] = +20;
  hymin[3] = -6;
  hymax[3] = +6;
  if (nPar_ == 5) {
    hymin[4] = -5;
    hymax[4] = +5;
  }

  hxmin = vector<double>( nPar_, -1 );
  hxmax = vector<double>( nPar_,  1 );

  hddMeasmin = vector<double>( 2, -1e-3 );
  hddMeasmax = vector<double>( 2,  1e-3 );

  hresmin = vector<double>( 2, -1e-2 );
  hresmax = vector<double>( 2,  1e-2 );

  hxaxtmin = vector<double>( nPar_, -1 );
  hxaxtmax = vector<double>( nPar_,  1 );

  hdxmin = vector<double>( nPar_, -1 );
  hdxmax = vector<double>( nPar_,  1 );

  hchi2min = 0; 
  hchi2max = 50; 

  maxNfitForDump_ = 10; 
  dump_ = false; 

  iLastPhiSec_ = 999;
  iLastEtaReg_ = 999;
}
  /*
L1track3D L1KalmanComb::singleStubSeed( const Stub* seed ) {
}

L1track3D L1KalmanComb::singleStubClusterSeed( const StubCluster* seed, const StubCluster* ){

  float qOverPt = seed->stubs()[0]->qOverPt();

  if ( seed->stubs().size() > 1 ) {
    for ( unsigned s = 1; s != seed->stubs().size(); s++ ) {
      if ( fabs(seed->stubs()[s]->qOverPt()) < qOverPt ) qOverPt = seed->stubs()[s]->qOverPt();
    }
  }

  float phi0 = ( seed->phi()+ seed->dphi() );
  float z0 = 0;
  float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));

  const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
  const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
  const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda

 
}
  */
L1fittedTrack L1KalmanComb::fitClusteredTrack( const L1track3D& l1track3D ){
  resetStates();
  deleteStubClusters();
  numUpdateCalls_ = 0;

  minStubLayersRed_ = Utility::numLayerCut("FIT", getSettings(), l1track3D.iPhiSec(), l1track3D.iEtaReg(), fabs(l1track3D.qOverPt()), l1track3D.eta());
  const TP* tpa(0);
  if( l1track3D.getMatchedTP() ){
    tpa = l1track3D.getMatchedTP();
  }
  tpa_ = tpa;

  const vector<const StubCluster*> cls = l1track3D.getStubClusters();

  //Kalman Filter
  std::vector<const kalmanState *> cands = doKF( l1track3D, cls, tpa );

  if( cands.size() ) {
    const kalmanState *cand = cands[0];

    //cout<<"Final KF candidate eta="<<cand->candidate().iEtaReg()<<" ns="<<cand->nSkippedLayers()<<" klid="<<cand->nextLayer()-1<<" n="<<cand->nStubLayers()<<endl;

    // Get track helix params.
    std::map<std::string, double> trackParams = getTrackParams(cand);

    L1fittedTrack returnTrk(getSettings(), l1track3D, cand->stubClusters(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], cand->chi2(), nPar_, true);

    bool consistentHLS = false;
    if (this->isHLS()) {
      unsigned int mBinHelixHLS, cBinHelixHLS;
      cand->getHLSextra(mBinHelixHLS, cBinHelixHLS, consistentHLS);
      if( getSettings()->kalmanDebugLevel() >= 3 ){
	// Check if (m,c) corresponding to helix params are correctly calculated by HLS code.
	bool HLS_OK = ((mBinHelixHLS == returnTrk.getCellLocationFit().first) && (cBinHelixHLS == returnTrk.getCellLocationFit().second));
	if (not HLS_OK) std::cout<<"WARNING HLS mBinHelix disagrees with C++:"
				 <<" (HLS,C++) m=("<<mBinHelixHLS<<","<<returnTrk.getCellLocationFit().first <<")"
				 <<" c=("<<cBinHelixHLS<<","<<returnTrk.getCellLocationFit().second<<")"<<endl;
      }
    }

    // Store supplementary info, specific to KF fitter.
    if(this->isHLS() && nPar_ == 4) {
      returnTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_, consistentHLS );
    } else {
      returnTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_ );
    }

    // If doing 5 parameter fit, optionally also calculate helix params & chi2 with beam-spot constraint applied,
    // and store inside L1fittedTrack object.
    if (getSettings()->kalmanAddBeamConstr()) {
      if (nPar_ == 5) {
	double chi2_bcon = 0.;
	std::map<std::string, double> trackParams_bcon = getTrackParams_BeamConstr(cand, chi2_bcon);
	returnTrk.setBeamConstr(trackParams_bcon["qOverPt"], trackParams_bcon["phi0"], chi2_bcon);
      }
    }

    // Currently not setup for full CKF
    // Fitted track params must lie in same sector as HT originally found track in.
    /*	if ( ! ( getSettings()->hybrid() ) || ( settings_->runFullKalman() ) ) { // consistentSector() function not yet working for Hybrid.
	if (! returnTrk.consistentSector()) {
	L1fittedTrack failedTrk(getSettings(), l1track3D, cand->stubClusters(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], cand->chi2(), nPar_, false);
	if(this->isHLS() && nPar_ == 4) {
	failedTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_, consistentHLS );
	} else {
	failedTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_ );
	}
	fittedTracks.push_back(failedTrk);
	}
	}
    */	

    //candidate dump
    if( getSettings()->kalmanDebugLevel() >= 3 ){
      cout << "------------------------------------" << endl;
      if( tpa && tpa->useForAlgEff() ){
	cout << "TP for eff. : index " << tpa->index() << endl;
      }
      cout << "Candidate : " << endl; 
      if( tpa && tpa->useForAlgEff() && returnTrk.getPurity() != 1 ){
	cout << "The candidate is not pure" << endl;
      }
      cand->dump( cout, tpa, true );
      cout << "------------------------------------" << endl;
    }
			
    //fill histograms for the selected state with TP for algEff
    if( getSettings()->kalmanFillInternalHists() ) fillCandHists( *cand, tpa );

    return returnTrk;

  }
  // If cands.size() < 1
  else {
    if (getSettings()->kalmanDebugLevel() >= 1) {
      bool goodTrack =  ( tpa && tpa->useForAlgEff() ); // Matches truth particle.
      if(goodTrack) {
	// Debug printout for Mark to understand why tracks are lost.

	int tpin=tpa->index();				
	cout<<"TRACK LOST: eta="<<l1track3D.iEtaReg()<<" pt="<<l1track3D.pt()<<" tp="<<tpin<<endl;
				
	for( auto stubCluster : cls ){
	  cout<<"    Stub: lay_red="<<stubCluster->layerIdReduced()<<" r="<<stubCluster->r()<<" z="<<stubCluster->z()<<"   assoc TPs =";
	  std::vector<const Stub *> stubs = stubCluster->stubs();
	  for( auto stub : stubs ){
	    for (const TP* tp_i : stub->assocTPs())  cout<<" "<<tp_i->index();
	    cout<<endl;
	    if(stub->assocTPs().size()==0) cout<<" none"<<endl;
	  }
	}
	cout<<"---------------------"<<endl;
	/*				
					for( it_last = last_states.begin(); it_last != last_states.end(); it_last++ ){
					const kalmanState *state = *it_last;
				
					//std::map<std::string, double> trackParams = getTrackParams(state);
					//L1fittedTrack returnTrk(getSettings(), l1track3D, state->stubs(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], state->chi2(), nPar_, true);
				
				
					std::vector<const Stub *> sstubs = state->stubs();
					for( auto stub : sstubs ){
				
					for (const TP* tp_i : stub->assocTPs()) {
					cout<<tp_i->index()<<endl;
					}
				
					cout<<stub->r()<<" "<<stub->z()<<" "<<state->nStubLayers()<<endl;
					}
				
					cout<<"---------------------"<<endl;
				
					}
	*/
	cout<<"====================="<<endl;
      }
    } // End debug level >=1
			
    //dump on the missed TP for efficiency calculation.
    if( getSettings()->kalmanDebugLevel() >= 3 ){
      if( tpa && tpa->useForAlgEff() ){
	cout << "TP for eff. missed addr. index : " << tpa << " " << tpa->index() << endl;
	//          printStubClusters( cout, cls );
	//	    printStubs( cout, stubs );
      }
    } // End debug level 3

    L1fittedTrack returnTrk(getSettings(), l1track3D, l1track3D.getStubClusters(), l1track3D.qOverPt(), 0, l1track3D.phi0(), l1track3D.z0(), l1track3D.tanLambda(), 9999, nPar_, false);
    returnTrk.setInfoKF( 0, numUpdateCalls_ );

    return returnTrk;
        
  }

}

L1fittedTrack L1KalmanComb::fit(const L1track3D& l1track3D){

  iLastPhiSec_ = iCurrentPhiSec_;
  iLastEtaReg_ = iCurrentEtaReg_;
  iCurrentPhiSec_ = l1track3D.iPhiSec();
  iCurrentEtaReg_ = l1track3D.iEtaReg();
  resetStates();
  deleteStubClusters();
  numUpdateCalls_ = 0;

  // Get cut on number of layers including variation due to dead sectors, pt dependence etc.
  minStubLayersRed_ = Utility::numLayerCut("FIT", getSettings(), l1track3D.iPhiSec(), l1track3D.iEtaReg(), fabs(l1track3D.qOverPt()), l1track3D.eta());

  //TP
  const TP* tpa(0);
  if( l1track3D.getMatchedTP() ){
    tpa = l1track3D.getMatchedTP();
  }
  tpa_ = tpa;

  //dump flag
  static unsigned nthFit(0);
  nthFit++;
  if( getSettings()->kalmanDebugLevel() >= 3 && nthFit <= maxNfitForDump_ ){
    if( tpa ) dump_ = true; 
    else dump_ = false;
  }
  else dump_ = false;

  //stub list from L1track3D, sorted in layer order - necessary for clustering only
  std::vector<const Stub*> stubs = l1track3D.getStubs();
		
  sort(stubs.begin(), stubs.end(), orderStubsByLayer); // Unnecessary?

#ifdef MERGE_STUBS
  // Eliminate identical duplicate stubs.
  for(unsigned i=0; i < stubs.size(); i++ ){
    const Stub *stub_a = stubs.at(i);
    for(unsigned j=i+1; j < stubs.size(); j++ ){
      const Stub *stub_b = stubs.at(j);
      if( stub_a->r() == stub_b->r() && stub_a->phi() == stub_b->phi() && stub_a->z() == stub_b->z() ){
	stubs.erase( stubs.begin() + j ); 
	if( getSettings()->kalmanFillInternalHists() ) 
	  hndupStub_->Fill(1);
	j--;
      }
    }
  }
#endif

  std::vector<const StubCluster *> stubcls;

  // Loop over each layer of stubs
  for( unsigned j_layer=0; j_layer < 16; j_layer++ ){

    // Create a vector of all the stubs on layer j_lauyer
    std::vector<const Stub *> layer_stubs;
    for(unsigned i=0; i < stubs.size(); i++ ){
      const Stub *stub = stubs.at(i);
      if( stub->layerId() == LayerId[j_layer] ){
	layer_stubs.push_back( stub );
      }
    }

#ifdef MERGE_STUBS
    if( LayerId[j_layer] < 10 ) 
      sort( layer_stubs.begin(), layer_stubs.end(), orderStubsByZ ); // barrel
    else
      sort( layer_stubs.begin(), layer_stubs.end(), orderStubsByR ); // endcap
#endif
    
    // If we are clustering all the stubs on the same layer ...
    if ( getSettings()->kalmanStubClustering() ) {

      std::vector<const Stub *> stubs_for_cls;

      // Loop over all the stubs in this individual layer
      for(unsigned i=0; i < layer_stubs.size(); i++ ){ // Stubs in single layer, ordered by z or r.

	stubs_for_cls.push_back(layer_stubs.at(i));

#ifdef MERGE_STUBS
	while( layer_stubs.at(i) != layer_stubs.back() ){
	  if( isOverlap( layer_stubs.at(i), layer_stubs.at(i+1), TYPE_NORMAL ) ){
	    stubs_for_cls.push_back( layer_stubs.at(i+1) );
	    if( getSettings()->kalmanFillInternalHists() ) 
	      hnmergeStub_->Fill(0);
	    i++;
	  }
	  else break;
	}
#endif
      } // new end of layers loop - new closing of brackets

      if( getSettings()->kalmanFillInternalHists() ) {

	if( tpa && tpa->useForAlgEff() ){

	  if( stubs_for_cls.size() > 1 ){
	    
	    std::set<const TP*> s_tps = stubs_for_cls.at(0)->assocTPs();
	    if( s_tps.find( tpa ) != s_tps.end() ){

	      const Stub *sa = stubs_for_cls.front();
	      const Stub *sb = stubs_for_cls.back();

	      double drphi = fabs( sa->r() * wrapRadian( sa->phi() - sectorPhi() ) - sb->r() * wrapRadian( sb->phi() - sectorPhi() ) ); 
	      double dz    = fabs( sa->z() - sb->z() );
	      double dr    = fabs( sa->r() - sb->r() );
	      TString hname;
	      if( LayerId[j_layer] < 10 ){

		hname = Form( "hBarrelStubMaxDistanceLayer%02d", LayerId[j_layer] );

		if( hBarrelStubMaxDistanceMap.find( hname ) == hBarrelStubMaxDistanceMap.end() ){
		  cout << hname << " does not exist." << endl;
		}
		else{
		  hBarrelStubMaxDistanceMap[hname]->Fill( drphi, dz );
		}
	      }
	      else{
		hname = Form( "hEndcapStubMaxDistanceRing%02d", sa->endcapRing()  );

		if( hEndcapStubMaxDistanceMap.find( hname ) == hEndcapStubMaxDistanceMap.end() ){
		  cout << hname << " does not exist." << endl;
		}
		else{
		  hEndcapStubMaxDistanceMap[hname]->Fill( drphi, dr );
		}
	      }
	    } // End if s.tps find
	  } // End if stubs for cluster size check
	} // End tpa && tpa->useForAlgEff
      } // End fill internal histos

      if( stubs_for_cls.size() < 1 ) continue;

      // dl error now disabled
      StubCluster *stbcl = new StubCluster( stubs_for_cls, sectorPhi(), 0 );
      stbcl_list_.push_back( stbcl );
      stubcls.push_back( stbcl );
      
      if( getSettings()->kalmanFillInternalHists() ) {
	if( !stbcl->barrel() ){
	  TString hname = Form( "hphiErrorRatioRing%d", stbcl->endcapRing() );
	  if( hphiErrorRatioMap.find(hname) == hphiErrorRatioMap.end() ){
	    cout << hname << " does not exist." << endl;
	  }
	  else{
	    hphiErrorRatioMap[hname]->Fill( fabs( stbcl->deltai() + 0.5 ), fabs( stbcl->dphi_dr() ) / stbcl->dphi_dl() );
	  }
	}
      }

      //  } // old closing of brackets for loop over layer stubs
    }

    // Else we are fitting each stub - i.e. one stub = one stub cluster
    else {
      for(unsigned i=0; i < layer_stubs.size(); i++ ){ // Stubs in single layer, ordered by z or r.

	std::vector<const Stub *> stubs_for_cls;
	stubs_for_cls.push_back(layer_stubs.at(i));

#ifdef MERGE_STUBS
	while( layer_stubs.at(i) != layer_stubs.back() ){
	  if( isOverlap( layer_stubs.at(i), layer_stubs.at(i+1), TYPE_NORMAL ) ){
	    stubs_for_cls.push_back( layer_stubs.at(i+1) );
	    if( getSettings()->kalmanFillInternalHists() ) 
	      hnmergeStub_->Fill(0);
	    i++;
	  }
	  else break;
	}
#endif

	if( getSettings()->kalmanFillInternalHists() ) {

	  if( tpa && tpa->useForAlgEff() ){

	    if( stubs_for_cls.size() > 1 ){

	      std::set<const TP*> s_tps = stubs_for_cls.at(0)->assocTPs();
	      if( s_tps.find( tpa ) != s_tps.end() ){

		const Stub *sa = stubs_for_cls.front();
		const Stub *sb = stubs_for_cls.back();

		double drphi = fabs( sa->r() * wrapRadian( sa->phi() - sectorPhi() ) - sb->r() * wrapRadian( sb->phi() - sectorPhi() ) ); 
		double dz    = fabs( sa->z() - sb->z() );
		double dr    = fabs( sa->r() - sb->r() );
		TString hname;
		if( LayerId[j_layer] < 10 ){

		  hname = Form( "hBarrelStubMaxDistanceLayer%02d", LayerId[j_layer] );

		  if( hBarrelStubMaxDistanceMap.find( hname ) == hBarrelStubMaxDistanceMap.end() ){
		    cout << hname << " does not exist." << endl;
		  }
		  else{
		    hBarrelStubMaxDistanceMap[hname]->Fill( drphi, dz );
		  }
		}
		else{
		  hname = Form( "hEndcapStubMaxDistanceRing%02d", sa->endcapRing()  );

		  if( hEndcapStubMaxDistanceMap.find( hname ) == hEndcapStubMaxDistanceMap.end() ){
		    cout << hname << " does not exist." << endl;
		  }
		  else{
		    hEndcapStubMaxDistanceMap[hname]->Fill( drphi, dr );
		  }
		}
	      }
	    }
	  }
	}

	// dl error now disabled
	StubCluster *stbcl = new StubCluster( stubs_for_cls, sectorPhi(), 0 );
	stbcl_list_.push_back( stbcl );
	stubcls.push_back( stbcl );

	if( getSettings()->kalmanFillInternalHists() ) {
	  if( !stbcl->barrel() ){
	    TString hname = Form( "hphiErrorRatioRing%d", stbcl->endcapRing() );
	    if( hphiErrorRatioMap.find(hname) == hphiErrorRatioMap.end() ){
	      cout << hname << " does not exist." << endl;
	    }
	    else{
	      hphiErrorRatioMap[hname]->Fill( fabs( stbcl->deltai() + 0.5 ), fabs( stbcl->dphi_dr() ) / stbcl->dphi_dl() );
	    }
	  }
	}
      }
    }
    if( getSettings()->kalmanFillInternalHists() ){ 
      if( tpa && tpa->useForAlgEff() ){
	hTrackEta_->Fill( tpa->eta() ); 
	static set<const TP *> set_tp;
	if( iCurrentPhiSec_ < iLastPhiSec_ && iCurrentEtaReg_ < iLastEtaReg_ ) set_tp.clear();
	if( set_tp.find( tpa ) == set_tp.end() ){
	  hUniqueTrackEta_->Fill( tpa->eta() );
	}
	set_tp.insert( tpa );
      }
    }
    
  } // end else

  if( getSettings()->kalmanFillInternalHists() ){ 
    if( tpa && tpa->useForAlgEff() ){
      hTrackEta_->Fill( tpa->eta() ); 
      static set<const TP *> set_tp;
      if( iCurrentPhiSec_ < iLastPhiSec_ && iCurrentEtaReg_ < iLastEtaReg_ ) set_tp.clear();
      if( set_tp.find( tpa ) == set_tp.end() ){
	hUniqueTrackEta_->Fill( tpa->eta() );
      }
      set_tp.insert( tpa );
    }
  }


  //track information dump
  if( getSettings()->kalmanDebugLevel() >= 1 ){

    std::cout << "===============================================================================" << endl;
    std::cout << "Input track cand: [phiSec,etaReg]=[" << l1track3D.iPhiSec() << "," << l1track3D.iEtaReg() << "]";
    std::cout <<" HT(m,c)=("<<l1track3D.getCellLocationHT().first << "," 
	                        <<l1track3D.getCellLocationHT().second << ") q/pt="
	      <<l1track3D.qOverPt()<<" tanL="<<l1track3D.tanLambda()<< " z0="<<l1track3D.z0()<< " phi0="<<l1track3D.phi0()
                                <<" nStubs="<<l1track3D.getNumStubs()<<std::endl;
    if (not getSettings()->hybrid()) printTP( cout, tpa );
    if( getSettings()->kalmanDebugLevel() >= 2 ){
      printStubLayers( cout, stubs );
      printStubClusters( cout, stubcls );
    }
  }

  //Kalman Filter
  std::vector<const kalmanState *> cands = doKF( l1track3D, stubcls, tpa );

 
  //return L1fittedTrk for the selected state (if KF produced one it was happy with).
  if( cands.size() ) {

    const kalmanState *cand = cands[0];

    //cout<<"Final KF candidate eta="<<cand->candidate().iEtaReg()<<" ns="<<cand->nSkippedLayers()<<" klid="<<cand->nextLayer()-1<<" n="<<cand->nStubLayers()<<endl;

    // Get track helix params.
    std::map<std::string, double> trackParams = getTrackParams(cand);

    L1fittedTrack returnTrk(getSettings(), l1track3D, cand->stubClusters(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], cand->chi2(), nPar_, true);

    bool consistentHLS = false;
    if (this->isHLS()) {
      unsigned int mBinHelixHLS, cBinHelixHLS;
      cand->getHLSextra(mBinHelixHLS, cBinHelixHLS, consistentHLS);
      if( getSettings()->kalmanDebugLevel() >= 3 ){
        // Check if (m,c) corresponding to helix params are correctly calculated by HLS code.
        bool HLS_OK = ((mBinHelixHLS == returnTrk.getCellLocationFit().first) && (cBinHelixHLS == returnTrk.getCellLocationFit().second));
        if (not HLS_OK) std::cout<<"WARNING HLS mBinHelix disagrees with C++:"
                                 <<" (HLS,C++) m=("<<mBinHelixHLS<<","<<returnTrk.getCellLocationFit().first <<")"
                                 <<" c=("<<cBinHelixHLS<<","<<returnTrk.getCellLocationFit().second<<")"<<endl;
      }
    }

    // Store supplementary info, specific to KF fitter.
    if(this->isHLS() && nPar_ == 4) {
      returnTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_, consistentHLS );
    } else {
      returnTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_ );
    }

    // If doing 5 parameter fit, optionally also calculate helix params & chi2 with beam-spot constraint applied,
    // and store inside L1fittedTrack object.
    if (getSettings()->kalmanAddBeamConstr()) {
      if (nPar_ == 5) {
	double chi2_bcon = 0.;
	std::map<std::string, double> trackParams_bcon = getTrackParams_BeamConstr(cand, chi2_bcon);
	returnTrk.setBeamConstr(trackParams_bcon["qOverPt"], trackParams_bcon["phi0"], chi2_bcon);
      }
    }

    // Fitted track params must lie in same sector as HT originally found track in.
    if (! getSettings()->hybrid() ) { // consistentSector() function not yet working for Hybrid.
      if (! returnTrk.consistentSector()) {
        L1fittedTrack failedTrk(getSettings(), l1track3D, cand->stubClusters(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], cand->chi2(), nPar_, false);
        if(this->isHLS() && nPar_ == 4) {
          failedTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_, consistentHLS );
        } else {
          failedTrk.setInfoKF( cand->nSkippedLayers(), numUpdateCalls_ );
        }
        return failedTrk;
      }
    }

    //candidate dump
    if( getSettings()->kalmanDebugLevel() >= 3 ){
      cout << "------------------------------------" << endl;
      if( tpa && tpa->useForAlgEff() ){
	cout << "TP for eff. : index " << tpa->index() << endl;
      }
      cout << "Candidate : " << endl; 
      if( tpa && tpa->useForAlgEff() && returnTrk.getPurity() != 1 ){
	cout << "The candidate is not pure" << endl;
      }
      cand->dump( cout, tpa, true );
      cout << "------------------------------------" << endl;
    }
			
    //fill histograms for the selected state with TP for algEff
    if( getSettings()->kalmanFillInternalHists() ) fillCandHists( *cand, tpa );
			
    return returnTrk;

  } else {

    if (getSettings()->kalmanDebugLevel() >= 1) {
      bool goodTrack =  ( tpa && tpa->useForAlgEff() ); // Matches truth particle.
      if(goodTrack) {
	// Debug printout for Mark to understand why tracks are lost.

	int tpin=tpa->index();				
	cout<<"TRACK LOST: eta="<<l1track3D.iEtaReg()<<" pt="<<l1track3D.pt()<<" tp="<<tpin<<endl;
				
	for( auto stubCluster : stubcls ){
	  cout<<"    Stub: lay_red="<<stubCluster->layerIdReduced()<<" r="<<stubCluster->r()<<" z="<<stubCluster->z()<<"   assoc TPs =";
	  std::vector<const Stub *> stubs = stubCluster->stubs();
	  for( auto stub : stubs ){
	    for (const TP* tp_i : stub->assocTPs())  cout<<" "<<tp_i->index();
	    cout<<endl;
	    if(stub->assocTPs().size()==0) cout<<" none"<<endl;
	  }
	}
	cout<<"---------------------"<<endl;
	/*				
					for( it_last = last_states.begin(); it_last != last_states.end(); it_last++ ){
					const kalmanState *state = *it_last;
				
					//std::map<std::string, double> trackParams = getTrackParams(state);
					//L1fittedTrack returnTrk(getSettings(), l1track3D, state->stubs(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], state->chi2(), nPar_, true);
				
				
					std::vector<const Stub *> sstubs = state->stubs();
					for( auto stub : sstubs ){
				
					for (const TP* tp_i : stub->assocTPs()) {
					cout<<tp_i->index()<<endl;
					}
				
					cout<<stub->r()<<" "<<stub->z()<<" "<<state->nStubLayers()<<endl;
					}
				
					cout<<"---------------------"<<endl;
				
					}
	*/
	cout<<"====================="<<endl;
      }
    }
			
    //dump on the missed TP for efficiency calculation.
    if( getSettings()->kalmanDebugLevel() >= 3 ){
      if( tpa && tpa->useForAlgEff() ){
	cout << "TP for eff. missed addr. index : " << tpa << " " << tpa->index() << endl;
	printStubClusters( cout, stubcls );
	printStubs( cout, stubs );
      }
    }

    L1fittedTrack returnTrk(getSettings(), l1track3D, l1track3D.getStubs(), l1track3D.qOverPt(), 0, l1track3D.phi0(), l1track3D.z0(), l1track3D.tanLambda(), 9999, nPar_, false);
    returnTrk.setInfoKF( 0, numUpdateCalls_ );
    return returnTrk;
  }
  
}

std::vector <L1fittedTrack> L1KalmanComb::findAndFit(const vector<const Stub*> inputStubs, const unsigned int iPhiSec, const unsigned int iEtaReg, 
		  const float etaMinSector, const float etaMaxSector, const float phiCentreSector ){

  iLastPhiSec_ = iCurrentPhiSec_;
  iLastEtaReg_ = iCurrentEtaReg_;
  iCurrentEtaReg_ = iEtaReg;
  iCurrentPhiSec_ = iPhiSec;
  resetStates();
  deleteStubClusters();
  numUpdateCalls_ = 0;

  float phiMinSector    = 2.*M_PI * (float(iCurrentPhiSec_)) / (float(getSettings()->numPhiSectors())) - M_PI;
  float phiMaxSector    = 2.*M_PI * (1.0 + float(iCurrentPhiSec_)) / (float(getSettings()->numPhiSectors())) - M_PI;

  // Calculate output opto-link ID from HT, assuming there is no MUX stage.
  // Required for the L1track3D constructor
  const unsigned numPhiSecPerOct =  settings_->numPhiSectors() / settings_->numPhiOctants();
  const unsigned optoLinkID = iCurrentEtaReg_ * numPhiSecPerOct + iCurrentPhiSec_;

  vector <L1track3D> trackCandidates; // Initialise vector of L1track3D candidates
  vector <L1fittedTrack> fittedTracks; // Initialise vector of L1fittedTracks to be returned

  // Time to work magic on input stubs

  const unsigned int seedingOption {settings_->kalmanSeedingOption()};
 
  if ( seedingOption == 10 ) {

    vector<const StubCluster*> seedClusters;
    vector<const StubCluster*> otherClusters;

    // number of phi and eta bins
    unsigned int nBinsKalmanSeedPhiAxis_  = settings_->kalmanSeedNbinsPhiAxis();
    unsigned int nBinsKalmanSeedEtaAxis_  = settings_->kalmanSeedNbinsEtaAxis();

    // init matrix
    matrix< vector<const Stub*> > kfStubArray_ (nBinsKalmanSeedPhiAxis_, nBinsKalmanSeedEtaAxis_);

    // binning increments
    float phiInc = (phiMaxSector - phiMinSector) / float(nBinsKalmanSeedPhiAxis_);
    float etaInc = (etaMaxSector - etaMinSector) / float(nBinsKalmanSeedEtaAxis_);

    // Fill kf seed stub array
    for ( auto stub : inputStubs ) {

      float stubPhi = stub->phi();
      float stubEta = stub->eta();

      unsigned int phiBin = std::ceil( (stubPhi-phiMinSector)/phiInc ) - 1;
      unsigned int etaBin = std::ceil( (stubEta-etaMinSector)/etaInc ) - 1;

      // overflow given uncert from poor bend resolution
      if ( phiBin >= nBinsKalmanSeedPhiAxis_ ) phiBin = nBinsKalmanSeedPhiAxis_ - 1;
      if ( etaBin >= nBinsKalmanSeedEtaAxis_ ) etaBin = nBinsKalmanSeedEtaAxis_ - 1;
      // underflow given uncert from poor bend resolution
//      if ( phiBin < 0 ) phiBin = 0;
//      if ( etaBin < 0 ) etaBin = 0;

      vector<const Stub*>& arrayStubs = kfStubArray_(phiBin,etaBin);
      arrayStubs.push_back(stub);
    }

    // Create stub clusters for each kf stub array cell
    std::vector<const StubCluster *> stubcls;

    for ( unsigned int phiBin = 0; phiBin != nBinsKalmanSeedPhiAxis_; phiBin++ ) {
      for ( unsigned int etaBin = 0; etaBin != nBinsKalmanSeedEtaAxis_; etaBin++ ) {
        // Access stubs in each KF array bin
        const vector<const Stub*>& arrayStubs = kfStubArray_(phiBin,etaBin);

        // create vector of stubs for each layer
        for ( unsigned j_layer=0; j_layer < 16; j_layer++ ){

          std::vector<const Stub *> layer_stubs;
          for(unsigned i=0; i < arrayStubs.size(); i++ ){
            const Stub *stub = arrayStubs.at(i);
            if( stub->layerId() == LayerId[j_layer] ){
              layer_stubs.push_back( stub );
            } // pushed back stub to layer vector
          } // end loop over stubs

	  if ( layer_stubs.size() == 0 ) continue; // if no stubs in layer, do not make stub cluster! 

          StubCluster *stbcl = new StubCluster( layer_stubs, sectorPhi(), 0);
          stubcls.push_back( stbcl );

        } // end layer loop
      } // end eta loop
    } // end phi loop

    // Create l1track3D seeds 
    for ( auto cls : stubcls ) {
      unsigned int layerId = cls->layerId();
      if ( layerId == 1 || layerId == 11 || layerId == 21 ) seedClusters.push_back(cls);
      else if ( layerId != 1 ) otherClusters.push_back(cls);
    }

    for ( auto seed : seedClusters ) {

      float qOverPt = seed->stubs()[0]->qOverPt();

      if ( seed->stubs().size() > 1 ) {
        for ( unsigned s = 1; s != seed->stubs().size(); s++ ) {
          if ( fabs(seed->stubs()[s]->qOverPt()) < qOverPt ) qOverPt = seed->stubs()[s]->qOverPt();
        }
      }


      float phi0 = ( seed->phi()+ seed->dphi() );
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));

      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda

      vector<const StubCluster*> cls {seed};
      cls.insert( cls.end(), otherClusters.begin(), otherClusters.end() );

      L1track3D l1track3D(getSettings(), cls, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);

      // Now cand is all setup ... do KF!
      fittedTracks.push_back( L1KalmanComb::fitClusteredTrack(l1track3D) );
    }

    // Fit each seed and return fitted track!
    return fittedTracks;
  }

  // seeding (no clustering) using layer 1+2+3
  else if ( seedingOption == 5 ) {

    vector<const Stub*> layer1Stubs;
    vector<const Stub*> layer2Stubs;
    vector<const Stub*> layer3Stubs;
    vector<const Stub*> otherStubs;

    // Default seeding option
    for ( auto stub : inputStubs ) {
      unsigned int reducedStubLayer = stub->layerIdReduced();
      if ( reducedStubLayer == 1 ) layer1Stubs.push_back(stub);
      else if ( reducedStubLayer == 2 ) layer2Stubs.push_back(stub);
      else if ( reducedStubLayer == 2 ) layer3Stubs.push_back(stub);
      else otherStubs.push_back(stub);
    }

    // Create seeds from layer 1 + all layer 2+ stubs
    for ( auto stub : layer1Stubs ) {
      float qOverPt = stub->qOverPt();
      float phi0 = stub->beta();
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));
      
      vector<const Stub*> stubs {stub};
      stubs.insert( stubs.end(), layer2Stubs.begin(), layer2Stubs.end() );
      stubs.insert( stubs.end(), layer3Stubs.begin(), layer3Stubs.end() );
      stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );
      
      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda
      
      L1track3D l1Trk3D(getSettings(), stubs, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);
      trackCandidates.push_back(l1Trk3D);
    }

    // Create seeds from layer 2 + all layer 3+ stubs
    for ( auto stub : layer2Stubs ) {
      float qOverPt = stub->qOverPt();
      float phi0 = stub->beta();
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));
      
      vector<const Stub*> stubs {stub};
      stubs.insert( stubs.end(), layer3Stubs.begin(), layer3Stubs.end() );
      stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );
      
      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda
      
      L1track3D l1Trk3D(getSettings(), stubs, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);
      trackCandidates.push_back(l1Trk3D);
    }

    // Create seeds from layer 3+ stubs
    for ( auto stub : layer3Stubs ) {
      float qOverPt = stub->qOverPt();
      float phi0 = stub->beta();
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));
      
      vector<const Stub*> stubs {stub};
      stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );
      
      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda
      
      L1track3D l1Trk3D(getSettings(), stubs, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);
      trackCandidates.push_back(l1Trk3D);
    }

    for ( auto trk : trackCandidates ) {
      fittedTracks.push_back(L1KalmanComb::fit(trk));
    }
    return fittedTracks;

  }

  // Seeding (no clustering) from layers 1+2  
  else if ( seedingOption == 1 ) {

    vector<const Stub*> seedStubs;

    vector<const Stub*> otherStubs;
    vector<const Stub*> otherStubs2;

    vector<const Stub*> posDiskStubs;
    vector<const Stub*> negDiskStubs;
    vector<const Stub*> posDiskStubs2;
    vector<const Stub*> negDiskStubs2;

    // Default seeding option
    for ( auto stub : inputStubs ) {

      unsigned int stubLayer = stub->layerId();

      if ( stubLayer == 1 || stubLayer == 11 || stubLayer == 21 ) seedStubs.push_back(stub);
      if ( stubLayer == 2 || stubLayer == 12 || stubLayer == 22 ) seedStubs.push_back(stub);
 
      if ( stubLayer != 1 ) otherStubs.push_back(stub);
      if ( stubLayer > 2 ) otherStubs2.push_back(stub);
      if ( stubLayer > 11 && stubLayer < 20 ) posDiskStubs.push_back(stub);
      if ( stubLayer > 12 && stubLayer < 20 ) posDiskStubs2.push_back(stub);
      if ( stubLayer > 21 ) negDiskStubs.push_back(stub);
      if ( stubLayer > 22 ) negDiskStubs2.push_back(stub);

   }

    // Create seeds 
    for ( auto stub : seedStubs ) {
      float qOverPt = stub->qOverPt();
      float phi0 = stub->beta();
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));
      
      vector<const Stub*> stubs {stub};
      const unsigned int layerId {stub->layerId()};

      if ( layerId == 1 )  stubs.insert( stubs.end(), otherStubs.begin(),    otherStubs.end() );
      if ( layerId == 11 ) stubs.insert( stubs.end(), posDiskStubs.begin(),  posDiskStubs.end() );
      if ( layerId == 21 ) stubs.insert( stubs.end(), negDiskStubs.begin(),  negDiskStubs.end() );
      if ( layerId == 2 )  stubs.insert( stubs.end(), otherStubs2.begin(),   otherStubs2.end() );
      if ( layerId == 12 ) stubs.insert( stubs.end(), posDiskStubs2.begin(), posDiskStubs2.end() );
      if ( layerId == 22 ) stubs.insert( stubs.end(), negDiskStubs2.begin(), negDiskStubs2.end() );
      
      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda
      
      L1track3D l1Trk3D(getSettings(), stubs, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);
      trackCandidates.push_back(l1Trk3D);
    }

    for ( auto trk : trackCandidates ) {
      fittedTracks.push_back(L1KalmanComb::fit(trk));
    }
    return fittedTracks;

  }
  
  // option 0 - seeding with layer 1 only and no clustering
  else {

    vector<const Stub*> seedStubs;
    vector<const Stub*> otherStubs;
    vector<const Stub*> posDiskStubs;
    vector<const Stub*> negDiskStubs;

    // Default seeding option
    for ( auto stub : inputStubs ) {

      const unsigned int stubLayer {stub->layerId()};
      if ( stubLayer == 1 || stubLayer == 11 || stubLayer == 21 ) seedStubs.push_back(stub);
      if ( stubLayer != 1 ) otherStubs.push_back(stub);
//      otherStubs.push_back(stub);
      if ( stubLayer > 11 && stubLayer < 20 ) posDiskStubs.push_back(stub);
      if ( stubLayer > 21 ) negDiskStubs.push_back(stub);

    }

    // create single stub seed

    for ( auto stub : seedStubs ) {
      float qOverPt = stub->qOverPt();
      float phi0 = stub->beta();
      float z0 = 0;
      float tan_lambda = 0.5*(1/tan(2*atan(exp(-etaMinSector))) + 1/tan(2*atan(exp(-etaMaxSector))));
      
      unsigned int stubLayer = stub->layerId();

      vector<const Stub*> stubs {stub};
      if ( stubLayer == 1 ) stubs.insert( stubs.end(), otherStubs.begin(), otherStubs.end() );
      if ( stubLayer == 11 ) stubs.insert( stubs.end(), posDiskStubs.begin(), posDiskStubs.end() );
      if ( stubLayer == 21 ) stubs.insert( stubs.end(), negDiskStubs.begin(), negDiskStubs.end() );

      const pair<unsigned int, unsigned int> cellLocation { make_pair(0,0) }; // No HT seed location - use dummy location
      const pair< float, float > helixParamsRphi { make_pair(qOverPt, phi0) }; // q/Pt + phi0
      const pair< float, float > helixParamsRz { make_pair(z0, tan_lambda) }; // z0, tan_lambda
      
      L1track3D l1Trk3D(getSettings(), stubs, cellLocation, helixParamsRphi, helixParamsRz, iCurrentPhiSec_, iCurrentEtaReg_, optoLinkID, false);
      trackCandidates.push_back(l1Trk3D);
    }
    
    for ( auto trk : trackCandidates ) {
      fittedTracks.push_back(L1KalmanComb::fit(trk));
    }
    return fittedTracks;
  }
}

std::vector<const kalmanState *> L1KalmanComb::doKF( const L1track3D& l1track3D, const std::vector<const StubCluster *> &stubClusters, const TP *tpa ){

#ifdef RECALC_DEBUG
  cout<<"FITTER new track: HT cell=("<<l1track3D.getCellLocationHT().first<<","<<l1track3D.getCellLocationHT().second<<")"<<endl;
#endif

  // output container (contains 0 or 1 states).
  std::vector<const kalmanState *> finished_states;

  std::map<unsigned int, const kalmanState *, std::greater<unsigned int> > best_state_by_nstubs; // Best state (if any) for each viable no. of stubs on track value. 
	
  // seed helix params & their covariance.
  std::vector<double> x0 = seedx(l1track3D);
  TMatrixD pxx0 = seedP(l1track3D);
  TMatrixD K( nPar_, 2 );
  TMatrixD dcov( 2, 2 );
	
  const kalmanState *state0 = mkState( l1track3D, 0, 0, 0, 0, x0, pxx0, K, dcov, 0, 0 );
	
  if( getSettings()->kalmanFillInternalHists() ) fillSeedHists( state0, tpa );
	
	
  // internal containers - i.e. the state FIFO. Contains estimate of helix params in last/next layer, with multiple entries if there were multiple stubs, yielding multiple states.
  std::vector<const kalmanState *> new_states;
  std::vector<const kalmanState *> prev_states;
  prev_states.push_back( state0 );
	

  // === Layer Mapping (i.e. layer order in which stubs should be processed) ===

  // index across is ian encoded layer id (where barrel layers=1,2,7,5,4,3 & endcap wheels=3,4,5,6,7 & 0 never occurs)
  // index down is eta reg
  // element is kalman layer where 7 is invalid
  // assumes we are in barrel, endcap adjustments later
  // should really be defined once in constructor 
  
  unsigned layerMap[18][8] = 
    { 
      { 7,  0,  7,  1,  2,  3,  4,  5 },
      { 7,  0,  7,  1,  2,  3,  4,  5 },
      { 7,  0,  1,  2,  3,  4,  5,  5 },
      { 7,  0,  1,  2,  3,  4,  5,  2 },
      { 7,  0,  1,  3,  4,  3,  6,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  5,  4,  3,  7,  2 },
      { 7,  0,  1,  3,  4,  3,  6,  2 },
      { 7,  0,  1,  2,  3,  4,  5,  2 },
      { 7,  0,  1,  2,  3,  4,  5,  5 },
      { 7,  0,  7,  1,  2,  3,  4,  5 },
      { 7,  0,  7,  1,  2,  3,  4,  5 },
    };
  
  // arrange stubs into Kalman layers according to eta region
  int etaReg = l1track3D.iEtaReg();
  std::map<int, std::vector<const StubCluster *> > layerStubs;

  // Get dead layers, if any.
  // They are assumed to be idetnical to those defined in StubKiller.cc
  bool remove2PSCut = getSettings()->kalmanRemove2PScut();
  set<unsigned> kalmanDeadLayers = getKalmanDeadLayers( layerMap, remove2PSCut );

  for( auto stubCluster : stubClusters ){
		
    int kalmanLayer = layerMap[etaReg][stubCluster->layerIdReduced()];

    if ( !stubCluster->barrel() ) {
			
      switch ( etaReg ) {
      case 3:
      case 14:
	if (stubCluster->layerIdReduced()==7) kalmanLayer = 6;
	break;
      case 4:
      case 13:
	if (stubCluster->layerIdReduced()==5) kalmanLayer = 5;
	break;
      case 5:
      case 12:
	if (stubCluster->layerIdReduced()==4) kalmanLayer = 5;
	break;
      default:
	break;
      }
			
    }
		
    if (layerStubs[kalmanLayer].size() < getSettings()->kalmanMaxStubsPerLayer()) {
      if (kalmanLayer != 7) layerStubs[kalmanLayer].push_back( stubCluster );
    }
  }

  // iterate using state->nextLayer() to determine next Kalman layer(s) to add stubs from
  const unsigned int maxIterations = 6;       // Increase if you want to allow 7 stubs per fitted track.
  for( unsigned iteration = 0; iteration < maxIterations; iteration++ ){   

    int combinations_per_iteration = 0;
		
    bool easy = (l1track3D.getNumStubs() < getSettings()->kalmanMaxStubsEasy());
    unsigned int kalmanMaxSkipLayers = easy ? getSettings()->kalmanMaxSkipLayersEasy() : getSettings()->kalmanMaxSkipLayersHard();
		
    // update each state from previous iteration (or seed) using stubs in next Kalman layer
    std::vector<const kalmanState *>::const_iterator i_state = prev_states.begin();
    for(; i_state != prev_states.end(); i_state++ ){ 
		
      const kalmanState *the_state = *i_state;
			

      unsigned layer = the_state->nextLayer();
      unsigned skipped = the_state->nSkippedLayers();

      // If this layer is known to be dead, skip to the next layer (layer+1)
      // The next_states_skipped will then look at layer+2
      // However, if there are stubs in this layer, then don't skip (e.g. our phi/eta boundaries might not line up exactly with a dead region)
      // Continue to skip until you reach a functioning layer (or a layer with stubs)
      unsigned nSkippedDeadLayers = 0;
      while ( kalmanDeadLayers.find(layer) != kalmanDeadLayers.end() && layerStubs[layer].size() == 0 ) {
	layer += 1;
	++nSkippedDeadLayers;
      }

      // containers for updated state+stub combinations
      std::vector<const kalmanState *> next_states;
      std::vector<const kalmanState *> next_states_skipped;

			
      // find stubs for this layer
      std::vector<const StubCluster *> stubs = layerStubs[layer]; // If layer > 6, this will return empty vector, so safe.

      // find stubs for next layer if we skip a layer, except when we are on the penultimate layer,
      // or we have exceeded the max skipped layers
      std::vector<const StubCluster *> next_stubs ;

      // If the next layer (layer+1) is a dead layer, then proceed to the layer after next (layer+2), if possible
      // Also note if we need to increase "skipped" by one more for these states
      unsigned nSkippedDeadLayers_nextStubs = 0;
      if ( skipped < kalmanMaxSkipLayers ) {
        if ( kalmanDeadLayers.find(layer+1) != kalmanDeadLayers.end()  && layerStubs[layer+1].size() == 0 ) {
	  next_stubs = layerStubs[layer+2];
	  nSkippedDeadLayers_nextStubs += 1;
        } else {
	  next_stubs = layerStubs[layer+1];
	}
      }

      // If track was not rejected by isGoodState() is previous iteration, failure here usually means the tracker ran out of layers to explore.
      // (Due to "kalmanLayer" not having unique ID for each layer within a given eta sector).
      if ( getSettings()->kalmanDebugLevel() >= 2 && best_state_by_nstubs.size() == 0 && stubs.size() == 0 && next_stubs.size() == 0) cout<<"State is lost by start of iteration "<<iteration<<" : #stubs="<<stubs.size()<<" #next_stubs="<<next_stubs.size()<<" layer="<<layer<<" eta="<<l1track3D.iEtaReg()<<endl;

      // If we skipped over a dead layer, only increment "skipped" after the stubs in next+1 layer have been obtained
      skipped += nSkippedDeadLayers;
		
      // check to guarantee no fewer than 2PS hits per state at iteration 1 (r<60cm)
      // iteration 0 will always include a PS hit, but iteration 1 could use 2S hits unless we include this
      if (iteration==1 && !remove2PSCut) {
	std::vector<const StubCluster *> temp_stubs;
	std::vector<const StubCluster *> temp_nextstubs;
	for (auto stub : stubs) {
	  if (stub->r()<60.0) temp_stubs.push_back(stub);
	}
	for (auto stub : next_stubs) {
	  if (stub->r()<60.0) temp_nextstubs.push_back(stub);
	}
	stubs = temp_stubs;
	next_stubs = temp_nextstubs;
      }

			
      combinations_per_iteration += stubs.size() + next_stubs.size();
			
			
      // loop over each stub in this layer and check for compatibility with this state
      for( unsigned i=0; i < stubs.size()  ; i++ ){
	
	const StubCluster * next_stubCluster = stubs[i];
				
	// Update helix params by adding this stub.
	const kalmanState * new_state = kalmanUpdate( skipped, layer+1, next_stubCluster, *the_state, tpa );
				
	if( getSettings()->kalmanFillInternalHists() ) fillStepHists( tpa, iteration, new_state );
				
	// Cut on track chi2, pt etc.
	if(isGoodState( *new_state ) ) next_states.push_back( new_state );
      }

      // loop over each stub in next layer if we skip, and check for compatibility with this state
      for( unsigned i=0; i < next_stubs.size()  ; i++ ){
	
	const StubCluster * next_stubCluster = next_stubs[i];
				
	const kalmanState * new_state = kalmanUpdate( skipped+1+nSkippedDeadLayers_nextStubs, layer+2+nSkippedDeadLayers_nextStubs, next_stubCluster, *the_state, tpa );
				
	if( getSettings()->kalmanFillInternalHists() ) fillStepHists( tpa, iteration, new_state );
				
	if(isGoodState( *new_state ) ) next_states_skipped.push_back( new_state );
      }		
			
      // post Kalman filter local sorting per state
      sort( next_states.begin(), next_states.end(), kalmanState::orderChi2);
      sort( next_states_skipped.begin(), next_states_skipped.end(), kalmanState::orderChi2);
			
			
      int i, max_states, max_states_skip;
			
      // If layer contained several stubs, so several states now exist, select only the best ones.
      // -- Disable this by setting to large values, as not used in latest KF firmware.
      // (But not too big as this wastes CPU).

      switch ( iteration ) {
      case 0:
	max_states = 15;
	max_states_skip = 15;
	break;
      case 1:
	max_states = 15;
	max_states_skip = 15;
	break;
      case 2:
	max_states = 15;
	max_states_skip = 15;
	break;
      case 3:
	max_states = 15;
	max_states_skip = 15;
	break;
      case 4:
	max_states = 15;
	max_states_skip = 15;
	break;
      case 5:
	max_states = 15;
	max_states_skip = 15;
	break;
      default:
	max_states = 15;
	max_states_skip = 15;
	break;
      }
			
			
      i = 0;
      for( auto state : next_states ){
					
	if( i < max_states ){
	  new_states.push_back( state );
	} else {
	  break;
	}
	i++;
	
      }
			
      i = 0; 
      for( auto state : next_states_skipped ){
	
	if( i < max_states_skip ){
	  new_states.push_back( state );
	} else {
	  break;
	}
	i++;
	
      }
			
    } //end of state loop


    if( getSettings()->kalmanFillInternalHists() ) {
      TString hname = Form( "hstubComb_itr%d", iteration );
      if( hstubCombMap.find(hname) == hstubCombMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else{
	hstubCombMap[hname]->Fill( combinations_per_iteration );
      }
    }
		
			 
    // copy new_states into prev_states for next iteration or end if we are on 
    // last iteration by clearing all states and making final state selection
		
    sort( new_states.begin(), new_states.end(), kalmanState::orderMinSkipChi2); // Sort by chi2*(skippedLayers+1)

    unsigned int nStubs = iteration + 1;
    // Success. We have at least one state that passes all cuts. Save best state found with this number of stubs.
    if (nStubs >= getSettings()->kalmanMinNumStubs() && new_states.size() > 0) best_state_by_nstubs[nStubs] = new_states[0]; 

    //if ( getSettings()->kalmanDebugLevel() >= 1 && best_state_by_nstubs.size() == 0 && new_states.size() == 0) cout<<"Track is lost by end iteration "<<iteration<<" : eta="<<l1track3D.iEtaReg()<<endl;

    if( nStubs == getSettings()->kalmanMaxNumStubs() ){ 
      // We're done.
      prev_states.clear();
      new_states.clear();
			
    } else {
			
      // Continue iterating.
      prev_states = new_states;
      new_states.clear(); 
			
    }
				
    /*
      int i = 0;
      bool found = false;
      for( auto best_state : best_states4 ){
			
      if( tpa && tpa->useForAlgEff() ) {
      std::map<std::string, double> trackParams = getTrackParams(best_state);
      L1fittedTrack returnTrk(getSettings(), l1track3D, best_state->stubs(), trackParams["qOverPt"], trackParams["d0"], trackParams["phi0"], trackParams["z0"], trackParams["t"], best_state->chi2(), nPar_, true);
      if (returnTrk.getNumMatchedLayers()>=4) {
      //temp_states.push_back(best_state);
      if(i==0) found = true;
      if (!found) cout<<"Lost this cand "<<i<<" "<<best_state->chi2()<<" "<<best_state->reducedChi2()<<" "<<best_state->path()<<" chose instead "<<best_states4[0]->chi2()<<" "<<best_states4[0]->reducedChi2()<<" "<<best_statesn4[0]->path()<<endl;
      }
      }*/
		
  }

  if (best_state_by_nstubs.size()) {
    // Select state with largest number of stubs.
    const kalmanState* stateFinal = best_state_by_nstubs.begin()->second; // First element has largest number of stubs.
    finished_states.push_back(stateFinal);
    if ( getSettings()->kalmanDebugLevel() >= 1 ) {
      cout<<"Track found! final state selection: nLay="<<stateFinal->nStubLayers()<<" etaReg="<<l1track3D.iEtaReg();
      std::map<std::string, double> y = getTrackParams( stateFinal );
      cout<<" q/pt="<<y["qOverPt"]<<" tanL="<<y["t"]<<" z0="<<y["z0"]<<" phi0="<<y["phi0"];
      cout<<" chosen from states:";
      for (const auto& p : best_state_by_nstubs) cout<<" "<<p.second->chi2()<<"/"<<p.second->nStubLayers();
      cout<<endl;
    }
  } else {
    if ( getSettings()->kalmanDebugLevel() >= 1 ) {
      cout<<"Track lost"<<endl;
    }
  }

  return finished_states;
}


//--- Update a helix state by adding a stub. 
//--- ("layer" is not the layer of the stub being added now, but rather the next layer that will be searched after this stub has been added).

const kalmanState *L1KalmanComb::kalmanUpdate( unsigned skipped, unsigned layer, const StubCluster *stubCluster, const kalmanState &state, const TP *tpa ){

  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "---------------" << endl;
    cout << "kalmanUpdate" << endl;
    cout << "---------------" << endl;
    printStubCluster( cout, stubCluster );
  }

  numUpdateCalls_++; // For monitoring, count calls to updator per track.

  // Helix params & their covariance.
  std::vector<double> xa     = state.xa();
  TMatrixD            cov_xa = state.pxxa(); 
  if( state.barrel() && !stubCluster->barrel() ){ 
    if( getSettings()->kalmanDebugLevel() >= 4 ) {
      cout << "STATE BARREL TO ENDCAP BEFORE " << endl;
      cout << "state : " << xa.at(0) << " " << xa.at(1) << " " << xa.at(2) << " " << xa.at(3) << endl;
      cout << "cov(x): " << endl; 
      cov_xa.Print();
    }
    barrelToEndcap( state.r(), stubCluster, xa, cov_xa );
    if( getSettings()->kalmanDebugLevel() >= 4 ){
      cout << "STATE BARREL TO ENDCAP AFTER " << endl;
      cout << "state : " << xa.at(0) << " " << xa.at(1) << " " << xa.at(2) << " " << xa.at(3) << endl;
      cout << "cov(x): " << endl; 
      cov_xa.Print();
    }
  }
  // Matrix to propagate helix params from one layer to next (=identity matrix).
  TMatrixD f = F(stubCluster, &state );
  TMatrixD ft(TMatrixD::kTransposed, f );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "f" << endl;
    f.Print();
    cout << "ft" << endl;
    ft.Print();
  }

  std::vector<double> fx = Fx( f, xa ); // Multiply matrices to get helix params at next layer.
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "fx = ["; 
    for( unsigned i = 0; i < nPar_; i++ ) cout << fx.at(i) << ", ";
    cout << "]" << endl;
  }

  std::vector<double> delta = residual(stubCluster, fx, state.candidate().qOverPt() );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "delta = " << delta[0] << ", " << delta[1] << endl;
  }

  // Derivative of predicted (phi,z) intercept with layer w.r.t. helix params.
  TMatrixD h = H(stubCluster);
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "h" << endl;
    h.Print();
  }


  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "previous state covariance" << endl;
    cov_xa.Print();
  }
  // Get contribution to helix parameter covariance from scattering (NOT USED).
  TMatrixD pxxm = PxxModel( &state, stubCluster );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "model xcov" << endl;
    pxxm.Print();
  }
  // Get covariance on helix parameters.
  TMatrixD pxcov = f * cov_xa * ft + pxxm;
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "forcast xcov + model xcov" << endl;
    pxcov.Print();
  }
  // Get hit position covariance matrix.
  TMatrixD dcov = PddMeas( stubCluster, &state );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "dcov" << endl;
    dcov.Print();
  }
  // Calculate Kalman Gain matrix.
  TMatrixD k = GetKalmanMatrix( h, pxcov, dcov );  
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "k" << endl;
    k.Print();
  }
	 
  std::vector<double> new_xa(nPar_);
  TMatrixD new_pxxa;
  GetAdjustedState( k, pxcov, fx, stubCluster, delta, new_xa, new_pxxa );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    if( nPar_ == 4 )
      cout << "adjusted x = " << new_xa[0] << ", " << new_xa[1] << ", " << new_xa[2] << ", " << new_xa[3] << endl;
    else if( nPar_ == 5 )
      cout << "adjusted x = " << new_xa[0] << ", " << new_xa[1] << ", " << new_xa[2] << ", " << new_xa[3] << ", " << new_xa[4] << endl;
    cout << "adjusted covx " << endl;
    new_pxxa.Print();
  }

  const kalmanState *new_state = mkState( state.candidate(), skipped, layer, stubCluster->layerId(), &state, new_xa, new_pxxa, k, dcov, stubCluster, 0 );
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "new state" << endl;
    new_state->dump( cout, tpa  );
  }


  return new_state;
}


double L1KalmanComb::calcChi2( const kalmanState &state )const{

  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "calcChi2 " << endl;
  }
  double chi2(0), chi2_p(0);

  if( state.last_state() ) {
    chi2 = state.last_state()->chi2();
			
    const StubCluster *stubCluster = state.stubCluster();
			
#ifdef RECALC_DEBUG
    unsigned int ID = (stubCluster != nullptr)  ?  stubCluster->stubs()[0]->index()  :  99999;
#endif

    if( stubCluster ){
	
      std::vector<double> delta = residual( stubCluster, state.last_state()->xa(), state.last_state()->candidate().qOverPt() );
      TMatrixD dcov = PddMeas( stubCluster, &state );
#ifdef RECALC_DEBUG
      cout<<"    FITTER SIGMA:      rphi="<<1000*sqrt(dcov(0,0))<<" rz="<<sqrt(dcov(1,1))<<" ID="<<ID<<endl;
#endif

      if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "dcov" << endl;
	dcov.Print();
	cout << "xcov" << endl;
	state.last_state()->pxxa().Print();
      }
      TMatrixD h = H(stubCluster);
      TMatrixD hxxh = HxxH( h, state.last_state()->pxxa() );
      if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "h" << endl;
	h.Print();
	cout << "hxcovh" << endl;
	hxxh.Print();
      }
      TMatrixD covR = dcov + hxxh;
      if( getSettings()->kalmanDebugLevel() >= 4 ){
	cout << "covR" << endl;
	covR.Print();
	cout << "---" << endl;
	cout << scientific << "delta = " << delta[0] << ", " << delta[1] << endl;
      }
      chi2_p = Chi2( covR, delta );  
	
    }
    chi2 += chi2_p;
#ifdef RECALC_DEBUG
    cout<<"  FITTER CHI2 UPDATE = "<<chi2<<" delta chi2="<<chi2_p<<" ID="<<ID<<endl;
#endif
  }

  return chi2;
}


double L1KalmanComb::Chi2( const TMatrixD &dcov, const std::vector<double> &delta, bool debug )const
{
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "dcov" << endl;
    dcov.Print();
  }

  if( dcov.Determinant() == 0 ) return 999;


  TMatrixD dcovi( dcov );
  dcovi.Invert();

  vector<double> tmp(2,0);
  for( int i=0; i < dcovi.GetNrows(); i++ ){ 
    for( int j=0; j < dcovi.GetNcols(); j++ ){ 
      tmp.at(j) += delta.at(i) * dcovi(i,j); 
    }
  }
  double chi2(0);
  for( int j=0; j < 2; j++ ){ 
    chi2 += tmp.at(j) * delta.at(j);
  }

#ifdef RECALC_DEBUG
  cout<<"    FITTER DELTA CHI2: rphi="<<dcovi(0,0)*delta.at(0)*delta.at(0)
      <<" rz="<<dcovi(1,1)*delta.at(1)*delta.at(1)<<endl;
#endif

  if( debug ){
    cout << "CHI SQUARE OUTPUT" << endl;
    cout << "cov" << endl;
    dcov.Print();
    cout << "cov inv" << endl;
    dcovi.Print();
    for( unsigned i=0; i < delta.size(); i++ ) cout << delta.at(i) << " ";
    cout << endl;
  }
  return chi2;
}


std::map<std::string, double> L1KalmanComb::getTrackParams( const L1KalmanComb *p, const kalmanState *state )
{
  return p->getTrackParams( state );
}


std::vector<double> L1KalmanComb::Hx( const TMatrixD &pH, const std::vector<double> &x )const
{
  std::vector<double> m( (unsigned) pH.GetNrows(), 0 );
  if( pH.GetNcols() != (int) x.size() ) { cerr << "Hx() : H and x have different dimensions" << endl; }
  else{

    for( int i=0; i < pH.GetNcols(); i++ ){ 
      for( int j=0; j < pH.GetNrows(); j++ ){ 
	m.at(j) += pH(j,i) * x.at(i);
      }
    }
  }
  return m;
}


std::vector<double> L1KalmanComb::Fx( const TMatrixD &pF, const std::vector<double> &x )const
{
  return Hx( pF, x );
}


TMatrixD L1KalmanComb::HxxH( const TMatrixD &pH, const TMatrixD &xx )const
{
  int nd = (unsigned) pH.GetNrows(); 
  TMatrixD tmp(nd,nPar_);
  TMatrixD mHxxH(nd,nd);
  if( pH.GetNcols() != xx.GetNcols() || pH.GetNcols() != xx.GetNrows() ) { cerr << "HxxH() : H and xx have different dimensions" << endl; }
  else{

    for( int i=0; i < pH.GetNrows(); i++ ){ 
      for( int j=0; j < xx.GetNrows(); j++ ){ 
	for( int k=0; k < xx.GetNcols(); k++ ){ 
	  tmp(i,k) += pH(i,j) * xx(j,k);
	}
      }
    }
    for( int i=0; i < tmp.GetNrows(); i++ ){ 
      for( int j=0; j < pH.GetNcols(); j++ ){ 
	for( int k=0; k < pH.GetNrows(); k++ ){ 
	  mHxxH(i,k) += tmp(i,j) * pH(k,j); 
	}
      }
    }
  }
  return mHxxH;

}


TMatrixD L1KalmanComb::GetKalmanMatrix( const TMatrixD &h, const TMatrixD &pxcov, const TMatrixD &dcov )const
{

  TMatrixD pxcovht(pxcov.GetNrows(),2);
  for( int i=0; i<pxcov.GetNrows(); i++ ){
    for( int j=0; j<pxcov.GetNcols(); j++ ){
      for( int k=0; k<h.GetNrows(); k++ ){
	pxcovht(i,k) += pxcov(i,j) * h(k,j);
      }
    }
  }
  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "pxcovht" << endl;
    pxcovht.Print();
  }

  TMatrixD tmp(dcov.GetNrows(), dcov.GetNcols() );
  TMatrixD hxxh = HxxH( h, pxcov );
  tmp = dcov + hxxh; 

  if( getSettings()->kalmanDebugLevel() >= 4 ){
    cout << "hxxh" << endl;
    hxxh.Print();
    cout << "dcov + hxxh " << endl;
    tmp.Print();
  }

  TMatrixD K( pxcovht.GetNrows(), tmp.GetNcols() );

  if(tmp.Determinant() == 0 ) return K; 
  tmp.Invert();

  for( int i=0; i<pxcovht.GetNrows(); i++ ){
    for( int j=0; j<pxcovht.GetNcols(); j++ ){
      for( int k=0; k<tmp.GetNcols(); k++ ){
	K(i,k)+=pxcovht(i,j)*tmp(j,k);
      }
    }
  }
  return K;
}


void L1KalmanComb::GetAdjustedState( const TMatrixD &K, const TMatrixD &pxcov, 
				     const std::vector<double> &x, const StubCluster *stubCluster, 
				     const std::vector<double>& delta,
				     std::vector<double> &new_x, TMatrixD &new_xcov )const
{
  TMatrixD h = H(stubCluster);

  for( int i=0; i < K.GetNrows(); i++ ){
    new_x.at(i) = x.at(i);  
    for( int j=0; j < K.GetNcols(); j++ ){
      new_x.at(i) += K(i,j) * delta.at(j);
    }
  }

  TMatrixD tmp(K.GetNrows(), h.GetNcols() );
  for( int i=0; i< K.GetNrows(); i++ ){
    tmp(i,i) = 1;
  }
  for( int i=0; i< K.GetNrows(); i++ ){
    for( int j=0; j< K.GetNcols(); j++ ){
      for( int k=0; k< h.GetNcols(); k++ ){
	tmp(i,k) += -1 * K(i,j) * h(j,k);
      }
    }
  }
  new_xcov.Clear();
  new_xcov.ResizeTo(pxcov.GetNrows(), pxcov.GetNcols());
  for( int i=0; i< tmp.GetNrows(); i++ ){
    for( int j=0; j< tmp.GetNcols(); j++ ){
      for( int k=0; k< pxcov.GetNcols(); k++ ){
	new_xcov(i,k) += tmp(i,j) * pxcov(j,k);
      }
    }
  }
}


void L1KalmanComb::resetStates()
{
  for( unsigned int i=0; i < state_list_.size(); i++ ){

    delete state_list_.at(i);
  }
  state_list_.clear();
}


const kalmanState *L1KalmanComb::mkState( const L1track3D &candidate, unsigned skipped, unsigned layer, unsigned layerId, const kalmanState *last_state, 
					  const std::vector<double> &x, const TMatrixD &pxx, const TMatrixD &K, const TMatrixD &dcov, const StubCluster* stubCluster, double chi2 )
{

  kalmanState *new_state = new kalmanState( candidate, skipped, layer, layerId, last_state, x, pxx, K, dcov, stubCluster, chi2, this, &getTrackParams );

  if( chi2 == 0 ){
    double new_state_chi2 = calcChi2( *new_state ); 
    new_state->setChi2( new_state_chi2 );
  }

  state_list_.push_back( new_state );
  return new_state;
}


std::vector<double> L1KalmanComb::residual(const StubCluster* stubCluster, const std::vector<double> &x, double candQoverPt )const{

  std::vector<double> vd = d(stubCluster); // Get (phi relative to sector, z) of hit.
  std::vector<double> hx = Hx( H(stubCluster), x ); // Ditto for intercept of helix with layer, in linear approximation.
  std::vector<double> delta(2);
  for( unsigned i=0; i<2; i++ ) delta.at(i) = vd.at(i) - hx.at(i);

  // Calculate higher order corrections to residuals.

  if (not getSettings()->kalmanHOdodgy()) {

    std::vector<double> correction = {0.,0.};

    float inv2R = (getSettings()->invPtToInvR()) * 0.5 * candQoverPt; // alternatively use x().at(0)
    float tanL = x.at(2);
    float z0 = x.at(3);

    float deltaS = 0.;
    if (getSettings()->kalmanHOhelixExp()) {
      // Higher order correction correction to circle expansion for improved accuracy at low Pt.
      double corr = stubCluster->r() * inv2R; 

      // N.B. In endcap 2S, this correction to correction[0] is exactly cancelled by the deltaS-dependent correction to it below.
      correction[0] += (1./6.)*pow(corr, 3); 

      deltaS = (1./6.)*(stubCluster->r())*pow(corr, 2);
      correction[1] -= deltaS * tanL;
    }

    if ( (not stubCluster->barrel()) && not (stubCluster->psModule())) {
      // These corrections rely on inside --> outside tracking, so r-z track params in 2S modules known.
      float rShift = (stubCluster->z() - z0)/tanL - stubCluster->r();

      // The above calc of rShift is approximate, so optionally check it with MC truth.
      // if (tpa_ != nullptr) rShift = (stubCluster->z() - tpa_->z0())/tpa_->tanLambda() - stubCluster->r();

      if (getSettings()->kalmanHOhelixExp()) rShift -= deltaS;

      if (getSettings()->kalmanHOprojZcorr() == 1) {
	// Add correlation term related to conversion of stub residuals from (r,phi) to (z,phi).
	correction[0] += inv2R * rShift; 
      }

      if (getSettings()->kalmanHOalpha()     == 1) {
	// Add alpha correction for non-radial 2S endcap strips..
	correction[0] += stubCluster->alpha() * rShift;
      }

      //cout<<"ENDCAP 2S STUB: (r,z)=("<<stubCluster->r()<<","<<stubCluster->z()<<") r*delta="<<stubCluster->r() * correction[0]<<" r*alphaCorr="<<stubCluster->r() * stubCluster->alpha() * rShift<<" rShift="<<rShift<<endl;
    }

    // Apply correction to residuals.
    delta[0] += correction[0];
    delta[1] += correction[1];
  }

  delta.at(0) = wrapRadian(delta.at(0));

  return delta;
}


void L1KalmanComb::bookHists(){

  if ( getSettings()->kalmanFillInternalHists() ) {

    edm::Service<TFileService> fs_;
    string dirName;
    if( fitterName_.compare("") == 0 ) dirName = "L1KalmanCombInternal";
    else dirName = fitterName_ + "Internal";

    TFileDirectory inputDir = fs_->mkdir(dirName.c_str());


    TString hname;
    hTrackEta_ = inputDir.make<TH1F>( "hTrackEta", "Track #eta; #eta", 50, -2.5, 2.5 );
    hUniqueTrackEta_ = inputDir.make<TH1F>( "hUniqueTrackEta", "Unique Track #eta; #eta", 50, -2.5, 2.5 );
    hndupStub_ = inputDir.make<TH1F>( "hndupStub", "# of duplicated stubs", 1, 0, 1 );
    hnmergeStub_ = inputDir.make<TH1F>( "hnmergeStub", "# of merged stubs", 1, 0, 1 );


    for( unsigned j_layer=0; j_layer < 6; j_layer++ ){
      hname = Form( "hBarrelStubMaxDistanceLayer%02d", LayerId[j_layer] );
      hBarrelStubMaxDistanceMap[hname] = inputDir.make<TH2F>( hname, Form( "max distance of stubs in barrel Layer %02d; dr#phi; dz", LayerId[j_layer] ), 
							      100, 0, 1., 100, 0, 10 );
    }

    for( unsigned j_ecring=1; j_ecring < 16; j_ecring++ ){
      hname = Form( "hEndcapStubMaxDistanceRing%02d", j_ecring );
      hEndcapStubMaxDistanceMap[hname] = inputDir.make<TH2F>( hname, Form( "max distance of stubs in endcap Ring %02d; dr#phi; dr", j_ecring ), 
							      100, 0, 1., 100, 0, 10 );
      hname = Form( "hphiErrorRatioRing%d", j_ecring );
      hphiErrorRatioMap[hname] = inputDir.make<TH2F>( hname, Form( "; fabs( strip id - 0.5 x nStrips + 0.5 ); #delta #phi_{r} / #delta #phi_{l}" ), 508, 0.0, 508.0 , 50, -0.5, 49.5 );
    }



    float nbins(2002);
    for( unsigned i=0; i < nPar_; i++ ){
      hname = Form( "hyt_%d", i );
      hytMap[hname] = inputDir.make<TH1F>( hname, Form( "; true track parameter values %d", i ), nbins, hymin[i], hymax[i] );
      hname = Form( "hy0_%d", i );
      hy0Map[hname] = inputDir.make<TH1F>( hname, Form( "; after HT track parameter values %d", i ), nbins, hymin[i], hymax[i] );
      hname = Form( "hyf_%d", i );
      hyfMap[hname] = inputDir.make<TH1F>( hname, Form( "; after KF track parameter values %d", i ), nbins, hymin[i], hymax[i] );
      hname = Form( "hx_%d", i );
      hxMap[hname] = inputDir.make<TH1F>( hname, Form( "; x values %d", i ), nbins, hxmin[i], hxmax[i] );
    }

 
    for( unsigned itr=0; itr<=5; itr++ ){
		
      hname = Form( "hstubComb_itr%d", itr );
      hstubCombMap[hname] = inputDir.make<TH1F>( hname, Form( "; #state+stub combinations, iteration %d ", itr ), 100, 0., 100.);

      for( unsigned i=0; i < nPar_; i++ ){
	for( unsigned j=0; j <= i; j++ ){
	
	  hname = Form( "hxcov_itr%d_%d_%d", itr, i, j );
	  hxcovMap[hname] = inputDir.make<TH1F>( hname, Form( "; state covariance adjusted values, iteration %d (%d,%d)", itr, i, j ), 
						 nbins, -1 * hdxmin[i]*hdxmin[j], hdxmax[i]*hdxmax[j] );
	}
      }
      for( unsigned i=0; i < nPar_; i++ ){
	for( unsigned j=0; j < nMeas_; j++ ){
	  hname = Form( "hk_itr%d_%d_%d", itr, i, j );
	  hkMap[hname] = inputDir.make<TH1F>( hname, Form( "; K(%d,%d), Iteration %d", i, j, itr ), 200, -1., 1. );
	}
      }
      for( unsigned i=0; i < nMeas_; i++ ){
	hname = Form( "hres_itr%d_%d", itr, i );
	hresMap[hname] = inputDir.make<TH1F>( hname, Form( "; residual values, iteration %d (%d)", itr, i ), 
					      nbins, hresmin[i], hresmax[i] );
	for( unsigned j=0; j <= i; j++ ){
	  hname = Form( "hmcov_itr%d_%d_%d", itr, i, j );
	  hmcovMap[hname] = inputDir.make<TH1F>( hname, Form( "; measurement covariance values, iteration %d (%d,%d)", itr, i, j ), 
						 nbins, -1 * hddMeasmin[i]*hddMeasmin[i], hddMeasmax[i]*hddMeasmax[j] );
	}
      }
    }
  }  
}


void L1KalmanComb::fillCandHists( const kalmanState &state, const TP *tpa )
{
  if( tpa && tpa->useForAlgEff() ){

    const kalmanState *the_state = &state;
    while( the_state ){
      if( the_state->stubCluster() ){
	std::vector<double> x = the_state->xa();
	for( unsigned i=0; i < nPar_; i++ ){
	  TString hname = Form( "hx_%d", i );
	  if( hxMap.find(hname) == hxMap.end() ){
	    cout << hname << " does not exist." << endl;
	  }
	  else hxMap[hname]->Fill(x.at(i));
	}
      }
      the_state = the_state->last_state();
    }


    std::map<std::string, double> mx = getTrackParams( &state );
    std::vector<double> vx(nPar_);
    vx[0] = mx["qOverPt"];
    vx[1] = mx["phi0"];
    vx[2] = mx["z0"];
    vx[3] = mx["t"];
    if( nPar_ == 5 ) vx[4] = mx["d0"];
    for( unsigned i=0; i < nPar_; i++ ){
      TString hname = Form( "hyf_%d", i );
      if( hyfMap.find(hname) == hyfMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hyfMap[hname]->Fill(vx[i]);
    }
  }
}


void L1KalmanComb::fillSeedHists( const kalmanState *state, const TP *tpa ){

  std::vector<double> x0   = state->xa();
  TMatrixD            pxx0 = state->pxxa();
  //Histogram Fill : seed pxxa 
  for( unsigned i=0; i < nPar_; i++ ){
    for( unsigned j=0; j <= i; j++ ){
      TString hname = Form( "hxcov_itr%d_%d_%d", 0, i, j );
      if( hxcovMap.find( hname ) == hxcovMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hxcovMap[hname]->Fill( pxx0(i,j) );
    }
  }

  if( tpa && tpa->useForAlgEff() ){
    std::vector<double> tpParams(nPar_);
    tpParams[0] = tpa->qOverPt();
    tpParams[1] = tpa->phi0();
    tpParams[2] = tpa->z0();
    tpParams[3] = tpa->tanLambda();
    if( nPar_ == 5 ) tpParams[4] = tpa->d0();
    for( unsigned i=0; i < nPar_; i++ ){
      TString hname = Form( "hyt_%d", i );
      if( hytMap.find(hname) == hytMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hytMap[hname]->Fill(tpParams[i]);
    }
    //Histogram Fill : Seed state 
    std::map<std::string, double> trackParams = getTrackParams( state );
    std::vector<double> trackParVec(nPar_);
    trackParVec[0] = trackParams["qOverPt"];
    trackParVec[1] = trackParams["phi0"];
    trackParVec[2] = trackParams["z0"];
    trackParVec[3] = trackParams["t"];
    if( nPar_ == 5 ) trackParVec[4] = trackParams["d0"];
    for( unsigned i=0; i < nPar_; i++ ){
      TString hname = Form( "hy0_%d", i );
      if( hy0Map.find(hname) == hy0Map.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hy0Map[hname]->Fill(trackParVec[i]);
    }
  }
}


void L1KalmanComb::fillStepHists( const TP *tpa, unsigned nItr, const kalmanState *new_state )
{
  unsigned path = 0;

  const std::vector<double> &xa = new_state->xa();
  const StubCluster *stubCluster = new_state->stubCluster();
  const TMatrixD &pxxa = new_state->pxxa();

  TString hname;

  for( unsigned i=0; i < nPar_; i++ ){

    for( unsigned j=0; j <= i; j++ ){
      hname = Form( "hxcov_itr%d_%d_%d", nItr, i, j );
      if( hxcovMap.find( hname ) == hxcovMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hxcovMap[hname]->Fill( pxxa(i,j) );
    }
  }
  for( unsigned i=0; i < nPar_; i++ ){
    for( int j=0; j < 2; j++ ){
      TString hname = Form( "hk_itr%d_%d_%d", nItr, i, j );
      if( hkMap.find( hname ) == hkMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hkMap[hname]->Fill( new_state->K()(i,j) );
    }
  }
  std::vector<double> delta_new = residual(stubCluster, xa, new_state->candidate().qOverPt() );
  for( unsigned int i=0; i < delta_new.size(); i++ ){
    TString hname = Form( "hres_itr%d_%d", nItr, i );
    if( hresMap.find(hname) == hresMap.end() ){
      cout << hname << " does not exist." << endl;
    }
    else hresMap[hname]->Fill( delta_new[i] );  
  }
  for( int i=0; i < 2; i++ ){
    for( int j=0; j < i; j++ ){
      TString hname = Form( "hmcov_itr%d_%d_%d", nItr, i, j );
      if( hmcovMap.find( hname ) == hmcovMap.end() ){
	cout << hname << " does not exist." << endl;
      }
      else hmcovMap[hname]->Fill( new_state->dcov()(i,j) );
    }
  }
}


void L1KalmanComb::deleteStubClusters()
{
  for( unsigned int i=0; i < stbcl_list_.size(); i++ ){
    delete stbcl_list_.at(i);
  }
  stbcl_list_.clear();
}


double L1KalmanComb::DeltaRphiForClustering( unsigned layerId, unsigned endcapRing )
{
  static double barrel_drphi[6] = { 0.05, 0.04, 0.05, 0.12, 0.13, 0.19 }; 
  if( layerId < 10 ) return barrel_drphi[layerId - 1];

  static double ec_drphi[16] =  
    { 0.04, 0.05, 0.04, 0.06, 0.06, 0.04, 0.06, 0.07, 0.15, 0.08, 0.27, 0.08, 0.27, 0.12, 0.09 };
  return ec_drphi[endcapRing - 1];
};


double L1KalmanComb::DeltaRForClustering( unsigned endcapRing )
{
  static double ec_dr[16] =  
    { 0.52, 0.56, 0.59, 0.86, 0.66, 0.47, 0.55, 0.72, 1.53, 1.10, 2.72, 0.91, 2.69, 0.67, 0.09 };
  return ec_dr[endcapRing - 1];

}


bool L1KalmanComb::isOverlap( const Stub* a, const Stub*b, OVERLAP_TYPE type ){

  std::set<const TP*> a_tps = a->assocTPs();
  std::set<const TP*> b_tps = b->assocTPs();
  double drphi = DeltaRphiForClustering( a->layerId(), a->endcapRing() );
  double dr(0);
  switch ( type ){

  case TYPE_NORMAL:
    if( a->layerId() != b->layerId() ) return false;

    if( a->layerId() < 7 ){
      if( fabs( b->z() - a->z() ) > 0.5 * b->stripLength() || fabs( wrapRadian( b->phi() - sectorPhi() ) * b->r() - wrapRadian( a->phi() - sectorPhi() ) * a->r() ) > 0.5 * b->stripPitch() ) return false;
    }
    else{
      dr = DeltaRForClustering( a->endcapRing() ); 
      if( fabs( b->r() - a->r() ) > 0.5 * b->stripLength() || fabs( wrapRadian( b->phi() - sectorPhi() ) * b->r() - wrapRadian( a->phi() - sectorPhi() ) * a->r() ) > 0.5 * b->stripPitch() ) return false;
    }
    return true;
  case TYPE_V2:
    if( a->layerId() != b->layerId() ) return false;

    if( a->layerId() < 7 ){
      if( fabs( b->z() - a->z() ) > 0.5 * b->stripLength() || fabs( wrapRadian( b->phi() - sectorPhi() ) * b->r() - wrapRadian( a->phi() - sectorPhi() ) * a->r() ) > drphi ) return false;
    }
    else{
      dr = DeltaRForClustering( a->endcapRing() ); 
      if( fabs( b->r() - a->r() ) > dr || fabs( wrapRadian( b->phi() - sectorPhi() ) * b->r() - wrapRadian( a->phi() - sectorPhi() ) * a->r() ) > drphi ) return false;
    }
    return true;

  case TYPE_NOCLUSTERING:
    return false;

  case TYPE_TP:
    for( auto a_tp : a_tps ) 
      if( b_tps.find( a_tp ) != b_tps.end() ) return true;
    return false;
  default:
    return false;
  }
}

set<unsigned> L1KalmanComb::getKalmanDeadLayers( unsigned layerMap[18][8], bool& remove2PSCut ) const {

  // By which Stress Test scenario (if any) are dead modules being emulated?
  const unsigned int killScenario = getSettings()->killScenario(); 
  // Should TMTT tracking be modified to reduce efficiency loss due to dead modules?
  const bool killRecover = getSettings()->killRecover();

  set<unsigned> deadLayers;

  if (killRecover) {
    if ( killScenario == 1 ) {
      deadLayers = {4};
      if ( iCurrentEtaReg_ < 5 || iCurrentEtaReg_ > 8 || iCurrentPhiSec_ < 8 || iCurrentPhiSec_ > 11 ) {
	deadLayers.clear();
      }

    }
    else if ( killScenario == 2 ) {
      deadLayers = {1};
      if ( iCurrentEtaReg_ > 8 || iCurrentPhiSec_ < 8 || iCurrentPhiSec_ > 11 ) {
	deadLayers.clear();
      }
      remove2PSCut = true;
    }
    else if ( killScenario == 3 ) {
      deadLayers = {1,2};
      if ( iCurrentEtaReg_ > 8 || iCurrentPhiSec_ < 8 || iCurrentPhiSec_ > 11 ) {
	deadLayers.clear();
      }
      else if ( iCurrentEtaReg_ < 1 ) {
	deadLayers = {0};
      }
      remove2PSCut = true;
    }
    else if ( killScenario == 4 ) {
      deadLayers = {1,2};
      if ( iCurrentEtaReg_ > 8 || iCurrentPhiSec_ < 8 || iCurrentPhiSec_ > 11 ) {
	deadLayers.clear();
      }
      else if ( iCurrentEtaReg_ > 3 ) {
	deadLayers = {0};
      }
      remove2PSCut = true;
    }
  }

  set<unsigned> kalmanDeadLayers;
  for ( auto layer : deadLayers ) {
    kalmanDeadLayers.insert( layerMap[iCurrentEtaReg_][layer] );
  }

  return kalmanDeadLayers;
}

//=== Function to calculate approximation for tilted barrel modules (aka B) copied from Stub class.

float L1KalmanComb::getApproxB(float z, float r) const {
  return getSettings()->bApprox_gradient() * fabs(z)/r + getSettings()->bApprox_intercept();
}

}

