///=== This is the Kalman Combinatorial Filter for 4 running parameters track fit algorithm.
 
 
#include "L1Trigger/TrackFindingTMTT/interface/KFRunningComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/kalmanState.h"
#include "L1Trigger/TrackFindingTMTT/interface/StubCluster.h"
#define CKF_DEBUG

namespace TMTT {

/*
// Scattering constants - HISTORIC NOT USED.

static unsigned nlayer_eta[25] = 
{ 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
7, 7, 7, 7, 6, 6, 6, 6, 6, 6};

static double matx_outer[25] = {
0.16, 0.17, 0.18, 0.19, 0.20, 
0.21, 0.26, 0.22, 0.26, 0.38,
0.41, 0.40, 0.44, 0.50, 0.54,
0.60, 0.44, 0.48, 0.60, 0.68,
0.50, 0.48, 0.64, 0.39, 0.20
};

static double matx_inner[25] = {
0.14, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 
0.12, 0.1, 0.1, 0.1, 0.15,
0.20, 0.25, 0.25, 0.3, 0.3,
0.35, 0.40, 0.40, 0.6, 0.6
};
*/

static double wrapRadian( double t ){

  if( t > 0 ){
    while( t > M_PI ) t-= 2*M_PI; 
  }
  else{
    while( t < - M_PI ) t+= 2*M_PI; 
  }
  return t;
}

KFRunningComb::KFRunningComb(const Settings* settings, const uint nPar, const string &fitterName ) : L1KalmanComb(settings, nPar, fitterName ){

  hdxmin[RPHI] = -1.;
  hdxmax[RPHI] = +1.;
  hdxmin[P] = -300.;
  hdxmax[P] = +300.;
  hdxmin[INV2R] = -0.3;
  hdxmax[INV2R] = +0.3;
  hdxmin[T] = -6.;
  hdxmax[T] = +6.;

  hxmin[RPHI] = -1;
  hxmax[RPHI] = +1;
  hxmin[P] = -300;
  hxmax[P] = +300;
  hxmin[INV2R] = -0.3;
  hxmax[INV2R] = +0.3;
  hxmin[T] = -6.;
  hxmax[T] = +6.;

  hddMeasmin[RPHI] = -1.e1;
  hddMeasmax[RPHI] = +1.e1;

  hresmin[RPHI] = -10.;
  hresmax[RPHI] = +10.;

  hxaxtmin[RPHI] = -1.e-1;
  hxaxtmax[RPHI] = +1.e-1;
  hxaxtmin[P] = -10.;
  hxaxtmax[P] = +10.;
  hxaxtmin[INV2R] = -1.e-3;
  hxaxtmax[INV2R] = +1.e-3;
  hxaxtmin[T] = -1.e-0;
  hxaxtmax[T] = +1.e-0;

}


std::map<std::string, double> KFRunningComb::getTrackParams(const kalmanState *state )const{

  std::vector<double> x = state->xa();
  std::map<std::string, double> y;
  if ( state->barrel() ) {
    y["qOverPt"] = 2. * x.at(INV2R) / getSettings()->invPtToInvR(); 
    y["phi0"] = x.at(RPHI)/state->r() - x.at(INV2R)*state->r();
    y["z0"] = x.at(P) - x.at(T)*state->r();
    y["t"] = x.at(T);
  }
  else { 
    y["qOverPt"] = 2. * x.at(INV2R) * x.at(T) / getSettings()->invPtToInvR();
    y["phi0"] = x.at(RPHI)/state->r() - x.at(T)*x.at(INV2R) * state->r();
    y["z0"] = state->z() - 1.0/x.at(T) * state->r();
    y["t"] = 1.0/x.at(T);
  }
  return y;
}

/* If using 5 param helix fit, get track params with beam-spot constraint & track fit chi2 from applying it. */
std::map<std::string, double> KFRunningComb::getTrackParams_BeamConstr( const kalmanState *state, double& chi2 ) const {
  return (this->getTrackParams(state));
}

 
TMatrixD KFRunningComb::H(const StubCluster* stubCluster)const{
  TMatrixD h(2, nPar_);
  h(RPHI_S, RPHI_S) = 1.0;
  h(P_S, P_S) = 1.0;

  return h;
}
 
/* Seed the state vector */
std::vector<double> KFRunningComb::seedx(const L1track3D& l1track3D)const{

  // Determine if innermost stub/cluster is barrel or not
  unsigned int innerLayer {99};
  bool isBarrel;
  if ( l1track3D.getStubClusters().size() != 0 ) {
    for ( auto cls : l1track3D.getStubClusters() ) {
      unsigned int clsLayer = Utility::layerMap(l1track3D.iEtaReg(),cls->layerIdReduced());
      if ( clsLayer < innerLayer ) isBarrel = cls->barrel();
    }
  }
  else {
    for ( auto s : l1track3D.getStubs() ) {
      unsigned int stubLayer = Utility::layerMap(l1track3D.iEtaReg(),s->layerIdReduced());
      if ( stubLayer < innerLayer ) isBarrel = s->barrel();
    }
  }

  std::vector<double> x(nPar_);
  if ( isBarrel ) {
    x[RPHI]  = 0.0;
    x[P]     = 0.0;
    x[INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
    x[T]     = l1track3D.tanLambda();
  }
  else {
    x[RPHI]  = 0.0;
    x[P]     = 0.0;
    x[INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/(2*l1track3D.tanLambda());
    x[T]     = 1.0/l1track3D.tanLambda();
  }
  return x;
}

/* Seed the covariance matrix */
TMatrixD KFRunningComb::seedP(const L1track3D& l1track3D, const bool seedPair)const{
  TMatrixD p(nPar_,nPar_);

  const double c = getSettings()->invPtToInvR() / 2; 

  // Assumed track seed (from HT) uncertainty in transverse impact parameter.
  const float d0Sigma = 1.0;

  if (getSettings()->hybrid()) {

    p(RPHI,RPHI) = 0.0051 * 0.0051 * 4; 
    p(P,P) = 5.0 * 5.0; // N.B. r-z seed uncertainties could be smaller for hybrid, except if seeded in 2S?
    p(INV2R,INV2R) = 0.0157 * 0.0157 * c * c * 4; 
    p(T,T) = 0.25 * 0.25 * 4;

  } else if ( getSettings()->runFullKalman() ) {

    p(RPHI,RPHI) = 0.0051 * 0.0051 * 4 * 4;
    p(P,P) = 5.0 * 5.0;
    p(INV2R,INV2R) = 0.0157 * 0.0157 * c * c * 4; // 250;
    p(T,T) = 0.25 * 0.25 * 4; // IRT: increased by factor 4, as was affecting fit chi2.

  } else {
    // optimised for 18x2 with additional error factor in pt/phi to avoid pulling towards wrong HT params
    p(RPHI,RPHI) = 0.0051 * 0.0051 * 4; // Based on HT cell size.
    p(P,P) = 5.0 * 5.0; 
    p(INV2R,INV2R) = 0.0157 * 0.0157 * c * c * 4;  // Base on HT cell size
    p(T,T) = 0.25 * 0.25 * 4; // IRT: increased by factor 4, as was affecting fit chi2.

    if ( getSettings()->numEtaRegions() <= 12 ) {    
      // Inflate eta errors
      p(T,T) = p(T,T) * 2 * 2;
    }
  }

  return p;
}

/* The forecast matrix
 * (does not equal identity matrix like the helix params version) */
TMatrixD KFRunningComb::F(const StubCluster* stubCluster, const kalmanState *state )const{
  TMatrixD F(nPar_,nPar_); 

  F(P, P) = 1.0;
  F(INV2R, INV2R) = 1.0;
  F(T, T) = 1.0;

  const bool isBarrel {stubCluster->barrel()};
  if ( isBarrel ) {
    F(RPHI, RPHI)   = stubCluster->r() / state->r();
    F(RPHI, INV2R) = -1.0 * stubCluster->r() / ( stubCluster->r() / state->r() );
    F(P, T) = stubCluster->r() - state->r();
  }
  else {
    F(RPHI, RPHI)   = stubCluster->r() / state->r();
    F(RPHI, INV2R) = -1.0 * stubCluster->r() / ( stubCluster->z() / state->z() );
    F(P, T) = stubCluster->z() - state->z();
  }
  return F;
}

/* the vector of measurements */
std::vector<double> KFRunningComb::d(const StubCluster* stubCluster )const{
  const bool isBarrel {stubCluster->barrel()};

  std::vector<double> meas;
  meas.resize(2);
  meas[RPHI_S] = stubCluster->r() * wrapRadian( stubCluster->phi() - sectorPhi() );
  if ( isBarrel ) meas[P_S] = stubCluster->z();
  else meas[P_S] = stubCluster->r();
  return meas;
}

// Assumed hit resolution in (phi,z)
TMatrixD KFRunningComb::PddMeas(const StubCluster* stubCluster, const kalmanState *state )const{

  double inv2R = (getSettings()->invPtToInvR()) * 0.5 * state->candidate().qOverPt();
  double inv2R2 = inv2R * inv2R;

  double tanl = state->xa().at(T);  // factor of 0.9 improves rejection
  double tanl2 = tanl * tanl; 

  TMatrixD p(2,2);

  double vphi(0);
  double vz(0);
  double vcorr(0);

  vphi = stubCluster->sigmaX();
  vz = stubCluster->sigmaZ();

  if ( !stubCluster->barrel() ) vcorr = stubCluster->dphi_dr();

  p(RPHI, RPHI) = vphi*vphi + vcorr*vcorr;
  p(P, P) = vz*vz;

  return p;

}

// State uncertainty due to scattering -- HISTORIC NOT USED
TMatrixD KFRunningComb::PxxModel( const kalmanState *state, const StubCluster *stubCluster )const
{
  TMatrixD p(nPar_,nPar_);
  return p;
}

bool KFRunningComb::isGoodState( const kalmanState &state, const bool seedPair )const
{
  const unsigned int option = getSettings()->kalmanSeedingOption();

  // Cut values. (Layer 0 entry here is dummy). -- todo : make configurable

  vector<float> z0Cut, ptTolerance, d0Cut, chi2Cut;
  //  Layer   =    0      1      2     3     4      5      6
  z0Cut       = { 999.,  999.,   15.,  15.,  15.,   15.,   15.};
  ptTolerance = { 999., 999., 0.3, 0.1, 0.05, 0.05, 0.05}; /// optimum for Full CKF 5 runnning
  d0Cut       = { 999.,  999.,  999.,  10.,   5.,    5.,   5.};  // Only used for 5 param helix fit
  chi2Cut     = { 999.,  999.,   999.,  999.,  999.,  120.,  160.};  // Consider reducing chi2 cut 2 to 7.

  unsigned nStubLayers = state.nStubLayers();
  bool goodState( true );

  std::map<std::string, double> y = getTrackParams( &state );
  double qOverPt = y["qOverPt"]; 
  double pt=fabs( 1/qOverPt ); 
  double z0=fabs( y["z0"] ); 

  // state parameter selections

  if (z0 > z0Cut[nStubLayers] ) goodState = false;
  if( pt < getSettings()->houghMinPt() - ptTolerance[nStubLayers] ) goodState = false;

  // chi2 selection

  if (getSettings()->kalmanMultiScattTerm() > 0.0001) {   // Scattering taken into account

    if (state.chi2() > chi2Cut[nStubLayers]) goodState=false; // No separate pT selection needed

  } else {  // scattering ignored - HISTORIC

    // N.B. Code below changed by Alexander Morton to allow tracking down to Pt = 2 GeV.
    if( nStubLayers == 2 ) {
      if (state.chi2() > 15.0) goodState=false; // No separate pT selection needed
    } else if ( nStubLayers == 3 ) {
      if (state.chi2() > 100.0 && pt > 2.7) goodState=false;
      if (state.chi2() > 120.0 && pt <= 2.7) goodState=false;
    } else if ( nStubLayers == 4 ) { 
      if (state.chi2() > 320.0 && pt > 2.7) goodState=false;
      if (state.chi2() > 1420.0 && pt <= 2.7) goodState=false;
    } else if ( nStubLayers == 5 ) {  // NEEDS TUNING FOR 5 OR 6 LAYERS !!!
      if (state.chi2() > 480.0 && pt > 2.7) goodState=false;
      if (state.chi2() > 2130.0 && pt <= 2.7) goodState=false;
    } else if ( nStubLayers >= 6 ) {  // NEEDS TUNING FOR 5 OR 6 LAYERS !!!
      if (state.chi2() > 640.0 && pt > 2.7) goodState=false;
      if (state.chi2() > 2840.0 && pt <= 2.7) goodState=false;
    }

  }

  const bool countUpdateCalls = false; // Print statement to count calls to Updator.

  if ( countUpdateCalls ||
       (getSettings()->kalmanDebugLevel() >= 2 && tpa_ != nullptr && tpa_->useForAlgEff() ) ||
       (getSettings()->kalmanDebugLevel() >= 2 && getSettings()->hybrid()) ) {
    if (not goodState) cout<<"State veto:"; 
    if (goodState)     cout<<"State kept:";
    cout<<" nlay="<<nStubLayers<<" nskip="<<state.nSkippedLayers()<<" chi2="<<state.chi2();
    cout<<" pt="<<pt<<" q/pt="<<qOverPt<<" tanL="<<y["t"]<<" z0="<<y["z0"]<<" phi0="<<y["phi0"];
    if (nPar_ == 5) cout<<" d0="<<y["D0"];
    cout<<" fake"<<(tpa_ == nullptr);
    if (tpa_ != nullptr) cout<<" pt(mc)="<<tpa_->pt();
    cout<<endl;
  }

  return goodState;
}

}

