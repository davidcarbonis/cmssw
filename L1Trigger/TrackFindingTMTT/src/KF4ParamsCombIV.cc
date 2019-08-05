///=== This is the Kalman Combinatorial Filter for 4 helix parameters track fit algorithm.
 
#include "L1Trigger/TrackFindingTMTT/interface/KF4ParamsCombIV.h"
#include "L1Trigger/TrackFindingTMTT/interface/kalmanState.h"
#include "L1Trigger/TrackFindingTMTT/interface/StubCluster.h"
#define CKF_DEBUG

namespace TMTT {
 
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

static double wrapRadian( double t ){

    if( t > 0 ){
	while( t > M_PI ) t-= 2*M_PI; 
    }
    else{
	while( t < - M_PI ) t+= 2*M_PI; 
    }
    return t;
}

KF4ParamsCombIV::KF4ParamsCombIV(const Settings* settings, const uint nPar, const string &fitterName ) : L1KalmanComb(settings, 4, fitterName ){

    hdxmin[0] = -0.0060;
    hdxmax[0] = +0.0060;
    hdxmin[1] = -4.1;
    hdxmax[1] = +4.1;
    hdxmin[2] = -1.1e-4;
    hdxmax[2] = +1.1e-4;
    hdxmin[3] = -1.;
    hdxmax[3] = +1.;

    hxmin[0] = -50;
    hxmax[0] = +50;
    hxmin[1] = -120;
    hxmax[1] = +120;
    hxmin[2] = -0.3 * 0.0057;
    hxmax[2] = +0.3 * 0.0057;
    hxmin[3] = -6.;
    hxmax[3] = +6.;

    hddMeasmin[1] = -1.e1;
    hddMeasmax[1] = +1.e1;

    hresmin[0] = -0.5;
    hresmax[0] = +0.5;
    hresmin[1] = -10.;
    hresmax[1] = +10.;


    hxaxtmin[0] = -1.e-3;
    hxaxtmax[0] = +1.e-3;
    hxaxtmin[1] = -1.e-1;
    hxaxtmax[1] = +1.e-1;
    hxaxtmin[2] = -10.;
    hxaxtmax[2] = +10.;
    hxaxtmin[3] = -6.e-0;
    hxaxtmax[3] = +6.e-0;
}

std::map<std::string, double> KF4ParamsCombIV::getTrackParams( kalmanState *state )const{

    std::map<std::string, double> y;
    if( state->barrel() ){
//	y["phi0"] = wrapRadian( state->xa().at(BP_RHOPHI) / state->r() + state->r() * state->xa().at(BP_INV2R) ) + sectorPhi() ;
	y["phi0"] = wrapRadian(state->xa().at(BP_RHOPHI)/state->r() + asin(state->r() * state->xa().at(BP_INV2R)));;
	y["z0"]   = state->xa().at(BP_Z)   - state->r() * state->xa().at(BP_T);
	y["qOverPt"] = state->xa().at(BP_INV2R) / getSettings()->invPtToInvR() * 2.; 
	y["t"] = state->xa().at(BP_T);

    }
    else{
	double t = 1./state->xa().at(EP_INVT);
	double inv2R = state->xa().at(EP_INV2RT) * t;
//	y["phi0"] = wrapRadian( state->xa().at(EP_RHOPHI) / state->xa().at(EP_RHO) + state->xa().at(EP_RHO) * inv2R + sectorPhi() );
	y["phi0"] = wrapRadian( state->xa().at(EP_RHOPHI) / state->xa().at(EP_RHO) + asin(state->xa().at(EP_RHO) * inv2R) );
	y["z0"]   = state->z() - state->xa().at(EP_RHO) * t;
	y["qOverPt"] = inv2R / getSettings()->invPtToInvR() * 2.; 
	y["t"] = t;
    }
    y["d0"] = 0;
    return y;
}
/* The Kalman measurement matrix
 * Here I always measure phi(r), and z(r) */
TMatrixD KF4ParamsCombIV::H(const StubCluster* stubCluster)const{
    TMatrixD h(2, 4);
    if( stubCluster->barrel() ){
	h(BM_RHOPHI,BP_RHOPHI) = 1;
	h(BM_Z,BP_Z)     = 1;
    }
    else{
	h(EM_RHOPHI,BP_RHOPHI) = 1;
	h(EM_RHO,EP_RHO) = 1;
    }
    return h;
}

/* Seed the state vector */
std::vector<double> KF4ParamsCombIV::seedx(const L1track3D& l1track3D)const{

    std::vector<double> x(nPar_);

    // placeholder as rphi and z need first stub considered info
    x[BP_RHOPHI] = 0.1 * wrapRadian( l1track3D.phi0() - sectorPhi() );
    x[BP_Z]   = l1track3D.z0();
    x[BP_INV2R] = getSettings()->invPtToInvR() * l1track3D.qOverPt()/2;
    x[BP_T]     = l1track3D.tanLambda();

    return x;
}

// Set the seed running params state vector
void KF4ParamsCombIV::setSeedX( kalmanState &state, const StubCluster *stubCluster )const{

  std::vector<double> x(nPar_);
  double r = stubCluster->r();
  double phi = stubCluster->phi();
  double rphi = r*phi;
  double inv2R = getSettings()->invPtToInvR() * state.candidate().qOverPt() / 2.;
  double t = state.candidate().tanLambda();

  if ( state.barrel() ) {
    x[BP_RHOPHI] = rphi;
    x[BP_Z]	 = stubCluster->z();
    x[BP_INV2R]  = inv2R;
    x[BP_T]	 = t;
  }
  else {
    x[EP_RHOPHI] = rphi;
    x[EP_RHO]    = r;
    if ( t ) {
      x[EP_INV2RT] = inv2R / t;
      x[EP_INVT]   = 1./t;
    }
    else {
      x[EP_INV2RT] = inv2R * t;
      x[EP_INVT]   = 1.*t;
    }
  }
  state.setSeed(x);
}

/* Seed the covariance matrix
 * Note: 1024 is an arbitrary 'large' value */
TMatrixD KF4ParamsCombIV::seedP(const L1track3D& l1track3D, const bool seedPair)const{
    TMatrixD p(4,4);
    double c = getSettings()->invPtToInvR() / 2; 
    //0.0057
    //
    p(BP_RHOPHI,BP_RHOPHI) = 0.1 * 0.006* 0.1 * 0.006; 
    p(BP_Z,BP_Z) = 4.0*4.0; 
//    p(BP_INV2R,BP_INV2R) = 0.020*0.020 * c * c; 
    p(BP_INV2R,BP_INV2R) = 0.0157 * 0.0157 * c * c * 4; 
    //0.020*c = 0.0001
    p(BP_T,BP_T) = 0.25 * 0.25 *4;

/*
    if( getSettings()->numEtaRegions() == 9 ) { 
	switch ( iCurrentEtaReg_ ){
	    case 0 :
	    case 8 :
		p(BP_T,BP_T) = 0.5*0.5; 
		break;
	    case 1 :
	    case 7 :
		p(BP_T,BP_T) = 0.4*0.4; 
		break;
	    case 2 :
	    case 6 :
		p(BP_T,BP_T) = 0.37*0.37; 
		break;
	    case 3 :
	    case 5 :
		p(BP_T,BP_T) = 0.27*0.27; 
		break;
	    case 4 :
		p(BP_T,BP_T) = 0.27*0.27; 
		break;
	}
    }
    else{
	p(BP_T,BP_T) = 0.5*0.5; 
    }
*/
    return p;
}
double KF4ParamsCombIV::getZ( kalmanState *state )const{

    double z(0);
    if( state->barrel() ){
	z = state->xa().at(BP_Z);
    }
    else{
	z = state->z();
    }
    return z;
}

double KF4ParamsCombIV::getZ0( kalmanState *state )const{

    double z0(0); 
    if( state->barrel() ){
	if( state->nStubLayers() == 0 ) z0 = state->xa().at(BP_Z);
	else
	    z0 = state->xa().at(BP_Z) - state->r() * state->xa().at(BP_T);
    }
    else{
	z0 = state->z() - state->xa().at(EP_RHO) / state->xa().at(EP_INVT);
    }
    return z0;
}
double KF4ParamsCombIV::getZVariance( kalmanState *state )const{

    double vz(0); 
    if( state->barrel() ){
	vz = state->pxxa()(BP_Z, BP_Z);
    }

    return vz;
}

double KF4ParamsCombIV::getZ0Variance( kalmanState *state )const{

    double vz0(0); 
    if( state->barrel() ){
	if( state->nStubLayers() == 0 ) state->pxxa()(BP_Z, BP_Z); 
	else
	    vz0 = state->pxxa()(BP_Z, BP_Z) + state->r() * state->r() * state->pxxa()(BP_T, BP_T);
    }
    else{
	double r = state->xa().at(EP_RHO);
	double vr = state->pxxa()(EP_RHO,EP_RHO);
	double invt = state->xa().at(EP_INVT);
	double vinvt = state->pxxa()(EP_INVT, EP_INVT);
	vz0 = r/invt * r/invt * ( vr/(r*r) + vinvt/(invt*invt) ); 
    }
    return vz0;
}
void KF4ParamsCombIV::barrelToEndcap( double r, const StubCluster *stubCluster, std::vector<double> &x )const{

    double rhophi = x[BP_RHOPHI];
    double Inv2R = x[BP_INV2R];
    double t  = x[BP_T];

    x[EP_RHOPHI]   = rhophi;
    x[EP_RHO]      = r;
    if( t ){
	x[EP_INV2RT]   = Inv2R / t; 
	x[EP_INVT]     = 1./t;
    }
    else{
	//	if( iCurrentEtaReg_ > 4  ) t = 1.;
	if( iCurrentEtaReg_ > getSettings()->numEtaRegions() / 2 ) t = 1.;
	else t = -1.;
	x[EP_INV2RT]   = Inv2R * t; 
	x[EP_INVT]     = 1. * t;
    }
}


/* The forecast matrix */
TMatrixD KF4ParamsCombIV::F(const StubCluster* stubCluster, kalmanState *state )const{

    TMatrixD F(4,4);
    for(int n = 0; n < 4; n++)
	F(n, n) = 1;

/*
    std::cout << __LINE__ << " : " << __FILE__ << std::endl;
    std::cout << "stubCluster->r(): " << stubCluster->r() << std::endl;
    std::cout << "stubCluster->phi(): " << stubCluster->phi() << std::endl;
    std::cout << "stubCluster layerIdReduced: " << stubCluster->layerIdReduced() << std::endl;;  

    std::cout << __LINE__ << " : " << __FILE__ << std::endl;
    std::cout << "state->nextLayer(): " << state->nextLayer() << std::endl;
    std::cout << "state->nStubLayers(): " << state->nStubLayers() << std::endl;
    std::cout << "state->r(): " << state->r() << std::endl;
    std::cout << "state->layerId: " << state->layerId() << std::endl;
*/

    if ( state->nextLayer() == 0 ) return F;

    if( stubCluster->barrel() ){ 
	double deltar = stubCluster->r() - state->r();
	F(BP_RHOPHI,BP_RHOPHI) = stubCluster->r()/state->r();
	F(BP_RHOPHI,BP_INV2R) = - deltar * stubCluster->r();
	F(BP_Z,BP_T)       =   deltar;
    }
    else{
	double statez = getZ(state);
	double deltaz = stubCluster->z() - statez;
	double z0 = getZ0( state );

	//	F(EP_RHOPHI,EP_RHOPHI) = 1 + z0 * ( 1./statez - 1./stub->z() );
	if( state->barrel() ){
	    double t = state->xa().at(BP_T);
	    if( state->nStubLayers() == 0 ){
		double deltar = stubCluster->r() - state->r();
		F(EP_RHOPHI,EP_RHOPHI) = stubCluster->r() / state->r();
		F(EP_RHOPHI,EP_INV2RT) = - stubCluster->r()*stubCluster->z();
		F(EP_RHO,EP_INVT)   =   stubCluster->z();
	    }
	    else{
		//	F(EP_RHOPHI,EP_RHOPHI) = 1 + z0 * ( 1./statez - 1./stub->z() );
		F(EP_RHOPHI,EP_RHOPHI) = stubCluster->r() / state->r(); 
		F(EP_RHOPHI,EP_INV2RT) = - deltaz * stubCluster->r();
		F(EP_RHO,EP_INVT)   =   deltaz;
	    }
	}
	else{
	    //F(EP_RHOPHI,EP_RHOPHI) = 1 + z0 * ( 1./statez - 1./stub->z() );
	    F(EP_RHOPHI,EP_RHOPHI) = stubCluster->r()/ state->r(); 
	    F(EP_RHOPHI,EP_INV2RT) = - deltaz * stubCluster->r(); 
	    F(EP_RHO,EP_INVT)   =   deltaz;
	}
    }

    return F;
}

/* the vector of measurements */
std::vector<double> KF4ParamsCombIV::d(const StubCluster* stubCluster )const{
    std::vector<double> meas;
    meas.resize(2);
    if( stubCluster->barrel() ){ 
//	meas[0] = stubCluster->r() * wrapRadian( stubCluster->phi() - sectorPhi() );
	meas[0] = stubCluster->r() * stubCluster->phi();
	meas[1] = stubCluster->z();
    }
    else{
//	meas[0] = stubCluster->r() * wrapRadian( stubCluster->phi() - sectorPhi() );
	meas[0] = stubCluster->r() * stubCluster->phi();
	meas[1] = stubCluster->r();
    }
    return meas;
}

TMatrixD KF4ParamsCombIV::PddMeas(const StubCluster* stubCluster, kalmanState *state )const{


    double rdphi = stubCluster->r() * stubCluster->dphi();
    TMatrixD p(2,2);
    p(0,0) = rdphi * rdphi; 
    p(1,1) = stubCluster->sigmaZ() * stubCluster->sigmaZ(); 

    return p;
}

/* State uncertainty */
TMatrixD KF4ParamsCombIV::PxxModel( kalmanState *state, const StubCluster *stubCluster )const
{

    TMatrixD p(4,4);
/*
    if( getSettings()->kalmanMultiScattFactor() ){
	unsigned i_eta = abs( stubCluster->eta() / 0.1 );
	if( i_eta > 24 ) i_eta = 24;
	double dl = matx_outer[i_eta] / nlayer_eta[i_eta];

	kalmanState * last_update_state = state->last_update_state();
	unsigned last_itr(1);
	if( last_update_state ) last_itr = last_update_state->nIterations();
	dl = ( stub_itr - last_itr ) * dl; 

	if( dl ){
	    std::map<std::string, double> y = getTrackParams( state );
	    double dtheta0 = 1./sqrt(3) * 0.0136 * fabs(y["qOverPt"]) * sqrt(dl)*( 1+0.038*log(dl) ); 
	    dtheta0 *= getSettings()->kalmanMultiScattFactor();
	    double drphi = stubCluster->r() * dtheta0;
	    p(0,0) = drphi * drphi; 
	}
    }
*/
    if( stubCluster->barrel() || state->nStubLayers() == 0 ) return p;
    //    if( stubCluster->barrel() ) return p;

    //error from z_{k-1}
    if( state->barrel() ){

	double vz = state->pxxa()( BP_Z, BP_Z );

	std::vector<double> x = state->xa();
	barrelToEndcap(state->r(), stubCluster, x);

	double dfx = stubCluster->r() * x.at(EP_INV2RT); 
	p(EP_RHOPHI, EP_RHOPHI) += vz * dfx * dfx;

	dfx = -1. * x.at(EP_INVT);
	p(EP_RHO, EP_RHO) += vz * dfx * dfx;

    }
    return p;
}

void KF4ParamsCombIV::barrelToEndcap( double r, const StubCluster *stubCluster, std::vector<double> &x, TMatrixD &cov_x )const
{

    double Inv2R = x[BP_INV2R];
    double t  = x[BP_T];

    barrelToEndcap( r, stubCluster, x );

    //    if( r != 0.1 ){
    TMatrixD G(4,4);
    /*
     * new_cov_x = y
     * y(0,0) = cov_x(0,0)^2 + cov_x(0,2)^2
     * y(0,2) = G(2,2) * ( cov_x(0,0) + cov_x(2,2) ) * cov_x(0,2)
     * y(2,0) = y(0,2)
     * y(2,2) = G(2,2)^2 * ( cov_x(0,2)^2 + cov_x(2,2)^2 ) + G(2,3) * ( cov_x(1,3)^2 + cov_x(3,3)^2 )
     * y(2,3) = G(2,3) * G(3,3) * ( cov_x(1,3)^2 + cov_x(3,3)^2 )
     * y(3,2) = y(2,3)
     * y(3,3) = G(3,3)^2 * ( cov_x(1,3)^2 + cov_x(3,3)^2 )
     */

    G(EP_RHOPHI,BP_RHOPHI)  = 1.;
    G(EP_RHOPHI,BP_Z)      = 0.;
    G(EP_RHOPHI,BP_INV2R)  = 0.; 
    G(EP_RHOPHI,BP_T)      = 0.; 

    G(EP_RHO,BP_RHOPHI)    = 0;
    G(EP_RHO,BP_Z)      = 0;
    G(EP_RHO,BP_INV2R)  = 0;
    G(EP_RHO,BP_T)      = 0;

    G(EP_INV2RT,BP_RHOPHI)    = 0.;
    G(EP_INV2RT,BP_Z)      = 0.;
    if( t ){
	G(EP_INV2RT,BP_INV2R)  = 1./t;
	G(EP_INV2RT,BP_T)      = Inv2R * -1./(t*t);
    }
    else{
	G(EP_INV2RT,BP_INV2R)  = 1.;
	G(EP_INV2RT,BP_T)      = Inv2R * -1.;
    }

    G(EP_INVT,BP_RHOPHI)    = 0;
    G(EP_INVT,BP_Z)      = 0;
    G(EP_INVT,BP_INV2R)  = 0;
    if( t ){
	G(EP_INVT,BP_T)      = -1./(t*t);
    }
    else{
	G(EP_INVT,BP_T)      = -1.;
    }

    TMatrixD org(cov_x);
    TMatrixD GT(TMatrixD::kTransposed, G );

    cov_x = G * org * GT;
    //   }
}


std::string KF4ParamsCombIV::getParams(){
    return "KF4ParamsCombIV";
}

/* Determine with a stub belongs/does not belong to a candidate
 * Decision based on hit and state uncertainty */
bool KF4ParamsCombIV::stubBelongs(const StubCluster* stubCluster, kalmanState& state, unsigned itr )const{

    return true;

}
double KF4ParamsCombIV::validationGateCutValue( const StubCluster *stubCluster, unsigned path )const
{
    double max_rchi2(100);
    if( path > 10 ) return max_rchi2;

    if( stubCluster->barrel() ){
	switch ( stubCluster->layerId() ){
	    case 3:
		max_rchi2 = 50;
	    case 4:
		max_rchi2 = 30;
	    case 5:
		max_rchi2 = 40;
	    case 6:
		max_rchi2 = 40;
	}
    }
    else{
	switch ( stubCluster->endcapRing() ){
	    case 2:
		max_rchi2 = 50;
	    case 3:
		max_rchi2 = 80;
	    case 4:
		max_rchi2 = 100;
	    case 5:
		max_rchi2 = 80;
	    case 6:
		max_rchi2 = 130;
	    case 7:
		max_rchi2 = 80;
	    case 8:
		max_rchi2 = 130;
	    case 9:
		max_rchi2 = 130;
	    case 10:
	    case 11:
	    case 12:
	    case 13:
	    case 14:
	    case 15:
		max_rchi2 = 50;
	}
    }

    return max_rchi2;
}

bool KF4ParamsCombIV::isGoodState( kalmanState &state )const
{

  vector<float> z0Cut, ptTolerance, d0Cut, chi2Cut;
  z0Cut       = { 999.,  999.,   15.,  15.,  15.,   15.,   15.};
  ptTolerance = { 999., 999., 0.3, 0.1, 0.05, 0.05, 0.05}; /// optimum for Full CKF 5 params
  d0Cut       = { 999.,  999.,  999.,  10.,   5.,    5.,   5.};  // Only used for 5 param helix fit
//  chi2Cut     = { 999.,  999.,   999.,  999.,  999.,  120.,  160.};  // Consider reducing chi2 cut 2 to 7.
  chi2Cut     = { 9999999.,  999999999.,   999999999.,  99999999999.,  9999999999.,  120.,  160.};  // Consider reducing chi2 cut 2 to 7.
  unsigned nStubLayers = state.nStubLayers();
  bool goodState( true );


  std::map<std::string, double> y = getTrackParams( &state );
  double qOverPt = y["qOverPt"]; 
  double pt=fabs( 1/qOverPt ); 
  double z0=fabs( y["z0"] ); 

  // state parameter selections

  if (z0 > z0Cut[nStubLayers] ) goodState = false;
  if( pt < getSettings()->houghMinPt() - ptTolerance[nStubLayers] ) goodState = false;

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

std::vector<double> KF4ParamsCombIV::residual(const StubCluster* stubCluster, const std::vector<double> &x )const{

    std::vector<double> vd = d(stubCluster);
    std::vector<double> hx = Hx( H(stubCluster), x ); 
    std::vector<double> delta(2);
    //    cout << "hx " << hx.at(0) << " " << hx.at(1) << endl; 
    //   cout << "vd " << vd.at(0) << " " << vd.at(1) << endl; 
    for( unsigned i=0; i<2; i++ ) delta.at(i) = vd.at(i) - hx.at(i);
    return delta;
}

}
