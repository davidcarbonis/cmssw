#ifndef PilotBladeStudy_h
#define PilotBladeStudy_h

/** \class PilotBladeStudy
 * Author: Tamas Almos Vami
 * Year: 2015
 * Desciption header file of the PilotBladeStudy.cc
 ************************************************************/

// ----------------------------------------------------------------------------------------

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
//#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TObject.h"
#include "TH1D.h"
#include "TFile.h"


#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonTrackMatcher.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/SiPixelRawData/interface/SiPixelRawDataError.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

//#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include <Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h>
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <DataFormats/Scalers/interface/Level1TriggerScalers.h>
#include <DataFormats/Common/interface/EDCollection.h>
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandPoissonQ.h"
// SimDataFormats
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// For ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStopwatch.h>

// STD
#include <memory>
#include <string>
#include <iostream>
#include <vector>

// For Resolution
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#define NOVAL_I -9999
#define NOVAL_F -9999.0

using namespace reco;
class TTree;
class TFile;
class RectangularPixelTopology;
class TObject;


class PilotBladeStudy : public edm::EDAnalyzer
{
 public:
  
  explicit PilotBladeStudy(const edm::ParameterSet&);
  virtual ~PilotBladeStudy();
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
 private:
  edm::ParameterSet iConfig_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  edm::ESHandle<MagneticField> magneticField_;
  edm::EDGetTokenT<edm::ConditionsInRunBlock> condInRunBlockToken_;
  edm::EDGetTokenT<edm::ConditionsInLumiBlock> condInLumiBlockToken_;
  
  edm::EDGetTokenT<reco::BeamSpot>                  		BSToken_;
  edm::EDGetTokenT<LumiSummary> lumiSummaryToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

  edm::EDGetTokenT< edm::DetSetVector<PixelDigi> >	siPixelDigisToken_;
  edm::EDGetTokenT< edm::DetSetVector<PixelDigi> >	PBDigisToken_;
  edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster> >	siPixelClustersToken_;
  edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster> >	PBClustersToken_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackCollToken_;
  
  edm::EDGetTokenT< edm::DetSetVector<SiPixelRawDataError> >	SiPixelRawDataErrorToken_;
  edm::EDGetTokenT<edm::EDCollection<DetId> > trackingErrorToken_;
  edm::EDGetTokenT<edm::EDCollection<DetId> > userErrorToken_;

  std::map<int,int> detIdFromFED40_;

  TTree* eventTree_;
  //TTree* lumiTree_;
  TTree* runTree_;
  TTree* trackTree_;
  TTree* digiTree_;
  TTree* clustTree_;
  TTree* trajTree_;
  TFile* outfile_;
  
  std::vector<std::string> triggerNames_; // Max 20 trigger names
  edm::InputTag triggerTag_;
  
  std::map<size_t,int> wbc;
  std::map<size_t,int> delay;
  
  bool usePixelCPE_;
  bool cosmicsCase_;
  bool isNewLS_;
  
  int nPixelHit = 0;
  int nStripHit = 0;
  int nPBHit    = 0;
    
 public:
  typedef std::vector< edm::ParameterSet > Parameters;

  struct PositionCorrection {
    float dx;
    float dy;
    PositionCorrection(float x=0.0, float y=0.0) { dx=x; dy=y; }
  };
  std::map<unsigned int,PositionCorrection> posCorr_;

  struct FiducialRegion {
    int rocX;
    int rocY;
    float marginX;
    float marginY;
    FiducialRegion(int rx, int ry, float x, float y) { rocX=rx; rocY=ry; marginX=x; marginY=y; }
    std::string print() { 
      std::stringstream ss;
      ss << "ROC(" << rocX << ", " << rocY << ") with margins ("<<marginX<<", "<<marginY<<")";
      return ss.str();
    };
  };
  std::map<unsigned int,std::vector<FiducialRegion> > fidReg_;


// ------------------------------- EventData info ------------------------------
  class EventData {
    public:
      int fill;
      int run;
      int ls;
      int orb;
      int bx;
      int evt;
      int trig;
      
      int wbc;
      int delay;

      int federrs_size;
      int federrs[16][2]; // [error index] [0:Nerror, 1:errorType]
    
      std::string list;

      EventData() { init(); };
      void init() {
        fill=NOVAL_I;
        run=NOVAL_I;
        ls=NOVAL_I;
        orb=NOVAL_I;
        bx=NOVAL_I;
        evt=NOVAL_I;
	trig=NOVAL_I;
	
	wbc=NOVAL_I;
	delay=NOVAL_I;

        federrs_size=0;
        for (size_t i=0; i<16; i++) federrs[i][0]=federrs[i][1]=NOVAL_I;

	list="fill/I:run/I:ls/I:orb/I:bx/I:evt/I:trig/I:wbc/I:delay/I:federrs_size/I:federrs[federrs_size][2]";
    }
  } evt_;
// ---------------------------- end of Event info ------------------------
  
// ------------------------------- RunData info ------------------------------
  class RunData {
   public:
    int fill;
    int run;

    std::string list;

    RunData() { init(); };
    void init() {
      fill=NOVAL_I;
      run=NOVAL_I;
      
      list="fill/I:run";
    }

  } run_;
// ---------------------------- end of Run info ----------------------------
  
// ------------------------------- TrackData info ------------------------------
  class TrackData {
   public:
    int nPixelHit;
    int nStripHit;
    float pt;
    float dxy;
    float dz;
    float eta;
    float phi;

    std::string list;
   
    TrackData() { init(); }
    void init() {
      nPixelHit=NOVAL_I;
      nStripHit=NOVAL_I;
      pt=NOVAL_F;
      dxy=NOVAL_F;
      dz=NOVAL_F;      
      eta=NOVAL_F;
      phi=NOVAL_F;

      list="nPixelHit/I:nStripHit/I:pt/F:dxy/F:dz/F:eta/F:phi/F";
    }
  };

  std::vector<TrackData> tracks_;
// --------------------------- end of Track info ---------------------------  
  
// ------------------------------- DigiData info ------------------------------
  class DigiData {
   public:
    int i;
    int row;
    int col;
    int adc;

    std::string list;

    DigiData() { init(); }
    void init() {
      i=NOVAL_I;
      row=NOVAL_I;
      col=NOVAL_I;
      adc=NOVAL_I;
      list="i/I:row:col:adc";
    }
  };
// ----------------------------- end of Digi info -----------------------------

// ------------------------------ ClusterData info ----------------------------
  class ClustData {
   public:
    float x;
    float y;
    float lx;
    float ly;
    float glx;
    float gly;
    float glz;
    
    int size;
    float charge;
    
    int nclu_mod;
    
    std::string list;

    ClustData() { init(); }
    void init() {
      x=NOVAL_F;
      y=NOVAL_F;
      lx=NOVAL_F;
      ly=NOVAL_F;
      glx=NOVAL_F;
      gly=NOVAL_F;
      glz=NOVAL_F;
      
      size=0;
      charge=NOVAL_F;
      
      nclu_mod=NOVAL_I;

      list="x/F:y/F:lx/F:ly/F:glx/F:gly/F:glz/F:size/I:charge/F:nclu_mod/I";
    }

    int getROCx() { return (x>=0 && x<160) ? int(x)/80-1+(int(x)/80) : -9999;  }
    int getROCy() { return (y>=0 && y<416) ? int(y)/52-4+(int(y)/208) : -9999; } 
  };
// --------------------------- end of ClusterData info ---------------------------- 
  
// ------------------------------ ModuleData info ---------------------------------
  class ModuleData {
   public:
    int det;
    int layer;
    int ladder;
    int half;
    int module;
    int disk;
    int blade;
    int panel;
    
    int federr;

    int side;
    int shl;

    int outer;

    unsigned int rawid;
    
    std::map<int, std::string> federrortypes;

    std::string list;

    ModuleData() { init(); }
    void init() {
      det=NOVAL_I;
      layer=NOVAL_I;
      ladder=NOVAL_I;
      half=NOVAL_I;
      module=NOVAL_I;
      disk=NOVAL_I;
      blade=NOVAL_I;
      panel=NOVAL_I;
      federr = NOVAL_I;
      side=NOVAL_I;
      shl=NOVAL_I;
      outer=NOVAL_I;
      rawid=abs(NOVAL_I);

      federrortypes.insert(std::pair<int, std::string>(25, "invalidROC"));
      federrortypes.insert(std::pair<int, std::string>(26, "gap word"));
      federrortypes.insert(std::pair<int, std::string>(27, "dummy word"));
      federrortypes.insert(std::pair<int, std::string>(28, "FIFO full error"));
      federrortypes.insert(std::pair<int, std::string>(29, "timeout error"));
      federrortypes.insert(std::pair<int, std::string>(30, "TBM error trailer"));
      federrortypes.insert(std::pair<int, std::string>(31, "event number error (TBM and FED event number mismatch)"));
      federrortypes.insert(std::pair<int, std::string>(32, "incorrectly formatted Slink Header"));
      federrortypes.insert(std::pair<int, std::string>(33, "incorrectly formatted Slink Trailer"));
      federrortypes.insert(std::pair<int, std::string>(34, "the event size encoded in the Slink Trailer is different than the size found at raw to digi conversion "));
      federrortypes.insert(std::pair<int, std::string>(35, "invalid FED channel number"));
      federrortypes.insert(std::pair<int, std::string>(36, "invalid ROC value "));
      federrortypes.insert(std::pair<int, std::string>(37, "invalid dcol or pixel value "));
      federrortypes.insert(std::pair<int, std::string>(38, "the pixels on a ROC weren't read out from lowest to highest row and dcol value"));
      federrortypes.insert(std::pair<int, std::string>(39, "CRC error"));
      federrortypes.insert(std::pair<int, std::string>(40, "Overflow error"));
      
      list="det/I:layer/I:ladder/I:half/I:module/I:disk/I:blade/I:panel/I:"
      "federr/I:side/I:shl/I:outer/I:rawid/i";
    }

    std::string shell() {
      std::ostringstream ss;
      if (det==0) {
	ss << "B" << ((module>0) ? "p" : "m") << ((ladder>0) ? "I" : "O");
      } else if (det==1) {
	ss << "B" << ((disk>0) ? "p" : "m") << ((blade>0) ? "I" : "O");
      }
      return ss.str();
    }

    int shell_num() {
      if (det==0)      return ((module>0) ? 0 : 2) + ((ladder>0) ? 0 : 1);
      else if (det==1) return ((disk>0) ? 0 : 2) + ((blade>0) ? 0 : 1);
      return -1;
    }
    
    std::string federr_name() {
      std::map<int, std::string>::const_iterator it=federrortypes.find(federr);
      return (it!=federrortypes.end()) ? it->second : "FED error not interpreted";
    }
    
  };
// --------------------------- end of ModuleData info ---------------------------------

// ------------------------------ TrajectoryData info ---------------------------------
  class TrajMeasData {
   public:
    int type;
    
    int nPixelHit;
    int nStripHit;
    int nPBHit;
    
    float lx;
    float ly;
    float lx_err;
    float ly_err;
    float glx;
    float gly;
    float glz;
    
    int onEdge;
    float alpha;
    float beta;
    
    float dx_cl;
    float dy_cl;
//  float dx_cl_corr;
//  float dy_cl_corr;
    float d_cl;

    float dx_hit;
    float dy_hit;
    float d_hit;


    //float norm_charge;

    std::string list;

    TrajMeasData() { init(); }
    void init() {
      type=NOVAL_I;
      
      nPixelHit = NOVAL_I;
      nStripHit = NOVAL_I;
      nPBHit    = NOVAL_I;
      
      lx=NOVAL_F;
      ly=NOVAL_F;
      lx_err=NOVAL_F;
      ly_err=NOVAL_F;
      glx=NOVAL_F;
      gly=NOVAL_F;
      glz=NOVAL_F;
      
      onEdge=NOVAL_I;
      alpha=NOVAL_F;
      beta=NOVAL_F;
            
      dx_cl=NOVAL_F;
      dy_cl=NOVAL_F;
//    dx_cl_corr=NOVAL_F;
//    dx_cl_corr=NOVAL_F;
      d_cl=NOVAL_F;
      
      dx_hit=NOVAL_F;
      dy_hit=NOVAL_F;
      d_hit=NOVAL_F;
//       norm_charge=NOVAL_F;

//       list="lx/F:ly/F:glx/F:gly/F:glz/F:res_dx/F:res_dz/F:hit_near/I:nPixelHit/I:nStripHit/I:nPBHits/I"
//       "dx_cl[2]/F:dx_cl_corr[2]/F:dy_cl[2]/F:dy_cl_corr[2]/F:dx_hit/F:dy_hit/F:"
//       "i/I:onEdge/I:type/I:lx_err/F:ly_err/F:d_cl[2]/F:alpha/F:beta/F:norm_charge/F";
      
      list="type/I:nPixelHit/I:nStripHit/I:nPBHit/I:lx/F:ly/F:lx_err/F:ly_err/F:glx/F:gly/F:glz/F:"
      "onEdge/I:alpha/F:beta/F:dx_cl/F:dy_cl/F:d_cl/F:dx_hit/F:dy_hit/F:d_hit/F";
    }

    int getROCx() { return (int(lx/0.8164)+int(lx/fabs(lx)));  }
    int getROCy() { return (int(ly/0.81)+int(ly/fabs(ly)));    }
    bool isWithinROCFiducial(float margX, float margY) { 
      return (fabs(fabs(fmod(fabs(ly),0.8164)-0.4082)-0.4082)>margX && fabs(fabs(fmod(fabs(ly),0.81)-0.405)-0.405)>margY);
    }
  };
// -------------------------- end of TrajectoryData info ------------------------------  
// ---------------------- end of *Data classes ----------------------
  
// ------------- from here we have the merged classes ---------------
// -------------------------- Digi ----------------------------------
  class Digi : public DigiData {
   public:

    ModuleData mod; // offline module number
    ModuleData mod_on; // online module number
    
    Digi() { mod.init(); mod_on.init(); }
    void init() {
      DigiData::init();
      mod.init();
      mod_on.init();
    }
  };

  std::vector<Digi> digis_;

  class Cluster : public ClustData {
   public:

    ModuleData mod; // offline module number
    ModuleData mod_on; // online module number

    Cluster() { mod.init(); mod_on.init();}
    void init() {
      ClustData::init();
      mod.init();
      mod_on.init();
    }
  };

  std::vector<Cluster> clust_;


  class TrajMeasurement : public TrajMeasData {
   public:
    ModuleData mod; // offline module number
    ModuleData mod_on; // online module number
    ClustData clu;
    TrackData trk;

    TrajMeasurement() { mod.init(); mod_on.init(); trk.init(); }
    void init() {
      TrajMeasData::init();
      mod.init();
      mod_on.init();
      trk.init();
    }
  };

  std::vector<std::vector<TrajMeasurement> > trajmeas_;
  
  void init_all() {
    evt_.init();
    tracks_.clear();
    trajmeas_.clear();
    clust_.clear();
    digis_.clear();
    nPixelHit=0;
    nStripHit=0;
    nPBHit=0;
  }

  
// ------------------------------ Used functions ------------------------------  
  ModuleData getModuleData(uint32_t rawId, const std::map<uint32_t, int>& federrors, std::string scheme="offline");

  void readFEDErrors(const edm::Event&, const edm::EventSetup&, 
		     edm::EDGetTokenT< edm::DetSetVector<SiPixelRawDataError> >, 
		     edm::EDGetTokenT< edm::EDCollection<DetId> >, std::map<uint32_t, int>&);
  
  void analyzeDigis(const edm::Event&, const edm::EventSetup&, 
		    edm::EDGetTokenT< edm::DetSetVector<PixelDigi> >, std::map<uint32_t, int>, 
		    int, bool);
  void analyzeClusters(const edm::Event&, const edm::EventSetup&, 
		       edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster> >, std::map<uint32_t, int>,
		       int, bool);
  void analyzeTracks(const edm::Event&, const edm::EventSetup&, 
		       edm::EDGetTokenT< TrajTrackAssociationCollection >, std::map<uint32_t, int>,
		       int, bool);
  void analyzeTrajs(const edm::Event&, const edm::EventSetup&, 
		       edm::EDGetTokenT< TrajTrackAssociationCollection >, std::map<uint32_t, int>,
		       int, bool);
  				
  void findClosestClusters(const edm::Event&, const edm::EventSetup&, uint32_t, 
			   float, float, float&, float&, 
			   edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster> >, ClustData&);
  
  void readExtraInfos();
			    

};

#endif
