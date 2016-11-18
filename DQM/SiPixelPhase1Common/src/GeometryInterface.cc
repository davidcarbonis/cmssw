// -*- C++ -*-
//
// Package:    SiPixelPhase1Common
// Class:      GeometryInterface
//
// Geometry depedence goes here.
//
// Original Author:  Marcel Schneider

#include "DQM/SiPixelPhase1Common/interface/GeometryInterface.h"

// edm stuff
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Tracker Geometry/Topology  stuff
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/GeometryObjects/interface/PTrackerParameters.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// C++ stuff
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iomanip>

// WTH is this needed? clang wants it for linking...
const GeometryInterface::Value GeometryInterface::UNDEFINED;

void GeometryInterface::load(edm::EventSetup const& iSetup) {
  //loadFromAlignment(iSetup, iConfig);
  loadFromTopology(iSetup, iConfig);
  loadTimebased(iSetup, iConfig);
  loadModuleLevel(iSetup, iConfig);
  loadFEDCabling(iSetup, iConfig);
  edm::LogInfo log("GeometryInterface");
  log << "Known colum names:\n";
  for (auto e : ids) log << "+++ column: " << e.first
    << " ok " << bool(extractors[e.second]) << " min " << min_value[e.second] << " max " << max_value[e.second] << "\n";
  is_loaded = true;
}

void GeometryInterface::loadFromTopology(edm::EventSetup const& iSetup, const edm::ParameterSet& iConfig) {
  // Get a Topology
  edm::ESHandle<TrackerTopology> trackerTopologyHandle;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopologyHandle);
  assert(trackerTopologyHandle.isValid());

  std::vector<ID> geomquantities;

  struct TTField {
    const TrackerTopology* tt;
    TrackerTopology::DetIdFields field;
    Value operator()(InterestingQuantities const& iq) {
      if (tt->hasField(iq.sourceModule, field))
        return tt->getField(iq.sourceModule, field);
      else
        return UNDEFINED;
    };
  };

  const TrackerTopology* tt = trackerTopologyHandle.operator->();


  std::vector<std::pair<std::string, TTField>> namedPartitions {
    {"PXEndcap" , {tt, TrackerTopology::PFSide}},

    {"PXLayer"  , {tt, TrackerTopology::PBLayer}},
    {"PXLadder" , {tt, TrackerTopology::PBLadder}},
    {"PXBModule", {tt, TrackerTopology::PBModule}},

    {"PXBlade"  , {tt, TrackerTopology::PFBlade}},
    {"PXDisk"   , {tt, TrackerTopology::PFDisk}},
    {"PXPanel"  , {tt, TrackerTopology::PFPanel}},
    {"PXFModule", {tt, TrackerTopology::PFModule}},
  };

  for (auto& e : namedPartitions) {
    geomquantities.push_back(intern(e.first));
    addExtractor(intern(e.first), e.second, UNDEFINED, UNDEFINED);
  }

  auto pxbarrel  = [] (InterestingQuantities const& iq) { return iq.sourceModule.subdetId() == PixelSubdetector::PixelBarrel ? 0 : UNDEFINED; };
  auto pxforward = [] (InterestingQuantities const& iq) { return iq.sourceModule.subdetId() == PixelSubdetector::PixelEndcap ? 0 : UNDEFINED; };
  addExtractor(intern("PXBarrel"),  pxbarrel,  0, 0);
  addExtractor(intern("PXForward"), pxforward, 0, 0);

  // Redefine the disk numbering to use the sign
  auto pxendcap = extractors[intern("PXEndcap")];
  auto diskid = intern("PXDisk");
  auto pxdisk = extractors[diskid];
  extractors[diskid] = [pxdisk, pxendcap] (InterestingQuantities const& iq) {
    auto disk = pxdisk(iq);
    if (disk == UNDEFINED) return UNDEFINED;
    auto endcap = pxendcap(iq);
    //if (disk != 3) return UNDEFINED; //HACK to hide no-PB, saves space.
    return endcap == 1 ? -disk : disk;
  };

  // Get a Geometry
  edm::ESHandle<TrackerGeometry> trackerGeometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometryHandle);
  assert(trackerGeometryHandle.isValid());

  // some parameters to record the ROCs here
  auto module_rows = iConfig.getParameter<int>("module_rows") - 1;
  auto module_cols = iConfig.getParameter<int>("module_cols") - 1;

  // We need to track some extra stuff here for the Shells later.
  auto pxlayer  = extractors[intern("PXLayer")];
  auto pxladder = extractors[intern("PXLadder")];
  auto pxmodule = extractors[intern("PXBModule")];
  auto pxblade  = extractors[intern("PXBlade")];
  std::vector<Value> maxladders;
  Value maxmodule = 0;
  Value innerring = iConfig.getParameter<int>("n_inner_ring_blades");
  Value outerring = 0;

  // Now travrse the detector and collect whatever we need.
  auto detids = trackerGeometryHandle->detIds();
  for (DetId id : detids) {
    if (id.subdetId() != PixelSubdetector::PixelBarrel && id.subdetId() != PixelSubdetector::PixelEndcap) continue;
    auto iq = InterestingQuantities{nullptr, id, 0, 0};
    auto layer = pxlayer(iq);
    if (layer != UNDEFINED) {
      if (layer >= Value(maxladders.size())) maxladders.resize(layer+1);
      auto ladder = pxladder(iq);
      if (ladder > maxladders[layer]) maxladders[layer] = ladder;
    }
    auto module = pxmodule(iq);
    if (module != UNDEFINED && module > maxmodule) maxmodule = module;
    auto blade = pxblade(iq);
    if (blade != UNDEFINED && blade > outerring) outerring = blade;

    // we record each module 4 times, one for each corner, so we also get ROCs
    // in booking (at least for the ranges)
    iq.row = 0; iq.col = 0;
    all_modules.push_back(iq);
    iq.row = module_rows; iq.col = 0;
    all_modules.push_back(iq);
    iq.row = 0; iq.col = module_cols;
    all_modules.push_back(iq);
    iq.row = module_rows; iq.col = module_cols;
    all_modules.push_back(iq);
  }

  outerring = outerring - innerring;

  // Shells are a concept that cannot be derived from bitmasks. 
  // Use hardcoded logic here.
  // This contains a lot more assumptions about general geometry than the rest
  // of the code, but it might work for Phase0 as well.
  addExtractor(intern("PXRing"), 
    [pxblade, innerring] (InterestingQuantities const& iq) {
      auto blade = pxblade(iq);
      if (blade == UNDEFINED) return UNDEFINED;
      if (blade <= innerring) return Value(1);
      else return Value(2);
    }
  );

  addExtractor(intern("HalfCylinder"),
    [pxendcap, pxblade, innerring, outerring] (InterestingQuantities const& iq) {
      auto ec = pxendcap(iq);
      if (ec == UNDEFINED) return UNDEFINED;
      auto blade = pxblade(iq);
      // blade 1 and 56 are at 3 o'clock. This is a mess.
      auto inring  = blade > innerring ? (innerring+outerring+1) - blade : blade;
      auto perring = blade > innerring ? outerring : innerring;

      // inring is now 1-based, 1 at 3 o'clock, upto perring.
      int frac = (int) ((inring-1) / float(perring) * 4); // floor semantics here
      if (frac == 0 || frac == 3) return 10*ec + 1; // inner half
      if (frac == 1 || frac == 2) return 10*ec + 2; // outer half
      assert(!"HalfCylinder logic problem");
      return UNDEFINED;
    }, 0, 0 // N/A
  );

  // For the '+-shape' (ladder vs. module) plots, we need online numbers with
  // (unused) 0-ladder/module at x=0/z=0. This means a lot of messing with the
  // ladder/shell numbering...
  addExtractor(intern("onlineLadder"),
    [pxbarrel, pxladder, pxlayer, maxladders] (InterestingQuantities const& iq) {
      if(pxbarrel(iq) == UNDEFINED) return UNDEFINED;
      auto layer  = pxlayer(iq);
      auto ladder = pxladder(iq);
      int frac = (int) ((ladder-1) / float(maxladders[layer]) * 4); // floor semantics
      Value quarter = maxladders[layer] / 4;
      if (frac == 0) return -ladder + quarter + 1; // top right - +1 for gap
      if (frac == 1) return -ladder + quarter; // top left -
      if (frac == 2) return -ladder + quarter; // bot left - same
      if (frac == 3) return -ladder  + 4*quarter + quarter + 1; // bot right - like top right but wrap around
      assert(!"Shell logic problem");
      return UNDEFINED;
    }
  );

  addExtractor(intern("onlineModule"),
    [pxmodule, maxmodule] (InterestingQuantities const& iq) {
      Value mod = pxmodule(iq);  // range 1..maxmodule
      if (mod == UNDEFINED) return UNDEFINED;
      mod -= (maxmodule/2 + 1); // range -(max_module/2)..-1, 0..
      if (mod >= 0) mod += 1;    // range -(max_module/2)..-1, 1..
      return mod;
    }
  );

  addExtractor(intern("onlineBlade"),
    [pxendcap, pxblade] (InterestingQuantities const& iq) {
      if(pxendcap(iq) == UNDEFINED) return UNDEFINED;
      auto blade = pxblade(iq);
      // Ring 1
      if ( 1 <= blade && blade < 6 ) return 6 - blade; // 5 on 1st quarter
      if ( 6 <= blade && blade < 17 ) return 5 - blade; // 11 on 2nd half
      if ( 17 <= blade && blade < 23 ) return 28 - blade; // 6 on 4th quarter
      // Ring 2
      if ( 23 <= blade && blade < 31 ) return 31 - blade; // 8 on 1st quarter
      if ( 31 <= blade && blade < 48 ) return 30 - blade; // 17 on 2nd half
      if ( 48 <= blade && blade < 57 ) return 65 - blade; // 9 on 4th quarter

      assert(!"Shell logic problem");
      return UNDEFINED;
    }
  );

  auto onlineladder = extractors[intern("onlineLadder")];
  auto onlinemodule = extractors[intern("onlineModule")];
  addExtractor(intern("Shell"),
    [onlineladder, onlinemodule] (InterestingQuantities const& iq) {
      auto sl = onlineladder(iq);
      auto sm = onlinemodule(iq);
      if (sl == UNDEFINED) return UNDEFINED;
      return Value((sm < 0 ? 10 : 20) + (sl < 0 ? 2 : 1)); // negative means outer shell!?
    }, 0, 0 // N/A
  );

  addExtractor(intern(""), // A dummy column. Not much special handling required.
    [] (InterestingQuantities const& iq) { return 0; },
    0, 0
  );

}

void GeometryInterface::loadTimebased(edm::EventSetup const& iSetup, const edm::ParameterSet& iConfig) {
  // extractors for quantities that are roughly time-based. We cannot book plots based on these; they have to
  // be grouped away in step1.
  addExtractor(intern("Lumisection"),
    [] (InterestingQuantities const& iq) {
      if(!iq.sourceEvent) return UNDEFINED;
      return Value(iq.sourceEvent->luminosityBlock());
    },
    1, iConfig.getParameter<int>("max_lumisection")
  );
  int onlineblock = iConfig.getParameter<int>("onlineblock");
  int n_onlineblocks = iConfig.getParameter<int>("n_onlineblocks");
  addExtractor(intern("OnlineBlock"),
    [onlineblock] (InterestingQuantities const& iq) {
      if(!iq.sourceEvent) return UNDEFINED;
      return Value(onlineblock + iq.sourceEvent->luminosityBlock() / onlineblock);
    },
    // note: this range is not visible anywhere (if the RenderPlugin does its job),
    // but the strange range allows the RenderPlugin to know the block size.
    onlineblock, onlineblock+n_onlineblocks-1
  );
  addExtractor(intern("BX"),
    [] (InterestingQuantities const& iq) {
      if(!iq.sourceEvent) return UNDEFINED;
      return Value(iq.sourceEvent->bunchCrossing());
    },
    1, iConfig.getParameter<int>("max_bunchcrossing")
  );
}

void GeometryInterface::loadModuleLevel(edm::EventSetup const& iSetup, const edm::ParameterSet& iConfig) {
  // stuff that is within modules. Might require some phase0/phase1/strip switching later
  addExtractor(intern("row"),
    [] (InterestingQuantities const& iq) {
      return Value(iq.row);
    },
    0, iConfig.getParameter<int>("module_rows") - 1
  );
  addExtractor(intern("col"),
    [] (InterestingQuantities const& iq) {
      return Value(iq.col);
    },
    0, iConfig.getParameter<int>("module_cols") - 1
  );

  int   n_rocs     = iConfig.getParameter<int>("n_rocs");
  float roc_cols   = iConfig.getParameter<int>("roc_cols");
  float roc_rows   = iConfig.getParameter<int>("roc_rows");
  auto  pxmodule   = extractors[intern("PXBModule")];
  auto  pxpanel    = extractors[intern("PXPanel")];
<<<<<<< HEAD
=======
  auto pxladder = extractors[intern("PXLadder")];
  auto pxlayer  = extractors[intern("PXLayer")];

  auto pxendcap = extractors[intern("PXEndcap")];
  auto pxblade  = extractors[intern("PXBlade")];
  Value innerring = iConfig.getParameter<int>("n_inner_ring_blades");
  Value outerring = 0;

  Value maxmodule = 0;
  std::vector<Value> maxladders;

  auto detids = trackerGeometryHandle->detIds();
  for (DetId id : detids) {
    auto iq = InterestingQuantities{.sourceModule = id };
    auto module = pxmodule(iq);
    if (module != UNDEFINED && module > maxmodule) maxmodule = module;

    auto blade = pxblade(iq);
    if (blade != UNDEFINED && blade > outerring) outerring = blade;

    if (id.subdetId() != PixelSubdetector::PixelBarrel && id.subdetId() != PixelSubdetector::PixelEndcap) continue;
    auto layer = pxlayer(iq);
    if (layer != UNDEFINED) {
      if (layer >= Value(maxladders.size())) maxladders.resize(layer+1);
      auto ladder = pxladder(iq);
      if (ladder > maxladders[layer]) maxladders[layer] = ladder;
    }
  }

  outerring = outerring - innerring;

  addExtractor(intern("ROC"),
    [n_rocs, roc_cols, roc_rows] (InterestingQuantities const& iq) {
      int fedrow = int(iq.row / roc_rows);
      int fedcol = int(iq.col / roc_cols);
      if (fedrow == 0) return Value(fedcol);
      if (fedrow == 1) return Value(n_rocs - 1 - fedcol);
      return UNDEFINED;
    }
  );

  // arbitrary per-ladder numbering (for inefficiencies)
  auto roc = extractors[intern("ROC")];
  addExtractor(intern("ROCinLadder"),
    [pxmodule, roc, n_rocs] (InterestingQuantities const& iq) {
      auto mod = pxmodule(iq);
      if (mod == UNDEFINED) return UNDEFINED;
      return Value(roc(iq) + n_rocs * (mod-1));
    }
  );
  addExtractor(intern("ROCinBlade"),
    [pxmodule, pxpanel, roc, n_rocs] (InterestingQuantities const& iq) {
      auto mod = pxpanel(iq);
      if (mod == UNDEFINED) return UNDEFINED;
      return Value(roc(iq) + n_rocs * (mod-1));
    }
  );

  addExtractor(intern("ROCinLayerRow"),
    [pxladder, pxlayer, maxladders, roc, roc_rows] (InterestingQuantities const& iq) {

      auto ladder = pxladder(iq);
      if (ladder == UNDEFINED) return UNDEFINED;

      int rocRow = int(iq.row / roc_rows);

      auto layer  = pxlayer(iq);
      int frac = (int) ((ladder-1) / float(maxladders[layer]) * 4); // floor semantics
      Value quarter = maxladders[layer] / 4;

      if (frac == 0) return Value( (-ladder + quarter + 1) * 2 + (rocRow) -1 );
      if (frac == 1) return Value( (-ladder + quarter) * 2 + (rocRow) );
      if (frac == 2) return Value( (-ladder + quarter) * 2 + (rocRow) );
      if (frac == 3) return Value( (-ladder  + 4*quarter + quarter + 1) * 2 + (rocRow) -1 );
      assert(!"Shell logic problem");
      return UNDEFINED;

    }
  );
  addExtractor(intern("ROCinLayerCol"),
    [pxmodule, maxmodule, roc, roc_cols] (InterestingQuantities const& iq) {
      auto mod = pxmodule(iq);
      if (mod == UNDEFINED) return UNDEFINED;

      int rocCol = int(iq.col / roc_cols);

      mod -= (maxmodule/2 + 1);
      if (mod >= 0) return Value ( -(mod + 1) * 8 + rocCol );    // range -(max_module/2)..-1, 1..
      else return Value ( -mod * 8 + rocCol -7 );    // range -(max_module/2)..-1, 0..

      assert(!"Shell logic problem");
      return UNDEFINED;
    }
  );

  addExtractor(intern("ROCinDiskRow"),
    [pxendcap, pxpanel, roc, roc_rows] (InterestingQuantities const& iq) {
      auto ec = pxendcap(iq);
      if (ec == UNDEFINED) return UNDEFINED;
      auto panel = pxpanel(iq);
      if (panel == UNDEFINED) return UNDEFINED;

      int rocRow = int(iq.row / roc_rows);

       return Value ( panel * 2 + rocRow );

      assert(!"Shell logic problem");
      return UNDEFINED;
    }
  );


  addExtractor(intern("ROCinDiskCol"),
    [pxendcap, pxblade, roc, roc_cols] (InterestingQuantities const& iq) {
      auto ec = pxendcap(iq);
      if (ec == UNDEFINED) return UNDEFINED;
      auto blade = pxblade(iq);
      if (blade == UNDEFINED) return UNDEFINED;

      int rocCol = int(iq.col / roc_cols);

      // Ring 1
      if ( 1 <= blade && blade < 6 ) return (6 - blade)*12 + rocCol -11; // 5 on 1st quarter
      if ( 6 <= blade && blade < 17 ) return (5 - blade)*12 + rocCol; // 11 on 2nd half
      if ( 17 <= blade && blade < 23 ) return (28 - blade)*12 + rocCol -11; // 6 on 4th quarter
      // Ring 2 
      if ( 23 <= blade && blade < 31 ) return (31 - blade)*12 + rocCol -11; // 8 on 1st quarter
      if ( 31 <= blade && blade < 48 ) return (30 - blade)*12 + rocCol; // 17 on 2nd half
      if ( 48 <= blade && blade < 57 ) return (65 - blade)*12 + rocCol -11; // 9 on 4th quarter

//      return Value ( blade * 8 + rocCol );
      assert(!"Blade logic problem.");
      return UNDEFINED;
    }
  );

  addExtractor(intern("DetId"),
    [] (InterestingQuantities const& iq) {
      uint32_t id = iq.sourceModule.rawId();
      return Value(id);
    }//,
//    0, 0 // No sane value possible here.
  );
}

void GeometryInterface::loadFEDCabling(edm::EventSetup const& iSetup, const edm::ParameterSet& iConfig) {
  auto cablingMapLabel = iConfig.getParameter<std::string>("CablingMapLabel");
  edm::ESHandle<SiPixelFedCablingMap> theCablingMap;
  iSetup.get<SiPixelFedCablingMapRcd>().get(cablingMapLabel, theCablingMap);
  std::map<DetId, Value> fedmap;
  uint32_t minFED = UNDEFINED, maxFED = 0;

  if (theCablingMap.isValid()) {
    auto map = theCablingMap.product();

    for(auto iq : all_modules) {
      std::vector<sipixelobjects::CablingPathToDetUnit> paths = map->pathToDetUnit(iq.sourceModule.rawId());
      for (auto p : paths) {
        //std::cout << "+++ cabling " << iq.sourceModule.rawId() << " " << p.fed << " " << p.link << " " << p.roc << "\n";
        fedmap[iq.sourceModule] = Value(p.fed);
        if (p.fed > maxFED) maxFED = p.fed;
        if (p.fed < minFED) minFED = p.fed;
      }
    }
  } else {
    edm::LogError("GeometryInterface") << "+++ No cabling map. Cannot extract FEDs.\n";
  }

  addExtractor(intern("FED"),
    [fedmap] (InterestingQuantities const& iq) {
      if (iq.sourceModule == 0xFFFFFFFF)
        return Value(iq.col); // hijacked for the raw data plugin
      auto it = fedmap.find(iq.sourceModule);
      if (it == fedmap.end()) return GeometryInterface::UNDEFINED;
      return it->second;
    }
  );
  addExtractor(intern("FEDChannel"),
    [] (InterestingQuantities const& iq) {
      // TODO: we also should be able to compute the channel from the ROC.
      // But for raw data, we only need this hack.
      //if (iq.sourceModule == 0xFFFFFFFF)
      return Value(iq.row); // hijacked for the raw data plugin
    },
    0, 39 // TODO: real range
  );
}
