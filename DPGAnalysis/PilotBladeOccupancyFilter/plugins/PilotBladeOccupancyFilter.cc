// -*- C++ -*-
//
// Package:    DPGAnalysis/PilotBladeOccupancyFilter
// Class:      PilotBladeOccupancyFilter
// 
/**\class PilotBladeOccupancyFilter PilotBladeOccupancyFilter.cc DPGAnalysis/PilotBladeOccupancyFilter/plugins/PilotBladeOccupancyFilter.cc

 Description: Filter based on number of clusters in Pixel close to the Pilot Blade

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marcel Andre Schneider
//         Created:  Tue, 15 Nov 2016 13:57:22 GMT
//
//


// system include files
#include <memory>
#include <stdint.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"


//
// class declaration
//

class PilotBladeOccupancyFilter : public edm::stream::EDFilter<> {
   public:
      explicit PilotBladeOccupancyFilter(const edm::ParameterSet&);
      ~PilotBladeOccupancyFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
      edm::EDGetTokenT<edm::DetSetVector<SiPixelCluster>> clusterToken_;
      int digisCut_;
      int clustersCut_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PilotBladeOccupancyFilter::PilotBladeOccupancyFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   digiToken_ = consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"));
   clusterToken_ = consumes<edm::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters"));
   digisCut_ = iConfig.getParameter<int>("mindigis");
   clustersCut_ = iConfig.getParameter<int>("minclusters");

}


PilotBladeOccupancyFilter::~PilotBladeOccupancyFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
PilotBladeOccupancyFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::DetSetVector<PixelDigi>> digis;
   iEvent.getByToken(digiToken_, digis);
   if (digisCut_ && !digis.isValid()) {
     std::cout << "+++ digis cut but no digis collection!\n";
     return false; 
   }
   edm::Handle<edm::DetSetVector<SiPixelCluster>> clusters;
   iEvent.getByToken(clusterToken_, clusters);
   if (clustersCut_ && !clusters.isValid()) {
     std::cout << "+++ clusters cut but no clusters collection!\n";
     return false; 
   }

   //Retrieve tracker topology and geometry
   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   const TrackerTopology* const tTopo = tTopoHandle.product();
   
   std::set<uint32_t> modules = {
                      344080644, 344080648, 344080652, 344080656, // 10 Panel 1
                      344080900, 344080904, 344080908,            // 10 Panel 2
                      344081668, 344081672, 344081676, 344081680, // 11 Panel 1
                      344081924, 344081928, 344081932             // 11 Panel 2
   };

   int ndigis = 0;
   int nclusters = 0;

   if (digisCut_) {
     edm::DetSetVector<PixelDigi>::const_iterator it;
     for (it = digis->begin(); it != digis->end(); ++it) {
       auto id = DetId(it->detId());
       if (id.subdetId() == PixelSubdetector::PixelEndcap 
         && tTopo->pxfDisk(id) == 2 && tTopo->pxfSide(id) == 1) {
          ndigis += it->size();
       }
     }
  }
  if (clustersCut_) {
     edm::DetSetVector<SiPixelCluster>::const_iterator cit;
     for (cit = clusters->begin(); cit != clusters->end(); ++cit) {
       auto id = DetId(cit->detId());
       if (id.subdetId() == PixelSubdetector::PixelEndcap 
         && tTopo->pxfDisk(id) == 2 && tTopo->pxfSide(id) == 1) {
         for (auto c : *cit) {
           if (c.size() > 1)
		nclusters++;
         }
       }
     }
   }
   return  ndigis >= digisCut_ && nclusters >= clustersCut_;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PilotBladeOccupancyFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PilotBladeOccupancyFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PilotBladeOccupancyFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PilotBladeOccupancyFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PilotBladeOccupancyFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PilotBladeOccupancyFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PilotBladeOccupancyFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PilotBladeOccupancyFilter);
