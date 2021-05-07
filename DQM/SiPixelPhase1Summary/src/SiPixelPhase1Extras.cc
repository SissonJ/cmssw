// -*- C++ -*-
//
// Package:    SiPixelPhase1Extras
// Class:      SiPixelPhase1Extras
//
/**\class 

 Description: Create the Phsae 1 pixel summary map

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Duncan Leggat
//         Created:  5th December 2016
//
//
#include "DQM/SiPixelPhase1Summary/interface/SiPixelPhase1Extras.h"
// Framework
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// DQM Framework
#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
#include "DQMServices/Core/interface/DQMStore.h"
// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
// DataFormats
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelNameUpgrade.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapNameUpgrade.h"
//
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace edm;

SiPixelPhase1Extras::SiPixelPhase1Extras(const edm::ParameterSet& iConfig) : conf_(iConfig) {
  LogInfo("PixelDQM") << "SiPixelPhase1Extras::SiPixelPhase1Extras: Hello!" << endl;
  effFolderName_ = conf_.getParameter<std::string>("EffFolderName");
  vtxFolderName_ = conf_.getParameter<std::string>("VtxFolderName");

}

SiPixelPhase1Extras::~SiPixelPhase1Extras() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  LogInfo("PixelDQM") << "SiPixelPhase1Extras::~SiPixelPhase1Extras: Destructor" << endl;
}

void SiPixelPhase1Extras::beginRun(edm::Run const& run, edm::EventSetup const& eSetup) {}


//------------------------------------------------------------------
// Method called for every event
//------------------------------------------------------------------
void SiPixelPhase1Extras::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {

  /// put in here the booking of histograms -- don't need a separate function I don't think!

  /// put in here the filling of those histograms
  /// can use "effFolderName_" and "vtxFolderName_" as variables
  //

  iBooker.setCurrentFolder(effFolderName_);

  //Book the new histos
  MonitorElement * eff_v_vtx_barrel = iBooker.book2D("hitefficiency_per_meanNvtx_per_PXLayer_PXBarrel", "hitefficiency_per_meanNvtx_per_PXLayer_PXBarrel; meanNvtx; PXLayer",101,-.5,100.5,3,.5,3.5);

  MonitorElement * eff_v_vtx_forward = iBooker.book2D("hitefficiency_per_meanNvtx_per_PXDisk_PXForward", "hitefficiency_per_meanNvtx_per_PXDisk_PXForward; meanNvtx; PXDisk",101,-.5,100.5,7,-3.5,3.5);
  
  //Get the existing histos
  MonitorElement * vtx_v_lumi = iGetter.get(vtxFolderName_ + "/NumberOfGoodPVtxVsLS_GenTk");
  
  MonitorElement * eff_v_lumi_barrel = iGetter.get(effFolderName_ + "/hitefficiency_per_Lumisection_per_PXLayer_PXBarrel");

  MonitorElement * eff_v_lumi_forward = iGetter.get(effFolderName_ + "/hitefficiency_per_Lumisection_per_PXDisk_PXForward");
  
  //initialize variables
  int numLumi = int(vtx_v_lumi->getTH1()->GetNbinsX());
  int nvtx = 0;
  double eff = 0.0;
  
  //For loop to loop through lumisections
  for(int iLumi = 1; iLumi<numLumi-1; iLumi++)
  {
    //get the meanNvtx for each lumi
    nvtx = int(vtx_v_lumi->getTH1()->GetBinContent(iLumi));
    if(nvtx !=0)
    {
      std::cout<<nvtx<<endl;
      std::cout<<"lumi: "<<iLumi<<endl;
      //loop through the layers
      for(int iLayer = 1; iLayer<8; iLayer++)
      {
        //get the eff at the lumisection and layer
        eff = eff_v_lumi_forward->getTProfile2D()->GetBinContent(iLumi-1,iLayer);
	std::cout<<"disk eff: "<<eff<<endl;
        //set the efficiency in the new histo
        eff_v_vtx_forward->getTH2F()->SetBinContent(nvtx+1, iLayer, eff);
      }

      //loop through the layers
      for(int iLayer = 1; iLayer<5; iLayer++)
      {
        //get the efficiency for each lumi at each layer
        eff = eff_v_lumi_barrel->getTProfile2D()->GetBinContent(iLumi-1, iLayer);
	std::cout<<"barrel eff: "<<eff<<endl;
        //set the efficiency
        eff_v_vtx_barrel->getTH2F()->SetBinContent(nvtx+1, iLayer, eff);
      }
    }
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelPhase1Extras);
