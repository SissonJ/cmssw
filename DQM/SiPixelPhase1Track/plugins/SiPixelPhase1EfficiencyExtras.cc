// -*- C++ -*-
//
// Package:    SiPixelPhase1EfficiencyExtras
// Class:      SiPixelPhase1EfficiencyExtras
//
/**\class 

 Description: Create the Phase 1 extra efficiency trend plots

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jack Sisson, Julie Hogan
//         Created:  7 July, 2021
//
//

// Framework
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DQM Framework
#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"
// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
// DataFormats
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelNameUpgrade.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapNameUpgrade.h"
//
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace edm;

class SiPixelPhase1EfficiencyExtras : public DQMEDHarvester {
public:
  explicit SiPixelPhase1EfficiencyExtras(const edm::ParameterSet& conf);
  ~SiPixelPhase1EfficiencyExtras() override;

  //       virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //         //void dqmBeginRun(const edm::Run&, edm::EventSetup const&) ;
  //           //virtual void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
protected:
  void beginRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

  void dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) override;

  std::string effFolderName_;
  std::string vtxFolderName_;
  std::string instLumiFolderName_;

private:
  edm::ParameterSet conf_;
};

SiPixelPhase1EfficiencyExtras::SiPixelPhase1EfficiencyExtras(const edm::ParameterSet& iConfig) : conf_(iConfig) {
  LogInfo("PixelDQM") << "SiPixelPhase1EfficiencyExtras::SiPixelPhase1EfficiencyExtras: Hello!" << endl;
  effFolderName_ = conf_.getParameter<std::string>("EffFolderName");
  vtxFolderName_ = conf_.getParameter<std::string>("VtxFolderName");
  instLumiFolderName_ = conf_.getParameter<std::string>("InstLumiFolderName");
}

SiPixelPhase1EfficiencyExtras::~SiPixelPhase1EfficiencyExtras() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  LogInfo("PixelDQM") << "SiPixelPhase1EfficiencyExtras::~SiPixelPhase1EfficiencyExtras: Destructor" << endl;
}

void SiPixelPhase1EfficiencyExtras::beginRun(edm::Run const& run, edm::EventSetup const& eSetup) {}

//------------------------------------------------------------------
// Method called for every event
//------------------------------------------------------------------
void SiPixelPhase1EfficiencyExtras::dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) {
  iBooker.setCurrentFolder(effFolderName_);

  //Get the existing histos
  MonitorElement* vtx_v_lumi = iGetter.get(vtxFolderName_ + "/NumberOfGoodPVtxVsLS_GenTk");

  MonitorElement* scalLumi_v_lumi = iGetter.get(instLumiFolderName_ + "/lumiVsLS");

  MonitorElement* eff_v_lumi_forward =
      iGetter.get(effFolderName_ + "/hitefficiency_per_Lumisection_per_PXDisk_PXForward");

  MonitorElement* eff_v_lumi_barrel =
      iGetter.get(effFolderName_ + "/hitefficiency_per_Lumisection_per_PXLayer_PXBarrel");

  if (!vtx_v_lumi) {
    edm::LogWarning("SiPixelPhase1EfficiencyExtras")
        << "no NumberOfGoodPVtxVsLS_GenTK ME is available in " << vtxFolderName_ << std::endl;
    return;
  } else if (!scalLumi_v_lumi) {
    edm::LogWarning("SiPixelPhase1EfficiencyExtras")
        << "no lumiVsLS ME is available in " << instLumiFolderName_ << std::endl;
    return;
  }

  //Get the max value of inst lumi for plot
  int yMax2 = scalLumi_v_lumi->getTProfile()->GetMaximum();
  yMax2 = yMax2 + yMax2 * .1;

  //Book new histos
  MonitorElement* eff_v_vtx_barrel =
      iBooker.book2D("hitefficiency_per_meanNvtx_per_PXLayer_PXBarrel",
                     "hitefficiency_per_meanNvtx_per_PXLayer_PXBarrel; meanNvtx; PXLayer",
                     500,
                     0,
                     100,
                     3,
                     .5,
                     3.5);

  MonitorElement* eff_v_vtx_forward =
      iBooker.book2D("hitefficiency_per_meanNvtx_per_PXDisk_PXForward",
                     "hitefficiency_per_meanNvtx_per_PXDisk_PXForward; meanNvtx; PXDisk",
                     500,
                     0,
                     100,
                     7,
                     -3.5,
                     3.5);

  MonitorElement* eff_v_scalLumi_barrel =
      iBooker.book2D("hitefficiency_per_scalLumi_per_PXLayer_PXBarrel",
                     "hitefficiency_per_scalLumi_per_PXLayer_PXBarrel; scal inst lumi E30; PXLayer",
                     500,
                     0,
                     yMax2,
                     3,
                     .5,
                     3.5);

  MonitorElement* eff_v_scalLumi_forward =
      iBooker.book2D("hitefficiency_per_scalLumi_per_PXDisk_PXForward",
                     "hitefficiency_per_scalLumi_per_PXDisk_PXForward; scal inst lumi E30; PXDisk",
                     500,
                     0,
                     yMax2,
                     7,
                     -3.5,
                     3.5);

  //initialize variables
  int numLumiNvtx = int(vtx_v_lumi->getTH1()->GetNbinsX());
  int numLumiScal = int(scalLumi_v_lumi->getTProfile()->GetNbinsX());
  double nvtx = 0.0;
  double scalLumi = 0.0;
  double eff = 0.0;
  int binNumVtx = 0;
  int binNumScal = 0;

  //For loop to loop through lumisections
  for (int iLumi = 1; iLumi < numLumiNvtx - 1; iLumi++) {
    //get the meanNvtx and inst lumi for each lumi
    nvtx = vtx_v_lumi->getTH1()->GetBinContent(iLumi);
    scalLumi = scalLumi_v_lumi->getTProfile()->GetBinContent(iLumi);

    //Filter out useless iterations
    if (nvtx != 0 || scalLumi != 0) {
      //Grab the bin number for the nvtx and inst lumi
      binNumVtx = eff_v_vtx_barrel->getTH2F()->FindBin(nvtx);
      binNumScal = eff_v_scalLumi_barrel->getTH2F()->FindBin(scalLumi);

      //loop through the layers
      for (int iLayer = 1; iLayer < 8; iLayer++) {
        //get the eff at the lumisection and layer
        eff = eff_v_lumi_forward->getTProfile2D()->GetBinContent(iLumi - 1, iLayer);

        //set the efficiency in the new histo
        eff_v_vtx_forward->getTH2F()->SetBinContent(binNumVtx, iLayer, eff);

        //Filter iLumi to be smaller than the max x value of inst lumi
        if (iLumi <= numLumiScal) {
          //set the efficiency in the new histo
          eff_v_scalLumi_forward->getTH2F()->SetBinContent(binNumScal, iLayer, eff);
        }
      }

      //loop through the layers
      for (int iLayer = 1; iLayer < 5; iLayer++) {
        //get the efficiency for each lumi at each layer
        eff = eff_v_lumi_barrel->getTProfile2D()->GetBinContent(iLumi - 1, iLayer);

        //set the efficiency
        eff_v_vtx_barrel->getTH2F()->SetBinContent(binNumVtx, iLayer, eff);

        //Filter iLumi to be smaller than the max x value of inst lumi
        if (iLumi <= numLumiScal) {
          //set the efficiency in the new histo
          eff_v_scalLumi_barrel->getTH2F()->SetBinContent(binNumScal, iLayer, eff);
        }
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelPhase1EfficiencyExtras);
