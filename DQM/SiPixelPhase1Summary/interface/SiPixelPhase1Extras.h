#ifndef SiPixelPhase1Extras_SiPixelPhase1Extras_h
#define SiPixelPhase1Extras_SiPixelPhase1Extras_h
// -*- C++ -*-
//
// Package:     SiPixelPhase1Extras
// Class  :     SiPixelPhase1Extras
//
/**

 Description: Extras map generation for the Phase 1 pixel

 Usage:
    <usage>

*/
//
// Original Author:  Duncan Leggat
//         Created:  2nd December 2016
//

//#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class SiPixelPhase1Extras : public DQMEDHarvester {
public:
  explicit SiPixelPhase1Extras(const edm::ParameterSet& conf);
  ~SiPixelPhase1Extras() override;

  //       virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //void dqmBeginRun(const edm::Run&, edm::EventSetup const&) ;
  //virtual void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
protected:
  void beginRun(edm::Run const& run, edm::EventSetup const& eSetup) override;

  void dqmEndJob(DQMStore::IBooker& iBooker, DQMStore::IGetter& iGetter) override;

  std::string effFolderName_;
  std::string vtxFolderName_;

private:
  edm::ParameterSet conf_;

};

#endif
