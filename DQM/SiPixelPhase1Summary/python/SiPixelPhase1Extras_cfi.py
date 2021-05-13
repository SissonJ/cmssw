import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester


SiPixelPhase1Extras = DQMEDHarvester("SiPixelPhase1Extras",
                                     EffFolderName = cms.string('PixelPhase1/Tracks/'),
                                     VtxFolderName = cms.string('Tracking/TrackParameters/generalTracks/GeneralProperties/'),
)
