#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/SiStripCluster/interface/SiStripApproximateCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetInfo.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDetInfoFileReader.h"
#include "RecoTracker/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include <vector>
#include <memory>

class SiStripClusters2ApproxClusters : public edm::stream::EDProducer<> {
public:
  explicit SiStripClusters2ApproxClusters(const edm::ParameterSet& conf);
  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  edm::InputTag inputClusters;
  edm::InputTag beamSpot; // member variable for BeamSpotTag
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > clusterToken;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken; // token for BeamSpot
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::FileInPath fileInPath;
  SiStripDetInfo detInfo;

  edm::ESGetToken<ClusterShapeHitFilter, CkfComponentsRecord> csfToken_;

  unsigned int maxNSat;
};

SiStripClusters2ApproxClusters::SiStripClusters2ApproxClusters(const edm::ParameterSet& conf) {
  inputClusters = conf.getParameter<edm::InputTag>("inputClusters");

  beamSpot = conf.getParameter<edm::InputTag>("beamSpot"); // initialising the new member variable

  maxNSat = conf.getParameter<unsigned int>("maxSaturatedStrips");

  clusterToken = consumes<edmNew::DetSetVector<SiStripCluster> >(inputClusters);
  beamSpotToken = consumes<reco::BeamSpot>(beamSpot); // initialising beamSpot token
  
  tkGeomToken_ = esConsumes(); 
  csfToken_ = esConsumes(edm::ESInputTag("", "ClusterShapeHitFilter"));

  fileInPath = edm::FileInPath(SiStripDetInfoFileReader::kDefaultFile);
  detInfo = SiStripDetInfoFileReader::read(fileInPath.fullPath());

  produces<edmNew::DetSetVector<SiStripApproximateCluster> >();
}

void SiStripClusters2ApproxClusters::produce(edm::Event& event, edm::EventSetup const& setup) {
  auto result = std::make_unique<edmNew::DetSetVector<SiStripApproximateCluster> >();
  const auto& clusterCollection = event.get(clusterToken);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByToken(beamSpotToken, beamSpotHandle); // retrive BeamSpot data
  reco::BeamSpot const* bs = nullptr;
  if (beamSpotHandle.isValid())
    bs = &(*beamSpotHandle);

  //const auto& theFilter = &setup.getData(csfToken_);
  const auto& tkGeom = &setup.getData(tkGeomToken_);

  unsigned short Nstrips;
  double stripLength;
  //float hitPredPos;
  //int hitStrips;

  for (const auto& detClusters : clusterCollection) {
    edmNew::DetSetVector<SiStripApproximateCluster>::FastFiller ff{*result, detClusters.id()};

    for (const auto& cluster : detClusters){
      const GeomDet* det = tkGeom->idToDet(detClusters.id());

      GlobalPoint beamspot(bs->position().x(), bs->position().y(), bs->position().z());
      const StripGeomDetUnit* stripDet = (const StripGeomDetUnit*)(&det);
      Nstrips = detInfo.getNumberOfApvsAndStripLength(detClusters.id()).first * 128;
      stripLength = detInfo.getNumberOfApvsAndStripLength(detClusters.id()).second;
      double y = cluster.barycenter() * stripLength / (Nstrips * 128.0);
      const LocalPoint& lp = det->surface().toLocal(GlobalPoint(cluster.barycenter(),y,stripDet->surface().position().z()));
      const GlobalPoint& gpos = det->surface().toGlobal(lp);
      std::cout << gpos << std::endl;

      //GlobalVector gdir = beamspot - gpos;
      //LocalVector ldir = det->toLocal(gdir);
      //LocalPoint lpos = det->toLocal(gpos);
      //std::cout << gdir << ldir << lpos << std::endl;

      //bool usable = theFilter->getSizes(detClusters.id(), cluster, lpos, ldir, hitStrips, hitPredPos);
      //std::cout << usable << std::endl;

      //ff.push_back(SiStripApproximateCluster(cluster, maxNSat, bs));
      ff.push_back(SiStripApproximateCluster(cluster,maxNSat));
    }
  }

  event.put(std::move(result));
}

void SiStripClusters2ApproxClusters::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputClusters", edm::InputTag("siStripClusters"));
  desc.add<unsigned int>("maxSaturatedStrips", 3);
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot")); // add BeamSpot tag
  descriptions.add("SiStripClusters2ApproxClusters", desc);
}

DEFINE_FWK_MODULE(SiStripClusters2ApproxClusters);

