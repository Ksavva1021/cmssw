#ifndef DataFormats_SiStripCluster_SiStripApproximateCluster_h
#define DataFormats_SiStripCluster_SiStripApproximateCluster_h

#include "FWCore/Utilities/interface/typedefs.h"
#include "RecoTracker/PixelLowPtUtilities/interface/SlidingPeakFinder.h"

class SiStripCluster;
class SiStripApproximateCluster {
public:
  SiStripApproximateCluster() {}

  explicit SiStripApproximateCluster(cms_uint16_t barycenter,
                                     cms_uint8_t width,
                                     cms_uint8_t avgCharge,
				     bool filter,
                                     bool isSaturated) {
                                     //bool trimFilter,
                                     //bool peakFilter) {
    barycenter_ = barycenter;
    width_ = width;
    avgCharge_ = avgCharge;
    filter_ = filter;
    isSaturated_ = isSaturated;
    //trimFilter_ = trimFilter;
    //peakFilter_ = peakFilter;
  }

  explicit SiStripApproximateCluster(const SiStripCluster& cluster, unsigned int maxNSat, float hitPredPos, bool peakFilter);

  cms_uint16_t barycenter() const { return barycenter_; }
  cms_uint8_t width() const { return width_; }
  cms_uint8_t avgCharge() const { return avgCharge_; }
  bool filter() const { return filter_; }
  bool isSaturated() const { return isSaturated_; }
  //bool trimFilter() const { return trimFilter_; }
  bool peakFilter() const { return peakFilter_; }

private:
  cms_uint16_t barycenter_ = 0;
  cms_uint8_t width_ = 0;
  cms_uint8_t avgCharge_ = 0;
  bool filter_ = false;
  bool isSaturated_ = false;
  //bool trimFilter_ = false;
  bool peakFilter_ = false;
};
#endif  // DataFormats_SiStripCluster_SiStripApproximateCluster_h
