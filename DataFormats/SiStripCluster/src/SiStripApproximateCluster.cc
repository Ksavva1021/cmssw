#include "DataFormats/SiStripCluster/interface/SiStripApproximateCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"

#include <algorithm>
#include <cmath>

//SiStripApproximateCluster::SiStripApproximateCluster(const SiStripCluster& cluster, unsigned int maxNSat,const reco::BeamSpot* bs) {
SiStripApproximateCluster::SiStripApproximateCluster(const SiStripCluster& cluster, unsigned int maxNSat, unsigned int detId, float hitPredPos, float mipnorm, const SiStripNoises* theNoise) {
  barycenter_ = std::round(cluster.barycenter() * 10);
  width_ = cluster.size();
  avgCharge_ = cluster.charge() / cluster.size();
  isSaturated_ = false;
  trimFilter_ = false;
  peakFilter_ = false;

  //mimicing the algorithm used in StripSubClusterShapeTrajectoryFilter...
  //Looks for 3 adjacent saturated strips (ADC>=254)
  const auto& ampls = cluster.amplitudes();
  unsigned int thisSat = (ampls[0] >= 254), maxSat = thisSat;
  for (unsigned int i = 1, n = ampls.size(); i < n; ++i) {
    if (ampls[i] >= 254) {
      thisSat++;
    } else if (thisSat > 0) {
      maxSat = std::max<int>(maxSat, thisSat);
      thisSat = 0;
    }
  }
  if (thisSat > 0) {
    maxSat = std::max<int>(maxSat, thisSat);
  }
  if (maxSat >= maxNSat) {
    isSaturated_ = true;
  }
 
  //trimming filter
  double trimMaxADC_ = 30.;
  double trimMaxFracTotal_ = .15;
  double trimMaxFracNeigh_ = .25;
  double maxTrimmedSizeDiffNeg_ = .7;
  double maxTrimmedSizeDiffPos_ = 1.;
  double subclusterWindow_ = .7;
  double seedCutMIPs_ = .35;
  double seedCutSN_ = 7.;
  double subclusterCutMIPs_ = .45;
  double subclusterCutSN_ = 12.;

  unsigned int hitStripsTrim = ampls.size();
  int sum = std::accumulate(ampls.begin(), ampls.end(), 0);
  uint8_t trimCut = std::min<uint8_t>(trimMaxADC_, std::floor(trimMaxFracTotal_ * sum));
  auto begin = ampls.begin();
  auto last = ampls.end() - 1;
  while (hitStripsTrim > 1 && (*begin < std::max<uint8_t>(trimCut, trimMaxFracNeigh_ * (*(begin + 1))))) {
    hitStripsTrim--;
    ++begin;
  }
  while (hitStripsTrim > 1 && (*last < std::max<uint8_t>(trimCut, trimMaxFracNeigh_ * (*(last - 1))))) {
    hitStripsTrim--;
    --last;
  }
  if (hitStripsTrim < std::floor(std::abs(hitPredPos) - maxTrimmedSizeDiffNeg_)) {
    trimFilter_ =  false;
  } else if (hitStripsTrim <= std::ceil(std::abs(hitPredPos) + maxTrimmedSizeDiffPos_)) {
    trimFilter_ = true;
  }
  
  //peakFinder
  SlidingPeakFinder pf(std::max<int>(2, std::ceil(std::abs(hitPredPos) + subclusterWindow_))); 
  //auto range = theNoise->getRange(detId);
  //std::cout << "[" << *range.first << ", " << *range.second << "]" << std::endl;
  std::cout << theNoise->getRange(detId);

  //PeakFinderTest test(mipnorm,detId,cluster.firstStrip(),theNoise,seedCutMIPs_,seedCutSN_,subclusterCutMIPs_,subclusterCutSN_);
  //if (pf.apply(cluster.amplitudes(), test)) {
  //  peakFilter_ = true;
  //} else {
  //  peakFilter_ = false;
  //}
}
