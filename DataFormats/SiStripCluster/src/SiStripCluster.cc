
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

SiStripCluster::SiStripCluster(const SiStripDigiRange& range) : firstStrip_(range.first->strip()), error_x(-99999.9) {
  std::vector<uint8_t> v;
  v.reserve(range.second - range.first);

  uint16_t lastStrip = 0;
  bool firstInloop = true;
  for (SiStripDigiIter i = range.first; i != range.second; i++) {
    /// check if digis consecutive
    if (!firstInloop && i->strip() != lastStrip + 1) {
      for (int j = 0; j < i->strip() - (lastStrip + 1); j++) {
        v.push_back(0);
      }
    }
    lastStrip = i->strip();
    firstInloop = false;

    v.push_back(i->adc());
  }
  amplitudes_ = v;
}

SiStripCluster::SiStripCluster(SiStripApproximateCluster cluster) : error_x(-99999.9) {
  barycenter_ = cluster.barycenter();
  charge_ = cluster.width()*cluster.avgCharge();
  amplitudes_.resize(cluster.width(), cluster.avgCharge()); //fill amplitudes_ with the average charge of the cluster
  firstStrip_ = cluster.barycenter() - cluster.width()/2; //initialize firstStrip_ to avoid bug
  isSaturated_ = cluster.isSaturated();
  //std::cout << "Saturated: " << cluster.isSaturated() << std::endl;
}

int SiStripCluster::charge() const {
  if ( barycenter_ > 0 ) return charge_;
  return std::accumulate(begin(), end(), int(0));
}

//bool SiStripCluster::isSaturated() const { 
//  if (barycenter_ > 0 ) return isSaturated_;
//  const auto& ampls = amplitudes_;
//  unsigned int thisSat = (ampls[0] >= 254), maxSat = thisSat;
//  for (unsigned int i = 1, n = ampls.size(); i < n; ++i) {
//    if (ampls[i] >= 254) {
//      thisSat++;
//    } else if (thisSat > 0) {
//      maxSat = std::max<int>(maxSat, thisSat);
//      thisSat = 0;
//    }
//  }
//  if (thisSat > 0) {
//    maxSat = std::max<int>(maxSat, thisSat);
//  }
//  if (maxSat >= 3) {
//    return true;
//  }
//  return false;
//}

float SiStripCluster::barycenter() const {
  if ( barycenter_ > 0 ) return barycenter_;

  int sumx = 0;
  int suma = 0;
  auto asize = size();
  for (auto i = 0U; i < asize; ++i) {
    sumx += i * amplitudes_[i];
    suma += amplitudes_[i];
  }

  // strip centers are offcet by half pitch w.r.t. strip numbers,
  // so one has to add 0.5 to get the correct barycenter position.
  // Need to mask off the high bit of firstStrip_, which contains the merged status.
  return float((firstStrip_ & stripIndexMask)) + float(sumx) / float(suma) + 0.5f;
}
