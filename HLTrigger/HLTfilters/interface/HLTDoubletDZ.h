#ifndef HLTDoubletDZ_h
#define HLTDoubletDZ_h

//
// Class imlements |dZ|<Max for a pair of two objects
//

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include<string>
#include<vector>
namespace trigger {
  class TriggerFilterObjectWithRefs;
}

//
// class declaration
//

template<typename T1, typename T2>
class HLTDoubletDZ : public HLTFilter {

   public:

      explicit HLTDoubletDZ(const edm::ParameterSet&);
      ~HLTDoubletDZ();
      static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
      virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct) const;

   private:
      // configuration
      const std::vector<edm::InputTag> originTag1_;  // input tag identifying originals 1st product
      const std::vector<edm::InputTag> originTag2_;  // input tag identifying originals 2nd product
      const edm::InputTag inputTag1_;   // input tag identifying filtered 1st product
      const edm::InputTag inputTag2_;   // input tag identifying filtered 2nd product
      const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> inputToken1_;
      const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> inputToken2_;
      const int triggerType1_;
      const int triggerType2_;
      const double minDR_;              // minimum dR between two objects to be considered a pair
      const double maxDZ_;              // number of pairs passing cuts required
      const int    min_N_;              // number of pairs passing cuts required
      const bool   checkSC_;            // make sure SC constituents are different
      const bool   same_;               // 1st and 2nd product are one and the same

      typedef std::vector<T1> T1Collection;
      typedef edm::Ref<T1Collection> T1Ref;
      typedef std::vector<T2> T2Collection;
      typedef edm::Ref<T2Collection> T2Ref;

};

#endif //HLTDoubletDZ_h
