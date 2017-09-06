#pragma once

#include <UHH2/core/include/fwd.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Utils.h>

#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>

#include <UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h>

#include <string>
#include <vector>
#include <unordered_map>

namespace uhh2examples {
    
/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */
  class DijetSelection: public uhh2::Selection {
  public:
    DijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float dphi_min, third_frac_max;
};

  //METCut
  class METCut : public Selection {

   public:
    explicit METCut(float, float max_met=infinity);
    virtual bool passes(const Event&) override;

   private:
    float min_met_, max_met_;
  };


}
