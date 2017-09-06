#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/KinFitTest/include/KinFitTestSetup.h>

namespace uhh2examples {
    
class KinFitTestFit{
 public:

  void fit(uhh2::Event & event, int bl, int bh, int h1, int h2, KinFitTestSetup::ObjectType leptonType, KinFitTestSetup::Param jetParam_, KinFitTestSetup::Param lepParam_, KinFitTestSetup::Param metParam_);
};
}
