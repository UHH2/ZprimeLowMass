#pragma once

#include "UHH2/core/include/Hists.h"
#include <UHH2/KinFitTest/include/HistsBASE.h>

namespace uhh2examples {

class KinFitTestHists: public HistsBASE {
public:
    KinFitTestHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;
    virtual ~KinFitTestHists();

 protected:
    virtual void init(); 

};

}
