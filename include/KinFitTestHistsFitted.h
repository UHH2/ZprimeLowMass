#pragma once

#include "UHH2/core/include/Hists.h"
#include <UHH2/KinFitTest/include/HistsBASE.h>

namespace uhh2examples {

class KinFitTestHistsFitted: public HistsBASE {
public:
    KinFitTestHistsFitted(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;
    virtual ~KinFitTestHistsFitted();

  uhh2::Event::Handle<TLorentzVector> rec_fit_bhad_v4;
  uhh2::Event::Handle<TLorentzVector> rec_fit_blep_v4;
  uhh2::Event::Handle<TLorentzVector> rec_fit_lepton_v4;
  uhh2::Event::Handle<TLorentzVector> rec_fit_neutrino_v4;
  uhh2::Event::Handle<TLorentzVector> rec_fit_jetq1_v4;
  uhh2::Event::Handle<TLorentzVector> rec_fit_jetq2_v4;

 protected:
    virtual void init(); 

};

}
