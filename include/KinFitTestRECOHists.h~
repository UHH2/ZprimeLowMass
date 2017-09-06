#include "TMatrixD.h"

#include <UHH2/core/include/Utils.h>

#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/KinFitTest/include/HistsBASE.h>

#include <TLorentzVector.h>

namespace uhh2examples {


class KinFitTestRECOHists: public HistsBASE {


 public:
  KinFitTestRECOHists(uhh2::Context&, const std::string&, const std::string& ttgen="ttbargen", const std::string& hyps="TTbarReconstruction", const std::string& disc="Chi2");
  virtual ~KinFitTestRECOHists();

  virtual void fill(const uhh2::Event&) override;
  virtual void fill(const uhh2::Event&, const ReconstructionHypothesis* hyp, const float, const TTbarGen* ttgen=0, const double wgt=1.);
  float deltaR_fit(const LorentzVector& p1, const LorentzVector p2);
  float deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2);

  protected:
  uhh2::Event::Handle<TTbarGen> h_ttbar_gen;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  std::string disc_name_;

  virtual void init() override;



  float delta_phi(const float, const float); 

  void boost_x1_to_x2CM(TLorentzVector&, const TLorentzVector&); 
  float cosThetaX (const LorentzVector&, const LorentzVector&, const LorentzVector&); 
  float cosThetaCS(const LorentzVector&, const LorentzVector&); 

};
}
