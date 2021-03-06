#include "TMatrixD.h"

#include <UHH2/core/include/Utils.h>

#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/KinFitTest/include/HistsBASE.h>

#include <TLorentzVector.h>

namespace uhh2examples {


class KinFitTestHistsMuonUnc: public HistsBASE {


 public:
  KinFitTestHistsMuonUnc(uhh2::Context&, const std::string&, const std::string& ttgen="ttbargen", const std::string& hyps="TTbarReconstruction", const std::string& disc="Chi2");
  virtual ~KinFitTestHistsMuonUnc();

  virtual void fill(const uhh2::Event&) override;
  virtual void fill(const uhh2::Event&, const ReconstructionHypothesis* hyp, const float, const TTbarGen* ttgen=0, const double wgt=1.);
  float deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2);

  protected:
  uhh2::Event::Handle<TTbarGen> h_ttbar_gen;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  std::string disc_name_;

  virtual void init() override;

  uhh2::Event::Handle<TLorentzVector> input_jetq1_v4;
  uhh2::Event::Handle<TLorentzVector> input_jetq2_v4;
  uhh2::Event::Handle<TLorentzVector> input_jetbh_v4;
  uhh2::Event::Handle<TLorentzVector> input_jetbl_v4;
  uhh2::Event::Handle<TLorentzVector> input_lepton_v4;
  uhh2::Event::Handle<TLorentzVector> input_neutrino_v4;

  uhh2::Event::Handle<TMatrixD> unc_had1;
  uhh2::Event::Handle<TMatrixD> unc_had2;
  uhh2::Event::Handle<TMatrixD> unc_hadb;
  uhh2::Event::Handle<TMatrixD> unc_lepb;
  uhh2::Event::Handle<TMatrixD> unc_lepton;
  uhh2::Event::Handle<TMatrixD> unc_neu;


  float delta_phi(const float, const float); 

};
}
