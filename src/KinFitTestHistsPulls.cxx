#include "UHH2/KinFitTest/include/KinFitTestHistsPulls.h"
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <iostream>


using namespace std;
using namespace uhh2;
using namespace uhh2examples;


KinFitTestHistsPulls::KinFitTestHistsPulls(uhh2::Context& ctx, const std::string& dir, const std::string& ttgen, const std::string& hyps, const std::string& disc):
  HistsBASE(ctx, dir) {

  h_ttbar_gen  = ctx.get_handle<TTbarGen>(ttgen);
  h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis> >(hyps);
  disc_name_ = disc;

  pull_had1 = ctx.get_handle<TMatrixD>("pull_had1");
  pull_had2 = ctx.get_handle<TMatrixD>("pull_had2");
  pull_hadb = ctx.get_handle<TMatrixD>("pull_hadb");
  pull_lepb = ctx.get_handle<TMatrixD>("pull_lepb");
  pull_lepton = ctx.get_handle<TMatrixD>("pull_lepton");
  pull_neu = ctx.get_handle<TMatrixD>("pull_neu");
  pull_gtr_had1_et = ctx.get_handle<float>("pull_gtr_had1_et");
  pull_gtr_had1_eta = ctx.get_handle<float>("pull_gtr_had1_eta");
  pull_gtr_had1_phi = ctx.get_handle<float>("pull_gtr_had1_phi");
  pull_gtr_had2_et = ctx.get_handle<float>("pull_gtr_had2_et");
  pull_gtr_had2_eta = ctx.get_handle<float>("pull_gtr_had2_eta");
  pull_gtr_had2_phi = ctx.get_handle<float>("pull_gtr_had2_phi");
  pull_gtr_hadb_et = ctx.get_handle<float>("pull_gtr_hadb_et");
  pull_gtr_hadb_eta = ctx.get_handle<float>("pull_gtr_hadb_eta");
  pull_gtr_hadb_phi = ctx.get_handle<float>("pull_gtr_hadb_phi");
  pull_gtr_lepb_et = ctx.get_handle<float>("pull_gtr_lepb_et");
  pull_gtr_lepb_eta = ctx.get_handle<float>("pull_gtr_lepb_eta");
  pull_gtr_lepb_phi = ctx.get_handle<float>("pull_gtr_lepb_phi");
  pull_gtr_lepton_et = ctx.get_handle<float>("pull_gtr_lepton_et");
  pull_gtr_lepton_eta = ctx.get_handle<float>("pull_gtr_lepton_eta");
  pull_gtr_lepton_phi = ctx.get_handle<float>("pull_gtr_lepton_phi");
  pull_gtr_neu_et = ctx.get_handle<float>("pull_gtr_neu_et");
  pull_gtr_neu_eta = ctx.get_handle<float>("pull_gtr_neu_eta");
  pull_gtr_neu_phi = ctx.get_handle<float>("pull_gtr_neu_phi");

  init();
}


void KinFitTestHistsPulls::init(){

  book_TH1F("pull_gm_had1_a"              , 100, -3., 3.);
  book_TH1F("pull_gm_had1_b"              , 100, -3., 3.);
  book_TH1F("pull_gm_had1_c"              , 100, -3., 3.);
  book_TH1F("pull_gm_had1_d"              , 100, -3., 3.);
  book_TH1F("pull_gm_had2_a"              , 100, -3., 3.);
  book_TH1F("pull_gm_had2_b"              , 100, -3., 3.);
  book_TH1F("pull_gm_had2_c"              , 100, -3., 3.);
  book_TH1F("pull_gm_had2_d"              , 100, -3., 3.);
  book_TH1F("pull_gm_hadb_a"              , 100, -3., 3.);
  book_TH1F("pull_gm_hadb_b"              , 100, -3., 3.);
  book_TH1F("pull_gm_hadb_c"              , 100, -3., 3.);
  book_TH1F("pull_gm_hadb_d"              , 100, -3., 3.);
  book_TH1F("pull_gm_lepb_a"              , 100, -3., 3.);
  book_TH1F("pull_gm_lepb_b"              , 100, -3., 3.);
  book_TH1F("pull_gm_lepb_c"              , 100, -3., 3.);
  book_TH1F("pull_gm_lepb_d"              , 100, -3., 3.);
  book_TH1F("pull_gm_lepton_a"            , 100, -3., 3.);
  book_TH1F("pull_gm_lepton_b"            , 100, -3., 3.);
  book_TH1F("pull_gm_lepton_c"            , 100, -3., 3.);
  book_TH1F("pull_gm_neu_a"               , 100, -3., 3.);
  book_TH1F("pull_gm_neu_b"               , 100, -3., 3.);
  book_TH1F("pull_gm_neu_c"               , 100, -3., 3.);
  book_TH1F("pull_gtr_had1_a"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had1_b"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had1_c"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had1_d"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had2_a"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had2_b"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had2_c"              , 100, -3., 3.);
  book_TH1F("pull_gtr_had2_d"              , 100, -3., 3.);
  book_TH1F("pull_gtr_hadb_a"              , 100, -3., 3.);
  book_TH1F("pull_gtr_hadb_b"              , 100, -3., 3.);
  book_TH1F("pull_gtr_hadb_c"              , 100, -3., 3.);
  book_TH1F("pull_gtr_hadb_d"              , 100, -3., 3.);
  book_TH1F("pull_gtr_lepb_a"              , 100, -3., 3.);
  book_TH1F("pull_gtr_lepb_b"              , 100, -3., 3.);
  book_TH1F("pull_gtr_lepb_c"              , 100, -3., 3.);
  book_TH1F("pull_gtr_lepb_d"              , 100, -3., 3.);
  book_TH1F("pull_gtr_lepton_a"            , 100, -3., 3.);
  book_TH1F("pull_gtr_lepton_b"            , 100, -3., 3.);
  book_TH1F("pull_gtr_lepton_c"            , 100, -3., 3.);
  book_TH1F("pull_gtr_neu_a"               , 100, -3., 3.);
  book_TH1F("pull_gtr_neu_b"               , 100, -3., 3.);
  book_TH1F("pull_gtr_neu_c"               , 100, -3., 3.);

}


void KinFitTestHistsPulls::fill(const uhh2::Event& event){

  const float weight(event.weight);
  const TTbarGen* ttgen(0);


  if(!event.isRealData){
    const auto& ttbargen = event.get(h_ttbar_gen);
    ttgen = &ttbargen;
   }

 const ReconstructionHypothesis* hyp(0);
  float hyp_val(0.);
  if(event.is_valid(h_ttbar_hyps)){

    const auto& ttbar_hyps = event.get(h_ttbar_hyps);

    hyp = get_best_hypothesis(ttbar_hyps, disc_name_);
    if(!hyp) throw std::runtime_error("EffyTTbarHists::fill -- null pointer to ReconstructionHypothesis object (from \"get_best_hypothesis\")");

    hyp_val = hyp->discriminator(disc_name_);
  }



  fill(event, hyp, hyp_val, ttgen, weight);



  return;
}

void  KinFitTestHistsPulls::fill(const uhh2::Event& event, const ReconstructionHypothesis* hyp, const float hyp_val, const TTbarGen* ttgen, const double weight){

   bool ttljets(false);

 if(ttgen){
   ttljets = (ttgen->DecayChannel() == TTbarGen::e_muhad || ttgen->DecayChannel() == TTbarGen::e_ehad);
  }
    int param = event.get(parametrisation_);
    TMatrixD pull_had1_ = event.get(pull_had1);
    TMatrixD pull_had2_ = event.get(pull_had2);
    TMatrixD pull_hadb_ = event.get(pull_hadb);
    TMatrixD pull_lepb_ = event.get(pull_lepb);
    TMatrixD pull_lepton_ = event.get(pull_lepton);
    TMatrixD pull_neu_ = event.get(pull_neu);
    float pull_gtr_had1_et_ = event.get(pull_gtr_had1_et);
    float pull_gtr_had1_eta_ = event.get(pull_gtr_had1_eta);
    float pull_gtr_had1_phi_ = event.get(pull_gtr_had1_phi);
    float pull_gtr_had2_et_ = event.get(pull_gtr_had2_et);
    float pull_gtr_had2_eta_ = event.get(pull_gtr_had2_eta);
    float pull_gtr_had2_phi_ = event.get(pull_gtr_had2_phi);
    float pull_gtr_hadb_et_ = event.get(pull_gtr_hadb_et);
    float pull_gtr_hadb_eta_ = event.get(pull_gtr_hadb_eta);
    float pull_gtr_hadb_phi_ = event.get(pull_gtr_hadb_phi);
    float pull_gtr_lepb_et_ = event.get(pull_gtr_lepb_et);
    float  pull_gtr_lepb_eta_ = event.get(pull_gtr_lepb_eta);
    float pull_gtr_lepb_phi_ = event.get(pull_gtr_lepb_phi);
    float pull_gtr_lepton_et_ = event.get(pull_gtr_lepton_et);
    float pull_gtr_lepton_eta_ = event.get(pull_gtr_lepton_eta);
    float pull_gtr_lepton_phi_ = event.get(pull_gtr_lepton_phi);
    float pull_gtr_neu_et_ = event.get(pull_gtr_neu_et);
    float pull_gtr_neu_eta_ = event.get(pull_gtr_neu_eta);
    float pull_gtr_neu_phi_ = event.get(pull_gtr_neu_phi);
    H1("pull_gm_had1_a")->Fill(pull_had1_(0,0));
    H1("pull_gm_had1_b")->Fill(pull_had1_(1,0));
    H1("pull_gm_had1_c")->Fill(pull_had1_(2,0));
    if(param==0)    H1("pull_gm_had1_d")->Fill(pull_had1_(3,0));
    H1("pull_gm_had2_a")->Fill(pull_had2_(0,0));
    H1("pull_gm_had2_b")->Fill(pull_had2_(1,0));
    H1("pull_gm_had2_c")->Fill(pull_had2_(2,0));
    if(param==0)    H1("pull_gm_had2_d")->Fill(pull_had2_(3,0));
    H1("pull_gm_hadb_a")->Fill(pull_hadb_(0,0));
    H1("pull_gm_hadb_b")->Fill(pull_hadb_(1,0));
    H1("pull_gm_hadb_c")->Fill(pull_hadb_(2,0));
    if(param==0)    H1("pull_gm_hadb_d")->Fill(pull_hadb_(3,0));
    H1("pull_gm_lepb_a")->Fill(pull_lepb_(0,0));
    H1("pull_gm_lepb_b")->Fill(pull_lepb_(1,0));
    H1("pull_gm_lepb_c")->Fill(pull_lepb_(2,0));
    if(param==0)    H1("pull_gm_lepb_d")->Fill(pull_lepb_(3,0));
    H1("pull_gm_lepton_a")->Fill(pull_lepton_(0,0));
    H1("pull_gm_lepton_b")->Fill(pull_lepton_(1,0));
    H1("pull_gm_lepton_c")->Fill(pull_lepton_(2,0));
    H1("pull_gm_neu_a")->Fill(pull_neu_(0,0));
    H1("pull_gm_neu_b")->Fill(pull_neu_(1,0));
    H1("pull_gm_neu_c")->Fill(pull_neu_(2,0));
    if(!(pull_gtr_had1_et_==-100.)){
      H1("pull_gtr_had1_a")->Fill(pull_gtr_had1_et_);
      H1("pull_gtr_had1_b")->Fill(pull_gtr_had1_eta_);
      H1("pull_gtr_had1_c")->Fill(pull_gtr_had1_phi_);
      //      if(param==0)    H1("pull_gtr_had1_d")->Fill(pull_gtr_had1_(3,0));
      H1("pull_gtr_had2_a")->Fill(pull_gtr_had2_et_);
      H1("pull_gtr_had2_b")->Fill(pull_gtr_had2_eta_);
      H1("pull_gtr_had2_c")->Fill(pull_gtr_had2_phi_);
      //      if(param==0)    H1("pull_gtr_had2_d")->Fill(pull_gtr_had2_(3,0));
      H1("pull_gtr_hadb_a")->Fill(pull_gtr_hadb_et_);
      H1("pull_gtr_hadb_b")->Fill(pull_gtr_hadb_eta_);
      H1("pull_gtr_hadb_c")->Fill(pull_gtr_hadb_phi_);
      //      if(param==0)    H1("pull_gtr_hadb_d")->Fill(pull_gtr_hadb_(3,0));
      H1("pull_gtr_lepb_a")->Fill(pull_gtr_lepb_et_);
      H1("pull_gtr_lepb_b")->Fill(pull_gtr_lepb_eta_);
      H1("pull_gtr_lepb_c")->Fill(pull_gtr_lepb_phi_);
      //      if(param==0)    H1("pull_gtr_lepb_d")->Fill(pull_gtr_lepb_(3,0));
      H1("pull_gtr_lepton_a")->Fill(pull_gtr_lepton_et_);
      H1("pull_gtr_lepton_b")->Fill(pull_gtr_lepton_eta_);
      H1("pull_gtr_lepton_c")->Fill(pull_gtr_lepton_phi_);
      H1("pull_gtr_neu_a")->Fill(pull_gtr_neu_et_);
      H1("pull_gtr_neu_b")->Fill(pull_gtr_neu_eta_);
      H1("pull_gtr_neu_c")->Fill(pull_gtr_neu_phi_);
     }

}


float KinFitTestHistsPulls::deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2){
    float deltaeta = eta1 - eta2;
    float dphi = delta_phi(phi1, phi2);
    return sqrt(deltaeta * deltaeta + dphi * dphi);
}

float KinFitTestHistsPulls::delta_phi(const float phi1, const float phi2){

  float a_dphi = fabs(phi1-phi2);
  if(a_dphi > M_PI) a_dphi = 2*M_PI - a_dphi;

  return a_dphi;
}



KinFitTestHistsPulls::~KinFitTestHistsPulls(){}

