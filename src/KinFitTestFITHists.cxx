#include "UHH2/KinFitTest/include/KinFitTestFITHists.h"
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <iostream>


using namespace std;
using namespace uhh2;
using namespace uhh2examples;


KinFitTestFITHists::KinFitTestFITHists(uhh2::Context& ctx, const std::string& dir, const std::string& ttgen, const std::string& hyps, const std::string& disc):
  HistsBASE(ctx, dir) {

  h_ttbar_gen  = ctx.get_handle<TTbarGen>(ttgen);
  h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis> >(hyps);
  disc_name_ = disc;
  parametrisation_ = ctx.get_handle<int>("parametrisation");
  rec_fit_chi2_ = ctx.get_handle<float>("rec_fit_chi2");
  fit_probab = ctx.get_handle<float>("fit_probab");
  rec_fit_status_ = ctx.get_handle<float>("rec_fit_status");
  rec_fit_fit_found_ = ctx.get_handle<bool>("rec_fit_fit_found");
  rec_fit_is_input_ = ctx.get_handle<bool>("rec_fit_is_input");
  rec_fit_whad_v4 = ctx.get_handle<TLorentzVector>("rec_fit_whad_v4");
  rec_fit_wlep_v4 = ctx.get_handle<TLorentzVector>("rec_fit_wlep_v4");
  rec_fit_thad_v4 = ctx.get_handle<TLorentzVector>("rec_fit_thad_v4");
  rec_fit_tlep_v4 = ctx.get_handle<TLorentzVector>("rec_fit_tlep_v4");
  rec_fit_bhad_v4 = ctx.get_handle<TLorentzVector>("rec_fit_bhad_v4");
  rec_fit_blep_v4 = ctx.get_handle<TLorentzVector>("rec_fit_blep_v4");
  rec_fit_lepton_v4 = ctx.get_handle<TLorentzVector>("rec_fit_lepton_v4");
  rec_fit_neutrino_v4 = ctx.get_handle<TLorentzVector>("rec_fit_neutrino_v4");
  rec_fit_jetq1_v4 = ctx.get_handle<TLorentzVector>("rec_fit_jetq1_v4");
  rec_fit_jetq2_v4 = ctx.get_handle<TLorentzVector>("rec_fit_jetq2_v4");
  rec_fit_dR_combined = ctx.get_handle<float>("rec_fit_dR_combined");
  input_jetq1_v4 = ctx.get_handle<TLorentzVector>("input_jetq1_v4");
  input_jetq2_v4 = ctx.get_handle<TLorentzVector>("input_jetq2_v4");
  input_jetbh_v4 = ctx.get_handle<TLorentzVector>("input_jetbh_v4");
  input_jetbl_v4 = ctx.get_handle<TLorentzVector>("input_jetbl_v4");
  input_lepton_v4 = ctx.get_handle<TLorentzVector>("input_lepton_v4");
  input_neutrino_v4 = ctx.get_handle<TLorentzVector>("input_neutrino_v4");
  jet_permutation_before_dm_cut = ctx.get_handle<int>("jet_permutation_before_dm_cut");
  jet_permutation_after_dm_cut = ctx.get_handle<int>("jet_permutation_after_dm_cut");
  jet_combi_1 = ctx.get_handle<int>("jet_combi_1");
  jet_combi_2 = ctx.get_handle<int>("jet_combi_2");
  jet_combi_3 = ctx.get_handle<int>("jet_combi_3");
  jet_combi_4 = ctx.get_handle<int>("jet_combi_4");

  unc_had1 = ctx.get_handle<TMatrixD>("unc_had1");
  unc_had2 = ctx.get_handle<TMatrixD>("unc_had2");
  unc_hadb = ctx.get_handle<TMatrixD>("unc_hadb");
  unc_lepb = ctx.get_handle<TMatrixD>("unc_lepb");
  unc_lepton = ctx.get_handle<TMatrixD>("unc_lepton");
  unc_neu = ctx.get_handle<TMatrixD>("unc_neu");
  init();
}


void KinFitTestFITHists::init(){

  book_TH1F("gen_ttbar__M" , 5000, 0, 5000);
  book_TH1F("gen_ttbar__pt", 100, 0, 1000);
  book_TH1F("gen_top__M"   , 200, 0, 400);
  book_TH1F("gen_thad__M" , 250, 50, 300);
  book_TH1F("gen_thad__pt", 80, 0, 800);
  book_TH1F("gen_tlep__M" , 250, 50, 300);
  book_TH1F("gen_tlep__pt", 80, 0, 800);
  book_TH1F("gen_tlep__px", 120, -1200, 1200);
  book_TH1F("gen_tlep__py", 120, -1200, 1200);
  book_TH1F("gen_tlep__pz", 120, -1200, 1200);
  book_TH1F("gen_tops__DR" , 60, 0, 6);
  book_TH1F("gen_tops__Dpt", 180, -450, 450);
  book_TH1F("gen_Wlep__M" , 90, 60, 150);
  book_TH1F("gen_Wlep__Mt", 180, 0, 180);
  book_TH1F("gen_Wlep__pt", 80, 0, 800);
  book_TH1F("gen_Wlep__px", 120, -1200, 1200);
  book_TH1F("gen_Wlep__py", 120, -1200, 1200);
  book_TH1F("gen_Wlep__pz", 120, -1200, 1200);
  book_TH1F("gen_Whad__M" , 90, 60, 150);
  book_TH1F("gen_Whad__pt", 80, 0, 800);
  book_TH1F("gen_Whad__px", 120, -1200, 1200);
  book_TH1F("gen_Whad__py", 120, -1200, 1200);
  book_TH1F("gen_Whad__pz", 120, -1200, 1200);
  book_TH1F("gen_lep__pt"       , 80, 0, 800);
  book_TH1F("gen_lep__eta"      , 60, -6, 6);
  book_TH1F("gen_lep__phi"      , 60, -3.15, 3.15);
  book_TH1F("gen_lep__cosThetaX", 40, -1, 1);
  book_TH1F("gen_neu__pt"       , 80, 0, 800);
  book_TH1F("gen_neu__phi"      , 60, -3.15, 3.15);
  book_TH1F("gen_neu__px"       , 120, -600, 600);
  book_TH1F("gen_neu__py"       , 120, -600, 600);
  book_TH1F("gen_neu__pz"       , 120, -600, 600);
  book_TH1F("gen_neu__cosThetaX", 40, -1, 1);
  book_TH1F("gen_blep__pt"       , 80, 0, 800);
  book_TH1F("gen_blep__eta"      , 60, -6, 6);
  book_TH1F("gen_blep__phi"      , 60, -3.15, 3.15);
  book_TH1F("gen_blep__px"       , 120, -1200, 1200);
  book_TH1F("gen_blep__py"       , 120, -1200, 1200);
  book_TH1F("gen_blep__pz"       , 120, -1200, 1200);
  book_TH1F("gen_blep__cosThetaX", 40, -1, 1);
  book_TH1F("gen_jet_l1__pt"            , 120, 0, 1200);
  book_TH1F("gen_jet_l1__eta"           , 60, -6, 6);
  book_TH1F("gen_jet_l1__phi"           , 60, -3.15, 3.15);
  book_TH1F("gen_jet_l2__pt"            , 120, 0, 1200);
  book_TH1F("gen_jet_l2__eta"           , 60, -6, 6);
  book_TH1F("gen_jet_l2__phi"           , 60, -3.15, 3.15);
  book_TH1F("gen_jet_l1__l2__Dpt"            , 120, 0, 1200);
  book_TH1F("gen_jet_l1__l2__Deta"           , 60, -6, 6);
  book_TH1F("gen_jet_l1__l2__Dphi"           , 60, 0, 3.15);
  book_TH1F("gen_jets__eta"      , 60, -6, 6);

  //Reconstruction Hypothesis
  book_TH1F("rec_chi2", 100, 0, 100);
  book_TH2F("rec_chi2__VS__rec_ttbar__M", 300, 0, 600, 500, 0, 5000);
  book_TH1F("rec_ttbar__M"          , 5000, 0, 5000);
  book_TH1F("rec_ttbar__pt"         , 100, 0, 1000);
  book_TH1F("rec_ttbar__gen_DM"     , 120, -600, 600);
  book_TH1F("rec_ttbar__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_ttbar__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_ttbar__gen_Dpt_pct", 120, -1.2, 1.2);
  book_TH2F("rec_ttbar_M__VS__rec_ttbar__gen_DM", 4000, 0, 4000, 120, -600, 600);
  book_TH1F("rec_top__M"          , 200, 100, 300);
  book_TH1F("rec_top__gen_DM"     , 120, -600, 600);
  book_TH1F("rec_top__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH2F("gen_top_M__VS__rec_top__gen_DM_pct", 300, 0, 600, 120, -1.2, 1.2);
  book_TH1F("rec_top__pt"         , 80, 0, 800);
  book_TH1F("rec_top__p"         , 80, 0, 800);
  book_TH1F("rec_thad__M"            , 250, 50, 300);
  book_TH1F("rec_thad__pt"           , 80, 0, 800);
  book_TH1F("rec_thad__jetN"         , 7, 0, 7);
  book_TH1F("rec_thad__gen_DM"       , 240, -120, 120);
  book_TH1F("rec_thad__gen_Dpt"      , 120, -600, 600);
  book_TH1F("rec_thad__gen_Deta"     , 120, -1.2, 1.2);
  book_TH1F("rec_thad__gen_Dphi"     , 60, 0, 5.);
  book_TH1F("rec_thad__gen_DR"       , 60, 0, 3);
  book_TH1F("rec_thad__gen_DM_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_thad__gen_Dpt_pct"  , 120, -1.2, 1.2);
  book_TH1F("rec_tlep__M"          , 300, 50, 350);
  book_TH1F("rec_tlep__pt"         , 80, 0, 800);
  book_TH1F("rec_tlep__px"         , 120, -1200, 1200);
  book_TH1F("rec_tlep__py"         , 120, -1200, 1200);
  book_TH1F("rec_tlep__pz"         , 120, -1200, 1200);
  book_TH1F("rec_tlep__jetN"       , 4, 0, 4);
  book_TH1F("rec_tlep__gen_DM"     , 240, -120, 120);
  book_TH1F("rec_tlep__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_tlep__gen_Dpx"    , 120, -600, 600);
  book_TH1F("rec_tlep__gen_Dpy"    , 120, -600, 600);
  book_TH1F("rec_tlep__gen_Dpz"    , 120, -600, 600);
  book_TH1F("rec_tlep__gen_Deta"   , 120, -1.2, 1.2);
  book_TH1F("rec_tlep__gen_Dphi"   , 60, 0, 5.);
  book_TH1F("rec_tlep__gen_DR"     , 60, 0, 3);
  book_TH1F("rec_tlep__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_tlep__gen_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_tlep__gen_Dpx_pct", 120, -1.2, 1.2);
  book_TH1F("rec_tlep__gen_Dpy_pct", 120, -1.2, 1.2);
  book_TH1F("rec_tlep__gen_Dpz_pct", 120, -1.2, 1.2);
  book_TH1F("rec_tops__DR" , 60, 0, 6);
  book_TH1F("rec_tops__Dpt", 180, -450, 450);
  book_TH1F("rec_Wlep__M"          , 90, 60, 150);
  book_TH1F("rec_Wlep__Mt"         , 360, 0, 360);
  book_TH1F("rec_Wlep__pt"         , 80, 0, 800);
  book_TH1F("rec_Wlep__px"         , 120, -1200, 1200);
  book_TH1F("rec_Wlep__py"         , 120, -1200, 1200);
  book_TH1F("rec_Wlep__pz"         , 120, -1200, 1200);
  book_TH1F("rec_Wlep__gen_DM"     , 240, -120, 120);
  book_TH1F("rec_Wlep__gen_DMt"    , 240, -120, 120);
  book_TH1F("rec_Wlep__gen_DR"     , 60, 0, 3);
  book_TH1F("rec_Wlep__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_Wlep__gen_Dpx"    , 120, -600, 600);
  book_TH1F("rec_Wlep__gen_Dpy"    , 120, -600, 600);
  book_TH1F("rec_Wlep__gen_Dpz"    , 120, -600, 600);
  book_TH1F("rec_Wlep__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_Wlep__gen_DMt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_Wlep__gen_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_Wlep__gen_Dpx_pct", 120, -1.2, 1.2);
  book_TH1F("rec_Wlep__gen_Dpy_pct", 120, -1.2, 1.2);
  book_TH1F("rec_Wlep__gen_Dpz_pct", 120, -1.2, 1.2);
  book_TH1F("rec_lep__pt"            , 60, 0, 600);
  book_TH1F("rec_lep__eta"           , 60, -6, 6);
  book_TH1F("rec_lep__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_lep__gen_DR"        , 60, 0, .06);
  book_TH1F("rec_lep__gen_Dpt"       , 60, -40, 40);
  book_TH1F("rec_lep__gen_Deta"      , 60, -.06, .06);
  book_TH1F("rec_lep__gen_Dphi"      , 60, 0., .06);
  book_TH1F("rec_lep__gen_DcosThetaX", 36, -1.8, 1.8);
  book_TH1F("rec_lep__gen_Dpt_pct"   , 60, -.6, .6);
  book_TH1F("rec_lep__gen_Deta_pct"  , 60, -.06, .06);
  book_TH1F("rec_neu__pt"            , 60, 0, 600);
  book_TH1F("rec_neu__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_neu__px"            , 120, -600, 600);
  book_TH1F("rec_neu__py"            , 120, -600, 600);
  book_TH1F("rec_neu__pz"            , 120, -600, 600);
  book_TH1F("rec_neu__cosThetaX"     , 40, -1, 1);
  book_TH1F("rec_neu__gen_DR"        , 60, 0, 6);
  book_TH1F("rec_neu__gen_Dpt"       , 120, -600, 600);
  book_TH1F("rec_neu__gen_Dphi"      , 60, 0., 5.);
  book_TH1F("rec_neu__gen_Dpx"       , 120, -600, 600);
  book_TH1F("rec_neu__gen_Dpy"       , 120, -600, 600);
  book_TH1F("rec_neu__gen_Dpz"       , 120, -600, 600);
  book_TH1F("rec_neu__gen_DcosThetaX", 36, -1.8, 1.8);
  book_TH1F("rec_neu__gen_Dpt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_neu__gen_Dpx_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_neu__gen_Dpy_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_neu__gen_Dpz_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_blep__pt"            , 60, 0, 600);
  book_TH1F("rec_blep__eta"           , 60, -6, 6);
  book_TH1F("rec_blep__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_blep__px"            , 120, -1200, 1200);
  book_TH1F("rec_blep__py"            , 120, -1200, 1200);
  book_TH1F("rec_blep__pz"            , 120, -1200, 1200);
  book_TH1F("rec_blep__CSV"           , 100, 0, 1);
  book_TH1F("rec_blep__cosThetaX"     , 40, -1, 1);
  book_TH1F("rec_blep__gen_DR"        , 60, 0, 3);
  book_TH1F("rec_blep__gen_Dpt"       , 60, -300, 300);
  book_TH1F("rec_blep__gen_Deta"      , 60, -.6, .6);
  book_TH1F("rec_blep__gen_Dphi"      , 60, 0., 5.);
  book_TH1F("rec_blep__gen_Dpx"       , 60, -300, 300);
  book_TH1F("rec_blep__gen_Dpy"       , 60, -300, 300);
  book_TH1F("rec_blep__gen_Dpz"       , 60, -300, 300);
  book_TH1F("rec_blep__gen_DcosThetaX", 36, -1.8, 1.8);
  book_TH1F("rec_blep__gen_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_blep__gen_Deta_pct"  , 60, -.3, .3);
  book_TH1F("rec_blep__gen_Dpx_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_blep__gen_Dpy_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_blep__gen_Dpz_pct"   , 60, -1.2, 1.2);
  book_TH1F("deltaR__gen_lep__gen_neu"   , 60, 0., 6.);
  book_TH1F("deltaPhi__gen_lep__gen_neu" , 60, 0., 3.15);
  book_TH1F("deltaR__gen_blep__gen_neu"  , 60, 0., 6.);
  book_TH1F("deltaPhi__gen_blep__gen_neu", 60, 0., 3.15);
  book_TH1F("deltaR__gen_lep__gen_blep"  , 60, 0., 6.);
  book_TH1F("deltaR__gen_lep__gen_thad"  , 60, 0., 6.);
  book_TH1F("deltaRsum__gen_thad_jets"   , 60, 0., 6.);
  book_TH1F("deltaR__rec_lep__rec_neu"   , 60, 0., 6.);
  book_TH1F("deltaPhi__rec_lep__rec_neu" , 60, 0., 3.15);
  book_TH1F("deltaR__rec_blep__rec_neu"  , 60, 0., 6.);
  book_TH1F("deltaPhi__rec_blep__rec_neu", 60, 0., 3.15);
  book_TH1F("deltaR__rec_lep__rec_blep"  , 60, 0., 6.);
  book_TH1F("deltaR__rec_lep__rec_thad"  , 60, 0., 6.);
  book_TH1F("deltaRsum__rec_thad_jets"   , 60, 0., 6.);

  //Fit
  book_TH1F("rec_fit_chi2", 100, 0, 100);
  book_TH1F("fit_probab", 100, 0, 1);
  book_TH2F("rec_fit_chi2__VS__rec_fit_probab",50,0,50, 100,0,1 );
  book_TH2F("rec_fit_probab__VS__rec_fit_ttbar__gen_DM",100,0,1, 120,-600,600 );
  book_TH2F("rec_fit_chi2__VS__rec_fit_ttbar__M", 300, 0, 100, 1500, -500, 2500);
  book_TH1F("rec_fit_status", 21, -10, 10);
  book_TH1F("rec_fit_fit_found", 5, 0, 4);
  book_TH2F("rec_fit_probab__VS__rec_fit_thad__input_DM",100,0,1, 300,-150,150 );
  book_TH2F("rec_fit_probab__VS__rec_fit_tlep__input_DM",100,0,1, 1000,-500,500 );
  book_TH2F("rec_fit_probab__VS__rec_fit_whad__input_DM",100,0,1, 300,-150,150 );
  book_TH2F("rec_fit_probab__VS__rec_fit_wlep__input_DM",100,0,1, 1000,-500,500 );
  book_TH1F("rec_fit_jet_permutations_before_dm_cut", 50, 0, 50);
  book_TH1F("rec_fit_jet_permutations_after_dm_cut", 50, 0, 50);

  book_TH1F("rec_fit_ttbar__M"           , 4000, 0, 4000);
  book_TH1F("rec_fit_ttbar__M_input"           , 4000, 0, 4000);
  book_TH1F("rec_fit_ttbar__pt"          , 120, 0, 600);
  book_TH1F("rec_fit_ttbar__gen_DM"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_Dpt"     , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_0_340"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_0_340"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_340_440"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_340_440"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_440_540"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_440_540"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_540_640"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_540_640"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_640_800"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_640_800"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_800_1000"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_800_1000"  , 120, -1.2, 1.2);
  book_TH1F("rec_fit_ttbar__gen_DM_genM_gr1000"      , 120, -600, 600);
  book_TH1F("rec_fit_ttbar__gen_DM_pct_genM_gr1000"  , 120, -1.2, 1.2);

  book_TH1F("rec_fit_ttbar__gen_Dpt_pct" , 120, -1.2, 1.2);
  book_TH2F("gen_ttbar_M__VS__rec_fit_ttbar__gen_DM_pct", 1000, 0, 5000, 120, -3.0, 3.0);
  book_TH2F("gen_ttbar_M__VS__rec_fit_ttbar__gen_DM", 5000, 0, 5000, 120, -600, 600);
  book_TH2F("rec_fit_ttbar_M__VS__rec_fit_ttbar__gen_DM", 4000, 0, 4000, 120, -600, 600);
  book_TH1F("rec_fit_top__M"          , 1000, 0, 500);
  book_TH1F("rec_fit_top__gen_DM"    , 120, -600, 600);
  book_TH1F("rec_fit_top__gen_DM_pct"  , 120, -2.0, 2.0);
  book_TH2F("gen_top_M__VS__rec_fit_top__gen_DM_pct", 300, 0, 600, 120, -4.0, 4.0);
  book_TH1F("rec_fit_top__pt"         , 80, 0, 800);
  book_TH1F("rec_fit_top__p"         , 80, 0, 800);
  book_TH1F("rec_fit_thad__M"           , 500, 0, 500);
  book_TH1F("rec_fit_thad__M_fine"           , 2000, 170, 176);
  book_TH1F("rec_fit_thad__M_input"           , 500, 0, 500);
  book_TH1F("rec_fit_thad__pt"          , 80, 0, 800);
  book_TH1F("rec_fit_thad__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_thad__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_thad__gen_DM"      , 300, -300, 300);
  book_TH1F("rec_fit_thad__gen_Dpt"      , 120, -600, 600);
  book_TH1F("rec_fit_thad__gen_Deta"     , 120, -2., 2.);
  book_TH1F("rec_fit_thad__gen_Dphi"     , 60, 0, 3.15);
  book_TH1F("rec_fit_thad__gen_DR"       , 60, 0, 6);
  book_TH1F("rec_fit_thad__gen_DM_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_thad__gen_Dpt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_thad__input_DM"      , 300, -300, 300);
  book_TH1F("rec_fit_thad__input_Dpt"      , 120, -600, 600);
  book_TH1F("rec_fit_thad__input_Dpt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_thad__input_DEt"      , 120, -600, 600);
  book_TH1F("rec_fit_thad__input_DEt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_thad__input_Deta"     , 120, -2., 2.);
  book_TH1F("rec_fit_thad__input_Dphi"     , 60, 0, 3.15);
  book_TH1F("rec_fit_thad__input_DR"       , 60, 0, 6);
  book_TH1F("rec_fit_tlep__M"          , 500, 0, 500);
  book_TH1F("rec_fit_tlep__M_fine"           , 2000, 170, 176);
  book_TH1F("rec_fit_tlep__M_input"          , 500, 0, 500);
  book_TH1F("rec_fit_tlep__pt"         , 60, 0, 600);
  book_TH1F("rec_fit_tlep__px"         , 120, -1200, 1200);
  book_TH1F("rec_fit_tlep__py"         , 120, -1200, 1200);
  book_TH1F("rec_fit_tlep__pz"         , 120, -1200, 1200);
  book_TH1F("rec_fit_tlep__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_tlep__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_tlep__gen_DM"     , 300, -300, 300);
  book_TH1F("rec_fit_tlep__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__gen_Dpx"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__gen_Dpy"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__gen_Dpz"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__gen_Deta"   , 120, -2., 2.);
  book_TH1F("rec_fit_tlep__gen_Dphi"   , 60, 0, 3.15);
  book_TH1F("rec_fit_tlep__gen_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_tlep__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__gen_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__gen_Dpx_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__gen_Dpy_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__gen_Dpz_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__input_DM"     , 300, -300, 300);
  book_TH1F("rec_fit_tlep__input_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__input_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__input_DEt"    , 120, -600, 600);
  book_TH1F("rec_fit_tlep__input_DEt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_tlep__input_Deta"   , 120, -2., 2.);
  book_TH1F("rec_fit_tlep__input_Dphi"   , 60, 0, 3.15);
  book_TH1F("rec_fit_tlep__input_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_tops__DR" , 60, 0, 6);
  book_TH1F("rec_fit_tops__Dpt", 180, -450, 450);

  book_TH1F("rec_fit_Wlep__M"          , 360, 0, 360);
  book_TH1F("rec_fit_Wlep__M_fine"          , 3000, 70, 90);
  book_TH1F("rec_fit_Wlep__M_input"          , 360, 0, 360);
  book_TH1F("rec_fit_Wlep__pt"         , 60, 0, 600);
  book_TH1F("rec_fit_Wlep__px"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Wlep__py"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Wlep__pz"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Wlep__gen_DM"     , 240, -120, 120);
  book_TH1F("rec_fit_Wlep__gen_DMt"    , 240, -120, 120);
  book_TH1F("rec_fit_Wlep__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__gen_Dpx"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__gen_Dpy"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__gen_Dpz"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__gen_DM_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_DMt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_Dpx_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_Dpy_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_Dpz_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__gen_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_Wlep__input_DM"     , 240, -120, 120);
  book_TH1F("rec_fit_Wlep__input_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__input_Dpt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__input_DEt"    , 120, -600, 600);
  book_TH1F("rec_fit_Wlep__input_DEt_pct", 120, -1.2, 1.2);
  book_TH1F("rec_fit_Wlep__input_Deta"   , 120, -2., 2.);
  book_TH1F("rec_fit_Wlep__input_Dphi"   , 60, 0, 3.15);
  book_TH1F("rec_fit_Wlep__input_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_Whad__gen_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_Whad__M"          , 360, 0, 360);
  book_TH1F("rec_fit_Whad__M_fine"          , 3000, 70, 90);
  book_TH1F("rec_fit_Whad__ht"         , 360, 0, 360);
  book_TH1F("rec_fit_Whad__pt"         , 60, 0, 600);
  book_TH1F("rec_fit_Whad__px"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Whad__py"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Whad__pz"         , 120, -1200, 1200);
  book_TH1F("rec_fit_Whad__gen_DM"     , 240, -120, 120);
  book_TH1F("rec_fit_Whad__gen_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__gen_Dpt_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_fit_Whad__gen_Dpx"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__gen_Dpy"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__gen_Dpz"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__M_input"          , 360, 0, 360);
  book_TH1F("rec_fit_Whad__input_DM"     , 240, -120, 120);
  book_TH1F("rec_fit_Whad__input_Dpt"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__input_Dpt_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_fit_Whad__input_DEt"    , 120, -600, 600);
  book_TH1F("rec_fit_Whad__input_DEt_pct" , 120, -1.2, 1.2);
  book_TH1F("rec_fit_Whad__input_Deta"   , 120, -2., 2.);
  book_TH1F("rec_fit_Whad__input_Dphi"   , 60, 0, 3.15);
  book_TH1F("rec_fit_Whad__input_DR"     , 60, 0, 6);
  book_TH1F("rec_fit_Whad__input_DE"     , 400, -800, 800);

  book_TH1F("rec_fit_lep__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_lep__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_lep__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_lep__px"            , 120, -600, 600);
  book_TH1F("rec_fit_lep__py"            , 120, -600, 600);
  book_TH1F("rec_fit_lep__pz"            , 120, -600, 600);
  book_TH1F("rec_fit_lep__gen_DR"        , 60, 0, .06);
  book_TH1F("rec_fit_lep__gen_Dpt"       , 60, -40, 40);
  book_TH2F("rec_fit_lep__gen_Dpt__vs__gen_pt"       , 60, -40, 40, 200, 0, 200);
  book_TH1F("rec_fit_lep__gen_Deta"      , 60, -.06, .06);
  book_TH1F("rec_fit_lep__gen_Dphi"      , 60, 0., .06);
  book_TH1F("rec_fit_lep__gen_Dpt_pct"   , 60, -.6, .6);
  book_TH1F("rec_fit_lep__gen_Deta_pct"  , 60, -.06, .06);
  book_TH1F("rec_fit_lep__input_DR"        , 60, 0, .06);
  book_TH1F("rec_fit_lep__input_Dpt"       , 60, -40, 40);
  book_TH2F("rec_fit_lep__input_Dpt__vs__input_pt"       , 60, -40, 40, 200, 0, 200);
  book_TH1F("rec_fit_lep__input_Dpt_pct"   , 60, -.6, .6);
  book_TH1F("rec_fit_lep__input_DEt"       , 60, -40, 40);
  book_TH1F("rec_fit_lep__input_DEt_pct"   , 60, -.6, .6);
  book_TH1F("rec_fit_lep__input_Deta"      , 60, -.06, .06);
  book_TH1F("rec_fit_lep__input_Dphi"      , 60, 0., .06);
  book_TH1F("input_lep__gen_Dpt"       , 60, -40, 40);
  book_TH2F("input_lep__gen_Dpt__vs__gen_pt"       , 60, -40, 40, 200, 0, 200);
  book_TH1F("gen_lep__input_Deta"      , 60, -.06, .06);
  book_TH1F("gen_lep__input_Dphi"      , 60, 0., .06);
  book_TH1F("rec_fit_neu__pt"            , 150, 0, 300);
  book_TH1F("rec_fit_neu__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_neu__px"            , 120, -600, 600);
  book_TH1F("rec_fit_neu__py"            , 120, -600, 600);
  book_TH1F("rec_fit_neu__pz"            , 120, -600, 600);
  book_TH1F("rec_fit_neu__gen_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_neu__gen_Dpt"       , 120, -600, 600);
  book_TH2F("rec_fit_neu__gen_Dpt__vs__gen_pt"       , 120, -600, 600, 300, 0, 300);
  book_TH1F("rec_fit_neu__gen_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_neu__gen_Dpx"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__gen_Dpy"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__gen_Dpz"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__gen_Dpt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__gen_Dpx_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__gen_Dpy_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__gen_Dpz_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__input_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_neu__input_Dpt"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__input_Dpt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__input_DEt"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__input_DEt_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_neu__input_Deta"      , 120, -2., 2.);
  book_TH1F("rec_fit_neu__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_neu__input_Dpx"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__input_Dpy"       , 120, -600, 600);
  book_TH1F("rec_fit_neu__input_Dpz"       , 120, -600, 600);
  book_TH1F("input_neu__gen_Dpt"       , 120, -600, 600);
  book_TH2F("input_neu__gen_Dpt__vs__gen_pt"       , 120, -600, 600, 300, 0, 300);
  book_TH1F("gen_neu__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("gen_neu__input_Dpx"       , 120, -600, 600);
  book_TH1F("gen_neu__input_Dpy"       , 120, -600, 600);
  book_TH1F("gen_neu__input_Dpz"       , 120, -600, 600);

  book_TH1F("rec_fit_blep__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_blep__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_blep__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_blep__px"            , 120, -1200, 1200);
  book_TH1F("rec_fit_blep__py"            , 120, -1200, 1200);
  book_TH1F("rec_fit_blep__pz"            , 120, -1200, 1200);
  book_TH1F("rec_fit_blep__gen_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_blep__gen_Dpt"       , 60, -300, 300);
  book_TH2F("rec_fit_blep__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("rec_fit_blep__gen_Deta"      , 60, -.6, .6);
  book_TH1F("rec_fit_blep__gen_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_blep__gen_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_blep__gen_Deta_pct"  , 60, -.3, .3);
  book_TH1F("rec_fit_blep__input_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_blep__input_Dpt"       , 60, -300, 300);
  book_TH1F("rec_fit_blep__input_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_blep__input_DEt"       , 60, -300, 300);
  book_TH1F("rec_fit_blep__input_DEt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_blep__input_Deta"      , 100, -1., 1.);
  book_TH1F("rec_fit_blep__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("input_blep__gen_Dpt"       , 60, -300, 300);
  book_TH2F("input_blep__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("gen_blep__input_Deta"      , 100, -1., 1.);
  book_TH1F("gen_blep__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_bhad__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_bhad__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_bhad__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_bhad__px"            , 120, -1200, 1200);
  book_TH1F("rec_fit_bhad__py"            , 120, -1200, 1200);
  book_TH1F("rec_fit_bhad__pz"            , 120, -1200, 1200);
  book_TH1F("rec_fit_bhad__gen_Dpt"       , 60, -300, 300);
  book_TH2F("rec_fit_bhad__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("rec_fit_bhad__gen_Deta"      , 60, -.6, .6);
  book_TH1F("rec_fit_bhad__gen_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_bhad__gen_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_bhad__gen_Deta_pct"  , 60, -.3, .3);
  book_TH1F("rec_fit_bhad__gen_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_bhad__input_Dpt"       , 60, -300, 300);
  book_TH1F("rec_fit_bhad__input_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_bhad__input_DEt"       , 60, -300, 300);
  book_TH1F("rec_fit_bhad__input_DEt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_bhad__input_Deta"      , 100, -1., 1.);
  book_TH1F("rec_fit_bhad__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_bhad__input_DR"        , 60, 0, 6);
  book_TH1F("input_bhad__gen_Dpt"       , 60, -300, 300);
  book_TH2F("input_bhad__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("gen_bhad__input_Deta"      , 100, -1., 1.);
  book_TH1F("gen_bhad__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_jet_l1__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_jet_l1__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_jet_l1__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_jet_l2__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_jet_l2__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_jet_l2__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_jet_l1__l2__Dpt"            , 120, 0, 1200);
  book_TH1F("rec_fit_jet_l1__l2__Deta"           , 60, -6, 6);
  book_TH1F("rec_fit_jet_l1__l2__Dphi"           , 60, 0, 3.15);
  book_TH1F("rec_fit_jet_l1__input__Dpt"        , 60, -300, 300);
  book_TH1F("rec_fit_jet_l1__input__Dpt_pct"    , 60, -1.2, 1.2);
  book_TH1F("rec_fit_jet_l1__input__Deta"           , 60, -6, 6);
  book_TH1F("rec_fit_jet_l1__input__Dphi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_jet_l2__input__Dpt"            , 60, -300, 300);
  book_TH1F("rec_fit_jet_l2__input__Dpt_pct"    , 60, -1.2, 1.2);
  book_TH1F("rec_fit_jet_l2__input__Deta"           , 60, -6, 6);
  book_TH1F("rec_fit_jet_l2__input__Dphi"           , 60, 0, 3.15);
  book_TH1F("rec_fit_no_b_jets__pt"            , 120, 0, 1200);
  book_TH1F("rec_fit_no_b_jets__eta"           , 60, -6, 6);
  book_TH1F("rec_fit_no_b_jets__phi"           , 60, -3.15, 3.15);
  book_TH1F("rec_fit_no_b_jets__gen_Dpt"       , 60, -300, 300);
  book_TH2F("rec_fit_no_b_jets__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("rec_fit_no_b_jets__gen_Deta"      , 60, -.6, .6);
  book_TH1F("rec_fit_no_b_jets__gen_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_no_b_jets__gen_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__gen_Deta_pct"  , 60, -.3, .3);
  book_TH1F("rec_fit_no_b_jets__gen_DR"        , 60, 0, 6);
  book_TH1F("rec_fit_no_b_jets__px"            , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__py"            , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__pz"            , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__gen_Dpx"       , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__gen_Dpy"       , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__gen_Dpz"       , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__gen_Dpx_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__gen_Dpy_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__gen_Dpz_pct"   , 120, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__input_p"            , 120, 0, 1200);
  book_TH1F("rec_fit_no_b_jets__input_Dp"       , 120, -600, 600);
  book_TH1F("rec_fit_no_b_jets__input_Dpt"       , 60, -300, 300);
  book_TH1F("rec_fit_no_b_jets__input_Dpt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__input_DEt"       , 60, -300, 300);
  book_TH1F("rec_fit_no_b_jets__input_DEt_pct"   , 60, -1.2, 1.2);
  book_TH1F("rec_fit_no_b_jets__input_Deta"      , 100, -1., 1.);
  book_TH1F("rec_fit_no_b_jets__input_Dphi"      , 60, 0., 3.15);
  book_TH1F("rec_fit_no_b_jets__input_DR"        , 60, 0, 6);
  book_TH1F("input_no_b_jets__gen_Dpt"       , 60, -300, 300);
  book_TH2F("input_no_b_jets__gen_Dpt__vs__gen_pt"       , 60, -300, 300, 300, 0, 300);
  book_TH1F("gen_no_b_jets__input_Deta"      , 100, -1., 1.);
  book_TH1F("gen_no_b_jets__input_Dphi"      , 60, 0., 3.15);

  book_TH1F("rec_fit_DR_combined"        , 60, 0, 6);

  book_TH1F("unc_file_jet1_et"           , 100, 0, 500);
  book_TH1F("unc_file_jet1_eta"          , 100, 0, 0.01);
  book_TH1F("unc_file_jet1_phi"          , 100, 0, 0.05);
  book_TH1F("unc_file_jet2_et"           , 100, 0, 500);
  book_TH1F("unc_file_jet2_eta"          , 100, 0, 0.01);
  book_TH1F("unc_file_jet2_phi"          , 100, 0, 0.05);
  book_TH1F("unc_file_jetbh_et"          , 100, 0, 500);
  book_TH1F("unc_file_jetbh_eta"         , 100, 0, 0.01);
  book_TH1F("unc_file_jetbh_phi"         , 100, 0, 0.05);
  book_TH1F("unc_file_jetbl_et"          , 100, 0, 500);
  book_TH1F("unc_file_jetbl_eta"         , 100, 0, 0.01);
  book_TH1F("unc_file_jetbl_phi"         , 100, 0, 0.05);
  book_TH1F("unc_file_lepton_et"         , 200, 0, 20);
  book_TH1F("unc_file_lepton_eta"        , 100, 0, 0.0000001);
  book_TH1F("unc_file_lepton_phi"        , 100, 0, 0.0000001);
  book_TH1F("unc_file_neu_et"            , 100, 0, 500);
  book_TH1F("unc_file_neu_eta"           , 100, 0, 50);
  book_TH1F("unc_file_neu_phi"           , 100, 0, 5);


}


void KinFitTestFITHists::fill(const uhh2::Event& event){

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
    if(!hyp) throw std::runtime_error("EffyTTbarFITHists::fill -- null pointer to ReconstructionHypothesis object (from \"get_best_hypothesis\")");

    hyp_val = hyp->discriminator(disc_name_);
  }



  fill(event, hyp, hyp_val, ttgen, weight);



  return;
}

void  KinFitTestFITHists::fill(const uhh2::Event& event, const ReconstructionHypothesis* hyp, const float hyp_val, const TTbarGen* ttgen, const double weight){

  bool ttljets(false);

  float gen_ttbar_M(-1.), gen_ttbar_pt(-1.);
  float gen_thad_M(-1.), gen_thad_pt(-1.);//, gen_thad_phi(-1.), gen_thad_eta(-1.);
  float gen_tlep_M(-1.), gen_tlep_pt(-1.), gen_tlep_px(-1.), gen_tlep_py(-1.), gen_tlep_pz(-1.);
  float gen_Wlep_M(-1.), gen_Wlep_pt(-1.), gen_Wlep_px(-1.), gen_Wlep_py(-1.), gen_Wlep_pz(-1.), gen_Wlep_Mt(-1.);
  float gen_Whad_M(-1.), gen_Whad_pt(-1.), gen_Whad_px(-1.), gen_Whad_py(-1.), gen_Whad_pz(-1.);
  float gen_lep_pt(-1.), gen_lep_eta(-1.), gen_lep_phi(-1.), gen_lep_cosThetaX(-10.);
  float gen_neu_pt(-1.), gen_neu_phi(-1.), gen_neu_px(-1.), gen_neu_py(-1.), gen_neu_pz(-1.), gen_neu_cosThetaX(-10.);
  float gen_blep_pt(-1.), gen_blep_eta(-1.), gen_blep_phi(-1.), gen_blep_px(-1.), gen_blep_py(-1.), gen_blep_pz(-1.), gen_blep_cosThetaX(-10.), gen_bhad_pt(-1.),gen_bhad_phi(-1.);
  float gen_bhad_eta(-1.), gen_hadq1_eta(-1.), gen_hadq2_eta(-1.);

if(ttgen){
  ttljets = (ttgen->DecayChannel() == TTbarGen::e_muhad || ttgen->DecayChannel() == TTbarGen::e_ehad);
  LorentzVector gen_ttbar = ttgen->Top().v4()+ttgen->Antitop().v4();
    gen_ttbar_M  = gen_ttbar.M();
    gen_ttbar_pt = gen_ttbar.Pt();
    hist("gen_ttbar__M")->Fill(gen_ttbar_M , weight);
    hist("gen_ttbar__pt")->Fill(gen_ttbar_pt, weight);

    if(ttljets){
      gen_thad_M  = ttgen->TopHad().v4().M();
      gen_thad_pt = ttgen->TopHad().v4().Pt();
      gen_tlep_M  = ttgen->TopLep().v4().M();
      gen_tlep_pt = ttgen->TopLep().v4().Pt();
      gen_tlep_px = ttgen->TopLep().v4().Px();
      gen_tlep_py = ttgen->TopLep().v4().Py();
      gen_tlep_pz = ttgen->TopLep().v4().Pz();
      gen_Wlep_M  = ttgen->WLep().v4().M();
      gen_Wlep_Mt = sqrt(2*ttgen->ChargedLepton().pt()*ttgen->Neutrino().pt()*(1.-cos(delta_phi(ttgen->ChargedLepton().phi(), ttgen->Neutrino().phi()))));
      gen_Wlep_pt = ttgen->WLep().v4().Pt();
      gen_Wlep_px = ttgen->WLep().v4().Px();
      gen_Wlep_py = ttgen->WLep().v4().Py();
      gen_Wlep_pz = ttgen->WLep().v4().Pz();
      gen_Whad_M  = ttgen->WHad().v4().M();
      gen_Whad_pt = ttgen->WHad().v4().Pt();
      gen_Whad_px = ttgen->WHad().v4().Px();
      gen_Whad_py = ttgen->WHad().v4().Py();
      gen_Whad_pz = ttgen->WHad().v4().Pz();
      gen_lep_pt        = ttgen->ChargedLepton().v4().Pt();
      gen_lep_eta       = ttgen->ChargedLepton().v4().Eta();
      gen_lep_phi       = ttgen->ChargedLepton().v4().Phi();
      gen_lep_cosThetaX = cosThetaX(ttgen->ChargedLepton().v4(), ttgen->TopLep().v4(), gen_ttbar);
      gen_neu_pt        = ttgen->Neutrino().v4().Pt();
      gen_neu_phi       = ttgen->Neutrino().v4().Phi();
      gen_neu_px        = ttgen->Neutrino().v4().Px();
      gen_neu_py        = ttgen->Neutrino().v4().Py();
      gen_neu_pz        = ttgen->Neutrino().v4().Pz();
      gen_neu_cosThetaX = cosThetaX(ttgen->Neutrino().v4(), ttgen->TopLep().v4(), gen_ttbar);
      gen_blep_pt        = ttgen->BLep().v4().Pt();
      gen_blep_eta       = ttgen->BLep().v4().Eta();
      gen_blep_phi       = ttgen->BLep().v4().Phi();
      gen_blep_px        = ttgen->BLep().v4().Px();
      gen_blep_py        = ttgen->BLep().v4().Py();
      gen_blep_pz        = ttgen->BLep().v4().Pz();
      gen_blep_cosThetaX = cosThetaX(ttgen->BLep().v4(), ttgen->TopLep().v4(), gen_ttbar);
      gen_bhad_pt        = ttgen->BHad().v4().Pt();
      gen_bhad_phi       = ttgen->BHad().v4().Phi();
      H1("gen_jets__eta") ->Fill(gen_blep_eta, weight);
      gen_bhad_eta       = ttgen->BHad().v4().Eta();
      H1("gen_jets__eta") ->Fill(gen_bhad_eta, weight);
      gen_hadq1_eta       = ttgen->Q1().v4().Eta();
      H1("gen_jets__eta") ->Fill(gen_hadq1_eta, weight);
      gen_hadq2_eta       = ttgen->Q2().v4().Eta();
      H1("gen_jets__eta") ->Fill(gen_hadq2_eta, weight);
      H1("gen_top__M") ->Fill(gen_thad_M , weight);
      H1("gen_top__M") ->Fill(gen_tlep_M , weight);
      H1("gen_thad__M") ->Fill(gen_thad_M , weight);
      H1("gen_thad__pt")->Fill(gen_thad_pt, weight);
      H1("gen_tlep__M") ->Fill(gen_tlep_M , weight);
      H1("gen_tlep__pt")->Fill(gen_tlep_pt, weight);
      H1("gen_tlep__px")->Fill(gen_tlep_px, weight);
      H1("gen_tlep__py")->Fill(gen_tlep_py, weight);
      H1("gen_tlep__pz")->Fill(gen_tlep_pz, weight);
      H1("gen_tops__Dpt")->Fill(gen_tlep_pt-gen_thad_pt, weight);
      H1("gen_tops__DR") ->Fill(uhh2::deltaR(ttgen->TopLep(), ttgen->TopHad()), weight);
      H1("gen_Wlep__M") ->Fill(gen_Wlep_M , weight);
      H1("gen_Wlep__Mt")->Fill(gen_Wlep_Mt, weight);
      H1("gen_Wlep__pt")->Fill(gen_Wlep_pt, weight);
      H1("gen_Wlep__px")->Fill(gen_Wlep_px, weight);
      H1("gen_Wlep__py")->Fill(gen_Wlep_py, weight);
      H1("gen_Wlep__pz")->Fill(gen_Wlep_pz, weight);
      H1("gen_Whad__M") ->Fill(gen_Whad_M , weight);
      H1("gen_Whad__pt")->Fill(gen_Whad_pt, weight);
      H1("gen_Whad__px")->Fill(gen_Whad_px, weight);
      H1("gen_Whad__py")->Fill(gen_Whad_py, weight);
      H1("gen_Whad__pz")->Fill(gen_Whad_pz, weight);
      H1("gen_lep__pt")       ->Fill(gen_lep_pt       , weight);
      H1("gen_lep__eta")      ->Fill(gen_lep_eta      , weight);
      H1("gen_lep__phi")      ->Fill(gen_lep_phi      , weight);
      H1("gen_lep__cosThetaX")->Fill(gen_lep_cosThetaX, weight);
      H1("gen_neu__pt")       ->Fill(gen_neu_pt       , weight);
      H1("gen_neu__phi")      ->Fill(gen_neu_phi      , weight);
      H1("gen_neu__px")       ->Fill(gen_neu_px       , weight);
      H1("gen_neu__py")       ->Fill(gen_neu_py       , weight);
      H1("gen_neu__pz")       ->Fill(gen_neu_pz       , weight);
      H1("gen_neu__cosThetaX")->Fill(gen_neu_cosThetaX, weight);
      H1("gen_blep__pt")       ->Fill(gen_blep_pt       , weight);
      H1("gen_blep__eta")      ->Fill(gen_blep_eta      , weight);
      H1("gen_blep__phi")      ->Fill(gen_blep_phi      , weight);
      H1("gen_blep__px")       ->Fill(gen_blep_px       , weight);
      H1("gen_blep__py")       ->Fill(gen_blep_py       , weight);
      H1("gen_blep__pz")       ->Fill(gen_blep_pz       , weight);
      H1("gen_blep__cosThetaX")->Fill(gen_blep_cosThetaX, weight);
      H1("gen_jet_l1__pt")       ->Fill(ttgen->Q1().v4().Pt()       , weight);
      H1("gen_jet_l1__eta")      ->Fill(ttgen->Q1().v4().Eta()      , weight);
      H1("gen_jet_l1__phi")      ->Fill(ttgen->Q1().v4().Phi()      , weight);
      H1("gen_jet_l2__pt")       ->Fill(ttgen->Q2().v4().Pt()       , weight);
      H1("gen_jet_l2__eta")      ->Fill(ttgen->Q2().v4().Eta()      , weight);
      H1("gen_jet_l2__phi")      ->Fill(ttgen->Q2().v4().Phi()      , weight);
      H1("gen_jet_l1__l2__Dpt")       ->Fill(ttgen->Q1().v4().Pt() - ttgen->Q2().v4().Pt()      , weight);
      H1("gen_jet_l1__l2__Deta")       ->Fill(ttgen->Q1().v4().Eta() - ttgen->Q2().v4().Eta()      , weight);
      H1("gen_jet_l1__l2__Dphi")       ->Fill(delta_phi(ttgen->Q1().v4().Phi(), ttgen->Q2().v4().Phi())      , weight);
    }
 }


//Reconstruction Hypothesis
  if(hyp){
    // ttbar
    const float rec_ttbar_M ((hyp->top_v4()+hyp->antitop_v4()).M());
    const float rec_ttbar_pt((hyp->top_v4()+hyp->antitop_v4()).Pt());
    H1("rec_chi2")->Fill(hyp_val, weight);
    H2("rec_chi2__VS__rec_ttbar__M")->Fill(hyp_val, rec_ttbar_M, weight);
    H1("rec_ttbar__M")      ->Fill(rec_ttbar_M , weight);
    H1("rec_ttbar__pt")     ->Fill(rec_ttbar_pt, weight);
    H1("rec_ttbar__gen_DM") ->Fill(rec_ttbar_M -gen_ttbar_M , weight);
    H1("rec_ttbar__gen_Dpt")->Fill(rec_ttbar_pt-gen_ttbar_pt, weight);
    if(gen_ttbar_M){
      H1("rec_ttbar__gen_DM_pct") ->Fill((rec_ttbar_M -gen_ttbar_M) /fabs(gen_ttbar_M) , weight);
    }
    H2("rec_ttbar_M__VS__rec_ttbar__gen_DM")->Fill(rec_ttbar_M, rec_ttbar_M -gen_ttbar_M, weight);
    if(gen_ttbar_pt) H1("rec_ttbar__gen_Dpt_pct")->Fill((rec_ttbar_pt-gen_ttbar_pt)/fabs(gen_ttbar_pt), weight);
    // thad
    float rec_thad_M (hyp->tophad_v4().M());
    float rec_thad_pt(hyp->tophad_v4().Pt());
    H1("rec_thad__M")   ->Fill(rec_thad_M , weight);
    H1("rec_thad__pt")  ->Fill(rec_thad_pt, weight);
    H1("rec_thad__jetN")->Fill(hyp->tophad_jets().size(), weight);
    if(ttljets){
      H1("rec_thad__gen_DM")  ->Fill(rec_thad_M -gen_thad_M , weight);
      H1("rec_thad__gen_Dpt") ->Fill(rec_thad_pt-gen_thad_pt, weight);
      H1("rec_thad__gen_Deta")->Fill(             hyp->tophad_v4().Eta()- ttgen->TopHad().v4().Eta() , weight);
      H1("rec_thad__gen_Dphi")->Fill(delta_phi(   hyp->tophad_v4().Phi(), ttgen->TopHad().v4().Phi()), weight);
      H1("rec_thad__gen_DR")  ->Fill(uhh2::deltaR(hyp->tophad_v4()      , ttgen->TopHad().v4())      , weight);
      if(gen_thad_M)  H1("rec_thad__gen_DM_pct") ->Fill((rec_thad_M -gen_thad_M) /fabs(gen_thad_M) , weight);
      if(gen_thad_pt) H1("rec_thad__gen_Dpt_pct")->Fill((rec_thad_pt-gen_thad_pt)/fabs(gen_thad_pt), weight);
    }
    const TopJet* tjet = hyp->tophad_topjet_ptr();
    if(tjet){
      LorentzVector tjet_subjet_sum;
      for(const auto& subj : tjet->subjets()) tjet_subjet_sum += subj.v4();
      const float rec_thad_Mgro(tjet_subjet_sum.M());
      const float rec_thad_Mpru(tjet->prunedmass());
      const float rec_thad_Msdp(tjet->softdropmass());
      H1("rec_thad__Mgro")->Fill(rec_thad_Mgro, weight);
      H1("rec_thad__Mpru")->Fill(rec_thad_Mpru, weight);
      H1("rec_thad__Msdp")->Fill(rec_thad_Msdp, weight);
      if(ttljets){
        H1("rec_thad__gen_DMgro")->Fill(rec_thad_Mgro-gen_thad_M, weight);
        H1("rec_thad__gen_DMpru")->Fill(rec_thad_Mpru-gen_thad_M, weight);
        H1("rec_thad__gen_DMsdp")->Fill(rec_thad_Msdp-gen_thad_M, weight);
        if(gen_thad_M) H1("rec_thad__gen_DMgro_pct")->Fill((rec_thad_Mgro-gen_thad_M)/fabs(gen_thad_M), weight);
        if(gen_thad_M) H1("rec_thad__gen_DMpru_pct")->Fill((rec_thad_Mpru-gen_thad_M)/fabs(gen_thad_M), weight);
        if(gen_thad_M) H1("rec_thad__gen_DMsdp_pct")->Fill((rec_thad_Msdp-gen_thad_M)/fabs(gen_thad_M), weight);
      }
    }
    // tlep
    float rec_tlep_M (hyp->toplep_v4().M());
    float rec_tlep_pt(hyp->toplep_v4().Pt());
    float rec_tlep_px(hyp->toplep_v4().Px());
    float rec_tlep_py(hyp->toplep_v4().Py());
    float rec_tlep_pz(hyp->toplep_v4().Pz());
    H1("rec_tlep__M")   ->Fill(rec_tlep_M               , weight);
    H1("rec_tlep__pt")  ->Fill(rec_tlep_pt              , weight);
    H1("rec_tlep__px")  ->Fill(rec_tlep_px              , weight);
    H1("rec_tlep__py")  ->Fill(rec_tlep_py              , weight);
    H1("rec_tlep__pz")  ->Fill(rec_tlep_pz              , weight);
    H1("rec_tlep__jetN")->Fill(hyp->toplep_jets().size(), weight);
    if(ttljets){
      H1("rec_tlep__gen_DM")  ->Fill(rec_tlep_M -gen_tlep_M , weight);
      H1("rec_tlep__gen_Dpt") ->Fill(rec_tlep_pt-gen_tlep_pt, weight);
      H1("rec_tlep__gen_Dpx") ->Fill(rec_tlep_px-gen_tlep_px, weight);
      H1("rec_tlep__gen_Dpy") ->Fill(rec_tlep_py-gen_tlep_py, weight);
      H1("rec_tlep__gen_Dpz") ->Fill(rec_tlep_pz-gen_tlep_pz, weight);
      H1("rec_tlep__gen_Deta")->Fill(             hyp->toplep_v4().Eta()- ttgen->TopLep().v4().Eta() , weight);
      H1("rec_tlep__gen_Dphi")->Fill(delta_phi(   hyp->toplep_v4().Phi(), ttgen->TopLep().v4().Phi()), weight);
      H1("rec_tlep__gen_DR")  ->Fill(uhh2::deltaR(hyp->toplep_v4()      , ttgen->TopLep().v4())      , weight);
      if(gen_tlep_M)  H1("rec_tlep__gen_DM_pct") ->Fill((rec_tlep_M -gen_tlep_M) /fabs(gen_tlep_M) , weight);
      if(gen_tlep_pt) H1("rec_tlep__gen_Dpt_pct")->Fill((rec_tlep_pt-gen_tlep_pt)/fabs(gen_tlep_pt), weight);
      if(gen_tlep_px) H1("rec_tlep__gen_Dpx_pct")->Fill((rec_tlep_px-gen_tlep_px)/fabs(gen_tlep_px), weight);
      if(gen_tlep_py) H1("rec_tlep__gen_Dpy_pct")->Fill((rec_tlep_py-gen_tlep_py)/fabs(gen_tlep_py), weight);
      if(gen_tlep_pz) H1("rec_tlep__gen_Dpz_pct")->Fill((rec_tlep_pz-gen_tlep_pz)/fabs(gen_tlep_pz), weight);
    }
    H1("rec_top__M")->Fill(rec_tlep_M, weight);
    H1("rec_top__M")->Fill(rec_thad_M, weight);
    if(gen_thad_M && gen_tlep_M){
      H1("rec_top__gen_DM")->Fill((rec_thad_M - gen_thad_M), weight);
      H1("rec_top__gen_DM")->Fill((rec_tlep_M - gen_tlep_M), weight);
    }
    if(gen_thad_M && gen_tlep_M) {
      H1("rec_top__gen_DM_pct")->Fill((rec_thad_M - gen_thad_M)/fabs(gen_thad_M), weight);
      H1("rec_top__gen_DM_pct")->Fill((rec_tlep_M - gen_tlep_M)/fabs(gen_tlep_M), weight);
    }
    if(gen_tlep_M && gen_thad_M) {
      H2("gen_top_M__VS__rec_top__gen_DM_pct")->Fill(gen_thad_M, (rec_thad_M - gen_thad_M)/fabs(gen_thad_M));
      H2("gen_top_M__VS__rec_top__gen_DM_pct")->Fill(gen_tlep_M, (rec_tlep_M - gen_tlep_M)/fabs(gen_tlep_M));
    }
    H1("rec_tops__Dpt")->Fill(rec_tlep_pt-rec_thad_pt, weight);
    H1("rec_tops__DR") ->Fill(uhh2::deltaR(hyp->toplep_v4(), hyp->tophad_v4()), weight);
    H1("rec_top__pt")->Fill(rec_tlep_pt, weight);
    H1("rec_top__pt")->Fill(rec_thad_pt, weight);
    H1("rec_top__p")->Fill(hyp->toplep_v4().P(), weight);
    H1("rec_top__p")->Fill(hyp->tophad_v4().P(), weight);
    // Wlep
    float rec_Wlep_M (hyp->wlep_v4().M());
    float rec_Wlep_Mt(sqrt(2*hyp->lepton().pt()*hyp->neutrino_v4().Pt()*(1.-cos(delta_phi(hyp->lepton().phi(), hyp->neutrino_v4().Phi())))));
    float rec_Wlep_pt(hyp->wlep_v4().Pt());
    float rec_Wlep_px(hyp->wlep_v4().Px());
    float rec_Wlep_py(hyp->wlep_v4().Py());
    float rec_Wlep_pz(hyp->wlep_v4().Pz());
    H1("rec_Wlep__M") ->Fill(rec_Wlep_M , weight);
    H1("rec_Wlep__Mt")->Fill(rec_Wlep_Mt, weight);
    H1("rec_Wlep__pt")->Fill(rec_Wlep_pt, weight);
    H1("rec_Wlep__px")->Fill(rec_Wlep_px, weight);
    H1("rec_Wlep__py")->Fill(rec_Wlep_py, weight);
    H1("rec_Wlep__pz")->Fill(rec_Wlep_pz, weight);
    if(ttljets){
      H1("rec_Wlep__gen_DM") ->Fill(rec_Wlep_M -gen_Wlep_M , weight);
      H1("rec_Wlep__gen_DMt")->Fill(rec_Wlep_Mt-gen_Wlep_Mt, weight);
      H1("rec_Wlep__gen_DR")  ->Fill(uhh2::deltaR(hyp->wlep_v4()      , ttgen->WLep().v4())      , weight);
      H1("rec_Wlep__gen_Dpt")->Fill(rec_Wlep_pt-gen_Wlep_pt, weight);
      H1("rec_Wlep__gen_Dpx")->Fill(rec_Wlep_px-gen_Wlep_px, weight);
      H1("rec_Wlep__gen_Dpy")->Fill(rec_Wlep_py-gen_Wlep_py, weight);
      H1("rec_Wlep__gen_Dpz")->Fill(rec_Wlep_pz-gen_Wlep_pz, weight);
      if(gen_Wlep_M)  H1("rec_Wlep__gen_DM_pct") ->Fill((rec_Wlep_M -gen_Wlep_M) /fabs(gen_Wlep_M) , weight);
      if(gen_Wlep_Mt) H1("rec_Wlep__gen_DMt_pct")->Fill((rec_Wlep_Mt-gen_Wlep_Mt)/fabs(gen_Wlep_Mt), weight);
      if(gen_Wlep_pt) H1("rec_Wlep__gen_Dpt_pct")->Fill((rec_Wlep_pt-gen_Wlep_pt)/fabs(gen_Wlep_pt), weight);
      if(gen_Wlep_px) H1("rec_Wlep__gen_Dpx_pct")->Fill((rec_Wlep_px-gen_Wlep_px)/fabs(gen_Wlep_px), weight);
      if(gen_Wlep_py) H1("rec_Wlep__gen_Dpy_pct")->Fill((rec_Wlep_py-gen_Wlep_py)/fabs(gen_Wlep_py), weight);
      if(gen_Wlep_pz) H1("rec_Wlep__gen_Dpz_pct")->Fill((rec_Wlep_pz-gen_Wlep_pz)/fabs(gen_Wlep_pz), weight);
    }
    // lepton
    const LorentzVector& rec_lep = hyp->lepton().v4();
    float rec_lep_pt (rec_lep.Pt());
    float rec_lep_eta(rec_lep.Eta());
    float rec_lep_phi(rec_lep.Phi());
    H1("rec_lep__pt")       ->Fill(rec_lep_pt       , weight);
    H1("rec_lep__eta")      ->Fill(rec_lep_eta      , weight);
    H1("rec_lep__phi")      ->Fill(rec_lep_phi      , weight);
    if(ttljets){
      H1("rec_lep__gen_DR")        ->Fill(uhh2::deltaR(rec_lep, ttgen->ChargedLepton()), weight);
      H1("rec_lep__gen_Dpt")       ->Fill(rec_lep_pt -gen_lep_pt                      , weight);
      H1("rec_lep__gen_Deta")      ->Fill(rec_lep_eta-gen_lep_eta                     , weight);
      H1("rec_lep__gen_Dphi")      ->Fill(delta_phi(rec_lep_phi, gen_lep_phi)         , weight);
      if(gen_lep_pt)  H1("rec_lep__gen_Dpt_pct") ->Fill((rec_lep_pt -gen_lep_pt) /fabs(gen_lep_pt) , weight);
      if(gen_lep_eta) H1("rec_lep__gen_Deta_pct")->Fill((rec_lep_eta-gen_lep_eta)/fabs(gen_lep_eta), weight);
    }
    // neutrino
    const LorentzVector& rec_neu = hyp->neutrino_v4();
    float rec_neu_pt (rec_neu.Pt());
    float rec_neu_phi(rec_neu.Phi());
    float rec_neu_px (rec_neu.Px());
    float rec_neu_py (rec_neu.Py());
    float rec_neu_pz (rec_neu.Pz());
    float rec_neu_cosThetaX(cosThetaX(rec_neu, hyp->toplep_v4(), (hyp->top_v4()+hyp->antitop_v4())));
    H1("rec_neu__pt")       ->Fill(rec_neu_pt       , weight);
    H1("rec_neu__phi")      ->Fill(rec_neu_phi      , weight);
    H1("rec_neu__px")       ->Fill(rec_neu_px       , weight);
    H1("rec_neu__py")       ->Fill(rec_neu_py       , weight);
    H1("rec_neu__pz")       ->Fill(rec_neu_pz       , weight);
    H1("rec_neu__cosThetaX")->Fill(rec_neu_cosThetaX, weight);
    if(ttljets){
      H1("rec_neu__gen_DR")        ->Fill(uhh2::deltaR(rec_neu, ttgen->Neutrino()), weight);
      H1("rec_neu__gen_Dpt")       ->Fill(rec_neu_pt-gen_neu_pt                  , weight);
      H1("rec_neu__gen_Dphi")      ->Fill(delta_phi(rec_neu_phi, gen_neu_phi)    , weight);
      H1("rec_neu__gen_Dpx")       ->Fill(rec_neu_px-gen_neu_px                  , weight);
      H1("rec_neu__gen_Dpy")       ->Fill(rec_neu_py-gen_neu_py                  , weight);
      H1("rec_neu__gen_Dpz")       ->Fill(rec_neu_pz-gen_neu_pz                  , weight);
      H1("rec_neu__gen_DcosThetaX")->Fill(rec_neu_cosThetaX-gen_neu_cosThetaX    , weight);
      if(gen_neu_pt) H1("rec_neu__gen_Dpt_pct")->Fill((rec_neu_pt-gen_neu_pt)/fabs(gen_neu_pt), weight);
      if(gen_neu_px) H1("rec_neu__gen_Dpx_pct")->Fill((rec_neu_px-gen_neu_px)/fabs(gen_neu_px), weight);
      if(gen_neu_py) H1("rec_neu__gen_Dpy_pct")->Fill((rec_neu_py-gen_neu_py)/fabs(gen_neu_py), weight);
      if(gen_neu_pz) H1("rec_neu__gen_Dpz_pct")->Fill((rec_neu_pz-gen_neu_pz)/fabs(gen_neu_pz), weight);
    }

    // blep
    const LorentzVector& rec_blep = (hyp->toplep_v4()-hyp->wlep_v4());
    float rec_blep_pt (rec_blep.Pt());
    float rec_blep_eta(rec_blep.Eta());
    float rec_blep_phi(rec_blep.Phi());
    float rec_blep_px (rec_blep.Px());
    float rec_blep_py (rec_blep.Py());
    float rec_blep_pz (rec_blep.Pz());
    float rec_blep_CSV((hyp->toplep_jets().size() == 1) ? hyp->toplep_jets().at(0).btag_combinedSecondaryVertex() : -1.);
    float rec_blep_cosThetaX(cosThetaX(rec_blep, hyp->toplep_v4(), (hyp->top_v4()+hyp->antitop_v4())));
    H1("rec_blep__pt")       ->Fill(rec_blep_pt       , weight);
    H1("rec_blep__eta")      ->Fill(rec_blep_eta      , weight);
    H1("rec_blep__phi")      ->Fill(rec_blep_phi      , weight);
    H1("rec_blep__px")       ->Fill(rec_blep_px       , weight);
    H1("rec_blep__py")       ->Fill(rec_blep_py       , weight);
    H1("rec_blep__pz")       ->Fill(rec_blep_pz       , weight);
    H1("rec_blep__CSV")      ->Fill(rec_blep_CSV      , weight);
    H1("rec_blep__cosThetaX")->Fill(rec_blep_cosThetaX, weight);
    if(ttljets){
      H1("rec_blep__gen_DR")        ->Fill(uhh2::deltaR(rec_blep, ttgen->BLep()) , weight);
      H1("rec_blep__gen_Dpt")       ->Fill(rec_blep_pt -gen_blep_pt             , weight);
      H1("rec_blep__gen_Deta")      ->Fill(rec_blep_eta-gen_blep_eta            , weight);
      H1("rec_blep__gen_Dphi")      ->Fill(delta_phi(rec_blep_phi, gen_blep_phi), weight);
      H1("rec_blep__gen_Dpx")       ->Fill(rec_blep_px -gen_blep_px             , weight);
      H1("rec_blep__gen_Dpy")       ->Fill(rec_blep_py -gen_blep_py             , weight);
      H1("rec_blep__gen_Dpz")       ->Fill(rec_blep_pz -gen_blep_pz             , weight);
      H1("rec_blep__gen_DcosThetaX")->Fill(rec_blep_cosThetaX-gen_blep_cosThetaX, weight);
      if(gen_blep_pt)  H1("rec_blep__gen_Dpt_pct") ->Fill((rec_blep_pt -gen_blep_pt) /fabs(gen_blep_pt) , weight);
      if(gen_blep_eta) H1("rec_blep__gen_Deta_pct")->Fill((rec_blep_eta-gen_blep_eta)/fabs(gen_blep_eta), weight);
      if(gen_blep_px)  H1("rec_blep__gen_Dpx_pct") ->Fill((rec_blep_px -gen_blep_px) /fabs(gen_blep_px) , weight);
      if(gen_blep_py)  H1("rec_blep__gen_Dpy_pct") ->Fill((rec_blep_py -gen_blep_py) /fabs(gen_blep_py) , weight);
      if(gen_blep_pz)  H1("rec_blep__gen_Dpz_pct") ->Fill((rec_blep_pz -gen_blep_pz) /fabs(gen_blep_pz) , weight);
    }

    if(ttljets){
      H1("deltaR__gen_lep__gen_neu") ->Fill(uhh2::deltaR(ttgen->ChargedLepton(), ttgen->Neutrino()), weight);
      H1("deltaR__gen_lep__gen_blep")->Fill(uhh2::deltaR(ttgen->ChargedLepton(), ttgen->BLep())    , weight);
      H1("deltaR__gen_lep__gen_thad")->Fill(uhh2::deltaR(ttgen->ChargedLepton(), ttgen->TopHad())  , weight);
      H1("deltaR__gen_blep__gen_neu")->Fill(uhh2::deltaR(ttgen->BLep()         , ttgen->Neutrino()), weight);
      H1("deltaPhi__gen_lep__gen_neu") ->Fill(delta_phi(ttgen->ChargedLepton().v4().Phi(), ttgen->Neutrino().v4().Phi()), weight);
      H1("deltaPhi__gen_blep__gen_neu")->Fill(delta_phi(ttgen->BLep()         .v4().Phi(), ttgen->Neutrino().v4().Phi()), weight);
    }
    H1("deltaR__rec_lep__rec_neu") ->Fill(uhh2::deltaR(rec_lep , rec_neu)         , weight);
    H1("deltaR__rec_lep__rec_blep")->Fill(uhh2::deltaR(rec_lep , rec_blep)        , weight);
    H1("deltaR__rec_lep__rec_thad")->Fill(uhh2::deltaR(rec_lep , hyp->tophad_v4()), weight);
    H1("deltaR__rec_blep__rec_neu")->Fill(uhh2::deltaR(rec_blep, rec_neu)         , weight);
    H1("deltaPhi__rec_lep__rec_neu") ->Fill(delta_phi(rec_lep .Phi(), rec_neu.Phi()), weight);
    H1("deltaPhi__rec_blep__rec_neu")->Fill(delta_phi(rec_blep.Phi(), rec_neu.Phi()), weight);
    float deltaRsum__rec_thad_jets(0.);
    for  (unsigned int i=0  ; i<hyp->tophad_jets().size(); ++i){
      for(unsigned int j=i+1; j<hyp->tophad_jets().size(); ++j){
        deltaRsum__rec_thad_jets += uhh2::deltaR(hyp->tophad_jets().at(i), hyp->tophad_jets().at(j));
      }
    }
    if(!deltaRsum__rec_thad_jets) deltaRsum__rec_thad_jets = -1.;
    H1("deltaRsum__rec_thad_jets")->Fill(deltaRsum__rec_thad_jets, weight);
  }

  float rec_fit_chi2__(-1.), rec_fit_ttbar_M_(-1.), rec_fit_ttbar_pt_(-1.), rec_fit_status(-100.);
  float rec_fit_thad_M_(-1.), rec_fit_thad_pt_(-1.),rec_fit_thad_Et_(-1.), rec_fit_thad_eta_(-1.), rec_fit_thad_phi_(-1.);
 float rec_fit_tlep_M_(-1.), rec_fit_tlep_pt_(-1.), rec_fit_tlep_Et_(-1.), rec_fit_tlep_px_(-1.), rec_fit_tlep_py_(-1.), rec_fit_tlep_pz_(-1.), rec_fit_tlep_eta_(-1.), rec_fit_tlep_phi_(-1.);
 float rec_fit_Wlep_M_(-1.), rec_fit_Wlep_Mt_(-1.), rec_fit_Wlep_pt_(-1.), rec_fit_Wlep_Et_(-1.), rec_fit_Wlep_px_(-1.), rec_fit_Wlep_py_(-1.), rec_fit_Wlep_pz_(-1.);
 float rec_fit_Whad_M_(-1.), rec_fit_Whad_ht_(-1.), rec_fit_Whad_pt_(-1.), rec_fit_Whad_Et_(-1.), rec_fit_Whad_px_(-1.), rec_fit_Whad_py_(-1.), rec_fit_Whad_pz_(-1.);
 float rec_fit_lep_pt_(-1.),rec_fit_lep_Et_(-1.), rec_fit_lep_eta_(-1.), rec_fit_lep_phi_(-1.), rec_fit_lep_px_(-1.), rec_fit_lep_py_(-1.), rec_fit_lep_pz_(-1.), rec_fit_neu_pt_(-1.),rec_fit_neu_Et_(-1.),rec_fit_neu_phi_(-1.), rec_fit_neu_px_(-1.),rec_fit_neu_py_(-1.) , rec_fit_neu_pz_(-1.) ;
 float  rec_fit_blep_pt_(-1.),rec_fit_blep_Et_(-1.), rec_fit_blep_eta_(-1.), rec_fit_blep_phi_(-1.), rec_fit_blep_px_(-1.), rec_fit_blep_py_(-1.), rec_fit_blep_pz_(-1.),  rec_fit_bhad_pt_(-1.), rec_fit_bhad_Et_(-1.), rec_fit_bhad_eta_(-1.), rec_fit_bhad_phi_(-1.), rec_fit_bhad_px_(-1.), rec_fit_bhad_py_(-1.), rec_fit_bhad_pz_(-1.);
 bool fit_found(false);
 TLorentzVector rec_fit_Thad = event.get(rec_fit_thad_v4);
 TLorentzVector rec_fit_Tlep = event.get(rec_fit_tlep_v4);
 TLorentzVector rec_fit_Wlep = event.get(rec_fit_wlep_v4);
 TLorentzVector rec_fit_Whad = event.get(rec_fit_whad_v4);
 TLorentzVector rec_fit_Bhad = event.get(rec_fit_bhad_v4);
 TLorentzVector rec_fit_Blep = event.get(rec_fit_blep_v4);
 TLorentzVector rec_fit_Lepton = event.get(rec_fit_lepton_v4);
 TLorentzVector rec_fit_Neu = event.get(rec_fit_neutrino_v4);
 TLorentzVector rec_fit_jetq1 = event.get(rec_fit_jetq1_v4);
 TLorentzVector rec_fit_jetq2 = event.get(rec_fit_jetq2_v4);
 TLorentzVector input_jetq1 = event.get(input_jetq1_v4);
 TLorentzVector input_jetq2 = event.get(input_jetq2_v4);
 TLorentzVector input_jetbh = event.get(input_jetbh_v4);
 TLorentzVector input_jetbl = event.get(input_jetbl_v4);
 TLorentzVector input_lepton = event.get(input_lepton_v4);
 TLorentzVector input_neutrino = event.get(input_neutrino_v4);
 TLorentzVector input_wlep = (input_lepton + input_neutrino);
 TLorentzVector input_whad = (input_jetq1 + input_jetq2);
 TLorentzVector input_tlep = (input_lepton + input_neutrino + input_jetbl);
 TLorentzVector input_thad = (input_jetq1 + input_jetq2 + input_jetbh);
 TLorentzVector input_ttbar = (input_tlep + input_thad);

 int jet_permutation_before = event.get(jet_permutation_before_dm_cut);
 int jet_permutation_after = event.get(jet_permutation_after_dm_cut);

 rec_fit_chi2__=event.get(rec_fit_chi2_);
 H1("rec_fit_chi2")->Fill(rec_fit_chi2__, weight);
 float rec_fit_probab_ = event.get(fit_probab);
 H1("fit_probab")->Fill(rec_fit_probab_, weight);
 H2("rec_fit_chi2__VS__rec_fit_probab")->Fill(rec_fit_chi2__, rec_fit_probab_, weight);
 rec_fit_status = event.get(rec_fit_status_);
 H1("rec_fit_status")->Fill(rec_fit_status);
 fit_found = event.get(rec_fit_fit_found_);
 H1("rec_fit_fit_found")->Fill(fit_found);
 H1("rec_fit_jet_permutations_before_dm_cut")->Fill(jet_permutation_before);
 H1("rec_fit_jet_permutations_after_dm_cut")->Fill(jet_permutation_after);
  //ttbar
  TLorentzVector rec_fit_ttbar(rec_fit_Thad + rec_fit_Tlep);
  rec_fit_ttbar_M_= rec_fit_ttbar.M();
  rec_fit_ttbar_pt_=rec_fit_ttbar.Pt();
  H2("rec_fit_chi2__VS__rec_fit_ttbar__M")->Fill(rec_fit_chi2__, rec_fit_ttbar_M_, weight);
  H2("rec_fit_probab__VS__rec_fit_ttbar__gen_DM")->Fill(rec_fit_probab_, rec_fit_ttbar_M_ - gen_ttbar_M, weight);
  H1("rec_fit_ttbar__M")->Fill(rec_fit_ttbar_M_, weight);
  H1("rec_fit_ttbar__M_input")->Fill(input_ttbar.M(), weight);
  H1("rec_fit_ttbar__pt")->Fill(rec_fit_ttbar_pt_, weight);
  H1("rec_fit_ttbar__gen_DM") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M<=340)  H1("rec_fit_ttbar__gen_DM_genM_0_340") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>340 && gen_ttbar_M<=440)  H1("rec_fit_ttbar__gen_DM_genM_340_440") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>440 && gen_ttbar_M<=540)  H1("rec_fit_ttbar__gen_DM_genM_440_540") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>540 && gen_ttbar_M<=640)  H1("rec_fit_ttbar__gen_DM_genM_540_640") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>640 && gen_ttbar_M<=800)  H1("rec_fit_ttbar__gen_DM_genM_640_800") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>800 && gen_ttbar_M<=1000)  H1("rec_fit_ttbar__gen_DM_genM_800_1000") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  if(gen_ttbar_M>1000)  H1("rec_fit_ttbar__gen_DM_genM_gr1000") ->Fill(rec_fit_ttbar_M_ - gen_ttbar_M , weight);
  H2("rec_fit_ttbar_M__VS__rec_fit_ttbar__gen_DM")->Fill(rec_fit_ttbar_M_, rec_fit_ttbar_M_ -gen_ttbar_M, weight);
  H1("rec_fit_ttbar__gen_Dpt")->Fill(rec_fit_ttbar_pt_ - gen_ttbar_pt, weight);
  if(gen_ttbar_M)  {
    H1("rec_fit_ttbar__gen_DM_pct") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M) /fabs(gen_ttbar_M) , weight);
    if(gen_ttbar_M<=340)  H1("rec_fit_ttbar__gen_DM_pct_genM_0_340") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M ) /fabs(gen_ttbar_M), weight);
    if(gen_ttbar_M>340 && gen_ttbar_M<=440)  H1("rec_fit_ttbar__gen_DM_pct_genM_340_440") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M) /fabs(gen_ttbar_M) , weight);
    if(gen_ttbar_M>440 && gen_ttbar_M<=540)  H1("rec_fit_ttbar__gen_DM_pct_genM_440_540") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M) /fabs(gen_ttbar_M)  , weight);
    if(gen_ttbar_M>540 && gen_ttbar_M<=640)  H1("rec_fit_ttbar__gen_DM_pct_genM_540_640") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M) /fabs(gen_ttbar_M) , weight);
    if(gen_ttbar_M>640 && gen_ttbar_M<=800)  H1("rec_fit_ttbar__gen_DM_pct_genM_640_800") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M ) /fabs(gen_ttbar_M), weight);
    if(gen_ttbar_M>800 && gen_ttbar_M<=1000)  H1("rec_fit_ttbar__gen_DM_pct_genM_800_1000") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M) /fabs(gen_ttbar_M) , weight);
    if(gen_ttbar_M>1000)  H1("rec_fit_ttbar__gen_DM_pct_genM_gr1000") ->Fill((rec_fit_ttbar_M_ - gen_ttbar_M ) /fabs(gen_ttbar_M), weight);
    H2("gen_ttbar_M__VS__rec_fit_ttbar__gen_DM_pct")->Fill(gen_ttbar_M, (rec_fit_ttbar_M_ -gen_ttbar_M) /fabs(gen_ttbar_M), weight);
    H2("gen_ttbar_M__VS__rec_fit_ttbar__gen_DM")->Fill(gen_ttbar_M, (rec_fit_ttbar_M_ -gen_ttbar_M), weight);
  }
  if(gen_ttbar_pt) H1("rec_fit_ttbar__gen_Dpt_pct")->Fill((rec_fit_ttbar_pt_ - gen_ttbar_pt)/fabs(gen_ttbar_pt), weight);
  //top hadronic side
  rec_fit_thad_M_=rec_fit_Thad.M();
  rec_fit_thad_pt_=rec_fit_Thad.Pt();
  rec_fit_thad_Et_=rec_fit_Thad.Et();
  rec_fit_thad_eta_=rec_fit_Thad.Eta();
  rec_fit_thad_phi_=rec_fit_Thad.Phi();
  H1("rec_fit_thad__M")->Fill(rec_fit_thad_M_ ,weight);
  H1("rec_fit_thad__M_fine")->Fill(rec_fit_thad_M_ ,weight);
  H1("rec_fit_thad__M_input")->Fill(input_thad.M() ,weight);
  H1("rec_fit_thad__pt")->Fill(rec_fit_thad_pt_ ,weight);
  H1("rec_fit_thad__eta")->Fill(rec_fit_thad_eta_ ,weight);
  H1("rec_fit_thad__phi")->Fill(rec_fit_thad_phi_ ,weight);
  if(ttljets){
    H1("rec_fit_thad__gen_DM")  ->Fill(rec_fit_thad_M_ - gen_thad_M , weight);
    H1("rec_fit_thad__gen_Dpt") ->Fill(rec_fit_thad_pt_ - gen_thad_pt, weight);
    H1("rec_fit_thad__gen_Deta")->Fill(rec_fit_thad_eta_ - ttgen->TopHad().v4().Eta() , weight);
    H1("rec_fit_thad__gen_Dphi")->Fill(delta_phi(rec_fit_thad_phi_, ttgen->TopHad().v4().Phi()), weight);
    //    H1("rec_fit_thad__gen_DR")  ->Fill(uhh2::deltaR(thad_v4, ttgen->TopHad().v4())      , weight);
    if(gen_thad_M)  H1("rec_fit_thad__gen_DM_pct") ->Fill((rec_fit_thad_M_ - gen_thad_M) /fabs(gen_thad_M) , weight);
    if(gen_thad_pt) H1("rec_fit_thad__gen_Dpt_pct")->Fill((rec_fit_thad_pt_ - gen_thad_pt)/fabs(gen_thad_pt), weight);
  }
  H1("rec_fit_thad__input_DM")  ->Fill(rec_fit_thad_M_ - input_thad.M() , weight);
  H1("rec_fit_thad__input_Dpt") ->Fill(rec_fit_thad_pt_ - input_thad.Pt(), weight);
  H1("rec_fit_thad__input_Dpt_pct") ->Fill((rec_fit_thad_pt_ - input_thad.Pt())/fabs(input_thad.Pt()), weight);
  H1("rec_fit_thad__input_DEt") ->Fill(rec_fit_thad_Et_ - input_thad.Et(), weight);
  H1("rec_fit_thad__input_DEt_pct") ->Fill((rec_fit_thad_Et_ - input_thad.Et())/fabs(input_thad.Et()), weight);
  H1("rec_fit_thad__input_Deta")->Fill(rec_fit_thad_eta_ - input_thad.Eta() , weight);
  H1("rec_fit_thad__input_Dphi")->Fill(delta_phi(rec_fit_thad_phi_, input_thad.Phi()), weight);
  H2("rec_fit_probab__VS__rec_fit_thad__input_DM")->Fill(rec_fit_probab_, rec_fit_thad_M_ - input_thad.M());
  //top leptonic side
  rec_fit_tlep_M_=rec_fit_Tlep.M();
  rec_fit_tlep_pt_=rec_fit_Tlep.Pt();
  rec_fit_tlep_Et_=rec_fit_Tlep.Et();
  rec_fit_tlep_px_=rec_fit_Tlep.Px();
  rec_fit_tlep_py_=rec_fit_Tlep.Py();
  rec_fit_tlep_pz_=rec_fit_Tlep.Pz();
  rec_fit_tlep_eta_=rec_fit_Tlep.Eta();
  rec_fit_tlep_phi_=rec_fit_Tlep.Phi();
  H1("rec_fit_tlep__M")->Fill(rec_fit_tlep_M_ ,weight);
  H1("rec_fit_tlep__M_fine")->Fill(rec_fit_tlep_M_ ,weight);
  H1("rec_fit_tlep__M_input")->Fill(input_tlep.M() ,weight);
  H1("rec_fit_tlep__pt")->Fill(rec_fit_tlep_pt_ ,weight);
  H1("rec_fit_tlep__px")->Fill(rec_fit_tlep_px_ ,weight);
  H1("rec_fit_tlep__py")->Fill(rec_fit_tlep_py_ ,weight);
  H1("rec_fit_tlep__pz")->Fill(rec_fit_tlep_pz_ ,weight);
  H1("rec_fit_tlep__eta")->Fill(rec_fit_tlep_eta_ ,weight);
  H1("rec_fit_tlep__phi")->Fill(rec_fit_tlep_phi_ ,weight);
  if(ttljets){
    H1("rec_fit_tlep__gen_DM")  ->Fill(rec_fit_tlep_M_ - gen_tlep_M , weight);
    H1("rec_fit_tlep__gen_Dpt") ->Fill(rec_fit_tlep_pt_ - gen_tlep_pt, weight);
    H1("rec_fit_tlep__gen_Dpx") ->Fill(rec_fit_tlep_px_ - gen_tlep_px, weight);
    H1("rec_fit_tlep__gen_Dpy") ->Fill(rec_fit_tlep_py_ - gen_tlep_py, weight);
    H1("rec_fit_tlep__gen_Dpz") ->Fill(rec_fit_tlep_pz_ - gen_tlep_pz, weight);
    H1("rec_fit_tlep__gen_Deta")->Fill(rec_fit_tlep_eta_ - ttgen->TopLep().v4().Eta() , weight);
    H1("rec_fit_tlep__gen_Dphi")->Fill(delta_phi(rec_fit_tlep_phi_, ttgen->TopLep().v4().Phi()), weight);
    if(gen_tlep_M)  H1("rec_fit_tlep__gen_DM_pct") ->Fill((rec_fit_tlep_M_ - gen_tlep_M) /fabs(gen_tlep_M) , weight);
    if(gen_tlep_pt) H1("rec_fit_tlep__gen_Dpt_pct")->Fill((rec_fit_tlep_pt_ - gen_tlep_pt)/fabs(gen_tlep_pt), weight);
    if(gen_tlep_px) H1("rec_fit_tlep__gen_Dpx_pct")->Fill((rec_fit_tlep_px_ - gen_tlep_px)/fabs(gen_tlep_px), weight);
    if(gen_tlep_py) H1("rec_fit_tlep__gen_Dpy_pct")->Fill((rec_fit_tlep_py_ - gen_tlep_py)/fabs(gen_tlep_py), weight);
    if(gen_tlep_pz) H1("rec_fit_tlep__gen_Dpz_pct")->Fill((rec_fit_tlep_pz_ - gen_tlep_pz)/fabs(gen_tlep_pz), weight);
  }
  H1("rec_fit_tlep__input_DM")  ->Fill(rec_fit_tlep_M_ - input_tlep.M() , weight);
  H1("rec_fit_tlep__input_Dpt") ->Fill(rec_fit_tlep_pt_ - input_tlep.Pt(), weight);
  H1("rec_fit_tlep__input_Dpt_pct") ->Fill((rec_fit_tlep_pt_ - input_tlep.Pt())/fabs(input_tlep.Pt()), weight);
  H1("rec_fit_tlep__input_DEt") ->Fill(rec_fit_tlep_Et_ - input_tlep.Et(), weight);
  H1("rec_fit_tlep__input_DEt_pct") ->Fill((rec_fit_tlep_Et_ - input_tlep.Et())/fabs(input_tlep.Et()), weight);
  H1("rec_fit_tlep__input_Deta")->Fill(rec_fit_tlep_eta_ - input_tlep.Eta() , weight);
  H1("rec_fit_tlep__input_Dphi")->Fill(delta_phi(rec_fit_tlep_phi_, input_tlep.Phi()), weight);
  H2("rec_fit_probab__VS__rec_fit_tlep__input_DM")->Fill(rec_fit_probab_, rec_fit_tlep_M_ - input_tlep.M());
  H1("rec_fit_top__M")->Fill(rec_fit_tlep_M_, weight);
  H1("rec_fit_top__M")->Fill(rec_fit_thad_M_, weight);
   if(gen_thad_M && gen_tlep_M){
      H1("rec_fit_top__gen_DM")->Fill((rec_fit_thad_M_ - gen_thad_M), weight);
      H1("rec_fit_top__gen_DM")->Fill((rec_fit_tlep_M_ - gen_tlep_M), weight);
    }
    if(gen_thad_M && gen_tlep_M) {
      H1("rec_fit_top__gen_DM_pct")->Fill((rec_fit_thad_M_ - gen_thad_M)/fabs(gen_thad_M), weight);
      H1("rec_fit_top__gen_DM_pct")->Fill((rec_fit_tlep_M_ - gen_tlep_M)/fabs(gen_tlep_M), weight);
    }
    if(gen_tlep_M && gen_thad_M) {
      H2("gen_top_M__VS__rec_fit_top__gen_DM_pct")->Fill(gen_thad_M, (rec_fit_thad_M_ - gen_thad_M)/fabs(gen_thad_M));
      H2("gen_top_M__VS__rec_fit_top__gen_DM_pct")->Fill(gen_tlep_M, (rec_fit_tlep_M_ - gen_tlep_M)/fabs(gen_tlep_M));
    }
  H1("rec_fit_tops__Dpt")->Fill(rec_fit_tlep_pt_ - rec_fit_thad_pt_, weight);
  H1("rec_fit_top__pt")->Fill(rec_fit_tlep_pt_, weight);
  H1("rec_fit_top__pt")->Fill(rec_fit_thad_pt_, weight);
  H1("rec_fit_top__p")->Fill(rec_fit_Tlep.P(), weight);
  H1("rec_fit_top__p")->Fill(rec_fit_Thad.P(), weight);
    // Wlep
  rec_fit_Wlep_M_ =rec_fit_Wlep.M();
  rec_fit_Wlep_pt_=rec_fit_Wlep.Pt();
  rec_fit_Wlep_Et_=rec_fit_Wlep.Et();
  rec_fit_Wlep_px_=rec_fit_Wlep.Px();
  rec_fit_Wlep_py_=rec_fit_Wlep.Py();
  rec_fit_Wlep_pz_=rec_fit_Wlep.Pz();
  H1("rec_fit_Wlep__M") ->Fill(rec_fit_Wlep_M_ , weight);
  H1("rec_fit_Wlep__M_fine") ->Fill(rec_fit_Wlep_M_ , weight);
  H1("rec_fit_Wlep__pt")->Fill(rec_fit_Wlep_pt_, weight);
  H1("rec_fit_Wlep__px")->Fill(rec_fit_Wlep_px_, weight);
  H1("rec_fit_Wlep__py")->Fill(rec_fit_Wlep_py_, weight);
  H1("rec_fit_Wlep__pz")->Fill(rec_fit_Wlep_pz_, weight);
  if(ttljets){
    H1("rec_fit_Wlep__gen_DM") ->Fill(rec_fit_Wlep_M_ - gen_Wlep_M , weight);
    H1("rec_fit_Wlep__gen_DMt")->Fill(rec_fit_Wlep_Mt_ - gen_Wlep_Mt, weight);
    H1("rec_fit_Wlep__gen_Dpt")->Fill(rec_fit_Wlep_pt_ - gen_Wlep_pt, weight);
    H1("rec_fit_Wlep__gen_Dpx")->Fill(rec_fit_Wlep_px_ - gen_Wlep_px, weight);
    H1("rec_fit_Wlep__gen_Dpy")->Fill(rec_fit_Wlep_py_ - gen_Wlep_py, weight);
    H1("rec_fit_Wlep__gen_Dpz")->Fill(rec_fit_Wlep_pz_ - gen_Wlep_pz, weight);
    if(gen_Wlep_M)  H1("rec_fit_Wlep__gen_DM_pct") ->Fill((rec_fit_Wlep_M_ -gen_Wlep_M) /fabs(gen_Wlep_M) , weight);
    if(gen_Wlep_Mt) H1("rec_fit_Wlep__gen_DMt_pct")->Fill((rec_fit_Wlep_Mt_-gen_Wlep_Mt)/fabs(gen_Wlep_Mt), weight);
    if(gen_Wlep_pt) H1("rec_fit_Wlep__gen_Dpt_pct")->Fill((rec_fit_Wlep_pt_-gen_Wlep_pt)/fabs(gen_Wlep_pt), weight);
    if(gen_Wlep_px) H1("rec_fit_Wlep__gen_Dpx_pct")->Fill((rec_fit_Wlep_px_-gen_Wlep_px)/fabs(gen_Wlep_px), weight);
    if(gen_Wlep_py) H1("rec_fit_Wlep__gen_Dpy_pct")->Fill((rec_fit_Wlep_py_-gen_Wlep_py)/fabs(gen_Wlep_py), weight);
    if(gen_Wlep_pz) H1("rec_fit_Wlep__gen_Dpz_pct")->Fill((rec_fit_Wlep_pz_-gen_Wlep_pz)/fabs(gen_Wlep_pz), weight);
  }
  H1("rec_fit_Wlep__M_input") ->Fill(input_wlep.M() , weight);
  H1("rec_fit_Wlep__input_DM") ->Fill(rec_fit_Wlep_M_ - input_wlep.M() , weight);
  H1("rec_fit_Wlep__input_Dpt")->Fill(rec_fit_Wlep_pt_ - input_wlep.Pt(), weight);
  H1("rec_fit_Wlep__input_Dpt_pct")->Fill((rec_fit_Wlep_pt_ - input_wlep.Pt())/fabs(input_wlep.Pt()), weight);
  H1("rec_fit_Wlep__input_DEt")->Fill(rec_fit_Wlep_Et_ - input_wlep.Et(), weight);
  H1("rec_fit_Wlep__input_DEt_pct")->Fill((rec_fit_Wlep_Et_ - input_wlep.Et())/fabs(input_wlep.Et()), weight);
  H1("rec_fit_Wlep__input_Deta")->Fill(rec_fit_Wlep.Eta() - input_wlep.Eta(), weight);
  H1("rec_fit_Wlep__input_Dphi")->Fill(delta_phi(rec_fit_Wlep.Phi(), input_wlep.Phi()), weight);

  H2("rec_fit_probab__VS__rec_fit_wlep__input_DM")->Fill(rec_fit_probab_, rec_fit_Wlep_M_ - input_wlep.M());

  rec_fit_Whad_M_ =rec_fit_Whad.M();
  rec_fit_Whad_ht_=(rec_fit_jetq1.Pt() + rec_fit_jetq2.Pt());
  rec_fit_Whad_pt_=rec_fit_Whad.Pt();
  rec_fit_Whad_Et_=rec_fit_Whad.Et();
  rec_fit_Whad_px_=rec_fit_Whad.Px();
  rec_fit_Whad_py_=rec_fit_Whad.Py();
  rec_fit_Whad_pz_=rec_fit_Whad.Pz();
  H1("rec_fit_Whad__M") ->Fill(rec_fit_Whad_M_ , weight);
  H1("rec_fit_Whad__M_fine") ->Fill(rec_fit_Whad_M_ , weight);
  H1("rec_fit_Whad__ht")->Fill(rec_fit_Whad_ht_, weight);
  H1("rec_fit_Whad__pt")->Fill(rec_fit_Whad_pt_, weight);
  H1("rec_fit_Whad__px")->Fill(rec_fit_Whad_px_, weight);
  H1("rec_fit_Whad__py")->Fill(rec_fit_Whad_py_, weight);
  H1("rec_fit_Whad__pz")->Fill(rec_fit_Whad_pz_, weight);
  H1("rec_fit_Whad__gen_DM") ->Fill(rec_fit_Whad_M_ - gen_Whad_M , weight);
  H1("rec_fit_Whad__gen_Dpt")->Fill(rec_fit_Whad_pt_ - gen_Whad_pt, weight);
  H1("rec_fit_Whad__gen_Dpx")->Fill(rec_fit_Whad_px_ - gen_Whad_px, weight);
  H1("rec_fit_Whad__gen_Dpy")->Fill(rec_fit_Whad_py_ - gen_Whad_py, weight);
  H1("rec_fit_Whad__gen_Dpz")->Fill(rec_fit_Whad_pz_ - gen_Whad_pz, weight);
  H1("rec_fit_Whad__M_input") ->Fill(input_whad.M() , weight);
  H1("rec_fit_Whad__input_DM") ->Fill(rec_fit_Whad_M_ - input_whad.M() , weight);
  H1("rec_fit_Whad__input_Dpt")->Fill(rec_fit_Whad_pt_ - input_whad.Pt(), weight);
  H1("rec_fit_Whad__input_Dpt_pct")->Fill((rec_fit_Whad_pt_ - input_whad.Pt())/fabs(input_whad.Pt()), weight);
  H1("rec_fit_Whad__input_DEt")->Fill(rec_fit_Whad_Et_ - input_whad.Et(), weight);
  H1("rec_fit_Whad__input_DEt_pct")->Fill((rec_fit_Whad_Et_ - input_whad.Et())/fabs(input_whad.Et()), weight);
  H1("rec_fit_Whad__input_Deta")->Fill(rec_fit_Whad.Eta() - input_whad.Eta(), weight);
  H1("rec_fit_Whad__input_Dphi")->Fill(delta_phi(rec_fit_Whad.Phi(), input_whad.Phi()), weight);
  H1("rec_fit_Whad__input_DE")->Fill(rec_fit_Whad.E() - input_whad.E(), weight);

  H2("rec_fit_probab__VS__rec_fit_whad__input_DM")->Fill(rec_fit_probab_,rec_fit_Whad_M_ - input_whad.M());
  if(ttljets){
    LorentzVector rec_fit_Whad_v4;
    rec_fit_Whad_v4.SetXYZT(rec_fit_Whad.X(), rec_fit_Whad.Y(), rec_fit_Whad.Z(), rec_fit_Whad.T());
    const LorentzVector gen_Whad_v4 = ttgen->WHad().v4();
    H1("rec_fit_Whad__gen_DR")        ->Fill(deltaR_fit(rec_fit_Whad_v4, gen_Whad_v4), weight);

    LorentzVector rec_fit_Wlep_v4;
    rec_fit_Wlep_v4.SetXYZT(rec_fit_Wlep.X(), rec_fit_Wlep.Y(), rec_fit_Wlep.Z(), rec_fit_Wlep.T());
    const LorentzVector gen_Wlep_v4 = ttgen->WLep().v4();
    H1("rec_fit_Wlep__gen_DR")        ->Fill(deltaR_fit(rec_fit_Wlep_v4, gen_Wlep_v4), weight);

    LorentzVector rec_fit_Thad_v4;
    rec_fit_Thad_v4.SetXYZT(rec_fit_Thad.X(), rec_fit_Thad.Y(), rec_fit_Thad.Z(), rec_fit_Thad.T());
    const LorentzVector gen_Thad_v4 = ttgen->TopHad().v4();
    H1("rec_fit_thad__gen_DR")        ->Fill(deltaR_fit(rec_fit_Thad_v4, gen_Thad_v4), weight);
  
    LorentzVector rec_fit_Tlep_v4;
    rec_fit_Tlep_v4.SetXYZT(rec_fit_Tlep.X(), rec_fit_Tlep.Y(), rec_fit_Tlep.Z(), rec_fit_Tlep.T());
    const LorentzVector gen_Tlep_v4 = ttgen->TopLep().v4();
    H1("rec_fit_tlep__gen_DR")        ->Fill(deltaR_fit(rec_fit_Tlep_v4, gen_Tlep_v4), weight);

    LorentzVector rec_fit_Bhad_v4;
    rec_fit_Bhad_v4.SetXYZT(rec_fit_Bhad.X(), rec_fit_Bhad.Y(), rec_fit_Bhad.Z(), rec_fit_Bhad.T());
    const LorentzVector gen_Bhad_v4 = ttgen->BHad().v4();
    H1("rec_fit_bhad__gen_DR")        ->Fill(deltaR_fit(rec_fit_Bhad_v4, gen_Bhad_v4), weight);

    LorentzVector rec_fit_Blep_v4;
    rec_fit_Blep_v4.SetXYZT(rec_fit_Blep.X(), rec_fit_Blep.Y(), rec_fit_Blep.Z(), rec_fit_Blep.T());
    const LorentzVector gen_Blep_v4 = ttgen->BLep().v4();
    H1("rec_fit_blep__gen_DR")        ->Fill(deltaR_fit(rec_fit_Blep_v4, gen_Blep_v4), weight);

    LorentzVector rec_fit_Lepton_v4;
    rec_fit_Lepton_v4.SetXYZT(rec_fit_Lepton.X(), rec_fit_Lepton.Y(), rec_fit_Lepton.Z(), rec_fit_Lepton.T());
    if(gen_lep_pt!=-1.){
      const LorentzVector gen_Lepton_v4 = ttgen->ChargedLepton().v4();
      H1("rec_fit_lep__gen_DR")        ->Fill(deltaR_fit(rec_fit_Lepton_v4, gen_Lepton_v4), weight);
    }
    LorentzVector rec_fit_Neu_v4;
    rec_fit_Neu_v4.SetXYZT(rec_fit_Neu.X(), rec_fit_Neu.Y(), rec_fit_Neu.Z(), rec_fit_Neu.T());
    const LorentzVector gen_Neu_v4 = ttgen->Neutrino().v4();
    H1("rec_fit_neu__gen_DR")        ->Fill(deltaR_fit(rec_fit_Neu_v4, gen_Neu_v4), weight);
    H1("rec_fit_tops__DR") ->Fill(deltaR_fit(rec_fit_Tlep_v4, rec_fit_Thad_v4), weight);

    LorentzVector input_thad_v4;
    input_thad_v4.SetXYZT(input_thad.X(), input_thad.Y(), input_thad.Z(), input_thad.T());
    H1("rec_fit_thad__input_DR")  ->Fill(deltaR_fit(rec_fit_Thad_v4, input_thad_v4)      , weight);

    LorentzVector input_tlep_v4;
    input_tlep_v4.SetXYZT(input_tlep.X(), input_tlep.Y(), input_tlep.Z(), input_tlep.T());
    H1("rec_fit_tlep__input_DR")  ->Fill(deltaR_fit(rec_fit_Tlep_v4, input_tlep_v4)      , weight);

    LorentzVector input_whad_v4;
    input_whad_v4.SetXYZT(input_whad.X(), input_whad.Y(), input_whad.Z(), input_whad.T());
    H1("rec_fit_Whad__input_DR")  ->Fill(deltaR_fit(rec_fit_Whad_v4, input_whad_v4)      , weight);

    LorentzVector input_wlep_v4;
    input_wlep_v4.SetXYZT(input_wlep.X(), input_wlep.Y(), input_wlep.Z(), input_wlep.T());
    H1("rec_fit_Wlep__input_DR")  ->Fill(deltaR_fit(rec_fit_Wlep_v4, input_wlep_v4)      , weight);

    LorentzVector input_lepton_v4;
    input_lepton_v4.SetXYZT(input_lepton.X(), input_lepton.Y(), input_lepton.Z(), input_lepton.T());
    H1("rec_fit_lep__input_DR")  ->Fill(deltaR_fit(rec_fit_Lepton_v4, input_lepton_v4)      , weight);

    LorentzVector input_neutrino_v4;
    input_neutrino_v4.SetXYZT(input_neutrino.X(), input_neutrino.Y(), input_neutrino.Z(), input_neutrino.T());
    H1("rec_fit_neu__input_DR")  ->Fill(deltaR_fit(rec_fit_Neu_v4, input_neutrino_v4)      , weight);

    LorentzVector input_bhad_v4;
    input_bhad_v4.SetXYZT(input_jetbh.X(), input_jetbh.Y(), input_jetbh.Z(), input_jetbh.T());
    H1("rec_fit_bhad__input_DR")        ->Fill(deltaR_fit(rec_fit_Bhad_v4, input_bhad_v4), weight);

    LorentzVector input_blep_v4;
    input_blep_v4.SetXYZT(input_jetbl.X(), input_jetbl.Y(), input_jetbl.Z(), input_jetbl.T());
    H1("rec_fit_blep__input_DR")        ->Fill(deltaR_fit(rec_fit_Blep_v4, input_blep_v4), weight);

    LorentzVector rec_fit_Jetq1_v4;
    rec_fit_Jetq1_v4.SetXYZT(rec_fit_jetq1.X(), rec_fit_jetq1.Y(), rec_fit_jetq1.Z(), rec_fit_jetq1.T());
    LorentzVector input_jetq1_v4;
    input_jetq1_v4.SetXYZT(input_jetq1.X(), input_jetq1.Y(), input_jetq1.Z(), input_jetq1.T());
    H1("rec_fit_no_b_jets__input_DR")->Fill(deltaR_fit(rec_fit_Jetq1_v4, input_jetq1_v4), weight);
    LorentzVector rec_fit_Jetq2_v4;
    rec_fit_Jetq2_v4.SetXYZT(rec_fit_jetq2.X(), rec_fit_jetq2.Y(), rec_fit_jetq2.Z(), rec_fit_jetq2.T());
    LorentzVector input_jetq2_v4;
    input_jetq2_v4.SetXYZT(input_jetq2.X(), input_jetq2.Y(), input_jetq2.Z(), input_jetq2.T());
    H1("rec_fit_no_b_jets__input_DR")->Fill(deltaR_fit(rec_fit_Jetq2_v4, input_jetq2_v4), weight);
   }
 
  rec_fit_lep_pt_ = rec_fit_Lepton.Pt();
  rec_fit_lep_Et_ = rec_fit_Lepton.Et();
  rec_fit_lep_eta_ = rec_fit_Lepton.Eta();
  rec_fit_lep_phi_ = rec_fit_Lepton.Phi();
  rec_fit_lep_px_ = rec_fit_Lepton.Px();
  rec_fit_lep_py_ = rec_fit_Lepton.Py();
  rec_fit_lep_pz_ = rec_fit_Lepton.Pz();
    H1("rec_fit_lep__pt")       ->Fill(rec_fit_lep_pt_       , weight);
    H1("rec_fit_lep__eta")      ->Fill(rec_fit_lep_eta_      , weight);
    H1("rec_fit_lep__phi")      ->Fill(rec_fit_lep_phi_      , weight);
    H1("rec_fit_lep__px")       ->Fill(rec_fit_lep_px_       , weight);
    H1("rec_fit_lep__py")       ->Fill(rec_fit_lep_py_       , weight);
    H1("rec_fit_lep__pz")       ->Fill(rec_fit_lep_pz_       , weight);
    if(ttljets){
  H1("rec_fit_lep__gen_Dpt")       ->Fill(rec_fit_lep_pt_ -gen_lep_pt                      , weight);
  H2("rec_fit_lep__gen_Dpt__vs__gen_pt")       ->Fill(rec_fit_lep_pt_ -gen_lep_pt, gen_lep_pt                      , weight);
  H1("rec_fit_lep__gen_Deta")      ->Fill(rec_fit_lep_eta_ - gen_lep_eta                     , weight);
  H1("rec_fit_lep__gen_Dphi")      ->Fill(delta_phi(rec_fit_lep_phi_, gen_lep_phi)         , weight);
    if(gen_lep_pt)  H1("rec_fit_lep__gen_Dpt_pct") ->Fill((rec_fit_lep_pt_ - gen_lep_pt) /fabs(gen_lep_pt) , weight);
    if(gen_lep_eta) H1("rec_fit_lep__gen_Deta_pct")->Fill((rec_fit_lep_eta_ - gen_lep_eta)/fabs(gen_lep_eta), weight);
    }
    H1("rec_fit_lep__input_Dpt")       ->Fill(rec_fit_lep_pt_ - input_lepton.Pt()                      , weight);
    H2("rec_fit_lep__input_Dpt__vs__input_pt")       ->Fill(rec_fit_lep_pt_ - input_lepton.Pt(), input_lepton.Pt()                      , weight);
    H1("rec_fit_lep__input_Dpt_pct")       ->Fill((rec_fit_lep_pt_ - input_lepton.Pt())/fabs(input_lepton.Pt())                      , weight);
    H1("rec_fit_lep__input_DEt")       ->Fill(rec_fit_lep_Et_ - input_lepton.Et()                      , weight);
    H1("rec_fit_lep__input_DEt_pct")       ->Fill((rec_fit_lep_Et_ - input_lepton.Et())/fabs(input_lepton.Et())                      , weight);
    H1("rec_fit_lep__input_Deta")      ->Fill(rec_fit_lep_eta_ - input_lepton.Eta()                     , weight);
    H1("rec_fit_lep__input_Dphi")      ->Fill(delta_phi(rec_fit_lep_phi_, input_lepton.Phi())         , weight);
   H1("input_lep__gen_Dpt")       ->Fill(input_lepton.Pt()  -   gen_lep_pt                 , weight);
   H2("input_lep__gen_Dpt__vs__gen_pt")       ->Fill( input_lepton.Pt() - gen_lep_pt,  gen_lep_pt                     , weight);
   H1("gen_lep__input_Deta")      ->Fill(gen_lep_eta - input_lepton.Eta()                     , weight);
   H1("gen_lep__input_Dphi")      ->Fill(delta_phi(gen_lep_phi, input_lepton.Phi())         , weight);
    float  dR_combined(0.);
    rec_fit_neu_pt_ = rec_fit_Neu.Pt();
    rec_fit_neu_Et_ = rec_fit_Neu.Et();
    rec_fit_neu_phi_ = rec_fit_Neu.Phi();
    rec_fit_neu_px_ = rec_fit_Neu.Px();
    rec_fit_neu_py_ = rec_fit_Neu.Py();
    rec_fit_neu_pz_ = rec_fit_Neu.Pz();
    H1("rec_fit_neu__pt")       ->Fill(rec_fit_neu_pt_       , weight);
    H1("rec_fit_neu__phi")      ->Fill(rec_fit_neu_phi_      , weight);
    H1("rec_fit_neu__px")       ->Fill(rec_fit_neu_px_       , weight);
    H1("rec_fit_neu__py")       ->Fill(rec_fit_neu_py_       , weight);
    H1("rec_fit_neu__pz")       ->Fill(rec_fit_neu_pz_       , weight);
    if(ttljets){
      H1("rec_fit_neu__gen_Dpt")       ->Fill(rec_fit_neu_pt_ - gen_neu_pt                  , weight);
      H2("rec_fit_neu__gen_Dpt__vs__gen_pt")       ->Fill(rec_fit_neu_pt_ - gen_neu_pt ,   gen_neu_pt             , weight);
      H1("rec_fit_neu__gen_Dphi")      ->Fill(delta_phi(rec_fit_neu_phi_, gen_neu_phi)    , weight);
      H1("rec_fit_neu__gen_Dpx")       ->Fill(rec_fit_neu_px_-gen_neu_px                  , weight);
      H1("rec_fit_neu__gen_Dpy")       ->Fill(rec_fit_neu_py_-gen_neu_py                  , weight);
      H1("rec_fit_neu__gen_Dpz")       ->Fill(rec_fit_neu_pz_-gen_neu_pz                  , weight);
      if(gen_neu_pt) H1("rec_fit_neu__gen_Dpt_pct")->Fill((rec_fit_neu_pt_-gen_neu_pt)/fabs(gen_neu_pt), weight);
      if(gen_neu_px) H1("rec_fit_neu__gen_Dpx_pct")->Fill((rec_fit_neu_px_-gen_neu_px)/fabs(gen_neu_px), weight);
      if(gen_neu_py) H1("rec_fit_neu__gen_Dpy_pct")->Fill((rec_fit_neu_py_-gen_neu_py)/fabs(gen_neu_py), weight);
      if(gen_neu_pz) H1("rec_fit_neu__gen_Dpz_pct")->Fill((rec_fit_neu_pz_-gen_neu_pz)/fabs(gen_neu_pz), weight);
    }
    H1("rec_fit_neu__input_Dpt")       ->Fill(rec_fit_neu_pt_ - input_neutrino.Pt()                  , weight);
    H1("rec_fit_neu__input_Dpt_pct")       ->Fill((rec_fit_neu_pt_ - input_neutrino.Pt())/fabs(input_neutrino.Pt())                  , weight);
    H1("rec_fit_neu__input_DEt")       ->Fill(rec_fit_neu_Et_ - input_neutrino.Et()                  , weight);
    H1("rec_fit_neu__input_DEt_pct")       ->Fill((rec_fit_neu_Et_ - input_neutrino.Et())/fabs(input_neutrino.Et())                  , weight);
    H1("rec_fit_neu__input_Dphi")      ->Fill(delta_phi(rec_fit_neu_phi_, input_neutrino.Phi())    , weight);
    H1("rec_fit_neu__input_Deta")      ->Fill(rec_fit_Neu.Eta() - input_neutrino.Eta()    , weight);
    H1("rec_fit_neu__input_Dpx")       ->Fill(rec_fit_neu_px_ - input_neutrino.Px()                  , weight);
    H1("rec_fit_neu__input_Dpy")       ->Fill(rec_fit_neu_py_-input_neutrino.Py()                  , weight);
    H1("rec_fit_neu__input_Dpz")       ->Fill(rec_fit_neu_pz_-input_neutrino.Pz()                  , weight);
    H1("input_neu__gen_Dpt")       ->Fill( input_neutrino.Pt() -  gen_neu_pt               , weight);
    H2("input_neu__gen_Dpt__vs__gen_pt")       ->Fill( input_neutrino.Pt() - gen_neu_pt,   gen_neu_pt             , weight);
    H1("gen_neu__input_Dphi")      ->Fill(delta_phi(gen_neu_phi, input_neutrino.Phi())    , weight);
    H1("gen_neu__input_Dpx")       ->Fill(gen_neu_px - input_neutrino.Px()                  , weight);
    H1("gen_neu__input_Dpy")       ->Fill(gen_neu_py - input_neutrino.Py()                  , weight);
    H1("gen_neu__input_Dpz")       ->Fill(gen_neu_pz - input_neutrino.Pz()                  , weight);
    // blep
    rec_fit_blep_pt_ =  rec_fit_Blep.Pt();
    rec_fit_blep_Et_ =  rec_fit_Blep.Et();
    rec_fit_blep_eta_ = rec_fit_Blep.Eta();
    rec_fit_blep_phi_ = rec_fit_Blep.Phi();
    rec_fit_blep_px_ =  rec_fit_Blep.Px();
    rec_fit_blep_py_ =  rec_fit_Blep.Py();
    rec_fit_blep_pz_ =  rec_fit_Blep.Pz();
    H1("rec_fit_blep__pt")       ->Fill(rec_fit_blep_pt_       , weight);
    H1("rec_fit_blep__eta")      ->Fill(rec_fit_blep_eta_      , weight);
    H1("rec_fit_blep__phi")      ->Fill(rec_fit_blep_phi_      , weight);
    H1("rec_fit_blep__px")       ->Fill(rec_fit_blep_px_       , weight);
    H1("rec_fit_blep__py")       ->Fill(rec_fit_blep_py_       , weight);
    H1("rec_fit_blep__pz")       ->Fill(rec_fit_blep_pz_       , weight);
    if(ttljets){
      H1("rec_fit_blep__gen_Dpt")       ->Fill(rec_fit_blep_pt_ - gen_blep_pt             , weight);
      H2("rec_fit_blep__gen_Dpt__vs__gen_pt")       ->Fill(rec_fit_blep_pt_ - gen_blep_pt,    gen_blep_pt          , weight);
      H1("rec_fit_blep__gen_Deta")      ->Fill(rec_fit_blep_eta_ - gen_blep_eta            , weight);
      H1("rec_fit_blep__gen_Dphi")      ->Fill(delta_phi(rec_fit_blep_phi_, gen_blep_phi), weight);
      if(gen_blep_pt)  H1("rec_fit_blep__gen_Dpt_pct") ->Fill((rec_fit_blep_pt_ -gen_blep_pt) /fabs(gen_blep_pt) , weight);
    }
    H1("rec_fit_blep__input_Dpt")       ->Fill(rec_fit_blep_pt_ - input_jetbl.Pt()             , weight);
    H1("rec_fit_blep__input_Dpt_pct")       ->Fill((rec_fit_blep_pt_ - input_jetbl.Pt())/fabs(input_jetbl.Pt())             , weight);
    H1("rec_fit_blep__input_DEt")       ->Fill(rec_fit_blep_Et_ - input_jetbl.Et()             , weight);
    H1("rec_fit_blep__input_DEt_pct")       ->Fill((rec_fit_blep_Et_ - input_jetbl.Et())/fabs(input_jetbl.Et())             , weight);
    H1("rec_fit_blep__input_Deta")      ->Fill(rec_fit_blep_eta_ - input_jetbl.Eta()            , weight);
    H1("rec_fit_blep__input_Dphi")      ->Fill(delta_phi(rec_fit_blep_phi_, input_jetbl.Phi()), weight);

    H1("input_blep__gen_Dpt")       ->Fill(input_jetbl.Pt() -   gen_blep_pt         , weight);
    H2("input_blep__gen_Dpt__vs__gen_pt")       ->Fill(input_jetbl.Pt() - gen_blep_pt,    gen_blep_pt          , weight);
    H1("gen_blep__input_Deta")      ->Fill(gen_blep_eta - input_jetbl.Eta()            , weight);
    H1("gen_blep__input_Dphi")      ->Fill(delta_phi(gen_blep_phi, input_jetbl.Phi()), weight);
   // bhad
    rec_fit_bhad_pt_ = rec_fit_Bhad.Pt();
    rec_fit_bhad_Et_ = rec_fit_Bhad.Et();
    rec_fit_bhad_eta_ = rec_fit_Bhad.Eta();
    rec_fit_bhad_phi_ = rec_fit_Bhad.Phi();
    rec_fit_bhad_px_ = rec_fit_Bhad.Px();
    rec_fit_bhad_py_ = rec_fit_Bhad.Py();
    rec_fit_bhad_pz_ = rec_fit_Bhad.Pz();
    H1("rec_fit_bhad__pt")       ->Fill(rec_fit_bhad_pt_       , weight);
    H1("rec_fit_bhad__eta")      ->Fill(rec_fit_bhad_eta_      , weight);
    H1("rec_fit_bhad__phi")      ->Fill(rec_fit_bhad_phi_      , weight);
    H1("rec_fit_bhad__px")       ->Fill(rec_fit_bhad_px_       , weight);
    H1("rec_fit_bhad__py")       ->Fill(rec_fit_bhad_py_       , weight);
    H1("rec_fit_bhad__pz")       ->Fill(rec_fit_bhad_pz_       , weight);
    if(ttljets){
      H1("rec_fit_bhad__gen_Dpt")       ->Fill(rec_fit_bhad_pt_ - gen_bhad_pt             , weight);
      H2("rec_fit_bhad__gen_Dpt__vs__gen_pt")       ->Fill(rec_fit_bhad_pt_ - gen_bhad_pt   ,    gen_bhad_pt      , weight);
      H1("rec_fit_bhad__gen_Deta")      ->Fill(rec_fit_bhad_eta_ - gen_bhad_eta            , weight);
      H1("rec_fit_bhad__gen_Dphi")      ->Fill(delta_phi(rec_fit_bhad_phi_, gen_bhad_phi), weight);
      if(gen_bhad_pt)  H1("rec_fit_bhad__gen_Dpt_pct") ->Fill((rec_fit_bhad_pt_ -gen_bhad_pt) /fabs(gen_bhad_pt) , weight);
      if(gen_bhad_eta) H1("rec_fit_bhad__gen_Deta_pct")->Fill((rec_fit_bhad_eta_ -gen_bhad_eta)/fabs(gen_bhad_eta), weight);
    }
    H1("rec_fit_bhad__input_Dpt")       ->Fill(rec_fit_bhad_pt_ - input_jetbh.Pt()             , weight);
    H1("rec_fit_bhad__input_Dpt_pct")       ->Fill((rec_fit_bhad_pt_ - input_jetbh.Pt())/fabs(input_jetbh.Pt())             , weight);
    H1("rec_fit_bhad__input_DEt")       ->Fill(rec_fit_bhad_Et_ - input_jetbh.Et()             , weight);
    H1("rec_fit_bhad__input_DEt_pct")       ->Fill((rec_fit_bhad_Et_ - input_jetbh.Et())/fabs(input_jetbh.Et())             , weight);
    H1("rec_fit_bhad__input_Deta")      ->Fill(rec_fit_bhad_eta_ - input_jetbh.Eta()            , weight);
    H1("rec_fit_bhad__input_Dphi")      ->Fill(delta_phi(rec_fit_bhad_phi_, input_jetbh.Phi()), weight);

    H1("input_bhad__gen_Dpt")       ->Fill(input_jetbh.Pt() - gen_bhad_pt            , weight);
    H2("input_bhad__gen_Dpt__vs__gen_pt")       ->Fill(input_jetbh.Pt() - gen_bhad_pt  ,    gen_bhad_pt      , weight);
    H1("gen_bhad__input_Deta")      ->Fill(gen_bhad_eta - input_jetbh.Eta()            , weight);
    H1("gen_bhad__input_Dphi")      ->Fill(delta_phi(gen_bhad_phi, input_jetbh.Phi()), weight);
   // no b jets
    H1("rec_fit_jet_l1__pt")          ->Fill(rec_fit_jetq1.Pt()       , weight);
    H1("rec_fit_jet_l1__eta")         ->Fill(rec_fit_jetq1.Eta()       , weight);
    H1("rec_fit_jet_l1__phi")         ->Fill(rec_fit_jetq1.Phi()       , weight);
    H1("rec_fit_jet_l2__pt")          ->Fill(rec_fit_jetq2.Pt()       , weight);
    H1("rec_fit_jet_l2__eta")         ->Fill(rec_fit_jetq2.Eta()       , weight);
    H1("rec_fit_jet_l2__phi")         ->Fill(rec_fit_jetq2.Phi()       , weight);
    H1("rec_fit_jet_l1__l2__Dpt")      ->Fill((rec_fit_jetq1.Pt() -  rec_fit_jetq2.Pt())        , weight);
    H1("rec_fit_jet_l1__l2__Deta")     ->Fill((rec_fit_jetq1.Eta() - rec_fit_jetq2.Eta())       , weight);
    H1("rec_fit_jet_l1__l2__Dphi")     ->Fill(delta_phi(rec_fit_jetq1.Phi(), rec_fit_jetq2.Phi())       , weight);
    H1("rec_fit_jet_l1__input__Dpt")      ->Fill((rec_fit_jetq1.Pt() -   input_jetq1.Pt())        , weight);
    H1("rec_fit_jet_l1__input__Dpt_pct")  ->Fill((rec_fit_jetq1.Pt() -   input_jetq1.Pt())/fabs(input_jetq1.Pt())        , weight);
    H1("rec_fit_jet_l1__input__Deta")     ->Fill((rec_fit_jetq1.Eta() -  input_jetq1.Eta())       , weight);
    H1("rec_fit_jet_l1__input__Dphi")     ->Fill(delta_phi(rec_fit_jetq1.Phi(),  input_jetq1.Phi())       , weight);
    H1("rec_fit_jet_l2__input__Dpt")      ->Fill((rec_fit_jetq2.Pt() -   input_jetq2.Pt())        , weight);
    H1("rec_fit_jet_l2__input__Dpt_pct")  ->Fill((rec_fit_jetq2.Pt() -   input_jetq2.Pt())/fabs(input_jetq2.Pt())        , weight);
    H1("rec_fit_jet_l2__input__Deta")     ->Fill((rec_fit_jetq2.Eta() -  input_jetq2.Eta())       , weight);
    H1("rec_fit_jet_l2__input__Dphi")     ->Fill(delta_phi(rec_fit_jetq2.Phi(),  input_jetq2.Phi())       , weight);
    H1("rec_fit_no_b_jets__pt")       ->Fill(rec_fit_jetq1.Pt()       , weight);
    H1("rec_fit_no_b_jets__eta")      ->Fill(rec_fit_jetq1.Eta()      , weight);
    H1("rec_fit_no_b_jets__phi")      ->Fill(rec_fit_jetq1.Phi()      , weight);
    H1("rec_fit_no_b_jets__px")       ->Fill(rec_fit_jetq1.Px()       , weight);
    H1("rec_fit_no_b_jets__py")       ->Fill(rec_fit_jetq1.Py()       , weight);
    H1("rec_fit_no_b_jets__pz")       ->Fill(rec_fit_jetq1.Pz()       , weight);
    H1("rec_fit_no_b_jets__pt")       ->Fill(rec_fit_jetq2.Pt()       , weight);
    H1("rec_fit_no_b_jets__eta")      ->Fill(rec_fit_jetq2.Eta()      , weight);
    H1("rec_fit_no_b_jets__phi")      ->Fill(rec_fit_jetq2.Phi()      , weight);
    H1("rec_fit_no_b_jets__px")       ->Fill(rec_fit_jetq2.Px()       , weight);
    H1("rec_fit_no_b_jets__py")       ->Fill(rec_fit_jetq2.Py()       , weight);
    H1("rec_fit_no_b_jets__pz")       ->Fill(rec_fit_jetq2.Pz()       , weight);
    if(ttljets){
      float dpt_q1_1 = rec_fit_jetq1.Pt() - ttgen->Q1().v4().Pt();
      float dpt_q1_2 = rec_fit_jetq1.Pt() - ttgen->Q2().v4().Pt();
      if(fabs(dpt_q1_1)<fabs(dpt_q1_2)){
	H1("rec_fit_no_b_jets__gen_Dpt")       ->Fill(dpt_q1_1             , weight);
	H2("rec_fit_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill(dpt_q1_1 ,  ttgen->Q1().v4().Pt()          , weight);
      	H1("rec_fit_no_b_jets__gen_Dpt")       ->Fill((rec_fit_jetq2.Pt() - ttgen->Q2().v4().Pt())             , weight);
      	H2("rec_fit_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill((rec_fit_jetq2.Pt() - ttgen->Q2().v4().Pt())  , ttgen->Q2().v4().Pt()          , weight);
	if(gen_bhad_pt)  H1("rec_fit_no_b_jets__gen_Dpt_pct") ->Fill((dpt_q1_1) /fabs(ttgen->Q1().v4().Pt()) , weight);
	if(gen_bhad_pt)  H1("rec_fit_no_b_jets__gen_Dpt_pct") ->Fill((rec_fit_jetq2.Pt() - ttgen->Q2().v4().Pt()) /fabs(ttgen->Q2().v4().Pt()) , weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_Dpt")       ->Fill(dpt_q1_2             , weight);
	H2("rec_fit_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill(dpt_q1_2 ,  ttgen->Q2().v4().Pt()          , weight);
      	H1("rec_fit_no_b_jets__gen_Dpt")       ->Fill((rec_fit_jetq2.Pt() - ttgen->Q1().v4().Pt())             , weight);
      	H2("rec_fit_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill((rec_fit_jetq2.Pt() - ttgen->Q1().v4().Pt())  , ttgen->Q1().v4().Pt()          , weight);
	if(ttgen->Q2().v4().Pt())  H1("rec_fit_no_b_jets__gen_Dpt_pct") ->Fill((dpt_q1_2) /fabs(ttgen->Q2().v4().Pt()) , weight);
	if(ttgen->Q1().v4().Pt())  H1("rec_fit_no_b_jets__gen_Dpt_pct") ->Fill((rec_fit_jetq2.Pt() - ttgen->Q1().v4().Pt()) /fabs(ttgen->Q1().v4().Pt()) , weight);
      }
      float deta_q1_1 = rec_fit_jetq1.Eta() - ttgen->Q1().v4().Eta();
      float deta_q1_2 = rec_fit_jetq1.Eta() - ttgen->Q2().v4().Eta();
      if(fabs(deta_q1_1)<fabs(deta_q1_2)){
	H1("rec_fit_no_b_jets__gen_Deta")      ->Fill(deta_q1_1            , weight);
	H1("rec_fit_no_b_jets__gen_Deta")      ->Fill(rec_fit_jetq2.Eta() - ttgen->Q2().v4().Eta()            , weight);
	if(ttgen->Q1().v4().Eta()) H1("rec_fit_no_b_jets__gen_Deta_pct")->Fill(deta_q1_1/fabs(ttgen->Q1().v4().Eta()), weight);
	if(ttgen->Q2().v4().Eta()) H1("rec_fit_no_b_jets__gen_Deta_pct")->Fill(rec_fit_jetq1.Eta() - ttgen->Q2().v4().Eta()/fabs(ttgen->Q2().v4().Eta()), weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_Deta")      ->Fill(deta_q1_2            , weight);
	H1("rec_fit_no_b_jets__gen_Deta")      ->Fill(rec_fit_jetq2.Eta() - ttgen->Q1().v4().Eta()            , weight);
	if(ttgen->Q2().v4().Eta()) H1("rec_fit_no_b_jets__gen_Deta_pct")->Fill(deta_q1_2/fabs(ttgen->Q2().v4().Eta()), weight);
	if(ttgen->Q1().v4().Eta()) H1("rec_fit_no_b_jets__gen_Deta_pct")->Fill(rec_fit_jetq2.Eta() - ttgen->Q1().v4().Eta()/fabs(ttgen->Q1().v4().Eta()), weight);
      } 
      float dphi_q1_1 = delta_phi(rec_fit_jetq1.Phi(), ttgen->Q1().v4().Phi());
      float dphi_q1_2 = delta_phi(rec_fit_jetq1.Phi(), ttgen->Q2().v4().Phi());
      if(fabs(dphi_q1_1)<fabs(dphi_q1_2)){
	H1("rec_fit_no_b_jets__gen_Dphi")      ->Fill(dphi_q1_1, weight);
	H1("rec_fit_no_b_jets__gen_Dphi")      ->Fill(delta_phi(rec_fit_jetq2.Phi(), ttgen->Q2().v4().Phi()), weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_Dphi")      ->Fill(dphi_q1_2, weight);
	H1("rec_fit_no_b_jets__gen_Dphi")      ->Fill(delta_phi(rec_fit_jetq2.Phi(), ttgen->Q1().v4().Phi()), weight);
      }
      float dpx_q1_1 = rec_fit_jetq1.Px() - ttgen->Q1().v4().Px();
      float dpx_q1_2 = rec_fit_jetq1.Px() - ttgen->Q2().v4().Px();
      if(fabs(dpx_q1_1)<fabs(dpx_q1_2)){
	H1("rec_fit_no_b_jets__gen_Dpx")       ->Fill(dpx_q1_1             , weight);
	H1("rec_fit_no_b_jets__gen_Dpx")       ->Fill( rec_fit_jetq2.Px() - ttgen->Q2().v4().Px()             , weight);
	if(ttgen->Q1().v4().Px())  H1("rec_fit_no_b_jets__gen_Dpx_pct") ->Fill(dpx_q1_1 /fabs(ttgen->Q1().v4().Px()) , weight);
	if(ttgen->Q2().v4().Px())  H1("rec_fit_no_b_jets__gen_Dpx_pct") ->Fill(rec_fit_jetq2.Px() - ttgen->Q2().v4().Px() /fabs(ttgen->Q2().v4().Px()) , weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_Dpx")       ->Fill(dpx_q1_2             , weight);
	H1("rec_fit_no_b_jets__gen_Dpx")       ->Fill( rec_fit_jetq2.Px() - ttgen->Q1().v4().Px()             , weight);
	if(ttgen->Q2().v4().Px())  H1("rec_fit_no_b_jets__gen_Dpx_pct") ->Fill(dpx_q1_2 /fabs(ttgen->Q2().v4().Px()) , weight);
	if(ttgen->Q1().v4().Px())  H1("rec_fit_no_b_jets__gen_Dpx_pct") ->Fill(rec_fit_jetq2.Px() - ttgen->Q1().v4().Px() /fabs(ttgen->Q1().v4().Px()) , weight);
      }
      float dpy_q1_1 = rec_fit_jetq1.Py() - ttgen->Q1().v4().Py();
      float dpy_q1_2 = rec_fit_jetq1.Py() - ttgen->Q2().v4().Py();
      if(fabs(dpy_q1_1)<fabs(dpy_q1_2)){
	H1("rec_fit_no_b_jets__gen_Dpy")       ->Fill(dpy_q1_1             , weight);
	H1("rec_fit_no_b_jets__gen_Dpy")       ->Fill( rec_fit_jetq2.Py() - ttgen->Q2().v4().Py()             , weight);
	if(ttgen->Q1().v4().Py())  H1("rec_fit_no_b_jets__gen_Dpy_pct") ->Fill(dpy_q1_1 /fabs(ttgen->Q1().v4().Py()) , weight);
	if(ttgen->Q2().v4().Py())  H1("rec_fit_no_b_jets__gen_Dpy_pct") ->Fill(rec_fit_jetq2.Py() - ttgen->Q2().v4().Py() /fabs(ttgen->Q2().v4().Py()) , weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_Dpy")       ->Fill(dpy_q1_2             , weight);
	H1("rec_fit_no_b_jets__gen_Dpy")       ->Fill( rec_fit_jetq2.Py() - ttgen->Q1().v4().Py()             , weight);
	if(ttgen->Q2().v4().Py())  H1("rec_fit_no_b_jets__gen_Dpy_pct") ->Fill(dpy_q1_2 /fabs(ttgen->Q2().v4().Py()) , weight);
	if(ttgen->Q1().v4().Py())  H1("rec_fit_no_b_jets__gen_Dpy_pct") ->Fill(rec_fit_jetq2.Py() - ttgen->Q1().v4().Py() /fabs(ttgen->Q1().v4().Py()) , weight);
      }
      float dpz_q1_1 = rec_fit_jetq1.Pz() - ttgen->Q1().v4().Pz();
      float dpz_q1_2 = rec_fit_jetq1.Pz() - ttgen->Q2().v4().Pz();
      if(fabs(dpz_q1_1)<fabs(dpz_q1_2)){
	H1("rec_fit_no_b_jets__gen_Dpz")       ->Fill(dpz_q1_1             , weight);
	H1("rec_fit_no_b_jets__gen_Dpz")       ->Fill( rec_fit_jetq2.Pz() - ttgen->Q2().v4().Pz()             , weight);
	if(ttgen->Q1().v4().Pz())  H1("rec_fit_no_b_jets__gen_Dpz_pct") ->Fill(dpz_q1_1 /fabs(ttgen->Q1().v4().Pz()) , weight);
	if(ttgen->Q2().v4().Pz())  H1("rec_fit_no_b_jets__gen_Dpz_pct") ->Fill(rec_fit_jetq2.Pz() - ttgen->Q2().v4().Pz() /fabs(ttgen->Q2().v4().Pz()) , weight);
      }
     else{
	H1("rec_fit_no_b_jets__gen_Dpz")       ->Fill(dpz_q1_2             , weight);
	H1("rec_fit_no_b_jets__gen_Dpz")       ->Fill( rec_fit_jetq2.Pz() - ttgen->Q1().v4().Pz()             , weight);
	if(ttgen->Q2().v4().Pz())  H1("rec_fit_no_b_jets__gen_Dpz_pct") ->Fill(dpz_q1_2 /fabs(ttgen->Q1().v4().Pz()) , weight);
	if(ttgen->Q1().v4().Pz())  H1("rec_fit_no_b_jets__gen_Dpz_pct") ->Fill(rec_fit_jetq2.Pz() - ttgen->Q1().v4().Pz() /fabs(ttgen->Q1().v4().Pz()) , weight);
      }
      LorentzVector rec_fit_jetq1_lv;
      rec_fit_jetq1_lv.SetXYZT(rec_fit_jetq1.X(), rec_fit_jetq1.Y(), rec_fit_jetq1.Z(), rec_fit_jetq1.T());
      LorentzVector rec_fit_jetq2_lv;
      rec_fit_jetq2_lv.SetXYZT(rec_fit_jetq2.X(), rec_fit_jetq2.Y(), rec_fit_jetq2.Z(), rec_fit_jetq2.T());

      float dR_q1_1 = deltaR_fit(rec_fit_jetq1_lv, ttgen->Q1().v4());
      float dR_q1_2 = deltaR_fit(rec_fit_jetq1_lv, ttgen->Q2().v4());
      if(fabs(dR_q1_1)<fabs(dR_q1_2)){
	H1("rec_fit_no_b_jets__gen_DR")        ->Fill(dR_q1_1, weight);
	H1("rec_fit_no_b_jets__gen_DR")        ->Fill(deltaR_fit(rec_fit_jetq2_lv, ttgen->Q2().v4()), weight);
      }
      else{
	H1("rec_fit_no_b_jets__gen_DR")        ->Fill(dR_q1_2, weight);
	H1("rec_fit_no_b_jets__gen_DR")        ->Fill(deltaR_fit(rec_fit_jetq2_lv, ttgen->Q1().v4()), weight);
      }
      H1("input_no_b_jets__gen_Dpt")       ->Fill(input_jetq1.Pt() - ttgen->Q1().v4().Pt()       , weight);
      H1("input_no_b_jets__gen_Dpt")       ->Fill( input_jetq2.Pt() -  ttgen->Q2().v4().Pt()     , weight);
      H2("input_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill( input_jetq1.Pt() - ttgen->Q1().v4().Pt(),  ttgen->Q1().v4().Pt()      , weight);
      H2("input_no_b_jets__gen_Dpt__vs__gen_pt")       ->Fill(input_jetq2.Pt() - ttgen->Q2().v4().Pt(),  ttgen->Q2().v4().Pt()      , weight);
      H1("gen_no_b_jets__input_Deta")       ->Fill(ttgen->Q1().v4().Eta() - input_jetq1.Eta()       , weight);
      H1("gen_no_b_jets__input_Deta")       ->Fill(ttgen->Q2().v4().Eta() - input_jetq2.Eta()       , weight);
    H1("gen_no_b_jets__input_Dphi")      ->Fill(delta_phi(ttgen->Q1().v4().Phi(), input_jetq1.Phi()), weight);
    H1("gen_no_b_jets__input_Dphi")      ->Fill(delta_phi(ttgen->Q2().v4().Phi(), input_jetq2.Phi()), weight);
    }
    H1("rec_fit_no_b_jets__input_p")->Fill(input_jetq1.P(), weight);
    H1("rec_fit_no_b_jets__input_p")->Fill(input_jetq2.P(), weight);
    H1("rec_fit_no_b_jets__input_Dp")       ->Fill((rec_fit_jetq1.P() - input_jetq1.P())       , weight);
    H1("rec_fit_no_b_jets__input_Dp")       ->Fill((rec_fit_jetq2.P() - input_jetq2.P())       , weight);

    H1("rec_fit_no_b_jets__input_Dpt")       ->Fill((rec_fit_jetq1.Pt() - input_jetq1.Pt())       , weight);
    H1("rec_fit_no_b_jets__input_Dpt")       ->Fill((rec_fit_jetq2.Pt() - input_jetq2.Pt())       , weight);
    H1("rec_fit_no_b_jets__input_DEt")       ->Fill((rec_fit_jetq1.Et() - input_jetq1.Et())       , weight);
    H1("rec_fit_no_b_jets__input_DEt")       ->Fill((rec_fit_jetq2.Et() - input_jetq2.Et())       , weight);
    H1("rec_fit_no_b_jets__input_Dpt_pct")   ->Fill((rec_fit_jetq1.Pt() - input_jetq1.Pt())/fabs(input_jetq1.Pt())       , weight);
    H1("rec_fit_no_b_jets__input_Dpt_pct")   ->Fill((rec_fit_jetq2.Pt() - input_jetq2.Pt())/fabs(input_jetq2.Pt())       , weight);
   H1("rec_fit_no_b_jets__input_DEt_pct")    ->Fill((rec_fit_jetq1.Et() - input_jetq1.Et())/fabs(input_jetq1.Et())       , weight);
    H1("rec_fit_no_b_jets__input_DEt_pct")   ->Fill((rec_fit_jetq2.Et() - input_jetq2.Et())/fabs(input_jetq2.Et())       , weight);
    H1("rec_fit_no_b_jets__input_Deta")      ->Fill(rec_fit_jetq1.Eta() - input_jetq1.Eta()   , weight);
    H1("rec_fit_no_b_jets__input_Deta")      ->Fill(rec_fit_jetq2.Eta() - input_jetq2.Eta()   , weight);
    H1("rec_fit_no_b_jets__input_Dphi")      ->Fill(delta_phi(rec_fit_jetq1.Phi(), input_jetq1.Phi()), weight);
    H1("rec_fit_no_b_jets__input_Dphi")      ->Fill(delta_phi(rec_fit_jetq2.Phi(), input_jetq2.Phi()), weight);
    H1("rec_fit_DR_combined")->Fill(dR_combined,weight);


     TMatrixD jet1_unc = event.get(unc_had1);
     TMatrixD jet2_unc = event.get(unc_had2);
     TMatrixD jetbh_unc = event.get(unc_hadb);
     TMatrixD jetbl_unc = event.get(unc_lepb);
     TMatrixD lepton_unc = event.get(unc_lepton);
     TMatrixD neu_unc = event.get(unc_neu);
     H1("unc_file_jet1_et")->Fill(jet1_unc(0,0));
     H1("unc_file_jet1_eta")->Fill(jet1_unc(1,1));
     H1("unc_file_jet1_phi")->Fill(jet1_unc(2,2));
     H1("unc_file_jet2_et")->Fill(jet2_unc(0,0));
     H1("unc_file_jet2_eta")->Fill(jet2_unc(1,1));
     H1("unc_file_jet2_phi")->Fill(jet2_unc(2,2));
     H1("unc_file_jetbh_et")->Fill(jetbh_unc(0,0));
     H1("unc_file_jetbh_eta")->Fill(jetbh_unc(1,1));
     H1("unc_file_jetbh_phi")->Fill(jetbh_unc(2,2));
     H1("unc_file_jetbl_et")->Fill(jetbl_unc(0,0));
     H1("unc_file_jetbl_eta")->Fill(jetbl_unc(1,1));
     H1("unc_file_jetbl_phi")->Fill(jetbl_unc(2,2));
     H1("unc_file_lepton_et")->Fill(lepton_unc(0,0));
     H1("unc_file_lepton_eta")->Fill(lepton_unc(1,1));
     H1("unc_file_lepton_phi")->Fill(lepton_unc(2,2));
     H1("unc_file_neu_et")->Fill(neu_unc(0,0));
     H1("unc_file_neu_eta")->Fill(neu_unc(1,1));
     H1("unc_file_neu_phi")->Fill(neu_unc(2,2));

}


float KinFitTestFITHists::deltaR_fit(const LorentzVector& p1, const LorentzVector p2){
    float deltaeta = p1.Eta() - p2.Eta();
    float dphi = deltaPhi(p1, p2);
    return sqrt(deltaeta * deltaeta + dphi * dphi);
}

float KinFitTestFITHists::deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2){
    float deltaeta = eta1 - eta2;
    float dphi = delta_phi(phi1, phi2);
    return sqrt(deltaeta * deltaeta + dphi * dphi);
}

float KinFitTestFITHists::delta_phi(const float phi1, const float phi2){

  float a_dphi = fabs(phi1-phi2);
  if(a_dphi > M_PI) a_dphi = 2*M_PI - a_dphi;

  return a_dphi;
}

void KinFitTestFITHists::boost_x1_to_x2CM(TLorentzVector& x1, const TLorentzVector& x2){

  if(!x2.E()) return;

  x1.Boost(-x2.Px()/x2.E(), -x2.Py()/x2.E(), -x2.Pz()/x2.E());

  return;
}

float KinFitTestFITHists::cosThetaX(const LorentzVector& tprod, const LorentzVector& top, const LorentzVector& ttbar0){

  TLorentzVector top_in_ttbarCM(   top.Px(),    top.Py(),    top.Pz(),    top.E());
  TLorentzVector tprod_in_topCM( tprod.Px(),  tprod.Py(),  tprod.Pz(),  tprod.E());
  TLorentzVector ttbar         (ttbar0.Px(), ttbar0.Py(), ttbar0.Pz(), ttbar0.E());

  boost_x1_to_x2CM(top_in_ttbarCM, ttbar);
  boost_x1_to_x2CM(tprod_in_topCM, ttbar);

  boost_x1_to_x2CM(tprod_in_topCM, top_in_ttbarCM);

  float thetaX = tprod_in_topCM.Angle(top_in_ttbarCM.Vect());

  return cos(thetaX);
}

float KinFitTestFITHists::cosThetaCS(const LorentzVector& top, const LorentzVector& ttbar0){

  const float sqrt_s = 13000.;

  TLorentzVector proton1       (         0.,          0.,   .5*sqrt_s,  .5*sqrt_s);
  TLorentzVector proton2       (         0.,          0.,  -.5*sqrt_s,  .5*sqrt_s);
  TLorentzVector top_in_ttbarCM(   top.Px(),    top.Py(),    top.Pz(),    top.E());
  TLorentzVector ttbar         (ttbar0.Px(), ttbar0.Py(), ttbar0.Pz(), ttbar0.E());

  boost_x1_to_x2CM(proton1       , ttbar);
  boost_x1_to_x2CM(proton2       , ttbar);
  boost_x1_to_x2CM(top_in_ttbarCM, ttbar);

  float thetaCS = top_in_ttbarCM.Angle((proton1+proton2).Vect());

  return cos(thetaCS);
}

KinFitTestFITHists::~KinFitTestFITHists(){}

