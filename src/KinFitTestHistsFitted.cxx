#include "UHH2/KinFitTest/include/KinFitTestHistsFitted.h"
#include "UHH2/core/include/Event.h"

#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>

#include "TH1F.h"
#include "TH2.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

KinFitTestHistsFitted::KinFitTestHistsFitted(Context & ctx, const string & dirname): HistsBASE(ctx, dirname){

  rec_fit_bhad_v4 = ctx.get_handle<TLorentzVector>("rec_fit_bhad_v4");
  rec_fit_blep_v4 = ctx.get_handle<TLorentzVector>("rec_fit_blep_v4");
  rec_fit_lepton_v4 = ctx.get_handle<TLorentzVector>("rec_fit_lepton_v4");
  rec_fit_neutrino_v4 = ctx.get_handle<TLorentzVector>("rec_fit_neutrino_v4");
  rec_fit_jetq1_v4 = ctx.get_handle<TLorentzVector>("rec_fit_jetq1_v4");
  rec_fit_jetq2_v4 = ctx.get_handle<TLorentzVector>("rec_fit_jetq2_v4");

  init();
}

void KinFitTestHistsFitted::init(){
  // primary vertices
  book_TH1F("wgt", 120, -6, 6);
  book_TH1F("pvN", 50, 0, 50);

  // muons
  book_TH1F("muoN", 10, 0, 10);
  //  book_TH1F("muo_charge", 5, -2, 3);
  book_TH1F("muo_pt", 40, 0, 200);
  book_TH1F("muo_eta", 40, -2.1, 2.1);
  // book_TH1F("muo_reliso", 40, 0, 0.5);

  // electrons
  book_TH1F("eleN", 10, 0, 10);
  //  book_TH1F("ele_charge"  , 5, -2, 3);
  book_TH1F("ele_pt", 40, 0, 200);
  book_TH1F("ele_eta",  40, -2.1, 2.1);
  //  book_TH1F("ele_reliso", 40, 0, 0.5);

  // jets
  book_TH1F("N_jets_Event", 20, 0, 20);  
  book_TH1F("jetN_CSVL_Event",10,0,10);
  book_TH1F("jetN_CSVM_Event",10,0,10);
  book_TH1F("jetN_CSVT_Event",10,0,10);

  book_TH1F("jetq1_pt", 90, 0, 900);
  book_TH1F("jetq1_eta", 40, -2.5, 2.5);
  book_TH1F("jetq2_pt",  90, 0, 900);
  book_TH1F("jetq2_eta", 40, -2.5, 2.5);
  book_TH1F("jet_bhad_pt",  90, 0, 900);
  book_TH1F("jet_bhad_eta",  40, -2.5, 2.5);
  book_TH1F("jet_blep_pt",  90, 0, 900);
  book_TH1F("jet_blep_eta",  40, -2.5, 2.5);
 
  // met
  book_TH1F("met_pt",90,0,900);
  book_TH1F("met_phi",60,-3.15,3.15);

  //Wlep
  // book_TH1F("wlep_ht", 90, 0, 900);
  // book_TH1F("wlep_pt", 90, 0, 900);

  //chi2 for fitted events
  // book_TH1F("chi2", 2000, 0, 10000);
  // book_TH1F("fitprob", 50, 0, 1);
}


void KinFitTestHistsFitted::fill(const Event & event){

  //weight
  double weight = event.weight;
  H1("wgt")->Fill(weight);

  int Npvs = event.pvs->size();
  H1("pvN")->Fill(Npvs, weight);

  TLorentzVector rec_fit_Bhad = event.get(rec_fit_bhad_v4);
  TLorentzVector rec_fit_Blep = event.get(rec_fit_blep_v4);
  TLorentzVector rec_fit_Lepton = event.get(rec_fit_lepton_v4);
  TLorentzVector rec_fit_Neu = event.get(rec_fit_neutrino_v4);
  TLorentzVector rec_fit_jetq1 = event.get(rec_fit_jetq1_v4);
  TLorentzVector rec_fit_jetq2 = event.get(rec_fit_jetq2_v4);

  //muons
  int Nmuons = event.muons->size();
  H1("muoN")->Fill(Nmuons, weight);
  if(Nmuons == 1){
    H1("muo_pt")->Fill(rec_fit_Lepton.Pt(), weight);
    //   H1("muo_charge")->Fill(rec_fit_Lepton.Charge(),weight);
    H1("muo_eta")->Fill(rec_fit_Lepton.Eta(), weight);
    //  H1("muo_reliso")->Fill(rec_fit_Lepton.relIso(), weight);
  }
 
  //electrons
  int Nele = event.electrons->size();
  H1("eleN")->Fill(Nele, weight);
  if(Nele == 1){
      H1("ele_pt")->Fill(rec_fit_Lepton.Pt(), weight);
      //      H1("ele_charge")->Fill(rec_fit_Lepton.charge(),weight);
      H1("ele_eta")->Fill(rec_fit_Lepton.Eta(), weight);
      //     H1("ele_reliso")->Fill(rec_fit_Lepton.relIso(), weight);
  } 


  //jets
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  H1("N_jets_Event")->Fill(Njets, weight);

  std::vector<float> jets_CSV;
  jets_CSV.reserve(event.jets->size());
  for(const auto& j : *event.jets) jets_CSV.push_back(j.btag_combinedSecondaryVertex());
  std::sort(jets_CSV.begin(), jets_CSV.end(), [](const float s1, const float s2){return s1 > s2;});

  int jetN__CSVL(0), jetN__CSVM(0), jetN__CSVT(0);
  for(unsigned int i=0; i<jets_CSV.size(); ++i){

    const float& csv = jets_CSV.at(i);

    if(csv > 0.605) ++jetN__CSVL;
    if(csv > 0.890) ++jetN__CSVM;
    if(csv > 0.970) ++jetN__CSVT;

    if(i > 2) continue;

  }

  H1("jetN_CSVL_Event")->Fill(jetN__CSVL, weight);
  H1("jetN_CSVM_Event")->Fill(jetN__CSVM, weight);
  H1("jetN_CSVT_Event")->Fill(jetN__CSVT, weight);
  

  H1("jetq1_pt")->Fill(rec_fit_jetq1.Pt(), weight);
  H1("jetq1_eta")->Fill(rec_fit_jetq1.Eta(), weight);
  H1("jetq2_pt")->Fill(rec_fit_jetq2.Pt(), weight);
  H1("jetq2_eta")->Fill(rec_fit_jetq2.Eta(), weight);
  H1("jet_bhad_pt")->Fill(rec_fit_Bhad.Pt(), weight);
  H1("jet_bhad_eta")->Fill(rec_fit_Bhad.Eta(), weight);
  H1("jet_blep_pt")->Fill(rec_fit_Blep.Pt(), weight);
  H1("jet_blep_eta")->Fill(rec_fit_Blep.Eta(), weight);
  
  //met
  H1("met_pt") ->Fill(rec_fit_Neu.Pt() , weight);
  H1("met_phi")->Fill(rec_fit_Neu.Phi(), weight);

  //Wlep
    // H1("wlep_ht")->Fill( rec_fit_Lepton.Pt()+rec_fit_Neu.Pt()      , weight);
    // H1("wlep_pt")->Fill((rec_fit_Lepton + rec_fit_Neu).Pt(), weight);



  // bool fit;
  // fit = event.get(fit_);
  // if(fit){
  //   float chi2(-1.);
  //   float fitprob(-1.);
  //   chi2 = event.get(chi2_allfitted_);
  //   fitprob = event.get(fitprob_allfitted_);
  //   H1("chi2")->Fill(chi2);
  //   H1("fitprob")->Fill(fitprob);
  // }

}



KinFitTestHistsFitted::~KinFitTestHistsFitted(){}
