#include "UHH2/KinFitTest/include/KinFitTestHists.h"
#include "UHH2/core/include/Event.h"

#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>

#include "TH1F.h"
#include "TH2.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

KinFitTestHists::KinFitTestHists(Context & ctx, const string & dirname): HistsBASE(ctx, dirname){
  init();
}

void KinFitTestHists::init(){
  // primary vertices
  book_TH1F("pvN", 50, 0, 50);
  book_TH1F("wgt", 120, -6, 6);

  // muons
  book_TH1F("muoN", 10, 0, 10);
  book_TH1F("muo_charge", 5, -2, 3);
  book_TH1F("muo_pt", 40, 0, 200);
  book_TH2F("muo_pt_vs_pt_error", 100, 0, 200, 20, 0, 10);
  book_TH1F("muo_eta", 40, -2.1, 2.1);
  book_TH1F("muo_reliso", 40, 0, 0.5);

  // electrons
  book_TH1F("eleN", 10, 0, 10);
  book_TH1F("ele_charge"  , 5, -2, 3);
  book_TH1F("ele_pt", 40, 0, 200);
  book_TH2F("ele_pt_vs_pt_error", 200, 0, 400, 200, 0, 400);
  book_TH1F("ele_eta",  40, -2.1, 2.1);
  book_TH2F("ele_eta_vs_eta_error", 20, 0, 3, 20, 0, 0.3);
  book_TH1F("ele_reliso", 40, 0, 0.5);

  // jets
  book_TH1F("N_jets", 20, 0, 20);  
  book_TH1F("jetN_CSVL",10,0,10);
  book_TH1F("jetN_CSVM",10,0,10);
  book_TH1F("jetN_CSVT",10,0,10);

  book_TH1F("jet1_pt", 90, 0, 900);
  book_TH1F("jet1_eta", 40, -2.5, 2.5);
  book_TH1F("jet2_pt",  90, 0, 900);
  book_TH1F("jet2_eta", 40, -2.5, 2.5);
  book_TH1F("jet3_pt",  90, 0, 900);
  book_TH1F("jet3_eta",  40, -2.5, 2.5);
  book_TH1F("jet4_pt",  90, 0, 900);
  book_TH1F("jet4_eta",  40, -2.5, 2.5);
 
  // met
  book_TH1F("met_pt",90,0,900);
  book_TH1F("met_phi",60,-3.15,3.15);

  //Wlep
  book_TH1F("wlep_ht", 90, 0, 900);
  book_TH1F("wlep_pt", 90, 0, 900);
  book_TH1F("wlep_Mt", 360, 0,  360);

}


void KinFitTestHists::fill(const Event & event){

  //weight
  double weight = event.weight;
  H1("wgt")->Fill(weight);

  //jets
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  H1("N_jets")->Fill(Njets, weight);

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

  H1("jetN_CSVL")->Fill(jetN__CSVL, weight);
  H1("jetN_CSVM")->Fill(jetN__CSVM, weight);
  H1("jetN_CSVT")->Fill(jetN__CSVT, weight);
  
  if(Njets>=1){
    H1("jet1_pt")->Fill(jets->at(0).pt(), weight);
    H1("jet1_eta")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    H1("jet2_pt")->Fill(jets->at(1).pt(), weight);
    H1("jet2_eta")->Fill(jets->at(1).eta(), weight);
  }
  if(Njets>=3){
    H1("jet3_pt")->Fill(jets->at(2).pt(), weight);
    H1("jet3_eta")->Fill(jets->at(2).eta(), weight);
  }
  if(Njets>=4){
    H1("jet4_pt")->Fill(jets->at(3).pt(), weight);
    H1("jet4_eta")->Fill(jets->at(3).eta(), weight);
  }

  //muons
  int Nmuons = event.muons->size();
  H1("muoN")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      H1("muo_pt")->Fill(thismu.pt(), weight);
      H1("muo_charge")->Fill(thismu.charge(),weight);
      H1("muo_eta")->Fill(thismu.eta(), weight);
      H1("muo_reliso")->Fill(thismu.relIso(), weight);
      double relres(-1.);
      if (thismu.pt()<100) relres = 0.01;
      else if (thismu.pt()<200) relres = 0.02;
      else if (thismu.pt()<400) relres = 0.035;
      else if (thismu.pt()<1000) relres = 0.07;
      double res = thismu.pt()*relres;
      H2("muo_pt_vs_pt_error")->Fill(double(thismu.pt()),double(res),weight);
  }

  //electrons
   int Nele = event.electrons->size();
  H1("eleN")->Fill(Nele, weight);
  for (const Electron & thisele : *event.electrons){
      H1("ele_pt")->Fill(thisele.pt(), weight);
      H1("ele_charge")->Fill(thisele.charge(),weight);
      H1("ele_eta")->Fill(thisele.eta(), weight);
      H1("ele_reliso")->Fill(thisele.relIso(), weight);
      H2("ele_pt_vs_pt_error")->Fill(double(thisele.pt()),double(thisele.ptError()),weight);
      H2("ele_eta_vs_eta_error")->Fill(double(thisele.eta()),double(thisele.etaError()),weight);
  } 
  int Npvs = event.pvs->size();
  H1("pvN")->Fill(Npvs, weight);

  //met
  H1("met_pt") ->Fill(event.met->pt() , weight);
  H1("met_phi")->Fill(event.met->phi(), weight);

  //Wlep
  const Particle* lep1(0);
  float max_lep_pt(0.);
  for(const auto& l : *event.muons)     if(l.pt() > max_lep_pt){ lep1 = &l; max_lep_pt = l.pt(); }
  for(const auto& l : *event.electrons) if(l.pt() > max_lep_pt){ lep1 = &l; max_lep_pt = l.pt(); }

  if(lep1){

    H1("wlep_ht")->Fill( event.met->pt()+lep1->pt()      , weight);
    H1("wlep_pt")->Fill((event.met->v4()+lep1->v4()).Pt(), weight);
    H1("wlep_Mt")->Fill(sqrt(2*event.met->pt()*lep1->pt()*(1.-cos(uhh2::deltaPhi(*event.met, *lep1)))), weight);
  }

}



KinFitTestHists::~KinFitTestHists(){}
