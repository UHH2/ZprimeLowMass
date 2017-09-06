#include "UHH2/KinFitTest/include/KinFitTestFit.h"
#include "UHH2/core/include/Event.h"

using namespace uhh2examples;
using namespace uhh2;

void KinFitTestFit::fit(uhh2::Event & event, int bl, int bh, int h1, int h2, KinFitTestSetup::ObjectType leptonType, KinFitTestSetup::Param jetParam_, KinFitTestSetup::Param lepParam_, KinFitTestSetup::Param metParam_){

    // dm_w_lep_all = -100.;    
    // dm_w_had_all = -100.;
    // dm_top_lep_all = -100.;    
    // dm_top_had_all = -100.;
  
    // p4Had1 = TLorentzVector(event.jets->at(h1).v4().Px(),event.jets->at(h1).v4().Py(),event.jets->at(h1).v4().Pz(),event.jets->at(h1).v4().E());
    // p4Had2 = TLorentzVector(event.jets->at(h2).v4().Px(),event.jets->at(h2).v4().Py(),event.jets->at(h2).v4().Pz(),event.jets->at(h2).v4().E());
    // p4HadB = TLorentzVector(event.jets->at(bh).v4().Px(),event.jets->at(bh).v4().Py(),event.jets->at(bh).v4().Pz(),event.jets->at(bh).v4().E());
    // p4LepB = TLorentzVector(event.jets->at(bl).v4().Px(),event.jets->at(bl).v4().Py(),event.jets->at(bl).v4().Pz(),event.jets->at(bl).v4().E());
   
    // switch(leptonType){
    // case KinFitTestSetup::kElectron :       p4Lep  = TLorentzVector(event.electrons->at(0).v4().Px(),event.electrons->at(0).v4().Py(),event.electrons->at(0).v4().Pz(),event.electrons->at(0).v4().E()); break;
    // case KinFitTestSetup::kMuon :           p4Lep  = TLorentzVector(event.muons->at(0).v4().Px(),event.muons->at(0).v4().Py(),event.muons->at(0).v4().Pz(),event.muons->at(0).v4().E()); break;
    // case KinFitTestSetup::kMet : break;
    // case KinFitTestSetup::kBJet : break;
    // case KinFitTestSetup::kUdscJet : break;
    // }

    // p4Neutrino  = TLorentzVector(event.met->v4().Px(),event.met->v4().Py(),0,sqrt(pow(event.met->v4().M(),2) + pow(event.met->v4().Pt(), 2)));  
 
    // // set covariance matrices of the objects to be fitted
    // KinFitTestSetup setup;
    // TMatrixD covHad1 = setup.setupMatrix(p4Had1, KinFitTestSetup::kUdscJet, jetParam_);
    // TMatrixD covHad2 = setup.setupMatrix(p4Had2, KinFitTestSetup::kUdscJet, jetParam_);
    // TMatrixD covHadB = setup.setupMatrix(p4HadB, KinFitTestSetup::kBJet, jetParam_);
    // TMatrixD covLepB = setup.setupMatrix(p4LepB, KinFitTestSetup::kBJet, jetParam_);
    // TMatrixD covLep  = setup.setupMatrix(p4Lep,  leptonType, lepParam_);
    // TMatrixD covMET  = setup.setupMatrix(p4Neutrino, KinFitTestSetup::kMet, metParam_);

    // event.set(unc_had1, covHad1);
    // event.set(unc_had2, covHad2);
    // event.set(unc_hadb, covHadB);
    // event.set(unc_lepb, covLepB);
    // event.set(unc_lepton, covLep);
    // event.set(unc_neu, covMET);

    // // set the kinematics of the objects to be fitted
    // jethad1_->setIni4Vec( &p4Had1 );
    // jethad2_->setIni4Vec( &p4Had2 );
    // jethadb_->setIni4Vec( &p4HadB );
    // jetlepb_->setIni4Vec( &p4LepB );
    // lepton_ ->setIni4Vec( &p4Lep  );
    // met_    ->setIni4Vec( &p4Neutrino  );


    // jethad1_->setCovMatrix( &covHad1 );
    // jethad2_->setCovMatrix( &covHad2 );
    // jethadb_->setCovMatrix( &covHadB );
    // jetlepb_->setCovMatrix( &covLepB );
    // lepton_ ->setCovMatrix( &covLep  );
    // met_    ->setCovMatrix( &covMET  );
    
    // if(constrainSumPt_){
    //   // setup Px and Py constraint for curent event configuration so that sum Pt will be conserved
    //   sumPxConstr_->setConstraint( p4Had1.Px() + p4Had2.Px() + p4HadB.Px() + p4LepB.Px() + p4Lep.Px() + p4Neutrino.Px() );
    //   sumPyConstr_->setConstraint( p4Had1.Py() + p4Had2.Py() + p4HadB.Py() + p4LepB.Py() + p4Lep.Py() + p4Neutrino.Py() );
    // }

    // //do the fit
    // myKinFitter->fit();

    // //get chi^2 of fit
    // chi2_prev = chi2_curr;
    // F_prev = F_curr;
    // chi2_curr = 1000.;
    // F_curr = 0.;

    // chi2_curr = myKinFitter->getS();
    // F_curr = myKinFitter->getF();

    // fitprob_curr = 0.;
    // fitprob_curr =  TMath::Prob(myKinFitter->getS(), myKinFitter->getNDF());
    // fit_status_each = myKinFitter->getStatus();

    // jet_combi_all.push_back(h1);
    // jet_combi_all.push_back(h2);
    // jet_combi_all.push_back(bh);
    // jet_combi_all.push_back(bl);

    // //TEST: Set values for all combinations
    // solution2.clear();
    // input_vectors_all.clear();

    // TLorentzVector sol_had1_2 = TLorentzVector(jethad1_->getCurr4Vec()->Px(), jethad1_->getCurr4Vec()->Py(), jethad1_->getCurr4Vec()->Pz(),jethad1_->getCurr4Vec()->E());
    // TLorentzVector sol_had2_2 = TLorentzVector(jethad2_->getCurr4Vec()->Px(), jethad2_->getCurr4Vec()->Py(),jethad2_->getCurr4Vec()->Pz(),jethad2_->getCurr4Vec()->E());
    // TLorentzVector sol_hadb_2 = TLorentzVector(jethadb_->getCurr4Vec()->Px(), jethadb_->getCurr4Vec()->Py(), jethadb_->getCurr4Vec()->Pz(), jethadb_->getCurr4Vec()->E());
    // TLorentzVector sol_lepb_2 = TLorentzVector(jetlepb_->getCurr4Vec()->Px(), jetlepb_->getCurr4Vec()->Py(), jetlepb_->getCurr4Vec()->Pz(), jetlepb_->getCurr4Vec()->E());
    // TLorentzVector sol_lepton_2 = TLorentzVector(lepton_->getCurr4Vec()->Px(), lepton_->getCurr4Vec()->Py(), lepton_->getCurr4Vec()->Pz(), lepton_->getCurr4Vec()->E());
    // TLorentzVector sol_met_2 =  TLorentzVector(met_->getCurr4Vec()->Px(), met_->getCurr4Vec()->Py(), met_->getCurr4Vec()->Pz(), met_->getCurr4Vec()->E());
 

    // solution2.push_back(sol_had1_2);
    // solution2.push_back(sol_had2_2);
    // solution2.push_back(sol_hadb_2);
    // solution2.push_back(sol_lepb_2);
    // solution2.push_back(sol_lepton_2);
    // solution2.push_back(sol_met_2);

    // input_vectors_all.push_back(p4Had1);
    // input_vectors_all.push_back(p4Had2);
    // input_vectors_all.push_back(p4HadB);
    // input_vectors_all.push_back(p4LepB);
    // input_vectors_all.push_back(p4Lep);
    // input_vectors_all.push_back(p4Neutrino);

    // dm_w_lep_all = (p4Lep+p4Neutrino).M() - (sol_lepton_2+sol_met_2).M();
    // dm_w_had_all = (p4Had1+p4Had2).M() - (sol_had1_2+sol_had2_2).M();
    // dm_top_lep_all = (p4Lep+p4Neutrino+p4LepB).M() - (sol_lepton_2+sol_met_2+sol_lepb_2).M();
    // dm_top_had_all = (p4Had1+p4Had2+p4HadB).M() - (sol_had1_2+sol_had2_2+sol_hadb_2).M();
    // if(debug){
    //   std::cout<<"wlep all "<<dm_w_lep_all<<std::endl;
    //   std::cout<<"whad all "<<dm_w_had_all<<std::endl;
    //   std::cout<<"top lep all "<<dm_top_lep_all<<std::endl;
    //   std::cout<<"top had all "<<dm_top_had_all<<std::endl;
    // }
    // //Status of Fitter: 0 converged, 1 not converged (reaches e.g. max iteration), -1 no fit performed, -10 aborted, 10 running
    // if(fit_status_each==0 && fitprob_curr>fit_prob){
    // 	fit_prob = fitprob_curr;
    // 	chi2_lowest = chi2_curr;
    // 	fit_status = fit_status_each;

    // 	jet_combi_best.clear();
    // 	jet_combi_best.push_back(h1);
    // 	jet_combi_best.push_back(h2);
    // 	jet_combi_best.push_back(bh);
    // 	jet_combi_best.push_back(bl);

    // 	solution.clear();
    // 	input_vectors_fit_found.clear();

    //     TLorentzVector sol_had1 = TLorentzVector(jethad1_->getCurr4Vec()->Px(), jethad1_->getCurr4Vec()->Py(), jethad1_->getCurr4Vec()->Pz(),jethad1_->getCurr4Vec()->E());
    // 	TLorentzVector sol_had2 = TLorentzVector(jethad2_->getCurr4Vec()->Px(), jethad2_->getCurr4Vec()->Py(),jethad2_->getCurr4Vec()->Pz(),jethad2_->getCurr4Vec()->E());
    // 	TLorentzVector sol_hadb = TLorentzVector(jethadb_->getCurr4Vec()->Px(), jethadb_->getCurr4Vec()->Py(), jethadb_->getCurr4Vec()->Pz(), jethadb_->getCurr4Vec()->E());
    // 	TLorentzVector sol_lepb = TLorentzVector(jetlepb_->getCurr4Vec()->Px(), jetlepb_->getCurr4Vec()->Py(), jetlepb_->getCurr4Vec()->Pz(), jetlepb_->getCurr4Vec()->E());
    // 	TLorentzVector sol_lepton = TLorentzVector(lepton_->getCurr4Vec()->Px(), lepton_->getCurr4Vec()->Py(), lepton_->getCurr4Vec()->Pz(), lepton_->getCurr4Vec()->E());
    // 	TLorentzVector sol_met =  TLorentzVector(met_->getCurr4Vec()->Px(), met_->getCurr4Vec()->Py(), met_->getCurr4Vec()->Pz(), met_->getCurr4Vec()->E());

    // 	solution.push_back(sol_had1);
    // 	solution.push_back(sol_had2);
    // 	solution.push_back(sol_hadb);
    // 	solution.push_back(sol_lepb);
    // 	solution.push_back(sol_lepton);
    // 	solution.push_back(sol_met);

    // 	input_vectors_fit_found.push_back(p4Had1);
    // 	input_vectors_fit_found.push_back(p4Had2);
    // 	input_vectors_fit_found.push_back(p4HadB);
    // 	input_vectors_fit_found.push_back(p4LepB);
    // 	input_vectors_fit_found.push_back(p4Lep);
    // 	input_vectors_fit_found.push_back(p4Neutrino);

    // 	int jet_variables_size(0);
    // 	switch(jetParam_){
    // 	case KinFitTestSetup::kEMom :         
    // 	  jet_variables_size = 4;
    // 	  event.set(parametrisation, 0);
    // 	  break;
    // 	case KinFitTestSetup::kEtEtaPhi :     
    // 	  jet_variables_size = 3;
    // 	  event.set(parametrisation, 1);
    // 	  break;
    // 	}
    // 	const TTbarGen* ttgen_test(0);
    // 	if(isMC_){
    // 	  const auto& ttbargen_test = event.get(h_ttbar_gen);
    // 	  ttgen_test = &ttbargen_test;
    // 	}
    // 	setPull(jet_variables_size, ttgen_test, 0, isMC_);

    // 	if(debug){
    // 	  std::cout<<"fit found"<<std::endl;
    // 	  std::cout<<"SOLUTION:"<<std::endl;
    // 	  for(unsigned int j=0; j<solution.size();j++){
    // 	    solution[j].Print();
    // 	    std::cout<<"ET "<<solution[j].Et()<<std::endl;
    // 	    std::cout<<"pT "<<solution[j].Pt()<<std::endl;
    // 	    std::cout<<"Mass "<<solution[j].M()<<std::endl;
    // 	  }
    // 	  std::cout<<"chi2 "<<chi2_curr<<std::endl;
    // 	  std::cout<<"Fit probability "<<fitprob_curr<<std::endl;
    // 	  std::cout<<"INPUT:"<<std::endl;
    // 	  for(unsigned int j=0; j<input_vectors_fit_found.size();j++){
    // 	    input_vectors_fit_found[j].Print();
    // 	    std::cout<<"Et "<<input_vectors_fit_found[j].Et()<<std::endl;
    // 	    std::cout<<"Pt "<<input_vectors_fit_found[j].Pt()<<std::endl;
    // 	    std::cout<<"Mass "<<input_vectors_fit_found[j].M()<<std::endl;
    // 	  }
    // 	}

    // 	fit_found = true;

    // 	dm_w_lep = (p4Lep+p4Neutrino).M() - (sol_lepton+sol_met).M();
    // 	dm_w_had = (p4Had1+p4Had2).M() - (sol_had1+sol_had2).M();
    // 	dm_top_lep = (p4Lep+p4Neutrino+p4LepB).M() - (sol_lepton+sol_met+sol_lepb).M();
    // 	dm_top_had = (p4Had1+p4Had2+p4HadB).M() - (sol_had1+sol_had2+sol_hadb).M();
    // if(debug){
    //   std::cout<<"%%%%%%%%%%% fit found"<<std::endl;
    //   std::cout<<"wlep  "<<dm_w_lep<<std::endl;
    //   std::cout<<"whad  "<<dm_w_had<<std::endl;
    //   std::cout<<"top lep  "<<dm_top_lep<<std::endl;
    //   std::cout<<"top had  "<<dm_top_had<<std::endl;
    // }
    // }
    

    // int jet_variables_size_all(0);
    // switch(jetParam_){
    // case KinFitTestSetup::kEMom :         
    //   jet_variables_size_all = 4; 	  
    //   event.set(parametrisation, 0);
    //  break;
    // case KinFitTestSetup::kEtEtaPhi :     
    //   jet_variables_size_all = 3; 
    //   event.set(parametrisation, 1);
    //   break;
    // }

    // const TTbarGen* ttgen_test1(0);
    // if(isMC_){
    //   const auto& ttbargen_test1 = event.get(h_ttbar_gen);
    //   ttgen_test1 = &ttbargen_test1;
    // }

    // setPull(jet_variables_size_all, ttgen_test1, 1, isMC_);
 
    // event.set(pull_had1, pull_had1_all);
    // event.set(pull_had2, pull_had2_all);
    // event.set(pull_hadb, pull_hadb_all);
    // event.set(pull_lepb, pull_lepb_all);
    // event.set(pull_lepton, pull_lepton_all);
    // event.set(pull_neu, pull_neu_all); 
    // event.set(pull_gtr_had1_et, -100.);
    // event.set(pull_gtr_had1_eta, -100.);
    // event.set(pull_gtr_had1_phi, -100.);
    // event.set(pull_gtr_had2_et, -100.);
    // event.set(pull_gtr_had2_eta, -100.);
    // event.set(pull_gtr_had2_phi, -100.);
    // event.set(pull_gtr_hadb_et, -100.);
    // event.set(pull_gtr_hadb_eta, -100.);
    // event.set(pull_gtr_hadb_phi, -100.);
    // event.set(pull_gtr_lepb_et, -100.);
    // event.set(pull_gtr_lepb_eta, -100.);
    // event.set(pull_gtr_lepb_phi, -100.);
    // event.set(pull_gtr_lepton_et, -100.);
    // event.set(pull_gtr_lepton_eta, -100.);
    // event.set(pull_gtr_lepton_phi, -100.);
    // event.set(pull_gtr_neu_et, -100.);
    // event.set(pull_gtr_neu_eta, -100.);
    // event.set(pull_gtr_neu_phi, -100.);

    // if((verbosity==2 || verbosity==3 ) ){
    //   std::cout<<"##############################"<<std::endl;
    //   std::cout<<"input: "<<std::endl;
    //   std::cout<<"lepton type "<<leptonType<<std::endl;
    //   std::cout<<"p4LepB = "<<std::endl;
    //   p4LepB.Print();
    //   std::cout<<"p4LepB pt = "<<p4LepB.Pt()<<std::endl;
    //   std::cout<<"p4LepB mass = "<<p4LepB.M()<<std::endl;
    //   std::cout<<"p4HadB = "<<std::endl;
    //   p4HadB.Print();
    //   std::cout<<"p4HadB pt = "<<p4HadB.Pt()<<std::endl;
    //   std::cout<<"p4HadB mass = "<<p4HadB.M()<<std::endl;
    //   std::cout<<"p4Had1 = "<<std::endl;
    //   p4Had1.Print();
    //   std::cout<<"p4Had1 pt = "<<p4Had1.Pt()<<std::endl;
    //   std::cout<<"p4Had1 mass = "<<p4Had1.M()<<std::endl;
    //   std::cout<<"p4Had2 = "<<std::endl;
    //   p4Had2.Print();
    //   std::cout<<"p4Had2 pt = "<<p4Had2.Pt()<<std::endl;
    //   std::cout<<"p4Had2 mass = "<<p4Had2.M()<<std::endl;
    //   std::cout<<"p4Lep  = "<<std::endl;
    //   p4Lep.Print();
    //   std::cout<<"p4Lep pt = "<<p4Lep.Pt()<<std::endl;
    //   std::cout<<"p4Lep mass = "<<p4Lep.M()<<std::endl;
    //   std::cout<<"p4Neutrino  = "<<std::endl;
    //   p4Neutrino.Print();
    //   std::cout<<"p4Neutrino pt = "<<p4Neutrino.Pt()<<std::endl;
    //   std::cout<<"p4Neutrino mass = "<<p4Neutrino.M()<<std::endl;
    //   std::cout<<"WHad mass = "<<(p4Had1+p4Had2).M()<<std::endl;
    //   std::cout<<"WLep mass = "<<(p4Lep+p4Neutrino).M()<<std::endl;
    //   std::cout<<"TLep mass = "<<(p4Lep+p4Neutrino+p4LepB).M()<<std::endl;
    //   std::cout<<"THad mass = "<<(p4Had1+p4Had2+p4HadB).M()<<std::endl;
    //   std::cout<<"##############################"<<std::endl;
    // }    

    // event.set(jet_permutation_before_dm_cut, 0);
    // event.set(jet_permutation_after_dm_cut, 0);


    // setBestFit(event, solution2, input_vectors_all, leptonType, chi2_curr, fitprob_curr, fit_status_each, fit_found, false, jet_combi_all);
    // h_fitted_allcomb->fill(event);
    // if(fit_status_each==0) h_fitted_allcomb_fit_converged->fill(event);

  }
