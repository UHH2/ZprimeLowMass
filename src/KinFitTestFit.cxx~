#include "UHH2/KinFitTest/include/KinFitTestSetup.h"
#include "UHH2/core/include/Event.h"

using namespace uhh2examples;
using namespace uhh2;


  TMatrixD KinFitTestSetup::setupMatrix(const TLorentzVector& object, const ObjectType objType, const Param param)
{
  TMatrixD CovM3 (3,3); CovM3.Zero();
  TMatrixD CovM4 (4,4); CovM4.Zero();
  const double pt  = object.Pt();
  const double et  = object.Et();
  const double eta = object.Eta();

  switch(objType) {
  case kUdscJet:
    // if object is a non-b jet
    {
      JetUncertainties jetUnc;
      HelperJet jetRes;
      switch(param) {
      case kEMom :
	CovM4(0,0) = pow(jetRes.a (pt, eta, HelperJet::kUds), 2);
	CovM4(1,1) = pow(jetRes.b (pt, eta, HelperJet::kUds), 2);
	CovM4(2,2) = pow(jetRes.c (pt, eta, HelperJet::kUds), 2);
	CovM4(3,3) = pow(jetRes.d (pt, eta, HelperJet::kUds), 2);
	return CovM4;
      case kEtEtaPhi : 
	//modified uncertainties from fall 2011
	CovM3(0,0) = pow(jetUnc.et (et, eta, JetUncertainties::kUds), 2);
	CovM3(1,1) = pow(jetUnc.eta(et, eta, JetUncertainties::kUds), 2);
	CovM3(2,2) = pow(jetUnc.phi(et, eta, JetUncertainties::kUds), 2);
	//default uncertainties
	// CovM3(0,0) = pow(jetRes.et (pt, eta, HelperJet::kUds), 2);
	// //   CovM3(0,0)*= pow(getEtaDependentScaleFactor(object)       , 2);
	// CovM3(1,1) = pow(jetRes.eta(pt, eta, HelperJet::kUds), 2);
	// CovM3(2,2) = pow(jetRes.phi(pt, eta, HelperJet::kUds), 2);
	return CovM3;
      }
    }
    break;
  case kBJet:
    // if object is a b jet
    {
      JetUncertainties jetUnc;
      HelperJet jetRes;
      switch(param) {
      case kEMom :
	CovM4(0,0) = pow(jetRes.a (pt, eta, HelperJet::kB), 2);
	CovM4(1,1) = pow(jetRes.b (pt, eta, HelperJet::kB), 2);
	CovM4(2,2) = pow(jetRes.c (pt, eta, HelperJet::kB), 2);
	CovM4(3,3) = pow(jetRes.d (pt, eta, HelperJet::kB), 2);
	return CovM4;
      case kEtEtaPhi : 
	//modified uncertainties from fall 2011
	CovM3(0,0) = pow(jetUnc.et (et, eta, JetUncertainties::kB), 2);
	CovM3(1,1) = pow(jetUnc.eta(et, eta, JetUncertainties::kB), 2);
	CovM3(2,2) = pow(jetUnc.phi(et, eta, JetUncertainties::kB), 2);
	//default uncertainties
	// CovM3(0,0) = pow(jetRes.et (pt, eta, HelperJet::kB), 2);
	// //   CovM3(0,0)*= pow(getEtaDependentScaleFactor(object)     , 2);
	// CovM3(1,1) = pow(jetRes.eta(pt, eta, HelperJet::kB), 2);
	// CovM3(2,2) = pow(jetRes.phi(pt, eta, HelperJet::kB), 2);
	return CovM3;
      }
    }
    break;
  case kMuon:
    // if object is a muon
    {
      MuonUncertainties muonUnc;
      HelperMuon muonRes;
      switch(param){
      case kEMom :
	CovM3(0,0) = pow(muonRes.a (pt, eta), 2);
	CovM3(1,1) = pow(muonRes.b (pt, eta), 2); 
	CovM3(2,2) = pow(muonRes.c (pt, eta), 2);
	return CovM3;
      case kEtEtaPhi :
	// uncertainties from fall 2011
	CovM3(0,0) = pow(muonUnc.et (et, eta), 2);
	CovM3(1,1) = pow(muonUnc.eta(et, eta), 2); 
	CovM3(2,2) = pow(muonUnc.phi(et, eta), 2);
	//default uncertainties
	// CovM3(0,0) = pow(muonRes.et (pt, eta), 2);
	// CovM3(1,1) = pow(muonRes.eta(pt, eta), 2); 
	// CovM3(2,2) = pow(muonRes.phi(pt, eta), 2);
	return CovM3;
      }
    }
    break;
  case kElectron:
    {
      // if object is an electron
      ElectronUncertainties elecUnc;
      HelperElectron elecRes;
      switch(param){
      case kEMom :
	CovM3(0,0) = pow(elecRes.a (pt, eta), 2);
	CovM3(1,1) = pow(elecRes.b (pt, eta), 2); 
	CovM3(2,2) = pow(elecRes.c (pt, eta), 2);
	return CovM3;
      case kEtEtaPhi :
	// uncertainties from fall 2011
	  CovM3(0,0) = pow(elecUnc.et (et, eta), 2);
	  CovM3(1,1) = pow(elecUnc.eta(et, eta), 2); 
	  CovM3(2,2) = pow(elecUnc.phi(et, eta), 2);
	  //default uncertainties
	  // CovM3(0,0) = pow(elecRes.et (pt, eta), 2);
	  // CovM3(1,1) = pow(elecRes.eta(pt, eta), 2); 
	  // CovM3(2,2) = pow(elecRes.phi(pt, eta), 2);
	return CovM3;
      }
    }
    break;
  case kMet:
    // if object is met
    {
      METUncertainties metUnc;
      HelperMET metRes;
      switch(param){
      case kEMom :
	CovM3(0,0) = pow(metRes.a(pt), 2);
	CovM3(1,1) = pow(metRes.b(pt), 2);
	CovM3(2,2) = pow(metRes.c(pt), 2);
	return CovM3;
      case kEtEtaPhi :
	// uncertainties from fall 2011
	CovM3(0,0) = pow(metUnc.et(et) , 2);
	CovM3(1,1) = pow(9999. , 2);
	//CovM3(1,1) = pow(metUnc.eta(et) , 2);
	CovM3(2,2) = pow(metRes.phi(et), 2);
	//default uncertainties
	// CovM3(0,0) = pow(metRes.et(pt) , 2);
	// CovM3(1,1) = pow(        9999. , 2);
	// CovM3(2,2) = pow(metRes.phi(pt), 2);
	return CovM3;
      }
    }
    break;
  }
  std::cout<<"Something went wrong while trying to setup a covariance matrix!"<<std::endl;
  return CovM4; //should never get here
}


//   void KinFitTestSetup::setupJets()
// {
//   TMatrixD empty3x3(3,3); 
//   TMatrixD empty4x4(4,4);
//   switch(jetParam_){ // setup jets according to parameterization
//   case kEMom :
//     jethad1_= new TFitParticleEMomDev   ("Jet1_pt", "Jet1_pt", 0, &empty4x4);
//     jethad2_= new TFitParticleEMomDev   ("Jet2_pt", "Jet2_pt", 0, &empty4x4);
//     jethadb_= new TFitParticleEMomDev   ("Jet3_pt", "Jet3_pt", 0, &empty4x4);
//     jetlepb_= new TFitParticleEMomDev   ("Jet4_pt", "Jet4_pt", 0, &empty4x4);
//      break;
//   case kEtEtaPhi :
//     jethad1_= new TFitParticleEtEtaPhi  ("Jet1_eta", "Jet1_eta", 0, &empty3x3);
//     jethad2_= new TFitParticleEtEtaPhi  ("Jet2_eta", "Jet2_eta", 0, &empty3x3);
//     jethadb_= new TFitParticleEtEtaPhi  ("Jet3_eta", "Jet3_eta", 0, &empty3x3);
//     jetlepb_= new TFitParticleEtEtaPhi  ("Jet4_eta", "Jet4_eta", 0, &empty3x3);
//     break;
//   }
// }

// void KinFitTestSetup::setupLeptons()
// {
//   TMatrixD empty3x3(3,3); 
//   switch(lepParam_){ // setup lepton according to parameterization
//   case kEMom       : 
//     lepton_  = new TFitParticleEScaledMomDev("Lepton",   "Lepton",   0, &empty3x3);
//     break;
//   case kEtEtaPhi   : 
//     lepton_  = new TFitParticleEtEtaPhi     ("Lepton",   "Lepton",   0, &empty3x3); 
//     break;
//   }
//   switch(metParam_){ // setup MET according to parameterization
//   case kEMom       : 
//     met_= new TFitParticleEScaledMomDev("MET", "MET", 0, &empty3x3);
//     break;
//   case kEtEtaPhi   : 
//     met_= new TFitParticleEtEtaPhi     ("MET", "MET", 0, &empty3x3);
//     break;
//   }



// void KinFitTestSetup::setupConstraints() 
// {
//   //Decay width top, approximately
//   // float decw_top = 0.175*pow((mTop_/mW_),3); //Result with mTop=173. and mW_=80.4 is 1.74 GeV
//   // float decw_w = 2.085;
  

//   massConstr_[kWHadMass] = new TFitConstraintM("WMassHad",      "WMassHad",      0, 0, mW_);
//   massConstr_[kWLepMass      ] = new TFitConstraintM("WMassLep",      "WMassLep",      0, 0, mW_);
//   massConstr_[kTopHadMass    ] = new TFitConstraintM("TopMassHad",    "TopMassHad",    0, 0, mTop_);
//   massConstr_[kTopLepMass    ] = new TFitConstraintM("TopMassLep",    "TopMassLep",    0, 0, mTop_);
//   //  massConstr_[kMET  ] = new TFitConstraintM("MET",  "MET",  0, 0,    0.);
//   //massConstr_[kEqualTopMasses] = new TFitConstraintM("EqualTopMasses","EqualTopMasses",0, 0,    0.);
//   // //  massConstr_[kEqualWMasses] = new TFitConstraintM("EqualWMasses","EqualWMasses",0, 0,    0.);
//   sumPxConstr_                 = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
//   sumPyConstr_                 = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);


//   massConstr_[kWHadMass      ]->addParticles1(jethad1_,   jethad2_    );
//   massConstr_[kWLepMass      ]->addParticles1(lepton_, met_);
//   massConstr_[kTopHadMass    ]->addParticles1(jethad1_, jethad2_, jethadb_);
//   massConstr_[kTopLepMass    ]->addParticles1(lepton_, met_, jetlepb_);
//   //  massConstr_[kMET           ]->addParticle1 (met_);
//   //massConstr_[kEqualTopMasses]->addParticles1(jethad1_, jethad2_, jethadb_);
//   //massConstr_[kEqualTopMasses]->addParticles2(lepton_, met_, jetlepb_);
//   // // massConstr_[kEqualWMasses]->addParticles1(jethad1_, jethad2_);
//   // // massConstr_[kEqualWMasses]->addParticles2(lepton_, met_);
//   sumPxConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);
//   sumPyConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);


//   // float decw_top = 1.7;
//   // float decw_w = 2.;
//   // float decw_top = 60.;
//   // float decw_w = 30.;

//   // massConstr_[kWHadMass] = new TFitConstraintMGaus("WMassHad",      "WMassHad",      0, 0, mW_  , decw_w);
//   // massConstr_[kWLepMass      ] = new TFitConstraintMGaus("WMassLep",      "WMassLep",      0, 0, mW_ , decw_w);
//   // massConstr_[kTopHadMass    ] = new TFitConstraintMGaus("TopMassHad",    "TopMassHad",    0, 0, mTop_, decw_top);
//   // massConstr_[kTopLepMass    ] = new TFitConstraintMGaus("TopMassLep",    "TopMassLep",    0, 0, mTop_, decw_top );
//   // // massConstr_[kEqualTopMasses] = new TFitConstraintMGaus("EqualTopMasses","EqualTopMasses",0, 0,    0.1, 1.);
//   // //massConstr_[kEqualWMasses] = new TFitConstraintMGaus("EqualWMasses","EqualWMasses",0, 0,    .1, 1.);
//   // sumPxConstr_                 = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
//   // sumPyConstr_                 = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);

//   // massConstr_[kWHadMass      ]->addParticles1(jethad1_,   jethad2_    );
//   // massConstr_[kWLepMass      ]->addParticles1(lepton_, met_);
//   // massConstr_[kTopHadMass    ]->addParticles1(jethad1_, jethad2_, jethadb_);
//   // massConstr_[kTopLepMass    ]->addParticles1(lepton_, met_, jetlepb_);
//   // // massConstr_[kEqualTopMasses]->addParticles1(jethad1_, jethad2_, jethadb_);
//   // // massConstr_[kEqualTopMasses]->addParticles2(lepton_, met_, jetlepb_);
//   // // massConstr_[kEqualWMasses]->addParticles1(jethad1_, jethad2_);
//   // // massConstr_[kEqualWMasses]->addParticles2(lepton_, met_);
//   // sumPxConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);
//   // sumPyConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);



//   // if(std::find(constrList_.begin(), constrList_.end(), kSumPt)!=constrList_.end())    constrainSumPt_ = true;
//   //else  constrainSumPt_ = false;
//   constrainSumPt_ = false;
//   //constrainSumPt_ = true;

//   }
