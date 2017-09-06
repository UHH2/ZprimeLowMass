#pragma once

#include <TLorentzVector.h>
#include "TMatrixD.h"

#include <UHH2/KinFitTest/include/KinFitTestUncertaintiesJets.h>
#include <UHH2/KinFitTest/include/KinFitTestUncertaintiesMuons.h>
#include <UHH2/KinFitTest/include/KinFitTestUncertaintiesMET.h>
#include <UHH2/KinFitTest/include/KinFitTestUncertaintiesElectrons.h>
#include <UHH2/KinFitTest/include/KinFitTestUncertaintiesDefault.h>


namespace uhh2examples {
    
class KinFitTestSetup{
 public:
  enum ObjectType{kUdscJet, kBJet, kMuon, kElectron, kMet};
  enum Param{ kEMom, kEtEtaPhi};

  TMatrixD setupMatrix(const TLorentzVector& object, const ObjectType objType, const Param param);
};
}
