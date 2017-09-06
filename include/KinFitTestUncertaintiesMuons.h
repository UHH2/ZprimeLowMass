/* Here the uncertainties for the covariance matrix are set.
   The uncertainties come in principle from:  cmssw/TopQuarkAnalysis/TopObjectResolutions/python/stringResolutions_etEtaPhi_Fall11_cff.py */



//Resolution of the Muons to set up covariance matrix
class MuonUncertainties {

 public:
  enum Flavor {kUds, kB};
    
  MuonUncertainties(){};
  ~MuonUncertainties(){};
  
  inline double et(double et, double eta);
  inline double eta(double et, double eta);
  inline double phi(double et, double eta);
};

inline double MuonUncertainties::et(double et, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.100)      res = et * (0.00517 + 0.000143 * et);
  else if(fabs(eta)<0.200) res = et * (0.00524 + 0.000143 * et);
  else if(fabs(eta)<0.300) res = et * (0.00585 + 0.000138 * et);
  else if(fabs(eta)<0.400) res = et * (0.0065 + 0.000133 * et);
  else if(fabs(eta)<0.500) res = et * (0.0071 + 0.000129 * et);
  else if(fabs(eta)<0.600) res = et * (0.00721 + 0.00013 * et);
  else if(fabs(eta)<0.700) res = et * (0.00757 + 0.000129 * et);
  else if(fabs(eta)<0.800) res = et * (0.0081 + 0.000127 * et);
  else if(fabs(eta)<0.900) res = et * (0.00916 + 0.000131 * et);
  else if(fabs(eta)<1.000) res = et * (0.0108 + 0.000151 * et);
  else if(fabs(eta)<1.100) res = et * (0.0115 + 0.000153 * et);
  else if(fabs(eta)<1.200) res = et * (0.013 + 0.000136 * et);
  else if(fabs(eta)<1.300) res = et * (0.0144 + 0.000131 * et);
  else if(fabs(eta)<1.400) res = et * (0.0149 + 0.000141 * et);
  else if(fabs(eta)<1.500) res = et * (0.014 + 0.000155 * et);
  else if(fabs(eta)<1.600) res = et * (0.0132 + 0.000169 * et);
  else if(fabs(eta)<1.700) res = et * (0.0129 + 0.0002 * et);
  else if(fabs(eta)<1.800) res = et * (0.0135 + 0.000264 * et);
  else if(fabs(eta)<1.900) res = et * (0.0144 + 0.00034 * et);
  else if(fabs(eta)<2.000) res = et * (0.0147 + 0.000441 * et);
  else if(fabs(eta)<2.100) res = et * (0.0154 + 0.000604 * et);
  else if(fabs(eta)<2.200) res = et * (0.0163 + 0.000764 * et);
  else if(fabs(eta)<2.300) res = et * (0.0173 + 0.000951 * et);
  else if(fabs(eta)<2.400) res = et * (0.0175 + 0.00126 * et);
  
  
  return res;
}
				       

inline double MuonUncertainties::eta(double et, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.100)      res = sqrt(pow(0.000433,2) + pow((0.000161/sqrt(et)),2) + pow((0.00334/et),2));
  else if(fabs(eta)<0.200) res = sqrt(pow(0.000381,2) + pow((0.000473/sqrt(et)),2) + pow((0.00259/et),2));
  else if(fabs(eta)<0.300) res = sqrt(pow(0.000337,2) + pow((0.000381/sqrt(et)),2) + pow((0.0023/et),2));
  else if(fabs(eta)<0.400) res = sqrt(pow(0.000308,2) + pow((0.000166/sqrt(et)),2) + pow((0.00249/et),2));
  else if(fabs(eta)<0.500) res = sqrt(pow(0.000289,2) + pow((5.37e-09/sqrt(et)),2) + pow((0.00243/et),2));
  else if(fabs(eta)<0.600) res = sqrt(pow(0.000279,2) + pow((0.000272/sqrt(et)),2) + pow((0.0026/et),2));
  else if(fabs(eta)<0.700) res = sqrt(pow(0.000282,2) + pow((3.63e-10/sqrt(et)),2) + pow((0.00288/et),2));
  else if(fabs(eta)<0.800) res = sqrt(pow(0.000265,2) + pow((0.000609/sqrt(et)),2) + pow((0.00212/et),2));
  else if(fabs(eta)<0.900) res = sqrt(pow(0.000241,2) + pow((0.000678/sqrt(et)),2) + pow((0.00221/et),2));
  else if(fabs(eta)<1.000) res = sqrt(pow(0.000228,2) + pow((0.000612/sqrt(et)),2) + pow((0.00245/et),2));
  else if(fabs(eta)<1.100) res = sqrt(pow(0.000217,2) + pow((0.000583/sqrt(et)),2) + pow((0.00307/et),2));
  else if(fabs(eta)<1.200) res = sqrt(pow(0.000195,2) + pow((0.000751/sqrt(et)),2) + pow((0.00282/et),2));
  else if(fabs(eta)<1.300) res = sqrt(pow(0.000183,2) + pow((0.000838/sqrt(et)),2) + pow((0.00227/et),2));
  else if(fabs(eta)<1.400) res = sqrt(pow(0.000196,2) + pow((0.000783/sqrt(et)),2) + pow((0.00274/et),2));
  else if(fabs(eta)<1.500) res = sqrt(pow(0.0002,2) + pow((0.000832/sqrt(et)),2) + pow((0.00254/et),2));
  else if(fabs(eta)<1.600) res = sqrt(pow(0.000205,2) + pow((0.0007/sqrt(et)),2) + pow((0.00304/et),2));
  else if(fabs(eta)<1.700) res = sqrt(pow(0.000214,2) + pow((0.000747/sqrt(et)),2) + pow((0.00319/et),2));
  else if(fabs(eta)<1.800) res = sqrt(pow(0.000238,2) + pow((0.000582/sqrt(et)),2) + pow((0.00343/et),2));
  else if(fabs(eta)<1.900) res = sqrt(pow(0.000263,2) + pow((0.000721/sqrt(et)),2) + pow((0.00322/et),2));
  else if(fabs(eta)<2.000) res = sqrt(pow(0.000284,2) + pow((0.000779/sqrt(et)),2) + pow((0.0031/et),2));
  else if(fabs(eta)<2.100) res = sqrt(pow(0.000316,2) + pow((0.000566/sqrt(et)),2) + pow((0.00384/et),2));
  else if(fabs(eta)<2.200) res = sqrt(pow(0.000353,2) + pow((0.000749/sqrt(et)),2) + pow((0.0038/et),2));
  else if(fabs(eta)<2.300) res = sqrt(pow(0.000412,2) + pow((0.00102/sqrt(et)),2) + pow((0.00351/et),2));
  else if(fabs(eta)<2.400) res = sqrt(pow(0.000506,2) + pow((0.000791/sqrt(et)),2) + pow((0.0045/et),2));

  return res;
}


inline double MuonUncertainties::phi(double et, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.100)      res = sqrt(pow(7.21e-05,2) + pow((7e-05/sqrt(et)),2) + pow((0.00296/et),2));
  else if(fabs(eta)<0.200) res = sqrt(pow(6.79e-05,2) + pow((0.000245/sqrt(et)),2) + pow((0.00274/et),2));
  else if(fabs(eta)<0.300) res = sqrt(pow(7.08e-05,2) + pow((6.75e-05/sqrt(et)),2) + pow((0.00307/et),2));
  else if(fabs(eta)<0.400) res = sqrt(pow(6.59e-05,2) + pow((0.000301/sqrt(et)),2) + pow((0.00281/et),2));
  else if(fabs(eta)<0.500) res = sqrt(pow(6.27e-05,2) + pow((0.000359/sqrt(et)),2) + pow((0.00278/et),2));
  else if(fabs(eta)<0.600) res = sqrt(pow(6.46e-05,2) + pow((0.00036/sqrt(et)),2) + pow((0.00285/et),2));
  else if(fabs(eta)<0.700) res = sqrt(pow(6.54e-05,2) + pow((0.000348/sqrt(et)),2) + pow((0.00301/et),2));
  else if(fabs(eta)<0.800) res = sqrt(pow(6.2e-05,2) + pow((0.000402/sqrt(et)),2) + pow((0.00304/et),2));
  else if(fabs(eta)<0.900) res = sqrt(pow(6.26e-05,2) + pow((0.000458/sqrt(et)),2) + pow((0.0031/et),2));
  else if(fabs(eta)<1.000) res = sqrt(pow(7.18e-05,2) + pow((0.000469/sqrt(et)),2) + pow((0.00331/et),2));
  else if(fabs(eta)<1.100) res = sqrt(pow(6.98e-05,2) + pow((0.000507/sqrt(et)),2) + pow((0.00338/et),2));
  else if(fabs(eta)<1.200) res = sqrt(pow(6.21e-05,2) + pow((0.000584/sqrt(et)),2) + pow((0.00345/et),2));
  else if(fabs(eta)<1.300) res = sqrt(pow(5.37e-05,2) + pow((0.000667/sqrt(et)),2) + pow((0.00352/et),2));
  else if(fabs(eta)<1.400) res = sqrt(pow(5.37e-05,2) + pow((0.000711/sqrt(et)),2) + pow((0.00358/et),2));
  else if(fabs(eta)<1.500) res = sqrt(pow(5.98e-05,2) + pow((0.000713/sqrt(et)),2) + pow((0.00362/et),2));
  else if(fabs(eta)<1.600) res = sqrt(pow(6.21e-05,2) + pow((0.000781/sqrt(et)),2) + pow((0.00348/et),2));
  else if(fabs(eta)<1.700) res = sqrt(pow(6.92e-05,2) + pow((0.000865/sqrt(et)),2) + pow((0.00337/et),2));
  else if(fabs(eta)<1.800) res = sqrt(pow(9.13e-05,2) + pow((0.000896/sqrt(et)),2) + pow((0.00348/et),2));
  else if(fabs(eta)<1.900) res = sqrt(pow(0.000102,2) + pow((0.000994/sqrt(et)),2) + pow((0.00337/et),2));
  else if(fabs(eta)<2.000) res = sqrt(pow(0.000123,2) + pow((0.00108/sqrt(et)),2) + pow((0.00315/et),2));
  else if(fabs(eta)<2.100) res = sqrt(pow(0.000169,2) + pow((0.000947/sqrt(et)),2) + pow((0.00422/et),2));
  else if(fabs(eta)<2.200) res = sqrt(pow(0.000176,2) + pow((0.00116/sqrt(et)),2) + pow((0.00423/et),2));
  else if(fabs(eta)<2.300) res = sqrt(pow(0.000207,2) + pow((0.00115/sqrt(et)),2) + pow((0.00469/et),2));
  else if(fabs(eta)<2.400) res = sqrt(pow(0.00027,2) + pow((0.00113/sqrt(et)),2) + pow((0.00528/et),2));

  return res;
}
