/* Here the uncertainties for the covariance matrix are set.
   The uncertainties come in principle from:  cmssw/TopQuarkAnalysis/TopObjectResolutions/python/stringResolutions_etEtaPhi_Fall11_cff.py */



//Resolution of the Electrons to set up covariance matrix
class ElectronUncertainties {

 public:
  enum Flavor {kUds, kB};
    
  ElectronUncertainties(){};
  ~ElectronUncertainties(){};
  
  inline double et(double et, double eta);
  inline double eta(double et, double eta);
  inline double phi(double et, double eta);
};

inline double ElectronUncertainties::et(double et, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.174)      res = et * (sqrt(pow(0.00534,2) + pow((0.079/sqrt(et)),2) + pow((0.163/et),2)));
  else if(fabs(eta)<0.261) res = et * (sqrt(pow(0.00518,2) + pow((0.0749/sqrt(et)),2) + pow((0.227/et),2)));
  else if(fabs(eta)<0.348) res = et * (sqrt(pow(0.00332,2) + pow((0.0879/sqrt(et)),2) + pow((0.12/et),2)));
  else if(fabs(eta)<0.435) res = et * (sqrt(pow(0.00445,2) + pow((0.0895/sqrt(et)),2) + pow((0.186/et),2)));
  else if(fabs(eta)<0.522) res = et * (sqrt(pow(0.00453,2) + pow((0.0893/sqrt(et)),2) + pow((0.21/et),2)));
  else if(fabs(eta)<0.609) res = et * (sqrt(pow(0.00308,2) + pow((0.0886/sqrt(et)),2) + pow((0.188/et),2)));
  else if(fabs(eta)<0.696) res = et * (sqrt(pow(0.00308,2) + pow((0.0914/sqrt(et)),2) + pow((0.182/et),2)));
  else if(fabs(eta)<0.783) res = et * (sqrt(pow(0.00442,2) + pow((0.0914/sqrt(et)),2) + pow((0.231/et),2)));
  else if(fabs(eta)<0.870) res = et * (sqrt(pow(0.00455,2) + pow((0.0949/sqrt(et)),2) + pow((0.335/et),2)));
  else if(fabs(eta)<0.957) res = et * (sqrt(pow(0.00181,2) + pow((0.102/sqrt(et)),2) + pow((0.333/et),2)));
  else if(fabs(eta)<1.044) res = et * (sqrt(pow(0.000764,2) + pow((0.108/sqrt(et)),2) + pow((0.42/et),2)));
  else if(fabs(eta)<1.131) res = et * (sqrt(pow(0.00114,2) + pow((0.128/sqrt(et)),2) + pow((0.55/et),2)));
  else if(fabs(eta)<1.218) res = et * (sqrt(pow(4.14e-09,2) + pow((0.155/sqrt(et)),2) + pow((0.674/et),2)));
  else if(fabs(eta)<1.305) res = et * (sqrt(pow(8.03e-09,2) + pow((0.144/sqrt(et)),2) + pow((0.8/et),2)));
  else if(fabs(eta)<1.392) res = et * (sqrt(pow(0.00842,2) + pow((0.118/sqrt(et)),2) + pow((0.951/et),2)));
  else if(fabs(eta)<1.479) res = et * (sqrt(pow(0.00684,2) + pow((0.144/sqrt(et)),2) + pow((0.892/et),2)));
  else if(fabs(eta)<1.653) res = et * (sqrt(pow(0.0245,2) + pow((0.196/sqrt(et)),2) + pow((0.555/et),2)));
  else if(fabs(eta)<1.740) res = et * (sqrt(pow(0.0174,2) + pow((0.127/sqrt(et)),2) + pow((0.894/et),2)));
  else if(fabs(eta)<1.830) res = et * (sqrt(pow(0.0144,2) + pow((0.133/sqrt(et)),2) + pow((0.708/et),2)));
  else if(fabs(eta)<1.930) res = et * (sqrt(pow(0.0149,2) + pow((0.126/sqrt(et)),2) + pow((0.596/et),2)));
  else if(fabs(eta)<2.043) res = et * (sqrt(pow(0.0143,2) + pow((0.12/sqrt(et)),2) + pow((0.504/et),2)));
  else if(fabs(eta)<2.172) res = et * (sqrt(pow(0.0162,2) + pow((0.0965/sqrt(et)),2) + pow((0.483/et),2)));
  else if(fabs(eta)<2.322) res = et * (sqrt(pow(0.0122,2) + pow((0.13/sqrt(et)),2) + pow((0.207/et),2)));
  else if(fabs(eta)<2.500) res = et * (sqrt(pow(0.0145,2) + pow((0.127/sqrt(et)),2) + pow((0.0782/et),2)));
  return res;
}


inline double ElectronUncertainties::eta(double et, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.174)      res = sqrt(pow(0.000452,2) + pow((0.000285/sqrt(et)),2) + pow((0.00376/et),2));
  else if(fabs(eta)<0.261) res = sqrt(pow(0.00038,2) + pow((0.000571/sqrt(et)),2) + pow((0.00276/et),2));
  else if(fabs(eta)<0.348) res = sqrt(pow(0.000351,2) + pow((1.36e-09/sqrt(et)),2) + pow((0.00324/et),2));
  else if(fabs(eta)<0.435) res = sqrt(pow(0.000319,2) + pow((0.00061/sqrt(et)),2) + pow((0.00182/et),2));
  else if(fabs(eta)<0.522) res = sqrt(pow(0.000301,2) + pow((0.000612/sqrt(et)),2) + pow((0.00146/et),2));
  else if(fabs(eta)<0.609) res = sqrt(pow(0.000297,2) + pow((0.000791/sqrt(et)),2) + pow((2.09e-08/et),2));
  else if(fabs(eta)<0.696) res = sqrt(pow(0.00032,2) + pow((0.000329/sqrt(et)),2) + pow((0.00325/et),2));
  else if(fabs(eta)<0.783) res = sqrt(pow(0.000309,2) + pow((0.000821/sqrt(et)),2) + pow((0.00119/et),2));
  else if(fabs(eta)<0.870) res = sqrt(pow(0.000293,2) + pow((0.000767/sqrt(et)),2) + pow((0.00211/et),2));
  else if(fabs(eta)<0.957) res = sqrt(pow(0.000275,2) + pow((0.000765/sqrt(et)),2) + pow((0.00227/et),2));
  else if(fabs(eta)<1.044) res = sqrt(pow(0.000274,2) + pow((0.000622/sqrt(et)),2) + pow((0.00299/et),2));
  else if(fabs(eta)<1.131) res = sqrt(pow(0.000269,2) + pow((0.000929/sqrt(et)),2) + pow((0.00183/et),2));
  else if(fabs(eta)<1.218) res = sqrt(pow(0.000268,2) + pow((0.000876/sqrt(et)),2) + pow((0.00234/et),2));
  else if(fabs(eta)<1.305) res = sqrt(pow(0.000258,2) + pow((0.000782/sqrt(et)),2) + pow((0.00246/et),2));
  else if(fabs(eta)<1.392) res = sqrt(pow(0.000269,2) + pow((0.000817/sqrt(et)),2) + pow((0.00278/et),2));
  else if(fabs(eta)<1.479) res = sqrt(pow(0.000267,2) + pow((0.000734/sqrt(et)),2) + pow((0.00327/et),2));
  else if(fabs(eta)<1.653) res = sqrt(pow(0.000268,2) + pow((0.000757/sqrt(et)),2) + pow((0.00295/et),2));
  else if(fabs(eta)<1.740) res = sqrt(pow(0.000274,2) + pow((1.77e-09/sqrt(et)),2) + pow((0.00435/et),2));
  else if(fabs(eta)<1.830) res = sqrt(pow(0.000274,2) + pow((0.00101/sqrt(et)),2) + pow((0.000982/et),2));
  else if(fabs(eta)<1.930) res = sqrt(pow(0.000299,2) + pow((0.000686/sqrt(et)),2) + pow((0.00341/et),2));
  else if(fabs(eta)<2.043) res = sqrt(pow(0.000329,2) + pow((3.05e-10/sqrt(et)),2) + pow((0.00439/et),2));
  else if(fabs(eta)<2.172) res = sqrt(pow(0.00037,2) + pow((1.32e-08/sqrt(et)),2) + pow((0.00447/et),2));
  else if(fabs(eta)<2.322) res = sqrt(pow(0.000442,2) + pow((4.03e-10/sqrt(et)),2) + pow((0.00544/et),2));
  else if(fabs(eta)<2.500) res = sqrt(pow(0.000577,2) + pow((0.000768/sqrt(et)),2) + pow((0.00331/et),2));
  return res;
}


inline double ElectronUncertainties::phi(double et, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.174)      res = sqrt(pow(0.000101,2) + pow((0.0011/sqrt(et)),2) + pow((0.00346/et),2));
  else if(fabs(eta)<0.261) res = sqrt(pow(9.3e-05,2) + pow((0.00115/sqrt(et)),2) + pow((0.0035/et),2));
  else if(fabs(eta)<0.348) res = sqrt(pow(0.000103,2) + pow((0.00117/sqrt(et)),2) + pow((0.00333/et),2));
  else if(fabs(eta)<0.435) res = sqrt(pow(0.00011,2) + pow((0.00115/sqrt(et)),2) + pow((0.00365/et),2));
  else if(fabs(eta)<0.522) res = sqrt(pow(0.000105,2) + pow((0.00122/sqrt(et)),2) + pow((0.00343/et),2));
  else if(fabs(eta)<0.609) res = sqrt(pow(0.000102,2) + pow((0.00129/sqrt(et)),2) + pow((0.00328/et),2));
  else if(fabs(eta)<0.696) res = sqrt(pow(0.000103,2) + pow((0.00139/sqrt(et)),2) + pow((0.00253/et),2));
  else if(fabs(eta)<0.783) res = sqrt(pow(0.000115,2) + pow((0.00139/sqrt(et)),2) + pow((0.00293/et),2));
  else if(fabs(eta)<0.870) res = sqrt(pow(0.000121,2) + pow((0.00158/sqrt(et)),2) + pow((0.00151/et),2));
  else if(fabs(eta)<0.957) res = sqrt(pow(0.000128,2) + pow((0.00169/sqrt(et)),2) + pow((1.93e-08/et),2));
  else if(fabs(eta)<1.044) res = sqrt(pow(0.000145,2) + pow((0.00179/sqrt(et)),2) + pow((1.69e-08/et),2));
  else if(fabs(eta)<1.131) res = sqrt(pow(0.000185,2) + pow((0.00182/sqrt(et)),2) + pow((2.99e-09/et),2));
  else if(fabs(eta)<1.218) res = sqrt(pow(0.000194,2) + pow((0.002/sqrt(et)),2) + pow((2.39e-08/et),2));
  else if(fabs(eta)<1.305) res = sqrt(pow(0.000226,2) + pow((0.00206/sqrt(et)),2) + pow((5.88e-08/et),2));
  else if(fabs(eta)<1.392) res = sqrt(pow(0.000247,2) + pow((0.00225/sqrt(et)),2) + pow((1.47e-09/et),2));
  else if(fabs(eta)<1.479) res = sqrt(pow(0.000234,2) + pow((0.00233/sqrt(et)),2) + pow((4.92e-09/et),2));
  else if(fabs(eta)<1.653) res = sqrt(pow(0.00025,2) + pow((0.00268/sqrt(et)),2) + pow((7.5e-09/et),2));
  else if(fabs(eta)<1.740) res = sqrt(pow(0.000284,2) + pow((0.00275/sqrt(et)),2) + pow((6.56e-09/et),2));
  else if(fabs(eta)<1.830) res = sqrt(pow(0.000356,2) + pow((0.00279/sqrt(et)),2) + pow((0.00261/et),2));
  else if(fabs(eta)<1.930) res = sqrt(pow(0.000347,2) + pow((0.00298/sqrt(et)),2) + pow((1.02e-08/et),2));
  else if(fabs(eta)<2.043) res = sqrt(pow(0.000302,2) + pow((0.00322/sqrt(et)),2) + pow((5.22e-08/et),2));
  else if(fabs(eta)<2.172) res = sqrt(pow(0.000287,2) + pow((0.00349/sqrt(et)),2) + pow((3e-11/et),2));
  else if(fabs(eta)<2.322) res = sqrt(pow(0.000214,2) + pow((0.00436/sqrt(et)),2) + pow((2.98e-09/et),2));
  else if(fabs(eta)<2.500) res = sqrt(pow(8.02e-05,2) + pow((0.00525/sqrt(et)),2) + pow((0.00581/et),2));
  return res;
}
