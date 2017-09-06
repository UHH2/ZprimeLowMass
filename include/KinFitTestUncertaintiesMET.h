/* Here the uncertainties for the covariance matrix are set.
   The uncertainties come in principle from:  cmssw/TopQuarkAnalysis/TopObjectResolutions/python/stringResolutions_etEtaPhi_Fall11_cff.py */



//Resolution of the MET to set up covariance matrix
class METUncertainties {

 public:
  enum Flavor {kUds, kB};
    
  METUncertainties(){};
  ~METUncertainties(){};
  
  inline double et(double et);
  inline double eta(double et);
  inline double phi(double et);
};

inline double METUncertainties::et(double et)
{
  double res = et * (sqrt(pow(0.0337,2) + pow((0.888/sqrt(et)),2) + pow((19.6/et),2)));
  return res;
}


inline double METUncertainties::eta(double et)
{
  double res = sqrt(pow(0,2) + pow((0/sqrt(et)),2) + pow((0/et),2));
  return res;
}


inline double METUncertainties::phi(double et)
{
 double res = sqrt(pow(1.28e-08,2) + pow((1.45/sqrt(et)),2) + pow((1.03/et),2));
 return res;
}
