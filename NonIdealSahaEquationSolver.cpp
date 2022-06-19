#include<NonIdealSahaEquationSolver.h>
#include<fstream>
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
using std::vector;
using std::map;
using std::pair;

extern double avogadro_number;

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::NonIdealSahaEquationSolver(IoData& iod_, VarFcnBase* vf_)
                          : SahaEquationSolver(Iod_, vf_),
                            eps0(iod_.ion.vacuum_permittivity), factor_deltaI(0.0)
{ }

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::NonIdealSahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, 
                                IoData& iod_, VarFcnBase* vf_, MPI_Comm* comm)
                          : SahaEquationSolver(iod_ion_mat_, iod_, vf_, comm),
                            eps0(iod_.ion.vacuum_permittivity)
{
  pi = 2.0*acos(0);
  factor_deltaI = e*e/(4.0*pi*eps0)
  factor_LambB  = h/sqrt(2.0*pi*me*kb);

  f.resize(elem.size(), vector<double>());
  alpha.resize(elem.size(), vector<double>());
  for(int j=0; j<f.size(); j++) {
    f[j].resize(elem[j].rmax+1, 0.0);
    alpha[j].resize(elem[j].rmax+1, 0.0);
  }
}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::~NonIdealSahaEquationSolver()
{ }

//--------------------------------------------------------------------------

double 
NonIdealSahaEquationSolver::ComputeDeltaI(int r, int j, double T, double one_over_lambD)
{
  double dI(0.0);
  if(iod_ion_mat->depression == MaterialIonizationModel::GRIEM) {
    dI = (r+1.0)*factor_deltaI*one_over_lambD;
  } else {//Ebeling
    if(one_over_lambD==0.0)
      dI = 0.0;
    else 
      dI = (r+1.0)*factor_deltaI/(1.0/one_over_lambD + 0.125*factor_LambB/sqrt(T));
  }
  assert(isfinite(deltaI));
  return dI;
}

//--------------------------------------------------------------------------

double 
NonIdealSahaEquationSolver::ComputeStateForElement(int j, double T, double nh, double zav, 
                                                   double one_over_lambD, bool compute_alpha)
{
  assert(zav>=0.0 && nh>0.0);

  if(zav==0.0) {
    double zej = elem[j].molar_fraction*elem[j].rmax;
    if(compute_alpha) {
      alpha[j][0] = elem[j].molar_fraction;
      for(int r=1; r<=rmax; r++)
        alpha[j][r] = 0.0;
    }
    return zej;
  }


  // Now, zav must be positive
  
  double kbT = kb*T;
  double fcore = pow( (2.0*pi*(me/h)*(kbT/h)), 1.5)/nh;

  int rmax = elem[j].rmax;

  // compute f_{r,j}, r = 1, ..., rmax
  double Ur0(0.0), Ur1(0.0), deltaI0(0.0), deltaI1(0.0);
  for(int r=0; r<=rmax; r++) {
    deltaI0 = deltaI1;
    Ur0     = Ur1;
    deltaI1 = ComputeDeltaI(r,j,T,one_over_lambD);
    Ur1     = elem[j].CalculatePartitionFunction(r, T, deltaI1);

    f[j][r] = r==0 ? 0.0 //f[j][0] is not used anyway
                   : 2.0*Ur1/Ur0*fcore*exp(-(elem[j].I[r-1]-deltaI0)/kbT);
  }

  vector<double> zav_power(rmax+1, 0.0);
  zav_power[0] = 1.0;
  for(int r=1; r<=rmax; r++)
    zav_power[r] = zav_power[r-1]*zav;

  double denominator = 0.0; 
  double numerator = 0.0;
  double fprod = 1.0;

  for(int r=1; r<=rmax; r++) {
    fprod *= f[j][r]*zav_power[rmax-r];
    numerator += (double)r*fprod; 
    denominator += fprod;
  }

  denominator += zav_power[rmax];

  //check denominator and compute zej
  double zej(0.0);
  if(denominator==0.0) //zav = 0 or extremely close to 0)
    zej = 0.0;
  else { 
    zej = numerator/denominator;
    if(!std::isfinite(zej)) {
      fprintf(stderr,"\033[0;31m*** Error: Non-Ideal Saha equation solver failed. "
                     "Z for element %d is not finite.\n\033[0m", j);
      exit(-1);
    }
    if(zej<0.0) {
      fprintf("\033[0;35mWarning: Found negative Z (%e) for element %d. Setting it to 0.\n\033[0m",
              zej, j);
      zej = 0.0;
    } else if (zej>rmax){
      fprintf("\033[0;35mWarning: Found Z greater than rmax (%e vs. %d) for element %d. Setting it to %d.\n\033[0m",
              zej, rmax, j, rmax);
      zej = (double)rmax;
    }
  }
  zej *= elem[j].molar_fraction;


  //compute molar fractions if needed
  if(compute_alpha) {
    assert(std::isfinite(numerator)); 
    if(numerator==0.0) {//zej must be 0
      assert(zej==0.0);
      alpha[j][0] = elem[j].molar_fraction;
      for(int r=1; r<=rmax; r++)
        alpha[j][r] = 0.0;
    }
    alpha[j][0] = std::min(zej*zav_power[rmax]/numerator, elem[j].molar_fraction);
    double ne = zav*nh;
    double summation = alpha[j][0];
    for(int r=1; r<=rmax; r++) {
      alpha[j][r] = alpha[j][r-1]/ne*f[j][r];
      summation += alpha[j][r];
      if(summation>elem[j].molar_fraction) {
        alpha[j][r] -= (summation - elem[j].molar_fraction)
        summation = elem[j].molar_fraction;
        for(rr=r+1; rr<=rmax; rr++)
          alpha[j][rr] = 0.0;
        break;
      }
    } 

    //make sure the sum of alphas is equal to molar_fraction
    double ratio = elem[j].molar_fraction/summation;
    if(ratio != 1.0) {
      for(auto&& alph : alpha[j]) 
        alph *= ratio;
    }
  }

  return zej;
}

//--------------------------------------------------------------------------


void
SahaEquationSolver::Solve(double* v, double& zav, double& nh, double& ne, 
                          map<int, vector<double> >& alpha_rj)
{

  if(!iod_ion_mat) { //dummy solver 
    zav = 0.0;
    ne = 0.0;
    nh = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = (it==alpha_rj.begin() && r==0) ? 1.0 : 0.0;
    }
    return;
  }


  double T = vf->GetTemperature(v[0], vf->GetInternalEnergyPerUnitMass(v[0], v[4]));
  nh = v[0]/molar_mass*avogadro_number;

  if(T<=Tmin) { //no ionization
    zav = 0.0;
    ne = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = (r==0) ? elem[it->first].molar_fraction : 0.0;
    }
    return;
  }

  // ------------------------------
  // Step 1: Solve for Zav 
  // ------------------------------
  ZavEquation fun(*this, T, nh);

  //Find initial bracketing interval (zav0, zav1)
  double zav0, zav1, f0, f1; 
  zav0 = 0.0;
  zav1 = max_atomic_number; //zav1>zav0
  f0 = fun(zav0);
  bool found_initial_interval = false;
  for(int i=0; i<iod_ion_mat->maxIts; i++) {
    f1 = fun(zav1);
    if(f0*f1<=0.0) {
      found_initial_interval = true;
      break;
    }
    zav1 /= 2.0;
  }
  if(!found_initial_interval) {
    fprintf(stderr,"\033[0;31m*** Error: Saha equation solver failed. "
            "Cannot find an initial bracketing interval. (p = %e, T = %e)\n\033[0m", v[4], T);
    exit(-1);
  }

  // Calling boost function for root-finding
  // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT  
  boost::uintmax_t maxit = iod_ion_mat->maxIts;
  double tol = iod_ion_mat->convergence_tol;
  if(f0==0.0) {
    zav = zav0; maxit = 0;
  } else if(f1==0.0) {
    zav = zav1; maxit = 0;
  } else {
    for(int trial = 0; trial < 2; trial++) {
      pair<double,double> sol; 
      sol = toms748_solve(fun, zav0, zav1, f0, f1,
                          [=](double r0, double r1){return r1-r0<std::min(tol,0.001*(zav1-zav0));},
                          maxit);
      zav = 0.5*(sol.first + sol.second);
      if(zav>=0) break;

      // fail-safe
      if(trial>0) break;

      if(!isfinite(sol.first) || sol.first<0)
        sol.first = 0.0;
      if(!isfinite(sol.second) || sol.second<sol.first)
        sol.second = zav1;
      zav0 = sol.first;
      zav1 = sol.second;
      f0 = fun(zav0);
      f1 = fun(zav1); 
      if(f0*f1>0)
        break;
    }
  }

  if(!(zav>0)) {
    zav = 0.0;
    ne = 0.0;
    for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
      vector<double> &alpha = it->second;
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = (r==0) ? elem[it->first].molar_fraction : 0.0;
    }
    return;
  }

  //*******************************************************************

#if DEBUG_SAHA_SOLVER == 1
  fprintf(stderr,"-- Saha equation solver converged in %d iterations, Zav = %.12e.\n", (int)maxit, zav);
#endif

  I AM HERE!!!

  //post-processing.
  ne = zav*nh;

  for(auto it = alpha_rj.begin(); it != alpha_rj.end(); it++) {
    int j = it->first; //element id
    vector<double> &alpha = it->second; //alpha_r

    if(j>=elem.size()) {//this material does not have element j
      for(int r=0; r<alpha.size(); r++)
        alpha[r] = 0.0;
      continue;
    }

    double zej = fun.GetZej(zav, j);
    double denom = 0.0;
    double zav_power = 1.0;
    for(int i=1; i<=elem[j].rmax; i++) {
      zav_power *= zav;
      denom += (double)i/zav_power*fun.GetFProd(i,j);
    }

    if(denom>0)
      alpha[0] = zej/denom;
    else {
      alpha[0] = 1.0;
      for(int r=1; r<alpha.size(); r++)
        alpha[r] = 0.0;
      return;
    }
 
    double fr(0.0);
    for(int r=1; r<alpha.size()-1; r++) {
      fr = (fun.GetFProd(r-1,j) == 0.0) ? 0.0 : fun.GetFProd(r,j)/fun.GetFProd(r-1,j);
      alpha[r] = (r<=elem[j].rmax) ? alpha[r-1]/zav*fr : 0.0;
    }

    int last_one = alpha.size()-1; 
    alpha[last_one] = elem[j].molar_fraction;
    for(int r=0; r<last_one; r++)
      alpha[last_one] -= alpha[r];

    //allow some roundoff error
    //assert(alpha[last_one]>=-1.0e-4);
    if(alpha[last_one]<0)
      alpha[last_one] = 0;
  }

}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::
LambDEquation::LambDEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_, double zav_)
             : saha(saha_), T(T_), nh(nh_), zav(zav_)
{
  assert(zav>=0.0 && nh>0.0); 

  factor_lambdaD = saha.e*saha.e/(saha.eps0*saha.kb*T);
}

//--------------------------------------------------------------------------

double NonIdealSahaEquationSolver::
LambDEquation::ComputeRHS(double one_over_lambD)
{
  double summation = zav; 
  for(int j=0; j<elem.size(); j++) {
    saha.ComputeStateForElement(j, T, nh, zav, one_over_lambD, true);  //we need alphas 
    for(int r=1; r<=elem[j].rmax; r++)
      summation += r*r*saha.alpha[j][r];
  }

  return factor_lambdaD*sqrt(nh*summation);
}

//--------------------------------------------------------------------------

NonIdealSahaEquationSolver::
ZavEquation::ZavEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_, double lambD_)
           : saha(saha_), T(T_), nh(nh_), lambD(lambD_)
{
  assert(zav>=0.0 && nh>0.0);
}

//--------------------------------------------------------------------------

double NonIdealSahaEquationSolver::
ZavEquation::ComputeRHS(double zav)
{
  double rhs = 0.0;
  for(int j=0; j<elem.size(); j++)
    rhs += saha.ComputeStateForElement(j, T, nh, zav, lambD, false); //we don't need alphas here.

  return rhs;
}

//--------------------------------------------------------------------------






//--------------------------------------------------------------------------

