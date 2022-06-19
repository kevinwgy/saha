#ifndef _SAHA_EQUATION_SOLVER_H_
#define _SAHA_EQUATION_SOLVER_H_

#include<SahaEquationSolver.h>
#include<NonIdealAtomicIonizationData.h>

//----------------------------------------------------------------
// Class NonIdealSahaEquationSolver is responsible for solving the 
// non-ideal Saha equation for one material (w/ a fixed id)
// Input: v and id 
// Output: Zav, ne, nh, alphas.
// (A dummy solver is defined for materials not undergoing ionization)
//----------------------------------------------------------------

class NonIdealSahaEquationSolver : public SahaEquationSolver {

  //! constants
  double eps0; //!< vacuum permittivity;
  double pi;
  double factor_deltaI, factor_LambB;

  //! variables for *temporary* use
  std::vector<std::vector<double> > f; //!< f_{r,j}, j=0,...,elem.size-1, r=0,...,rmax(j), (f[0][j]=0, not used)
  std::vector<std::vector<double> > alpha; //!< alpha_{r,j}, same dimensions as f

public:

  NonIdealSahaEquationSolver(IoData& iod, VarFcnBase* vf_); //!< creates a dummy solver

  NonIdealSahaEquationSolver(MaterialIonizationModel& iod_ion_mat_, IoData& iod_, VarFcnBase* vf_, MPI_Comm* comm);

  ~NonIdealSahaEquationSolver();

  void Solve(double* v, double& zav, double& nh, double& ne, std::map<int, std::vector<double> >& alpha_rj);

protected:

  // computes the depression of ionization energy, for a given Debye length lambD and temperature T
  double ComputeDeltaI(int r, int j, double T, double one_over_lambD);

  // returns Zej (for a given j), fills "f" if zav!=0. Also fills "alpha" if compute_alpha == true
  double ComputeStateForElement(int j, double T, double nh, double zav, 
                                double one_over_lambD, bool compute_alpha);


  //! nested class / functor: nonlinear equation for lambD (Debye length), given Zav
  //! note that the independent variable is actually 1/lambD, not lambD (which can be inf)
  class LambDEquation { 
    double factor_lambD;
    double T, nh, zav;
    NonIdealSahaEquationSolver& saha;
  public:
    LambDEquation(NonIdealSahaEquationSolver& saha_, double T_, double nh_, double zav_);
    ~LambDEquation() {}
    double operator() (double one_over_lambD) {return one_over_lambD - ComputeRHS(one_over_lambD);}
  private:
    double ComputeRHS(double one_over_lambD);
  };


  //! nested class / functor: nonlinear equation for Zav, given lambD
  class ZavEquation {
    double T, nh, lambD;
    NonIdealSahaEquationSolver& saha;
  public:
    ZavEquation(NonIdealSahaEquationSolver &saha_, double T_, double nh_, double lambD_);
    ~ZavEquation() {}
    double operator() (double zav) {return zav - ComputeRHS(zav);}
  private:
    double ComputeRHS(double zav); //!< compute the right-hand-side of the Zav equation
  };


};



#endif
