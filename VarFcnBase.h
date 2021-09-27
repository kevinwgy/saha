#ifndef _VAR_FCN_BASE_H_
#define _VAR_FCN_BASE_H_

#include <IoData.h>
#include <Vector3D.h>
#include <cmath>
#include <iostream>

/****************************************************************************
 * This class is the base class for the VarFcn classes where EOS can be
 * an arbitrary convex equation of state.
 * Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state, since it is assumed that the EOS that must be used at this point 
 * is known.
 * As a base class, it is also assumed that the Euler equations are 
 * considered, meaning there are a priori only 5 variables that
 * define each state. Thus for flows with more than 5 variables,
 * some transformation operators must be overwritten in
 * the appropriate class.
 *
 * lay-out of the base class is:
 *  - 1 -  Transformation Operators
 *  - 2 -  General Functions
 *  - 3 -  Equations of State Parameters
 *  - 4 -  EOS related functions
 ***************************************************************************/

extern int verbose;

class VarFcnBase {

public:
  
  enum Type{ STIFFENED_GAS = 0, MIE_GRUNEISEN = 1, JWL = 2} type;

  double rhomin,pmin;

  double failsafe_density;

  VarFcnBase(MaterialModelData &data) {
    rhomin = data.rhomin;
    pmin = data.pmin;
    failsafe_density = data.failsafe_density;
  }
  virtual ~VarFcnBase() {}
 
  //----- EOS-Specific Functions -----//
  //! get pressure from density (rho) and internal energy per unit mass (e)
  virtual double GetPressure(double rho, double e) const{
    print_error("*** Error:  GetPressure Function not defined\n");
    exit(-1); return 0.0;}

  //! get e (internal energy per unit mass) from density (rho) and pressure (p)
  virtual double GetInternalEnergyPerUnitMass(double rho, double p) const{
    print_error("*** Error:  GetInternalEnergyPerUnitMass Function not defined\n");
    exit(-1); return 0.0;}

  //! get rho (density) from p (pressure) and p (internal energy per unit mass)
  virtual double GetDensity(double p, double e) const{
    print_error("*** Error:  GetDensity Function not defined\n");
    exit(-1); return 0.0;}

  //! dpdrho = \frac{\partial p(\rho,e)}{\partial \rho}
  virtual double GetDpdrho(double rho, double e) const{
    print_error("*** Error:  GetDpdrho Function not defined\n");
    exit(-1); return 0.0;}

  //! BigGamma = 1/rho*(\frac{\partial p(\rho,e)}{\partial e})
  //  It is called "BigGamma" to distinguish it from the small "gamma" in perfect and stiffened EOS.
  virtual double GetBigGamma(double rho, double e) const{
    print_error("*** Error:  GetBigGamma Function not defined\n");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetTemperature(double rho, double e) const{
    print_error("*** Error:  GetTemperature Function not defined\n");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetReferenceTemperature() const{
    print_error("*** Error:  GetReferenceTemperature Function not defined\n");
    exit(-1); return 0.0;}

  //! temperature law, defined separately for each EOS
  virtual double GetInternalEnergyPerUnitMassFromTemperature(double rho, double T) const{
    print_error("*** Error:  GetInternalEnergyPerUnitMassFromTemperature Function not defined\n");
    exit(-1); return 0.0;}

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double rho, double p) const{
    if(m2c_isnan(rho) || m2c_isnan(p)) {
      fprintf(stderr, "*** Error: CheckState failed. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    if(rho <= 0.0) {
      if(verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    double e = GetInternalEnergyPerUnitMass(rho,p);
    double c2 = GetDpdrho(rho, e) + p/rho*GetBigGamma(rho, e);
    if(c2<=0){
      if(verbose>1)
        fprintf(stdout, "Warning: Negative density or violation of hyperbolicity. rho = %e, p = %e.\n", rho, p);
      return true;
    }
    return false;
  }

  //checks that the Euler equations are still hyperbolic
  virtual bool CheckState(double *V) const{
    if(m2c_isnan(V[0]) || m2c_isnan(V[1]) || m2c_isnan(V[2]) || m2c_isnan(V[3]) || m2c_isnan(V[4])) {
      fprintf(stderr, "*** Error: CheckState failed. V = %e %e %e %e %e\n", V[0], V[1], V[2], V[3], V[4]);
      return true;
    }
    return CheckState(V[0], V[4]); 
  }
 
  //check for phase transitions
  virtual bool CheckPhaseTransition(int id/*id of the other phase*/) const{
    return false; //by default, phase transition is not allowed/considered
  }

  //----- Transformation Operators -----//
  inline void ConservativeToPrimitive(double *U, double *V); 
  inline void PrimitiveToConservative(double *V, double *U);

  //----- General Functions -----//
  inline int GetType() const{ return type; }

  inline double ComputeSoundSpeed(double rho, double e);
  inline double ComputeSoundSpeedSquare(double rho, double e); //!< this one does not crash on negative c^2
  inline double ComputeMachNumber(double *V);
  inline double ComputeEnthalpyPerUnitMass(double rho, double p); //!< h = e + p/rho
  inline double ComputeTotalEnthalpyPerUnitMass(double *V); //!< H = 1/rho*(E + p)

  // Clipping
  inline bool ClipDensityAndPressure(double *V, double *U = 0);
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline
void VarFcnBase::ConservativeToPrimitive(double *U, double *V)
{
  V[0] = U[0];

  double invRho = 1.0 / U[0];

  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;

  double e = (U[4] - 0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])) * invRho;
  V[4] = GetPressure(V[0], e);
}

//------------------------------------------------------------------------------

inline
void VarFcnBase::PrimitiveToConservative(double *V, double *U)
{
  U[0] = V[0];

  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];

  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); //pass in rho and p
  U[4] = V[0]*(e + 0.5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]));
}

//------------------------------------------------------------------------------

inline 
double VarFcnBase::ComputeSoundSpeed(double rho, double e)
{
  double c2 = GetDpdrho(rho, e) + GetPressure(rho,e)/rho*GetBigGamma(rho, e);
  if(c2<=0) {
    fprintf(stderr,"*** Error: Cannot calculate speed of sound (Square-root of a negative number): rho = %e, e = %e.\n",
            rho, e);
    exit(-1);
  }
  return sqrt(c2);
}

//------------------------------------------------------------------------------

inline 
double VarFcnBase::ComputeSoundSpeedSquare(double rho, double e)
{
  return GetDpdrho(rho, e) + GetPressure(rho,e)/rho*GetBigGamma(rho, e);
}

//------------------------------------------------------------------------------

inline 
double VarFcnBase::ComputeMachNumber(double *V)
{
  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); 
  double c = ComputeSoundSpeedSquare(V[0], e);

  if(c<0) {
    fprintf(stderr,"*** Error: c^2 (square of sound speed) = %e in ComputeMachNumber. V = %e, %e, %e, %e, %e.\n",
            c, V[0], V[1], V[2], V[3], V[4]);
    exit(-1);
  } else
    c = sqrt(c);

  return sqrt(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])/c;
}

//------------------------------------------------------------------------------

inline
double VarFcnBase::ComputeEnthalpyPerUnitMass(double rho, double p)
{
  return GetInternalEnergyPerUnitMass(rho,p) + p/rho;
}

//------------------------------------------------------------------------------

inline 
double VarFcnBase::ComputeTotalEnthalpyPerUnitMass(double *V) //!< H = 1/rho*(E + p)
{
  double e = GetInternalEnergyPerUnitMass(V[0],V[4]); 
  return e + 0.5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]) + V[4]/V[0];
}

//------------------------------------------------------------------------------

inline
bool VarFcnBase::ClipDensityAndPressure(double *V, double *U)
{
//verification of density and pressure value
//if pressure/density < pmin/rhomin, set pressure/density to pmin/rhomin
//and rewrite V and U!!
  bool clip = false;

  if(V[0]<rhomin){
//    if(verbose)
//      fprintf(stderr,"clip density from %e to %e.\n", V[0], rhomin);
    V[0] = rhomin;
    clip = true;
  }

  if(V[4]<pmin){
//    if(verbose)
//      fprintf(stdout, "clip pressure from %e to %e\n", V[4], pmin);
    V[4] = pmin;
    clip = true;
  }

  if(clip && U) //also modify U
    PrimitiveToConservative(V,U);

  return clip;
}

//------------------------------------------------------------------------------



#endif
