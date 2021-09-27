#ifndef _VAR_FCN_MG_H
#define _VAR_FCN_MG_H

#include <VarFcnBase.h>
#include <fstream>

/********************************************************************************
 * This class is the VarFcn class for the Mie-Gruneisen EOS in Euler
 * Equations. Only elementary functions are declared and/or defined here.
 * All arguments must be pertinent to only a single grid node or a single
 * state.
 *
 * EOS: Pressure = rho0*c0*c0*eta*(1-Gamma0/2*eta)/(1-s*eta)^2 + rho0*Gamma0*(e-e0)
 *
 * where
 *
 *   e      : internal energy per unit mass.
 *   eta    : 1 - rho0/rho
 *   rho    : density
 *
 *   rho0   : ref. density
 *   c0     : bulk speed of sound
 *   Gamma0 : Gruneisen coefficient at ref. state
 *   s      : slope of shock Hugoniot
 *   e0     : internal energy at ref. state
 *
 * References: Shafquat Islam's report (01/2021), Allen Robinson's technical report (2019)
 *
 *   Note: A temperature law can be defined as e = cv*(T-T0), where cv is a constant
 *         specific heat (at constant volume), and T0 a ref. temperature
 *
 ********************************************************************************/

class VarFcnMG : public VarFcnBase {

private:
  double rho0;
  double c0;
  double Gamma0;
  double s;
  double e0;

  double cv;

  double rho0_c0_c0;    //!< rho0*c0*c0
  double Gamma0_over_2; //!< Gamma0/2
  double Gamma0_rho0;   //!< Gamma0*rho0

public:
  VarFcnMG(MaterialModelData &data);
  ~VarFcnMG() {}

  //! ----- EOS-Specific Functions -----
  inline double GetPressure(double rho, double e) const {
    double eta = 1.0 - rho0/rho;
    return rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta)) + Gamma0_rho0*(e-e0);}

  inline double GetInternalEnergyPerUnitMass(double rho, double p) const {
    double eta = 1.0 - rho0/rho;
    return (p - rho0_c0_c0*eta*(1.0 - Gamma0_over_2*eta)/((1.0-s*eta)*(1.0-s*eta)))/Gamma0_rho0 + e0;
  }

  inline double GetDensity(double p, double e) const;

  inline double GetDpdrho(double rho, double e) const {
    double eta = 1.0 - rho0/rho;
    double S = 1.0 - s*eta;
    return rho0_c0_c0*(1.0 + (s - Gamma0)*eta)/(S*S*S)*rho0/(rho*rho);
  }

  inline double GetBigGamma(double rho, double e) const {return Gamma0_rho0/rho;}

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnMG::VarFcnMG(MaterialModelData &data) : VarFcnBase(data) {

  if(data.eos != MaterialModelData::MIE_GRUNEISEN){
    fprintf(stderr, "*** Error: MaterialModelData is not of type Mie-Gruneisen.\n");
    exit(-1);
  }

  type   = MIE_GRUNEISEN;

  rho0   = data.mgModel.rho0;
  c0     = data.mgModel.c0;
  Gamma0 = data.mgModel.Gamma0;
  s      = data.mgModel.s;
  e0     = data.mgModel.e0;

  cv     = data.mgModel.cv;

  rho0_c0_c0 = rho0*c0*c0;
  Gamma0_over_2 = 0.5*Gamma0;
  Gamma0_rho0 = Gamma0*rho0;
}

//------------------------------------------------------------------------------

inline
double VarFcnMG::GetDensity(double p, double e) const {

  //solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
  double c = p - Gamma0_rho0*(e - e0);
  double a = c*s*s + Gamma0_over_2*rho0_c0_c0;
  double b = -(2.0*c*s + rho0_c0_c0);

  if(a==0) {//linear equation: eta = -c/b, rho = rho0/(1-eta)

    return rho0/(1 + c/b);

  } else {

    double b2m4ac = b*b - 4.0*a*c;
    if(b2m4ac<0) {
      fprintf(stderr, "*** Error: The M-G EOS is invalid for the given p(%e) and e(%e) --- unable to solve it for rho.\n",
              p, e);
      exit(-1);
    }

    b2m4ac = sqrt(b2m4ac);
    double rho1 = rho0/(1.0 - (-b + b2m4ac)/(2.0*a));
    double rho2 = rho0/(1.0 - (-b - b2m4ac)/(2.0*a));

    if(rho1>rho2)
      std::swap(rho1,rho2); //rho2 should be the bigger one

    if(rho1>0 || rho2<0) { //both are negative, or both are positive
      fprintf(stderr, "*** Error: Cannot determine the solution (rho) of the M-G EOS (rho1 = %e, rho2 = %e). \n", 
              rho1, rho2);
      exit(-1);
    }
    
    return rho2; 

  }
}


//------------------------------------------------------------------------------

#endif
