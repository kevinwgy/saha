//Units: m, kg, s, K, A

under NavierStokesEquations {
  under Material[0] { //He,Ne,Ar mixture
    EquationOfState = StiffenedGas;
    under StiffenedGasModel {
      SpecificHeatRatio = 1.6666667;
      PressureConstant = 0;
      SpecificHeatAtConstantVolume = 1.18389e3; 
    }
  }
//  under Material[1] { //Soda-lime glass
//    EquationOfState = StiffenedGas;
//    under StiffenedGasModel {
//      SpecificHeatRatio = 3.9;
//      PressureConstant = 2.6208e9; //Pa
//    }
//  }
}

under Ionization {

  PlanckConstant = 6.62607004e-34;  // J.s
  ElectronCharge = 1.60217662e-19;  // A.s  (Coulumbs)
  ElectronMass   = 9.10938356e-31;  // kg
  BoltzmannConstant = 1.38064852e-23;  // J/K
  VacuumPermittivity = 8.85418781762e-12; // A^2 s^4 kg^-1 m^-3

  under Material[0] {
    Type = NonIdealSahaEquation; //IdealSahaEquation;
    DepressionModel = Ebeling; //Griem; //Ebeling; //None; 
    PartitionFunctionEvaluation = OnTheFly; //CubicSplineInterpolation; //LinearInterpolation;
    MaxIts = 200;
    ConvergenceTolerance = 1.0e-5;

    under Element[0] { //He
      MolarFraction = 0.3;
      MolarMass = 4.002602e-3; //kg/mol
      AtomicNumber = 2;
      MaxChargeNumber = 2;
      IonizationEnergyFile = "AtomicData_kg_m_s_K_A/I_2.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_kg_m_s_K_A/E_2_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_kg_m_s_K_A/g_2_";
      DegeneracyFilesSuffix  = ".txt";
    }

    under Element[1] { //Ne
      MolarFraction = 0.1;
      MolarMass = 20.1797e-3; //kg/mol
      AtomicNumber = 10;
      MaxChargeNumber = 10;
      IonizationEnergyFile = "AtomicData_kg_m_s_K_A/I_10.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_kg_m_s_K_A/E_10_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_kg_m_s_K_A/g_10_";
      DegeneracyFilesSuffix  = ".txt";
    }

    under Element[2] { //Ar
      MolarFraction = 0.6;
      MolarMass = 39.948e-3; //kg/mol
      AtomicNumber = 18;
      MaxChargeNumber = 18;
      IonizationEnergyFile = "AtomicData_kg_m_s_K_A/I_18.txt";
      ExcitationEnergyFilesPrefix = "AtomicData_kg_m_s_K_A/E_18_";
      ExcitationEnergyFilesSuffix = ".txt";
      DegeneracyFilesPrefix  = "AtomicData_kg_m_s_K_A/g_18_";
      DegeneracyFilesSuffix  = ".txt";
    }
  }

}

under BoundaryConditions {
  under Inlet {
    MaterialID = 0;
    Density    = 1e-6; 
    Temperature = 5.803e4; 
//    Pressure   = 5e7; //1.42843877e+07;
  }
}

under Output {
  MaxChargeNumber = 9;
  VerboseScreenOutput = High;
  MolarFractionsElement0 = On;
  MolarFractionsElement1 = On;
  MolarFractionsElement2 = On;
}








