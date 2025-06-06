#include <time.h>
#include <Utils.h>
#include <IoData.h>
#include <VarFcnSG.h>
#include <VarFcnNASG.h>
#include <VarFcnMG.h>
#include <VarFcnJWL.h>
#include <VarFcnANEOSEx1.h>
#include <VarFcnDummy.h>
#include <SahaEquationSolver.h>
#include <NonIdealSahaEquationSolver.h>
#include <set>
using std::cout;
using std::endl;

//#include<chrono> //for timing only
//using namespace std::chrono;


/*************************************
 * Main Function
 ************************************/
int verbose = 0;
MPI_Comm m2c_comm;
int main(int argc, char* argv[])
{
  MPI_Init(NULL, NULL);
  MPI_Comm comm = MPI_COMM_WORLD;

  m2c_comm = MPI_COMM_WORLD;

  clock_t start_time = clock(); //for timing purpose only

  //! Initialize PETSc and MPI 
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m                 START                    \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
  print("\n");


  //! Read user's input file
  IoData iod(argc, argv);

  verbose = iod.output.verbose;

  iod.finalize();

  //! Initialize VarFcn (EOS, etc.) 

  std::vector<VarFcnBase *> vf;
  for(int i=0; i<(int)iod.eqs.materials.dataMap.size(); i++)
    vf.push_back(NULL); //allocate memory for the VarFcn pointers

  for(auto it = iod.eqs.materials.dataMap.begin(); it != iod.eqs.materials.dataMap.end(); it++) {
    int matid = it->first;
    if(matid < 0 || matid >= (int)vf.size()) {
      print_error("*** Error: Detected error in the specification of material indices (id = %d).\n", matid);
      exit_mpi();
    }
    if(it->second->eos == MaterialModelData::STIFFENED_GAS)
      vf[matid] = new VarFcnSG(*it->second);
    else if(it->second->eos == MaterialModelData::NOBLE_ABEL_STIFFENED_GAS)
      vf[matid] = new VarFcnNASG(*it->second);
    else if(it->second->eos == MaterialModelData::MIE_GRUNEISEN)
      vf[matid] = new VarFcnMG(*it->second);
    else if(it->second->eos == MaterialModelData::JWL)
      vf[matid] = new VarFcnJWL(*it->second);
    else if(it->second->eos == MaterialModelData::ANEOS_BIRCH_MURNAGHAN_DEBYE)
      vf[matid] = new VarFcnANEOSEx1(*it->second);
    else {
      print_error("*** Error: Unable to initialize variable functions (VarFcn) for the specified material model.\n");
      exit_mpi();
    }
  }


  // Create Saha solvers
  std::vector<SahaEquationSolver*> saha; //one "saha" for each materialid

  if(iod.ion.materialMap.dataMap.size() == 0) {
    print_error("*** Error: Unable to create ionization operator. No models specified in the input file.\n");
    exit_mpi();
  }
  int nMat = vf.size();
  saha.resize(nMat, NULL);
  for(auto it = iod.ion.materialMap.dataMap.begin(); it != iod.ion.materialMap.dataMap.end(); it++) {
    if(it->first<0 || it->first>=nMat) {
      print_error("*** Error: Ionization model specified for an unknown material id (%d).\n", it->first);
      exit_mpi();
    }
    switch (it->second->type) {
      case MaterialIonizationModel::SAHA_IDEAL :
        print("- Initializing ideal Saha Equation solver for material %d.\n", it->first);
        saha[it->first] = new SahaEquationSolver(*(it->second), iod, vf[it->first], &comm);
        break;
      case MaterialIonizationModel::SAHA_NONIDEAL :
        print("- Initializing non-ideal Saha Equation solver for material %d.\n", it->first);
        saha[it->first] = new NonIdealSahaEquationSolver(*(it->second), iod, vf[it->first], &comm);
        break;
      default :
        print_error("*** Error: Cannot initialize ionization solver for material %d. Unknown type.\n",
                    it->first);
        exit_mpi(); 
    }
  }
  for(int i=0; i<(int)saha.size(); i++) { //create dummy solvers for materials w/o ionization model
    if(saha[i] == NULL)
      saha[i] = new SahaEquationSolver(iod, vf[i]);
  }

  FILE *sols=fopen("sols.txt", "w+");
  FILE *alphas=fopen("alphas.txt", "w+");

  if (!sols || !alphas){
    fprintf(stderr, "output file could not be opened \n");
    exit(-1);
  }

  fprintf(sols, "Density Temperature Nh Zav Ne DebyeLength \n");
  fprintf(alphas, "Density Temperature j r alpha \n");

//NOTE: hack for rho/T
  double rho = iod.bc.inlet.density; //make input density as rho0
  double rhof = 1.0e2;
  double npoints = 100.0;
  double drho =  pow(rhof/rho, 1.0/npoints);
  while (rho<=rhof*drho) { //START RHO LOOP
    int id = iod.bc.inlet.materialid;

    double V[5];
    V[0] = rho;
    V[1] = iod.bc.inlet.velocity_x;
    V[2] = iod.bc.inlet.velocity_y;
    V[3] = iod.bc.inlet.velocity_z;

    // Read temperature from input, then compute pressure 
    // (converted back to temperature in the solver)
    double T = iod.bc.inlet.temperature;
    double e = vf[id]->GetInternalEnergyPerUnitMassFromTemperature(V[0], T);
    V[4] = vf[id]->GetPressure(V[0], e);

    //Verification
    double ecal = vf[id]->GetInternalEnergyPerUnitMass(V[0],V[4]);
    double Tcal = vf[id]->GetTemperature(V[0], ecal);
    if(fabs(T-Tcal)/fabs(T)>1.0e-6) {
      fprintf(stderr,"Error: Calculated temperature (%e) different from input (%e). Check EOS.\n",
		      Tcal, T);
      exit(-1);
    }

    print("Solving the Saha Ionization Equation...\n");
    print("Input: density = %e, pressure = %e (MaterialID: %d).\n",V[0],V[4],id);
    e = vf[id]->GetInternalEnergyPerUnitMass(V[0], V[4]);
    print("By EOS: e = %e.\n", e);

    int max_charge_in_output = iod.output.max_charge_number;

    std::map<int, vector<double> > nodal_alphas;
    for(int iSpecies=0; iSpecies<OutputData::MAXSPECIES; iSpecies++) 
      if(iod.output.molar_fractions[iSpecies] == OutputData::ON) 
        nodal_alphas[iSpecies] = vector<double>(max_charge_in_output+2);

    double zav, nh, ne, lambda;
    saha[id]->Solve(V, zav, nh, ne, nodal_alphas, &lambda);

    print("\n");
    print("Solution: Zav = %e, Nh = %e, Ne = %e.\n", zav, nh, ne);
    fprintf(sols, "%e %e %e %e %e %e \n", rho, T, nh, zav,  ne, lambda);

    for(auto it = nodal_alphas.begin(); it != nodal_alphas.end(); it++) {
//     print("  - Species %d:\n", it->first);
      for(int i=0; i<(int)it->second.size()-1; i++){
//        print("    o Molar fraction of charged state %d: %e.\n", i, it->second[i]);
        fprintf(alphas, "%e %e %d %d %e\n", rho, T, it->first, i, it->second[i]);
      }
//      print("    o Molar fraction of charged state %d+: %e.\n", it->second.size()-1, it->second[it->second.size()-1]);
      fprintf(alphas, "%e %e %d %ld %e\n", rho, T, it->first, it->second.size()-1, it->second[it->second.size()-1]);
    }

    rho=rho*drho;

  } // end rho loop


  fclose(sols);
  fclose(alphas);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m           NORMAL TERMINATION             \033[0m\n"); 
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");

  for(int i=0; i<(int)vf.size(); i++)
    delete vf[i];

  for(auto it = saha.begin(); it != saha.end(); it++)
    if(*it)
      delete *it;

  return 0;
}

//--------------------------------------------------------------
