/* ***************************************************************************** 
Class for input parameters
***************************************************************************** */
 
#include <string.h>
#include <iostream>
#include <fstream>
#include "Param.h"

using namespace std;

double PI_val=3.1415926535897932384626433832795;

void Param::readFile(string _file)
{
  filename = _file;

  ifstream testfile;
  // check if file exists
  testfile.open(_file.c_str());
  if(!testfile.good()) {
    cerr << " !!bio_abm_fem - ERROR!!: file " << _file << " not found." << endl;
    exit(1);
  }
   
  GetPot ifile(_file.c_str());

  // read input
  cout << " --- read file (GetPot) " << _file << endl;

  // [coupling]
  fileCells2FEM = ifile("coupling/fileCells2FEM","cells.txt");
  fileCellsDensity2FEM = ifile("coupling/fileCellsDensity2FEM","cell_density.txt");
  fileFEM2Cells = ifile("coupling/fileFEM2Cells","concentration_O2.txt");
  
 
  // [cells] --- cell parameters
  /// @todo change name of initial cell array size
  n_initial_cells  = ifile("cells/n_initial_cells",1);

   cout << " HERE " << endl;
  n_phenotypes = ifile("cells/n_phenotypes",1);
  
  n_steps = ifile("cells/n_steps",0); 
  //_space_scale = ifile("cells/space_scale",1.);
  time_step = ifile("cells/time_step",1.);
  time_death = ifile("cells/time_death",2880.);
  cout << " HERE " << endl;
  radius = ifile("cells/radius", 5.0);
  // parameters for more than one phenotype
  growth_rate.resize(n_phenotypes);
  cout << " HERE " << endl;
  alpha_birthrate.resize(n_phenotypes);
  alpha_YoungM.resize(n_phenotypes);
  alpha_PoissonNo.resize(n_phenotypes);
  alpha_gcm.resize(n_phenotypes);
  adhesion_value.resize(n_phenotypes);

  Gcm = ifile("cells/Gcm",0.01);
  // alternative growth rate (two phenotypes)
  for (int k=0; k<n_phenotypes; k++) {
    growth_rate[k] = ifile("cells/growth_rate",0.1,k);
    alpha_birthrate[k] = ifile("cells/alpha_birthrate",1e-3,k); 
    alpha_YoungM[k] = ifile("cells/alpha_YoungM",1e-3,k);
    alpha_PoissonNo[k] = ifile("cells/alpha_PoissonNo",.5,k);
    alpha_gcm[k] = ifile("cells/alpha_gcm",1.,k);
    adhesion_value[k] = ifile("cells/adhesion_value",3.72e-4,k);
  }

  compressibility = ifile("cells/compressibility", 6.0);
  contact_inhibition = ifile("cells/contact_inhibition",16.);
  // birth energy parameters
  hypoxic_birth = ifile("cells/hypoxic_birth",6.0e-5); // TOMMASO
  normoxic_birth = ifile("cells/normoxic_birth",6.0e-4); // TOMMASO
  death = ifile("cells/death",0.0); // TOMMASO
  hypoxic_friction = ifile("cells/hypoxic_friction",1.0); // TOMMASO
  be_displacement = ifile("cells/be_displacement",1.5);
  be_multiplier = ifile("cells/be_multiplier",12.0);
  variance_motion = ifile("cells/variance_motion",4e-3);
  variance_phenotype = ifile("cells/variance_phenotype",0.);
  variance_adhesion = ifile("cells/variance_adhesion",0.);
  polarity_x.resize(n_initial_cells);
  polarity_y.resize(n_initial_cells);
  polarity_z.resize(n_initial_cells);
  ic_phenotype.resize(n_initial_cells);
  ic_follower_leader.resize(n_initial_cells);
 
  for (int i=0; i<n_initial_cells; i++) {
    polarity_x[i] = ifile("cells/polarity_x",0.,i);
    polarity_y[i] = ifile("cells/polarity_y",0.,i);
    if (dimension == 3) {
      polarity_z[i] = ifile("cells/polarity_z",0.,i);
    } else {
      polarity_z[i] = 0.;
    }
    ic_phenotype[i] = ifile("cells/ic_phenotype",0,i);
    ic_follower_leader[i] = ifile("cells/ic_follower_leader",1,i);
  }
  follower_force = ifile("cells/follower_force",0.);
  follower_denominator = ifile("cells/follower_denominator",0.);

  // [mutations]
  initial_phenotype = ifile("mutations/initial_phenotype",0.0);
  mutation_amount = ifile("mutations/mutation_amount",0.05);
  mutation_probability = ifile("mutations/mutation_probability",0.01);
  alpha_s = ifile("mutations/alpha_s",1.5);

  // [oxygen]
  oxygen_response = ifile("oxygen/oxygen_response",.0);
  threshold_death = ifile("oxygen/threshold_death",.7);
  threshold_hypo = ifile("oxygen/threshold_hypo",7.);
  
  // [fem] --- pde (oxygen)
  femSolverType = ifile("fem/femSolverType",0);
  meshdir = ifile("fem/meshdir","./");
  meshname = ifile("fem/meshname","rectangle_7x4x1.5_4Knodes.mesh");	
  FreeFemCall = ifile("fem/FreeFemCall","FreeFem++");
  FreeFemFile = ifile("fem/FreeFemFile","diffusion3d.edp");
  ic_file_cells = ifile("cells/ic_file_cells","nil");

  // [geo] --- boxes and geometry
  dimension = ifile("geo/dimension",3);
  boxesx = ifile("geo/boxesx",-1000);
  boxesy = ifile("geo/boxesy",-1000);
  lattice_length_x = ifile("geo/lattice_length_x",100.1);
  lattice_length_y = ifile("geo/lattice_length_y",100.1);
  if (dimension==3) {
    boxesz = ifile("geo/boxesz",-1000);
    lattice_length_z = ifile("geo/lattice_length_z",100.1);
  } else {
    lattice_length_z = 1.1;
    boxesz = 1;
  }
  ic_cell_x.resize(n_initial_cells);
  ic_cell_y.resize(n_initial_cells);
  ic_cell_z.resize(n_initial_cells);
  for (int i=0; i<n_initial_cells; i++) {
    ic_cell_x[i] = ifile("geo/ic_cell_x",lattice_length_x/2.,i);
    ic_cell_y[i] = ifile("geo/ic_cell_y",lattice_length_y/2.,i);
    if (dimension == 3) {
      ic_cell_z[i] = ifile("geo/ic_cell_z",lattice_length_z/2.,i);
    } else {
      ic_cell_z[i] = 0.;
    }
  }
  max_cell = ifile("geo/max_cell",100000.);
 

  // [fibres] --- fibre parameters
  ///@todo support general number of subdomains
  n_sub_domain = 1;
  ///@todo make this parameters more general
  x_start.resize(n_sub_domain);
  x_end.resize(n_sub_domain);

  x_start[0] = 0.0;
  x_end[0] = lattice_length_x;
  
  if (n_sub_domain == 2) {
    x_start[1] = lattice_length_x/2.0;
    x_end[0] = lattice_length_x/2.0;
    x_end[1] = lattice_length_x;
  }
      
  n_initial_fibres.resize(2);
  fibre_orientation_distribution.resize(2);
  fibre_orientation_mean_phi.resize(2);
  fibre_orientation_variance_phi.resize(2);
  fibre_orientation_mean_theta.resize(2);
  fibre_orientation_variance_theta.resize(2);
  fibre_length_mean.resize(2);
  fibre_length_variance.resize(2);
  for (unsigned int i=0; i<2; i++) {
    n_initial_fibres[i] =  ifile("fibres/n_initial_fibres",0,i);
    fibre_orientation_distribution[i] =
      ifile("fibres/fibre_orientation_distribution",0,i);
    fibre_orientation_mean_phi[i] =
      ifile("fibres/fibre_orientation_mean_phi",PI_val,i);
    fibre_orientation_variance_phi[i] =
      ifile("fibres/fibre_orientation_variance_phi",1.,i);
    fibre_orientation_mean_theta[i] =
      ifile("fibres/fibre_orientation_mean_theta",PI_val,i);    
    fibre_orientation_variance_theta[i] =
      ifile("fibres/fibre_orientation_variance_theta",1.,i);
    if (dimension==2) {
      fibre_orientation_mean_theta[i] = PI_val/2.;
      fibre_orientation_variance_theta[i] = 0.;
    }
    fibre_length_mean[i] =  ifile("fibres/fibre_length_mean",50.,i);
    fibre_length_variance[i] =  ifile("fibres/fibre_length_variance",2.,i);
  }
  
  fibre_radius = ifile("fibres/fibre_radius",0.2);
  fibre_make_gap = ifile("fibres/fibre_make_gap",0);
  vel_adhesion = ifile("fibres/vel_adhesion",0.);
  vel_contact = ifile("fibres/vel_contact",0.);
  fib_deg = ifile("fibres/fib_deg",1);
  prob_cell_deg = ifile("fibres/prob_cell_deg",1e-3);
  prob_diff_deg = ifile("fibres/prob_diff_deg",1e-6);
  cout <<  " -- " << endl;


  // [vessels] --- vessel parameters
  n_initial_vessels = ifile("vessels/n_initial_vessels",0);
  vessel_length.resize(n_initial_vessels);
  vessel_radius.resize(n_initial_vessels);
  vessel_startx.resize(n_initial_vessels);
  vessel_starty.resize(n_initial_vessels);
  vessel_startz.resize(n_initial_vessels);
  vessel_directionx.resize(n_initial_vessels);
  vessel_directiony.resize(n_initial_vessels);
  vessel_directionz.resize(n_initial_vessels);
  for (int i=0; i<n_initial_vessels; i++) {
    vessel_length[i] =  ifile("vessels/vessel_length",50.,i);
    vessel_radius[i] =  ifile("vessels/vessel_radius",5.,i);
    vessel_startx[i] =  ifile("vessels/vessel_startx",300.,i);
    vessel_starty[i] =  ifile("vessels/vessel_starty",300.,i);
    vessel_startz[i] =  ifile("vessels/vessel_startz",300.,i);
    vessel_directionx[i] =  ifile("vessels/vessel_directionx",1.,i);
    vessel_directiony[i] =  ifile("vessels/vessel_directiony",0.,i);
    vessel_directionz[i] =  ifile("vessels/vessel_directionz",0.,i);
  }
  vessel_PoissonNo = ifile("vessels/vessel_PoissonNo",0.5);
  vessel_YoungM = ifile("vessels/vessel_YoungM",1e-3);
  vessel_adhesion = ifile("vessels/vessel_adhesion",3.72e-4); 


  // postprocessing
  verbose = ifile("postprocessing/verbose",1);
  // default: write only cells
  writeVtkCells = ifile("postprocessing/writeVtkCells",1);
  writeCellList = ifile("postprocessing/writeCellList",0);
  writeVtkFibres = ifile("postprocessing/writeVtkFibres",0);
  writeVtkVessels = ifile("postprocessing/writeVtkVessels",0);
  getGenealogy = ifile("postprocessing/getGenealogy",0);
  outputDirectory = ifile("postprocessing/outputDirectory","./");
  testcase = ifile("postprocessing/testcase","test");
  fileCellsVisualization =
    ifile("postprocessing/fileCellsVisualization","celulas.txt");
  cellTracking = ifile("postprocessing/cellTracking",0);
  fileCellsTracking = ifile("postprocessing/fileCellsTracking","track.txt");
  fileCells = ifile("postprocessing/fileCells","all_cells.txt");
  writeStatistics = ifile("postprocessing/writeStatistics",0);
  casename = ifile("postprocessing/casename","case");
  casedirectory = ifile("postprocessing/casedirectory","./");
  
  cout << " --- ... parameters read. " << endl;

  // print parameter database
  print();

}

void Param::print()
{
  cout << endl;
  cout << " # ======================= " << endl;
  cout << " # INPUT FILE (GetPot): " << filename << endl;
  cout << " # ======================= " << endl;
  cout << endl;

  cout << "[coupling]" << endl;
  //cout << "fileConcCells = " << fileConcCells << endl;
  cout << "fileCells2FEM = " << fileCells2FEM << endl;
  cout << "fileCellsDensity2FEM = " << fileCellsDensity2FEM << endl;
  cout << "fileFEM2Cells = " << fileFEM2Cells << endl;
  cout << endl;

  cout << "[fem]" << endl;
  cout << "meshdir = " << meshdir << endl;
  cout << "meshname = " << meshname << endl;
  cout << "FreeFemFile = " << FreeFemFile << endl;
  cout << endl;

  cout << "[cells]" << endl;
  cout << "n_initial_cells = " << n_initial_cells << endl;
  cout << "n_phenotypes = " << n_phenotypes << endl;
  cout << "n_steps = " << n_steps << endl;
  cout << "time_step = " << time_step << endl;
  cout << "time_death = " << time_death << endl;
  cout << "radius = " << radius << endl;
  cout << "compressibility = " << compressibility << endl;
  cout << "contact_inhibition = " << contact_inhibition << endl;
  cout << "growth_rate = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout << growth_rate[i] << " ";
  }
  cout << "'" << endl;
  cout << "alpha_birthrate = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout << alpha_birthrate[i] << " ";
  }
  cout << "'" << endl;
  cout << "alpha_YoungM = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout << alpha_YoungM[i] << " ";
  }
  cout << "'" << endl;
  cout << "alpha_PoissonNo = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout <<alpha_PoissonNo[i] << " ";
  }
  cout << "'" << endl;
  cout << "Gcm = " << Gcm << endl;
  cout << "alpha_gcm = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout << alpha_gcm[i] << " ";
  }
  cout << "'" << endl;
  cout << "adhesion_value = '";
  for (int i=0; i<n_phenotypes; i++) {
    cout << adhesion_value[i] << " ";
  }
  cout << "'" << endl;
  cout << "variance_motion = " << variance_motion << endl;
  cout << "variance_phenotype = " << variance_phenotype << endl;
  cout << "variance_adhesion = " << variance_adhesion << endl;
  cout << "ic_phenotype = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_phenotype[i] << " ";
  }
  cout << "'" << endl;
  cout << "polarity_x = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << polarity_x[i] << " ";
  }
  cout << "'" << endl;
  cout << "polarity_y = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << polarity_y[i] << " ";
  }
  cout << "'" << endl;
  cout << "polarity_z = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << polarity_z[i] << " ";
  }
  cout << "'" << endl;
  cout << "ic_follower_leader = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_follower_leader[i] << " ";
  }
  cout << "'" << endl;
  cout << "hypoxic_birth = " << hypoxic_birth << endl;
  cout << "normoxic_birth = " << normoxic_birth << endl;
  cout << "death = " << death << endl;
  cout << "hypoxic_fiction = " << hypoxic_friction << endl;
  cout << "follower_force = " << follower_force << endl;
  cout << "follower_denominator = " << follower_denominator << endl;
  cout << "be_multiplier = " << be_multiplier << endl;
  cout << "be_displacement = " << be_displacement << endl;
  cout << endl;

  cout << "[mutations]" << endl; // TOMMASO
  cout << "initial_phenotype = " << initial_phenotype << endl;
  cout << "mutation_amount = " << mutation_amount << endl;
  cout << "mutation_probability = " << mutation_probability << endl;
  cout << "alpha_s = " << alpha_s << endl;
  cout << endl;
  
  cout << "[oxygen]" << endl;
  cout << "oxygen_response = " << oxygen_response << endl;
  cout << "threshold_death = " << threshold_death << endl;
  cout << "threshold_hypo = " << threshold_hypo << endl;
  cout << endl;

  cout << "[fibres]" << endl;
  cout << "n_initial_fibres = '";
  for (unsigned int i=0; i<n_initial_fibres.size(); i++) {
    cout << n_initial_fibres[i] << " ";
  }
  cout << "'" << endl;
  
  cout << "fibre_orientation_distribution = '" << fibre_orientation_distribution[0] << " " << fibre_orientation_distribution[1] << "'" << endl;
  cout << "fibre_orientation_mean_phi = '" << fibre_orientation_mean_phi[0] << " " << fibre_orientation_mean_phi[1] << "'" << endl;
  cout << "fibre_orientation_variance_phi = '" << fibre_orientation_variance_phi[0] << " " << fibre_orientation_variance_phi[1] << "'" << endl;
  cout << "fibre_orientation_mean_theta = '" << fibre_orientation_mean_theta[0] << " " << fibre_orientation_mean_theta[1] << "'" << endl;
  cout << "fibre_orientation_variance_theta = '" << fibre_orientation_variance_theta[0] << " " << fibre_orientation_variance_theta[1]  << "'" << endl;
  cout << "fibre_length_mean = '" << fibre_length_mean[0] << " " << fibre_length_mean[1] << "'" << endl;
  cout << "fibre_length_variance = '" << fibre_length_variance[0] << " " << fibre_length_variance[1] << "'" << endl;
  cout << "fibre_radius = " << fibre_radius << endl;
  cout << "fibre_make_gap = " << fibre_make_gap << endl;
  cout << "vel_adhesion = " << vel_adhesion << endl;
  cout << "vel_contact = " << vel_contact << endl;
  cout << "fib_deg = " << fib_deg << endl;
  cout << endl;

  cout << "[vessels]" << endl;
  cout << "n_initial_vessels = " << n_initial_vessels << endl;
  cout << "vessel_length = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_length[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_radius = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_radius[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_startx = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_startx[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_starty = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_starty[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_startz = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_startz[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directionx = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directionx[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directiony = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directiony[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_directionz = '";
  for (int i=0; i<n_initial_vessels; i++) {
    cout << vessel_directionz[i] << " ";
  }
  cout << "'" << endl;
  cout << "vessel_PoissonNo = " << vessel_PoissonNo << endl;
  cout << "vessel_YoungM = " << vessel_YoungM << endl;
  cout << "vessel_adhesion = " << vessel_adhesion << endl;
  cout << endl;

  cout << "[geo]" << endl;
  cout << "dimension = " << dimension << endl;
  cout << "boxesx = " << boxesx << endl;
  cout << "boxesy = " << boxesy << endl;
  cout << "boxesz = " << boxesz << endl;
  cout << "lattice_length_x = " << lattice_length_x << endl;
  cout << "lattice_length_y = " << lattice_length_y << endl;
  cout << "lattice_length_z = " << lattice_length_z << endl;
  cout << "box dimensions = " << lattice_length_x/boxesx << " by " << lattice_length_y/boxesy << " by " <<  lattice_length_z/boxesz << endl;
  if (ic_file_cells!="nil") {
    cout << "ic_file_cells = '" << ic_file_cells << "'" << endl;
  }
  cout << "ic_cell_x = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_x[i] << " ";
  }
  cout << "'" << endl;
  cout << "ic_cell_y = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_y[i] << " ";
  }
  cout << "'" << endl;
  cout << "ic_cell_z = '";
  for (int i=0; i<n_initial_cells; i++) {
    cout << ic_cell_z[i] << " ";
  }
  cout << "'" << endl;
  cout << "max_cell = " << max_cell << endl;
  cout << endl;
  
  cout << "[postprocessing]" << endl;
  cout << "verbose = " << verbose << endl;
  cout << "testcase = " << testcase << endl;
  cout << "casename = " << casename << endl;
  cout << "outputDirectory = " << outputDirectory << endl;
  cout << "casedirectory = " << casedirectory << endl;
  cout << "writeVtkCells = " << writeVtkCells << endl;
  cout << "writeVtkFibres = " << writeVtkFibres << endl;
  cout << "writeVtkVessels = " << writeVtkVessels << endl;
  cout << "getGenealogy = " << getGenealogy << endl;
  cout << "writeCellList = " << writeCellList << endl;
  cout << "fileCells = " << fileCells << endl;
  cout << "writeStatistics = " << writeStatistics << endl;

  
  cout << endl;
}


