/** ************************************************************************ 
*
* @class      Param
* @brief      manage the input parameters
* 
* * 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       12.2014
*
****************************************************************************/
#include <string.h>
#include "GetPot.h"

#ifndef __PARAM_H
#define __PARAM_H

using namespace std;

class Param
{
 
 public:
  Param(){};
  ~Param(){};

  string filename;
  string fileConcCells;
  string fileCells2FEM;
  string fileCellsDensity2FEM;
  string fileFEM2Cells;

  // [cells] --- cell parameters
  int n_phenotypes;
  int n_initial_cells;
  unsigned int n_steps;
  /// @brief cell time step
  double time_step;
  /// @brief maximum time allowed in hypoxic state
  double time_death;
  /// @brief cell initial radius [micron]
  double radius;
  /// @brief growth_rate;
  vector<double> growth_rate;
  /// @brief birth_rate;
  //double birth_rate;
  /// @brief alternative Birth Rate vector for array of different birth rates
  vector<double> alpha_birthrate;
  double hypoxic_birth;
  double normoxic_birth;
  double death;
  double hypoxic_friction;
  /// @brief cell compressibility
  double compressibility;
  /// @brief contact inhibition
  double contact_inhibition;
  /// @brief Young modulus [kPa*muN/mu m^2]
  //double YoungM;
  /// @brief alternative YoungM vector for array of different Young's moduli
  vector<double> alpha_YoungM;
  /// @brief Poisson modulus
  //double PoissonNo;
  /// @brief alternative PoissonNo vector for array of different poisson numbers
  vector<double> alpha_PoissonNo;
  /// @brief tissue viscosity
  double Gcm;
  /// @brief factor to multiply GCM
  vector<double> alpha_gcm;
  /// @brief adhesion values
  vector<double> adhesion_value;
  double variance_motion;
  double variance_phenotype;
  double variance_adhesion;
  /// @brief birth energy function
  double be_displacement;
  double be_multiplier;
  /// @brief initial value of the phenotype
  vector<int> ic_phenotype;
  vector<double> polarity_x;
  vector<double> polarity_y;
  vector<double> polarity_z;
  // @brief whether the cell is initially a follower or a leader
  vector<int> ic_follower_leader;
  /// @brief follower force
  double follower_force;
  double follower_denominator;

  // [mutations] TOMMASO
  double initial_phenotype;
  double mutation_amount;
  double mutation_probability;
  double alpha_s;
    
  // [oxygen] --- pde (oxygen)
  double oxygen_response;
  double threshold_hypo;
  /// @brief threshold hypo->death
  float threshold_death;

  // [fibres] --- fibre parameters
  unsigned int n_sub_domain;
  vector<double> x_start, x_end;
  /// @brief number of fibres to be created within the domain
  vector<int> n_initial_fibres;
  /// @brief number of sub_domains
  //int sub_domains;
  /// @brief type of fibre distibution  0=uniform 1=normal
  vector<int> fibre_orientation_distribution;
  /// @brief mean azimuth angle
  vector<double> fibre_orientation_mean_phi;
  /// @brief variance of azimuth angle
  vector<double> fibre_orientation_variance_phi;
  /// @brief mean elevation angle
  vector<double> fibre_orientation_mean_theta;
  /// @brief variance of elevation angle
  vector<double> fibre_orientation_variance_theta;
  /// @brief mean length
  vector<double> fibre_length_mean;
  /// @brief variance of length
  vector<double> fibre_length_variance;
  /// @brief fibre radius
  double fibre_radius;
  /// @brief do we want a gap around cell  0=no 1=yes
  int fibre_make_gap;
  /// @brief do we want fibres in different parts of the domain to have different properties 0=n 1=yes
  //int fibre_split_domain;
  /// @brief coefficient multiplier for velocity update due to adhesion
  double vel_adhesion;
  /// @brief coefficient multiplier for velocity update due to contact
  double vel_contact;
  /// @brief fibre_degradation swtich 0=off 1=on
  int fib_deg;
  /// @brief probability of fibre degradation due to cell
  int prob_cell_deg;
  /// @brief probability of fibre degradation due to diffusion
  int prob_diff_deg;



  // [vessels] --- vessel parameters
  /// @brief number of vessels
  int n_initial_vessels;
  /// @brief vessel radius
  vector<double> vessel_radius;
  /// @brief vessel length
  vector<double> vessel_length;
  /// @brief vessel start point 
  vector<double> vessel_startx, vessel_starty, vessel_startz;
  /// @brief vessel direction
  vector<double> vessel_directionx, vessel_directiony, vessel_directionz;
  /// @brief vessel Poisson number
  double vessel_PoissonNo;
  /// @brief vessel Young's modulus
  double vessel_YoungM;
  /// @brief vessel adhesion coefficient
  double vessel_adhesion;

  // [fem] --- pde oxygen
  ///@brief type of diffusion solver (0: no solver, 1: FreeFem)
  int femSolverType;
  string meshdir;
  string meshname;
  string FreeFemFile;
  string FreeFemCall;

   // [geo] --- boxes and geometry
  unsigned int dimension;
  unsigned int boxesx, boxesy, boxesz;
  double lattice_length_x;
  double lattice_length_y;
  double lattice_length_z;
  vector<double> ic_cell_x;
  vector<double> ic_cell_y;
  vector<double> ic_cell_z;
  string ic_file_cells;
  /// @brief maximum number of cells
  double max_cell;

  // postprocessing
  int verbose;
  string testcase;
  string casename;
  string outputDirectory;
  string casedirectory;
  string fileCellsVisualization; // obsolete
  int cellTracking;
  string fileCellsTracking;
  int writeVtkCells;
  int writeVtkFibres; 
  int writeVtkVessels; 
  int getGenealogy;
  int writeCellList;
  string fileCells;
  int writeStatistics;
  
  void readFile(string f);
 
  void print();

};


#endif
