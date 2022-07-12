#ifndef COUPLEDMODEL_H
#define COUPLEDMODEL_H

#include "PDE.h"
#include "Cell.h"
#include "Fibre.h"
#include "Vessel.h"
#include <vector>
#include <string>
#include <list>
/** ************************************************************************ 
*
* @class      CoupledModel
* @brief      manage the coupled simulation IBM-FEM
* 
* * 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       12.2014
*
****************************************************************************/

struct Box {
  /// @brief cells belonging to box
  std::vector<Cell> cells;
  /// @brief (pointer to) fibres belonging to box
  std::vector<Fibre*> p_fibres;
  // @ brief mesh elements contained in the box
  std::vector<int>  v_triangles;
};

/// @brief object to keep track of genealogy
struct Birthday {
  int cell;
  int mother_cell;
  int time;  
};
  
class CoupledModel
{
 public:
  CoupledModel(){};
  ~CoupledModel(){};

  // ****************
  // class member

  ///@brief id of simulation (when running multiple simulations)
  int simulation_id;
  
  /// @brief input parameters
  Param params;
  
  std::string input_file_name;

  ///@brief the list of fibres in the model
  std::vector<Fibre> fibres;

  // ******** vessel related
  ///@brief the list of vessels in the model
  std::vector<Vessel> vessels;
  /// @brief set initial vessels in the system 
  void set_ic_vessels();
  /// @brief determines vessels in contact with cells
  void compute_cell_vessels_contact(const int u,
		      const int v,
		      const int w);
  /// @brief distance between a cell and a vessel
  double DISTANCE(const Cell& cell, const Vessel& vessel);
  /// @brief point on vessel with shortest distance to cell
  void VESSELPOINT(const Cell& cell, const Vessel& vessel);
  float vessel_point[3];
  // ****************************************

  /// @brief finite element solver
  PDE oxy_diff; 


  /// @brief boxes in the domain
  Box ***boxes_A,***boxes_new_A;

  int maxx,maxy,maxz;
  int minx,miny,minz;
  int new_maxx,new_maxy,new_maxz;
  int new_minx,new_miny,new_minz;

  /// @brief physical size of boxes
  double box_sizex, box_sizey, box_sizez;

  /// @brief total number of cells
  unsigned int total_no_of_cells;
   /// @brief total number of dead/removed cells // TOMMASO 12/9/19
  unsigned int total_no_of_removed_cells;
  /// @brief max cell number allowed
  unsigned int max_cell;


  // cell counters (per type)
  std::vector<unsigned int> totNorm,totHypo,totDead;
  // counter per phenotype at each time step
  double phenotype1_count, phenotype2_count;


  // NEW 17/0/17
  double Area_forcex;
  double Area_forcey;
  double Area_forcez;
  
  // ------------------------
  ///@brief maximum number of cell allowed in a box
  ///@todo to be removed when allocating box element dynamically
  unsigned int max_cell_in_box;

  /// @brief current time step
  unsigned int reloj;

  /// @brief interaction distance
  double epsilon;
  /// @brief cell initial radius [micron]
  double Ro;
  double max_radius_cell;

  unsigned int movez;
  // further model parameter (put in param class)
  ///@brief variance of random motion 
  double vf;

  ///@brief displacement of cells after birth
  double birth_step;

  /// @brief energy bound for birth
  double birth_energy;
  
  // to check
  int initial_cells_in_box;
  int cells_counter;
  unsigned int daughter1;
  int dauther_birth;
 
  // ****************
  // methods
  /// @brief initialize models (using input file)
  void init(std::string f);
  
  /// @brief display no. of cell for each box
  void count_cells_per_box();
  /// @brief store the number of cells of different status
  void count_cells_per_type();

  // ******************
  // random generators
  /// @brief random between 0 and 1
  double aleatorio(); 
  /// @brief random between a and b
  double aleatorio(const double a, const double b); 
  /// @brief sample from Normal distribution N(m,s)
  double box_muller(const double m, const double s); 
  /// @brief sample from Normal distribution N(m,s) + bounded 
  double box_muller(const double m, const double s, 
		    const double M, const double k); 
  
  // ******************
  // cell functions
  /// @brief distance between 2 cells
  double DISTANCE(const Cell& c1, const Cell& c2);
  /// @brief distance between a cell and a fiber
  double DISTANCE(const Cell& cell, const Fibre& fibre);

  /// @brief change the status of cell according to O2 concentration
  void oxy_in_cell(Cell& cell);
  /// @brief change the status of cell according to phenotype
  void phenotype_of_cell(Cell& cell); //TOMMASO
  /// @brief mutate the cell within phenotypic range
  void cell_mutation(Cell& cell); //TOMMASO
  /// @brief search a cell in a given box (-1: cell not found)
  int search_in_box(const vector<int>& box_number,
		    const unsigned int cell_name);// In solve_system
  // @brief revert cell phenotype hypoxic->normoxic (stochastic)
  void reverse_phenotype(Cell& cell);

  // ******************
  // fibre functions
  /// @brief cross-product function
  void CrossProduct(vector<double> vector_A, vector<double> vector_B, double C_P[]);
    /// @brief dot-product function
  double DotProduct(vector<double> vector_A, vector<double> vector_B);
  double DotProduct(vector<double> vector_A, double vector_B[]);
  double DotProduct(double vector_A[], double vector_B[]);
  /// @brief function to determine if a pair of fibres cross link
  void CROSSLINKS(const Fibre& this_fibre, const Fibre& other_fibre);
  /// @brief determines fibres in contact with cells
  void compute_fibre_crosslinks(const int j);
  //void compute_fibre_crosslinks(const int u, const int v, const int w, const int j);

  double crosslink_count = 0;
  
  // ******************
  // box functions
  /// @brief allocate memories for position and boxes
  void allocate_compare_box(); 
  /// @brief cleans the elements of the box and update the new cells
  void update_box();
  /// @brief updates the maximum value in boxes
  void update_maximum();
  
  
  // ******************
  /// @brief set initial fibres in the system 
  void set_ic_fibres();

  ///@brief assign different ID to different parts of the domain
  int get_sub_domain_id(double x,double y,double z);
  /// @brief set initial cells in the system (from input file)
  void set_ic_cells();
  /// @brief set initial cells in the system (from input file)
  void set_ic_cells(string filename);
  
  /// @brief caculates the new cells in the system 
  void  cell_birth(Cell& cell);
  /// @brief removes the  cells which have died
  void  cell_death(Cell& cell);


  /// @brief compute new cell position
  void movement(const Cell& cell, 
		const int u, const int v, const int w, 
		const unsigned int cont_cell);
  void fibre_induced_movement(Cell& cell); // - CICELY NEW

  /// @brief calculate contact forces for cells in a given box
  void contact_forces(const int u,
		      const int v,
		      const int w,
		      const unsigned int cells_number);

   /// @brief determines fibres in contact with cells
  void compute_cell_fibres_contact(const int u,
		      const int v,
		      const int w);

  /// @brief determines cells in box [u,v,w] in contact with fibre j
  void compute_cell_fibres_contact(const int u,
				   const int v,
				   const int w,
				   const int j);

  /// @brief determines contact with the cells in box u,v,w
  void compute_cell_cell_contact(const int u, const int v,const int w);

  /// @brief compute force acting on the cells in the box u,v,w
  void compute_all_forces(const int u, const int v,const int w);

  void cell_cell_interaction(Cell& cell);
  void cell_vessel_interaction(Cell& cell);

  void update_cell_velocity(Cell& cell);
  
  /// @brief compute cell-cell and cell-vessel forces
  void hertz(Cell& cell);

  ///@brief cell velocity correction due to surrounding fibres
  void cell_fibres_interaction(Cell& cell);
  void cell_fibres_interaction(Cell& cell, Fibre& fibre); 

  ///@brief set the fibre_exists = 0 if the fibre is damaged
  void fibre_degradation(Cell& cell, Fibre& fibre, double cell_fibre_min_dist);
  
  // cell-fem coupling
  /// @brief set tetra per box
  /// @todo check this function
  void setElementsInBox(const Mesh& m);
  /// @brief check no. of element assigned to each box
  void checkElementsInBoxes();
  
  void compare_elements(int u, int v, int w,Mesh& _mesh);
  
  
  /// @brief time loop of coupled model
  void loop();

  /// @brief print all cells
  void printAllCells();

  /// @brief collect children of a given cell
  void collect_all_children(unsigned int mother_id,
					 std::vector<int>& children,
					 std::vector<int>& children_times);

  /// @brief loop through all cells and gather all births info
  void gather_all_births();
  /**
     @brief write cell list in a .vtk file
     @param[in] onlyCoord = 1 : write only the coordinates and radius
  */
  void writeVtk(std::string f,unsigned int onlyCoord=0);
  void writeBoxesVtk(string f);
  void writeFibresVtk(string f);
  void writeVesselsVtk(string f);
  void writeParameterList();

  /// @brief perform final operations
  void end();

};
  
#endif
  
 
