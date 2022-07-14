#include "CoupledModel.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <list>

/* ****************************************************************************/
using namespace std;

double PIG=3.1415926535897932384626433832795;

/* *****************************************************************************
   ALEATORIO: Generates random numbers between 0 and 1              
   ***************************************************************************** */
double CoupledModel::aleatorio()
{
  double random_value = rand();
  return (random_value/RAND_MAX);
}

/* *****************************************************************************
   ALEATORIO_moves: Generates random numbers between a and b        
   ***************************************************************************** */
double CoupledModel::aleatorio(const double a, const double b)
{
  double random_value = 0.;
 
  if (a<=b) {
    random_value = b - (b-a)*this->aleatorio();
  } else {
    random_value = 1. - 2.*this->aleatorio();
  }
  return(random_value);
}

/* *****************************************************************************
   Routine for generating gaussian random numbers N(m,s)
   ***************************************************************************** */
double CoupledModel::box_muller(const double m, const double s)	
{				        
  double w, y1;
  double x1, x2;
  static double y2;
  static int use_last = 0;
  if (use_last)	{	        /* use value from previous call */
    y1 = y2;
    use_last = 0;
  } else {
    do {
      x1 = 2.0 * this->aleatorio() - 1.0;
      x2 = 2.0 * this->aleatorio() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  
  return( m + y1 * s );
}

/* ****************************************************************************
   Routine for generating gaussian random numbers N(m,s) + bound
   ***************************************************************************** */
double CoupledModel::box_muller(const double m, const double s, 
				const double MAX, const double MIN)	
{	
  double val = this->box_muller(m, s);
			        
  if (val < MIN)
    return MIN;
  else if (val > MAX)
    return MAX;
  else  
    return val;
}

/* ***************************************************************************
   Routine for determining distance between cells
   *************************************************************************** */
double CoupledModel::DISTANCE(const Cell& c1, const Cell& c2)
{
  double dist = 
    pow(c1.position[0] - c2.position[0],2)+
    pow(c1.position[1] - c2.position[1],2)+
    pow(c1.position[2] - c2.position[2],2);

  return(sqrt(dist));
}

/* *****************************************************************************
   Routine for determining the cross-product of two vectors (needed below)
   *************************************************************************** */
/*
   @brief
   simple routine to determine the cross product of two vectors
   vector_A is crossed with vector_B resulting in vector C_P
*/
void CoupledModel::CrossProduct(vector<double> vector_A, vector<double> vector_B, double C_P[]) {
    C_P[0] = vector_A[1] * vector_B[2] - vector_A[2] * vector_B[1];
    C_P[1] = -(vector_A[0] * vector_B[2] - vector_A[2] * vector_B[0]);
    C_P[2] = vector_A[0] * vector_B[1] - vector_A[1] * vector_B[0];
    return;
}

/* *****************************************************************************
   Routines for determining the dot product of two vectors (needed below)
   *************************************************************************** */
/*
   @brief
   simple routines to determine the dot product of two vectors
   vector_A is dotted with vector_B resulting in dotproduct
   different variants could do with being standardised
*/
double CoupledModel::DotProduct(vector<double> vector_A, vector<double> vector_B) {
   double dotproduct = 0;
   for (int i = 0; i < 3; i++){
       dotproduct = dotproduct + vector_A[i] * vector_B[i];
   }
   return dotproduct;
}

double CoupledModel::DotProduct(double vector_A[], double vector_B[]) {
    double dotproduct = 0;
    for (int i = 0; i < 3; i++){
        dotproduct = dotproduct + vector_A[i] * vector_B[i];
    }
    return dotproduct;
}

double CoupledModel::DotProduct(vector<double> vector_A, double vector_B[]) {
    double dotproduct = 0;
    for (int i = 0; i < 3; i++){
        dotproduct = dotproduct + vector_A[i] * vector_B[i];
    }
    return dotproduct;
}

/* *****************************************************************************
Routine for determining if a pair of vectors intersect 3D
*************************************************************************** */
/*
   @brief
   DESCRIPTION TO GO HERE
*/
void CoupledModel::CROSSLINKS(const Fibre& this_fibre,const Fibre& other_fibre)
{
    //cout << "fibre " << this_fibre.fibre_name << " tested against fibre " << other_fibre.fibre_name << std::endl;

    // determine all relevant vectors:
    vector<double> p1_to_p2(3, 0.0);
    vector<double> q1_to_q2(3, 0.0);
    vector<double> p1_to_q1(3, 0.0);
    for (int i = 0; i < 3; i++) {
        p1_to_p2[i] = this_fibre.length * this_fibre.direction[i];
        q1_to_q2[i] = other_fibre.length * other_fibre.direction[i];
        p1_to_q1[i] = other_fibre.start[i] - this_fibre.start[i];
    }

    double FCP[3];
    CrossProduct(p1_to_p2, q1_to_q2, FCP);
    double test1 = DotProduct(FCP,FCP);
    double test2 = DotProduct(p1_to_q1,FCP);
    if (test1 != 0 && abs(test2) < 0.1)
    {
        // calculate the required coefficients (including dot products):
        double a = DotProduct(p1_to_p2,p1_to_q1)/DotProduct(p1_to_p2,p1_to_p2);
        double b = DotProduct(p1_to_p2,q1_to_q2)/DotProduct(p1_to_p2,p1_to_p2);
        // form relevant additional vectors
        vector<double> c(3, 0.0);
        vector<double> n(3, 0.0);
        for (int i = 0; i < 3; i++) {
            c[i] = b*p1_to_p2[i]-q1_to_q2[i];
            n[i] = p1_to_q1[i] - a*p1_to_p2[i];
        }
        // calculate the t values for each fibre
        double tq = DotProduct(c,n)/DotProduct(c,c);
        double tp = a+b*tq;

        if (0 <= tq && tq <= 1 && 0 <= tp && tp <= 1) {
            /*cout << "fibres " << this_fibre.fibre_name << " and " << other_fibre.fibre_name << " cross-link " << endl;
            cout << "     fibre p is {" << this_fibre.start[0] << "," << this_fibre.start[1] << "," << this_fibre.start[2] << "}"
                 << " + t {" << p1_to_p2[0] << "," << p1_to_p2[1] << "," << p1_to_p2[2] << "}" << endl;
            cout << "     fibre q is {" << other_fibre.start[0] << "," << other_fibre.start[1] << "," << other_fibre.start[2] << "}"
                 << " + t {" << q1_to_q2[0] << "," << q1_to_q2[1] << "," << q1_to_q2[2] << "}" << endl;
            cout << "     vector p to q is {" << p1_to_q1[0] << "," << p1_to_q1[1] << "," << p1_to_q1[2] << "}" << endl;

            cout << "     tp is " << tp << " and tq is " << tq << endl;

            cout << "           FCP is " << FCP[0] << " " << FCP[1] << " " << FCP[2] << endl;
            cout << "           test2 is " << test2 << endl;*/
            crosslink_count++;
        }
    }
}


/* ************************************************************************
   compute fibre crosslinks
   ************************************************************************ */
/*void CoupledModel::compute_fibre_crosslinks(const int u,
                                               const int v,const int w,
                                               const int j)*/
/*
   @brief
   probably need to only check fibres within a certain radius for large numbers
*/
void CoupledModel::compute_fibre_crosslinks(const int j)
{
    //unsigned int total_fibres_in_box=this->boxes_A[u][v][w].p_fibres.size();
    unsigned int total_fibres_in_domain = params.n_initial_fibres[0]+params.n_initial_fibres[1];

    for(unsigned int i=0; i<total_fibres_in_domain; i++){
        CROSSLINKS(this->fibres[j],this->fibres[i]);
    }
    this->fibres[j].fibre_crosslinks = crosslink_count;
    if (1 <= this->fibres[j].fibre_crosslinks){
        cout << "fibre " << this->fibres[j].fibre_name << " has " << this->fibres[j].fibre_crosslinks << " crosslinks" << endl;
    }
}


/* *****************************************************************************
   Routine for determining distance between cell and fibre
   *************************************************************************** */
/*
   @brief 
   The distance between cell and fibre is computed using the
   scalar product Y = <c_to_f,fibre_vector> where
   c_to_f = cell_center - fibre_start
   fibre_vector = fibre.length*fibre.direction
   -> if Y<0, then the cell is closer to fibre.start
   -> if Y>fibre.length^2, then the cell is closer to fibre.end
   -> otherwise, we take the projection of cell center onto the fibre
*/
double CoupledModel::DISTANCE(const Cell& cell, const Fibre& fibre)
{
  double c_to_f[3]; // cell-to-fibre vector
  double c_to_f_length_2; // length of c_to_f squared
  double c_to_f_dot_v_f; // scalar product c_to_f * (f_end-f_start)
  
  c_to_f_length_2 = 0.;
  c_to_f_dot_v_f = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_to_f[i] = cell.position[i]-fibre.start[i];
    c_to_f_length_2 = c_to_f_length_2 + c_to_f[i]*c_to_f[i];
    c_to_f_dot_v_f = c_to_f_dot_v_f + c_to_f[i]*fibre.length*fibre.direction[i];
  }
  
  if (c_to_f_dot_v_f < 0.) {
    // in this case the cell is closest to the fibre start point
    //cout << "the cell is closest to the start of the fibre" << endl;
    return sqrt(c_to_f_length_2);
    
    } else if (c_to_f_dot_v_f > fibre.length*fibre.length) {
    // in this case the cell is closest to the fibre end point
    //cout << "the cell is closest to the end of the fibre" << endl;
    double distance = 0.;
    for (unsigned int i=0; i<3; i++) {
      distance = distance +
	(cell.position[i]-(fibre.start[i]+fibre.length*fibre.direction[i]))*
	(cell.position[i]-(fibre.start[i]+fibre.length*fibre.direction[i]));
    }
    return sqrt(distance);
    
    } else {
    // in this case the cell is closest to a point along the fibre
    //cout << "the cell is closest a point along the fibre" << endl;
    double c_to_f_length_cos_alpha_2 = c_to_f_dot_v_f*c_to_f_dot_v_f/
      (fibre.length*fibre.length);
    return sqrt( c_to_f_length_2 - c_to_f_length_cos_alpha_2);
    }
}

/* ***************************************************************************
   Routine for determining distance between cell and vessel
   ***************************************************************************** */
/*
   @brief 
   The distance between cell and vessel is computed using the
   scalar product Y = <c_to_v,vessel_vector> where
   c_to_v = cell_center - vessel_start
   vessel_vector = vessel.length*vessel.direction
   -> if Y<0, then the cell is closer to vessel.start
   -> if Y>vessel.length^2, then the cell is closer to vessel.end
   -> otherwise, we take the projection of cell center onto the vessel
*/
double CoupledModel::DISTANCE(const Cell& cell, const Vessel& vessel)
{
  double c_to_v[3]; // cell-to-vessel vector
  double c_to_v_length_2; // length of c_to_v squared
  double c_to_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_to_v_length_2 = 0.;
  c_to_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_to_v[i] = cell.position[i]-vessel.ves_start[i];
    c_to_v_length_2 = c_to_v_length_2 + c_to_v[i]*c_to_v[i];
    c_to_v_dot_v_v = c_to_v_dot_v_v + c_to_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_to_v_dot_v_v < 0.) {
    return sqrt(c_to_v_length_2);    
    } else if (c_to_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    double distance = 0.;
    for (unsigned int i=0; i<3; i++) {
      distance = distance +
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]))*
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]));
    }
    return sqrt(distance);
  } else {
    double c_to_v_length_cos_alpha_2 = c_to_v_dot_v_v*c_to_v_dot_v_v/
      (vessel.ves_length*vessel.ves_length);
    return sqrt( c_to_v_length_2 - c_to_v_length_cos_alpha_2);
    }
} 

/* **************************************************************************
   Routine for determining point on vessel which provides minimum distance
   ***************************************************************************** */
void CoupledModel::VESSELPOINT(const Cell& cell, const Vessel& vessel)
{
  double c_v[3]; // cell-to-vessel vector
  double c_v_length_2; // length of c_to_v squared
  double c_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_v_length_2 = 0.;
  c_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_v[i] = cell.position[i]-vessel.ves_start[i];
    c_v_length_2 = c_v_length_2 + c_v[i]*c_v[i];
    c_v_dot_v_v = c_v_dot_v_v + c_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_v_dot_v_v < 0.) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i];
    }
    //cout << " case 1 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else if (c_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i];
    }
    //cout << " case 2 " << vessel_point[0] << " " << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else {
    double dist_to_vessel_point = c_v_dot_v_v/(vessel.ves_length);
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+dist_to_vessel_point*vessel.ves_direction[i];
    }
    //cout << " case 3 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } 
}

/* **************************************************************************
   INITIALIZE the class (read input file, set parameters)
   ***************************************************************************** */
void CoupledModel::init(string f)
{
  // read input parameters
  this->input_file_name = f;
  //read an input file with GetPot
  params.readFile(f);

  // initialize simulation
  this->reloj = 0;
  // set box size
  this->box_sizex=params.lattice_length_x/(params.boxesx+0.0);
  this->box_sizey=params.lattice_length_y/(params.boxesy+0.0);
  this->box_sizez=params.lattice_length_z/(params.boxesz+0.0);

  //cout << box_sizex << " " << box_sizey << " " << box_sizez << endl;

  // random seed
  srand(static_cast<unsigned>(time(NULL)));

  this->vf = params.variance_motion;

  //allocate memory
  this->allocate_compare_box();
  
  // initial values for range of occupied boxes
  this->maxx=0;
  this->maxy=0;
  this->maxz=0;
  this->minx=10000;
  this->miny=10000;
  this->minz=10000;
  
  this->new_maxx=0;
  this->new_maxy=0;
  this->new_maxz=0;  
  this->new_minx=10000;
  this->new_miny=10000;
  if (params.dimension == 3)  this->new_minz=10000;
  else this->new_minz=0;
  
  /// @brief: max cell per box: depends on cell volume and box size
  this->Ro = params.radius;
  this->epsilon = 2*this->Ro; // CICELY QU.: is this actually used
  this->max_radius_cell = this->Ro;
  double cell_estimated_volume = 4./3.*PIG*pow(this->max_radius_cell,3.);
  double box_volume =  this->box_sizex*this->box_sizey*this->box_sizez;
  if (params.dimension == 2) {
    cell_estimated_volume = 2.*PIG*pow(this->max_radius_cell,2.); // CICELYQU. - SHOULD THIS BE 1.*PIG*pow(this->max_radius_cell,2.)
    box_volume =  this->box_sizex*this->box_sizey;
  }

  /// @brief: cell compressibility: determined from ...
  double cell_compressibility = params.compressibility;
  this->max_cell_in_box = box_volume/cell_estimated_volume*cell_compressibility;
  if (params.dimension == 2) this->max_cell_in_box *= 2;
  cout << " === max. number of cell per box allowed: " <<  this->max_cell_in_box << endl;
  this->max_cell = params.max_cell;

  cout << " set initial conditions for fibres" << endl;
  this->set_ic_fibres();

  cout << " set initial conditions for vessels" << endl;
  this->set_ic_vessels();

  cout << " set initial conditions for cells" << endl;
  // set IC for cells
  if ( this->params.ic_file_cells == "nil") {
    this->set_ic_cells();
  } else {
    this->set_ic_cells(this->params.ic_file_cells);
  }
  // updates the maximum value in boxes
  this->update_maximum(); 
  this->update_box();
  

  // print initial configuration on screen
  if ( this->params.verbose > 0) {
    cout << " === number of boxes (" 
	 << params.dimension << "D): " 
	 << params.boxesx << " x "
	 << params.boxesy << " x "
	 << params.boxesz << endl;

    cout << " === initial configuration: " << endl;
    for(unsigned int k=0; k<params.boxesx;k++) {
      for(unsigned int l=0; l<params.boxesy;l++) {
	for(unsigned int n=0; n<params.boxesz;n++) {
	  if (this->boxes_A[k][l][n].cells.size()>0) {
	    cout << " -> box " << k << "," << l << "," << n
		 << ": " << this->boxes_A[k][l][n].cells.size()
		 << " cell(s) " << endl;
	    cout << "type: ";
	  for (unsigned int j=0; j<this->boxes_A[k][l][n].cells.size(); j++) {
	    cout << this->boxes_A[k][l][n].cells[j].type << " ";
	  }
	  cout << endl;
	  }
	}
      }
    }	
    cout << endl;
  }

  ///@todo check the birth step: too small? depends on dimension?
  this->birth_step=.9/pow(2,1./3.) * this->Ro; // cell displacement due to newbirth
  
  if (params.dimension==2){
    this->movez = 0;
  } else if (params.dimension==3) {
    this->movez = 1;
  } else {
    cout << " ** ERROR: dimension = " << params.dimension
	 << " not supported. " << endl;
    exit(1);
  }

  // to check
  initial_cells_in_box = 0;
  cells_counter = 0;
  daughter1 = 1;
}

/* ***************************************************************************
   allocate_compare_box : allocates memory for name, distancia and posiciones
   *************************************************************************** */
void CoupledModel::allocate_compare_box()
{

  // dynamic allocated memory for cells in each box
  this->boxes_A = new Box** [params.boxesx];
  this->boxes_new_A = new Box** [params.boxesx];
  
  for(unsigned int i=0; i<params.boxesx; i++) {
    this->boxes_A[i] = new Box* [params.boxesy]; 
    this->boxes_new_A[i] = new Box* [params.boxesy]; 
  } 
  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      this->boxes_A[i][j] = new Box [params.boxesz]; 
      this->boxes_new_A[i][j] = new Box [params.boxesz]; 
    } 
  }

  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      for(unsigned int k=0; k<params.boxesz; k++) {
	this->boxes_A[i][j][k].cells.resize(0);
	this->boxes_new_A[i][j][k].cells.resize(0);
	this->boxes_A[i][j][k].v_triangles.resize(0);
	this->boxes_new_A[i][j][k].v_triangles.resize(0);
	this->boxes_A[i][j][k].p_fibres.resize(0);
	this->boxes_new_A[i][j][k].p_fibres.resize(0);
      }
    } 
  }
}

/* **************************************************************************
   place initial cells in the system - [read from file]
   ***************************************************************************** */
/*
  @brief file format
  
  N of cells
  ...
  x_i[3],polarity_i[3]{0},interaction_phenotype {0},follow/lead{1},adhesion{3.72e-4}
  ...
*/
void CoupledModel::set_ic_cells(string filename)
{

  ifstream ic_all_cells;
  ic_all_cells.open(filename.c_str(),ios::in);
	
  // read number of cells
  ic_all_cells >> this->total_no_of_cells;

  this->total_no_of_removed_cells = 0; // TOMMAS0 PROJ. - used for dead cells.

  cout << " CoupledModel::set_ic_cells reading from " << filename
       << ": " << this->total_no_of_cells << " cells " << endl;
  
  for (unsigned int l = 0; l<this->total_no_of_cells; l++) {
    
    Cell cell;
    // read position
    for (unsigned int j=0; j<3; j++){
      ic_all_cells >> cell.position[j];
    }
    
    // read polarity
    for (unsigned int j=0; j<3; j++){
      ic_all_cells >> cell.polarity[j];
    }
    // initial velocity = 0
    for (unsigned int j=0; j<3; j++){
      cell.vel[j]=0.;
    }

    // default parameters
    cell.name=l;  // ID
    cell.contacts=0; // no contacts (computed in a separate routine)
    cell.mother_name = -1; // (no mother)


    cell.birthday = this->reloj; // birthdate
    cell.radius=this->Ro;
    
    // read additional parameters 
    cell.type=1;
    cell.cont_pheno = params.initial_phenotype; //ADDED 25/6/19 TOMMASO
    cell.phenotype=params.threshold_hypo; 
    cell.phenotype_counter = 0; 
    cell.polarised = 0; 
    cell.hypoxic_count = 0;
    ic_all_cells >> cell.interaction_phenotype;
    ic_all_cells >> cell.is_follower;
    ic_all_cells >> cell.adhesion;
    
    // (default) oxygen concentration
    cell.O2 = 100.; // !! TODO: this should be given by the diffusion solver !!
    cell.dxO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dyO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dzO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    
    
    if (cell.interaction_phenotype==0) this->phenotype1_count++;
    else this->phenotype2_count++;

    //Cell box  
    int u=(int)(floor(cell.position[0]/this->box_sizex));
    int v=(int)(floor(cell.position[1]/this->box_sizey));
    int w=(int)(floor(cell.position[2]/this->box_sizez));
    
    cell.box[0]=u;
    cell.box[1]=v;
    cell.box[2]=w;
    
    cell.new_box[0]=u;
    cell.new_box[1]=v;
    cell.new_box[2]=w;
    cout << " cell " << l << " placed in box: "
	 << u << " " << v << " " << w << endl;
    cout << " position: " << cell.position[0] << " " 
	 << cell.position[1] << " " << cell.position[2] <<
      " phenotype " << cell.interaction_phenotype << endl;


    this->boxes_A[u][v][w].cells.push_back(cell);
    this->boxes_new_A[u][v][w].cells.push_back(cell);

    // compute the extrema of occupied boxes
    if(this->maxx<u && u<(int) params.boxesx) {
      this->new_maxx=u;	
    }
    if(this->maxy<v && v<(int) params.boxesy) {
      this->new_maxy=v;	
    }
    
    if(this->maxz<w && w<(int) params.boxesz) {
      this->new_maxz=w;	
    }	
    
    if(this->minx>u && u>0) {
      this->new_minx=u;	
    }
    
    if(this->miny>v && v>0) {
      this->new_miny=v;	
    }
    
    if(this->minz>w && w>0) {
      this->new_minz=w;	
    }	
    
    /// @todo can we set newbox=[u,v,w]?
    /// @todo what is the differnce between maxx and new_maxx?
    if(this->maxx<u && u<(int) params.boxesx) {
      this->maxx=u;	
    }
    if(this->maxy<v && v<(int) params.boxesy) {
      this->maxy=v;
    }
    if(this->maxz<w && w<(int) params.boxesz) {
      this->maxz=w;
    }	
    
    if(this->minx>u && u>=0 ) {
      this->minx=u;	
    }
    if(this->miny>v && v>=0 ) {
      this->miny=v;
    }
    if(this->minz>w && w>=0 ) {
      this->minz=w;	
    }	
  }  
}

/* **************************************************************************
   place initial cells in the system
   ***************************************************************************** */
void CoupledModel::set_ic_cells()
{
  Cell cell;
  int newbox[3];
  this->total_no_of_cells = params.n_initial_cells;

  this->total_no_of_removed_cells = 0; // TOMMAS0 PROJ. - used for dead cells.

  this->phenotype1_count = 0.0;
  this->phenotype2_count = 0.0;

  //Initialise the boxes cells number:
  for(unsigned int l=0; l < this->total_no_of_cells ; l++) {

    //Position
    cell.position[0]=params.ic_cell_x[l];
    cell.position[1]=params.ic_cell_y[l];
    cell.position[2]=params.ic_cell_z[l];

    cell.position_old[0]=params.ic_cell_x[l];
    cell.position_old[1]=params.ic_cell_y[l];
    cell.position_old[2]=params.ic_cell_z[l];

    //velocity
    cell.polarity[0]=params.polarity_x[l];
    cell.polarity[1]=params.polarity_y[l];
    cell.polarity[2]=params.polarity_z[l];
    
    for (unsigned int j=0; j<3; j++){
      cell.vel[j]=0.;
      // POSSIBLE CHANGE ??
      //cell.vel[j]=params.Gcm*cell.polarity[j];
    }

    //Cell characteristics
    cell.name=l;  
    cell.contacts=0;
    cell.mother_name = -1;

    cell.birthday = this->reloj;
    cell.radius=this->Ro;
    // initialising the cell contact area
    //cell.contact_area_old=0.0;
    //cell.variation_area=0.0;
    
    cell.type=1;
    cell.cont_pheno = params.initial_phenotype; /* ADDED 25/6/19 TOMMASO */
    cell.phenotype=params.threshold_hypo; 
    cell.phenotype_counter = 0; 

    cell.polarised = 0;

    // oxygen concentration
    cell.O2 = 100.; // !! TODO: this should be given by the diffusion solver !!
    cell.dxO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dyO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dzO2 = 0.; // !! TODO: this should be given by the diffusion solver !!

    // -----
    cell.hypoxic_count = 0;

    // set the interaction phenotype
    cell.interaction_phenotype = params.ic_phenotype[l];

    // set whether the cell is a follower or leader
    cell.is_follower = params.ic_follower_leader[l];

    // CICELY NEW - Change the adhesion_value 
    cell.adhesion=params.adhesion_value[cell.interaction_phenotype];
    
    if (cell.interaction_phenotype==0) this->phenotype1_count++;
    else this->phenotype2_count++;
    
    //Cell box  
    int u=(int)(floor(cell.position[0]/this->box_sizex));
    int v=(int)(floor(cell.position[1]/this->box_sizey));
    int w=(int)(floor(cell.position[2]/this->box_sizez));
    
    cell.box[0]=u;
    cell.box[1]=v;
    cell.box[2]=w;

    cell.new_box[0]=u;
    cell.new_box[1]=v;
    cell.new_box[2]=w;
    cout << " cell " << l << " placed in box: "
	 << u << " " << v << " " << w << endl;
    cout << " position: " << cell.position[0] << " " 
	 << cell.position[1] << " " << cell.position[2] <<
      " phenotype " << cell.interaction_phenotype << endl;

    this->boxes_A[u][v][w].cells.push_back(cell);
    this->boxes_new_A[u][v][w].cells.push_back(cell);

    // box of the new cell
    newbox[0]=u; 
    newbox[1]=v;
    newbox[2]=w;
   
    // compute the extrema of occupied boxes
    if(this->maxx<u && u<(int) params.boxesx) {
      this->new_maxx=u;	
    }
    if(this->maxy<v && v<(int) params.boxesy) {
      this->new_maxy=v;	
    }  
    if(this->maxz<w && w<(int) params.boxesz) {
      this->new_maxz=w;	
    }	
    
    if(this->minx>u && u>0) {
      this->new_minx=u;	
    }
    if(this->miny>v && v>0) {
      this->new_miny=v;	
    }
    if(this->minz>w && w>0) {
      this->new_minz=w;	
    }	
    
    /// @todo can we set newbox=[u,v,w]?
    /// @todo what is the differnce between maxx and new_maxx?
    if(this->maxx<newbox[0] && newbox[0]<(int) params.boxesx) {
      this->maxx=newbox[0];	
    }
    if(this->maxy<newbox[1]&& newbox[1]<(int) params.boxesy) {
      this->maxy=newbox[1];
    }
    if(this->maxz<newbox[2]&& newbox[2]<(int) params.boxesz) {
      this->maxz=newbox[2];
    }	
    
    if(this->minx>newbox[0]&& newbox[0]>=0 ) {
      this->minx=newbox[0];	
    }
    if(this->miny>newbox[1]&& newbox[1]>=0 ) {
      this->miny=newbox[1]	;
    }
    if(this->minz>newbox[2]&& newbox[2]>=0 ) {
      this->minz=newbox[2];	
    }	
  }

  cout << " phenotype 1 - " << this->phenotype1_count << " cells " << endl;
  cout << " phenotype 2 - " << this->phenotype2_count << " cells " << endl;
  
  //cout << " -- occupied box with smallest coordinate (min box):";
  //cout << this->minx << " " << this->miny << " " << this->minz << endl;
  //cout << " -- occupied box with largest coordinate (max box):";
  //cout << this->maxx << " " << this->maxy << " " << this->maxz << endl; 
}

/* ***************************************************************************
   Place fibres in the system 
   *************************************************************************** */
void CoupledModel::set_ic_fibres() {
    Fibre fibre;
    double fname_update = 0.0;

    /// @brief we divide the domain in subdomains
    for (unsigned int sub_domain = 0; sub_domain < params.n_sub_domain; sub_domain++) {

        for (int ll = 0; ll < this->params.n_initial_fibres[sub_domain]; ll++) {
            fibre.fibre_name = ll + fname_update;
            fibre.fibre_exists = 1.;

            fibre.fradius = params.fibre_radius;

            // for each fibre randomly position one end point within the domain
            /// @brief x co-ordinate for each fibre
            fibre.start[0] = aleatorio(params.x_start[sub_domain], params.x_end[sub_domain]);
            /// @brief y co-ordinate for each fibre
            fibre.start[1] = aleatorio(0.0, params.lattice_length_y);
            /// @brief z co-ordinate for each fibre
            if (params.dimension == 3) {
                fibre.start[2] = aleatorio(0.0, params.lattice_length_z);
            } else {
                fibre.start[2] = 0.;
            }

            /// @brief for each fibre randomly assign angles for direction
            if (params.fibre_orientation_distribution[sub_domain] == 0) {
                fibre.phi = aleatorio(0.0, 2. * PIG);
                if (params.dimension == 3) fibre.theta = aleatorio(0., PIG);
                else fibre.theta = PIG / 2.;
            }

            /// @brief create a normally distributed set of angles
            if (params.fibre_orientation_distribution[sub_domain] == 1) {
                fibre.phi = box_muller(params.fibre_orientation_mean_phi[sub_domain],
                                       params.fibre_orientation_variance_phi[sub_domain]);
                fibre.theta = box_muller(params.fibre_orientation_mean_theta[sub_domain],
                                         params.fibre_orientation_variance_theta[sub_domain]);
            }

            /// @brief for each fibre assign length from a normal distribution
            fibre.length = box_muller(params.fibre_length_mean[sub_domain],
                                      params.fibre_length_variance[sub_domain]);

            /// @brief for each fibre randomly align throughout the domain
            /// @brief x direction for each fibre
            fibre.direction[0] = sin(fibre.theta) * cos(fibre.phi);
            /// @brief y direction for each fibre
            fibre.direction[1] = sin(fibre.theta) * sin(fibre.phi);
            /// @brief z direction for each fibre
            fibre.direction[2] = cos(fibre.theta);
            if (params.dimension == 2) {
                fibre.direction[2] = 0.;
            }

            fibre.end[0] = fibre.start[0] + fibre.direction[0] * fibre.length;
            fibre.end[1] = fibre.start[1] + fibre.direction[1] * fibre.length;
            fibre.end[2] = fibre.start[2] + fibre.direction[2] * fibre.length;

            //cout << "fibre start " << fibre.start[0] << " " << fibre.start[1] << " " << fibre.start[2] << endl;
            //cout << "fibre end " << fibre.end[0] << " " << fibre.end[1] << " " << fibre.end[2] << endl;
            //cout << "fibre direction " << fibre.direction[0] << " " << fibre.direction[1] << " " << fibre.direction[2] << endl;

            /// @brief Determine whether end point of fibre lies within Domain - if not update fibre
            double latice_length[3];
            latice_length[0] = params.lattice_length_x;
            latice_length[1] = params.lattice_length_y;
            latice_length[2] = params.lattice_length_z;

            while (fibre.start[0] + fibre.length * fibre.direction[0] > latice_length[0] ||
                   fibre.start[0] + fibre.length * fibre.direction[0] < 0) {

                fibre.start[0] = aleatorio(params.x_start[sub_domain], params.x_end[sub_domain]);

                if (params.fibre_orientation_distribution[sub_domain] == 0) {
                    fibre.phi = aleatorio(0.0, 2. * PIG);
                    if (params.dimension == 3) fibre.theta = aleatorio(0., PIG);
                    else fibre.theta = PIG / 2.;
                }
                if (params.fibre_orientation_distribution[sub_domain] == 1) {
                    fibre.phi = box_muller(params.fibre_orientation_mean_phi[sub_domain],
                                           params.fibre_orientation_variance_phi[sub_domain]);
                    fibre.theta = box_muller(params.fibre_orientation_mean_theta[sub_domain],
                                             params.fibre_orientation_variance_theta[sub_domain]);
                }

                fibre.direction[0] = sin(fibre.theta) * cos(fibre.phi);
                fibre.end[0] = fibre.start[0] + fibre.direction[0] * fibre.length;
            }

            while (fibre.start[1] + fibre.length * fibre.direction[1] > latice_length[1] ||
                   fibre.start[1] + fibre.length * fibre.direction[1] < 0) {

                fibre.start[1] = aleatorio(0.0, params.lattice_length_y);

                if (params.fibre_orientation_distribution[sub_domain] == 0) {
                    fibre.phi = aleatorio(0.0, 2. * PIG);
                    if (params.dimension == 3) fibre.theta = aleatorio(0., PIG);
                    else fibre.theta = PIG / 2.;
                }
                if (params.fibre_orientation_distribution[sub_domain] == 1) {
                    fibre.phi = box_muller(params.fibre_orientation_mean_phi[sub_domain],
                                           params.fibre_orientation_variance_phi[sub_domain]);
                    fibre.theta = box_muller(params.fibre_orientation_mean_theta[sub_domain],
                                             params.fibre_orientation_variance_theta[sub_domain]);
                }

                fibre.direction[1] = sin(fibre.theta) * sin(fibre.phi);
                fibre.end[1] = fibre.start[1] + fibre.direction[1] * fibre.length;
            }

            if (params.dimension == 3) {

                while (fibre.start[2] + fibre.length * fibre.direction[2] > latice_length[2] ||
                       fibre.start[2] + fibre.length * fibre.direction[2] < 0) {

                    fibre.start[2] = aleatorio(0.0, params.lattice_length_z);

                    if (params.fibre_orientation_distribution[sub_domain] == 0) {
                        fibre.phi = aleatorio(0.0, 2. * PIG);
                        if (params.dimension == 3) fibre.theta = aleatorio(0., PIG);
                        else fibre.theta = PIG / 2.;
                    }
                    if (params.fibre_orientation_distribution[sub_domain] == 1) {
                        fibre.phi = box_muller(params.fibre_orientation_mean_phi[sub_domain],
                                               params.fibre_orientation_variance_phi[sub_domain]);
                        fibre.theta = box_muller(params.fibre_orientation_mean_theta[sub_domain],
                                                 params.fibre_orientation_variance_theta[sub_domain]);
                    }

                    fibre.direction[2] = cos(fibre.phi);
                    fibre.end[2] = fibre.start[2] + fibre.direction[2] * fibre.length;
                }
            }
            fibres.push_back(fibre);
        }
        fname_update += this->params.n_initial_fibres[sub_domain];
    }
}

/* ***************************************************************************
   Place vessels in the system 
   *************************************************************************** */
void CoupledModel::set_ic_vessels()
{
  Vessel vessel;

  for (int v=0; v < this->params.n_initial_vessels ; v++){
    vessel.vessel_name = v;
    vessel.ves_radius = params.vessel_radius[v];

    /// @brief vessel start position
    vessel.ves_start[0] = params.vessel_startx[v];
    vessel.ves_start[1] = params.vessel_starty[v];
    vessel.ves_start[2] = params.vessel_startz[v];
    
    /// @brief vessel length
    vessel.ves_length = params.vessel_length[v];
   
    /// @brief vessel direction
    vessel.ves_direction[0] = params.vessel_directionx[v];
    vessel.ves_direction[1] = params.vessel_directiony[v];
    vessel.ves_direction[2] = params.vessel_directionz[v];
   
    vessels.push_back(vessel);
  }
}

/* ***************************************************************************
   update minima and maxima 
   *************************************************************************** */
void CoupledModel::update_maximum()
{
  bool stopRun = false;
  if (this->new_maxx < this->new_minx) {
    cout << " ** WARNING: new_maxx " << new_maxx << ", new_minx " << new_minx << endl;
    cout << " ** WARNING: maxx " << maxx << ", minx " << minx << endl;
    stopRun = true;
  }
  if (this->new_maxy < this->new_miny) {
    cout << " ** WARNING: new_maxy " << new_maxy << ", new_miny " << new_miny << endl;
    cout << " ** WARNING: maxy " << maxy << ", miny " << miny << endl;
    stopRun = true;
  }
  if (params.dimension == 3) {
    if (this->new_maxz < this->new_minz) {
      cout << " ** WARNING: new_maxz " << new_maxz << ", new_minz " << new_minz << endl;
      cout << " ** WARNING: maxz " << maxz << ", minz " << minz << endl;
      stopRun = true;
    }
  }
  if (stopRun) exit(1);
  this->maxx = this->new_maxx;
  this->maxy = this->new_maxy;
  this->maxz = this->new_maxz;
  
  this->minx = this->new_minx;
  this->miny = this->new_miny;
  this->minz = this->new_minz;
}

/* ***************************************************************************
   Cleans the elements of the box and update the new cells 
   *************************************************************************** */
void CoupledModel::update_box()
{
  this->total_no_of_cells = 0;
  
  for(int u=this->minx; u<=this->maxx; u++) {
    for(int v=this->miny; v<=this->maxy; v++) {
      for(int w=this->minz; w<=this->maxz; w++) {
        // clear cell contacts
	    for(unsigned int j=0; j<this->boxes_new_A[u][v][w].cells.size(); j++)   {
	        this->boxes_new_A[u][v][w].cells[j].clear_contacts();
	    }
	
	    // boxes.cells = boxes_new.cells
	    this->boxes_A[u][v][w].cells = this->boxes_new_A[u][v][w].cells;

	    this->total_no_of_cells += this->boxes_new_A[u][v][w].cells.size();
	    //for(unsigned int j=0; j<this->boxes_A[u][v][w].cells.size(); j++)   {
	        //  clear cell contacts
	        //  this->boxes_A[u][v][w].cells[j].clear_contacts();
	    //}

	    // clear new box arrays
	    this->boxes_new_A[u][v][w].cells.clear();
	    // new: free also the used memory
	    vector<Cell> swap(this->boxes_new_A[u][v][w].cells);
	    //this->boxes_new_A[u][v][w].cells.resize(0)
      }  
    }  
  }
}

/* ***************************************************************************
   For each box, find the elements with barycenter inside the box
   **************************************************************************** */
void CoupledModel::setElementsInBox(const Mesh& _mesh)
{
  if (_mesh.dim==2) {

    // triangular mesh
    for(unsigned int l=0; l<_mesh.nTria; l++) {
      // barycenter of triangle
      double xT = (_mesh.xp[_mesh.tria[3*l]-1] +  _mesh.xp[_mesh.tria[3*l+1]-1] +  
		   _mesh.xp[_mesh.tria[3*l+2]-1])/3.;
      double yT = (_mesh.yp[_mesh.tria[3*l]-1] +  _mesh.yp[_mesh.tria[3*l+1]-1] +  
		   _mesh.yp[_mesh.tria[3*l+2]-1])/3.;
      double zT = params.lattice_length_z/2.;

      /// @todo floor() not needed here
      int u= (int)(floor( xT/this->box_sizex ));
      int v= (int)(floor( yT/this->box_sizey ));
      int w= (int)(floor( zT/this->box_sizez ));

      this->boxes_A[u][v][w].v_triangles.push_back(l);  
    }   
  } else {

        // tetrahedral mesh
        for(unsigned int l=0; l<_mesh.nTetra; l++) {
            double xT = (_mesh.xp[_mesh.tetra[4*l]-1] +  _mesh.xp[_mesh.tetra[4*l+1]-1] +
            _mesh.xp[_mesh.tetra[4*l+2]-1] + _mesh.xp[_mesh.tetra[4*l+3]-1])/4.;
            double yT = (_mesh.yp[_mesh.tetra[4*l]-1] +  _mesh.yp[_mesh.tetra[4*l+1]-1] +  
		    _mesh.yp[_mesh.tetra[4*l+2]-1] + _mesh.yp[_mesh.tetra[4*l+3]-1])/4.;
            double zT = (_mesh.zp[_mesh.tetra[4*l]-1] +  _mesh.zp[_mesh.tetra[4*l+1]-1] +  
		    _mesh.zp[_mesh.tetra[4*l+2]-1] + _mesh.zp[_mesh.tetra[4*l+3]-1])/4.;

            /// @todo floor() not needed here
            int u= (int)(floor( xT/this->box_sizex ));
            int v= (int)(floor( yT/this->box_sizey ));
            int w= (int)(floor( zT/this->box_sizez ));
      
            this->boxes_A[u][v][w].v_triangles.push_back(l);
        }
  }
}

/* ***************************************************************************
   check how many elements have been assigned to a box (good for debugging)
   **************************************************************************** */
void CoupledModel::checkElementsInBoxes()
{
  for(unsigned int k=0; k<params.boxesx; k++) {
    for(unsigned int l=0; l<params.boxesy; l++) {
      for(unsigned int n=0; n<params.boxesz; n++) {
	    cout << "box [" << k << " , "  << l << " , "  << n << "].tetra = ";
	    for (unsigned int ij=0; ij< this->boxes_A[k][l][n].v_triangles.size(); ij++) {
	        cout << this->boxes_A[k][l][n].v_triangles[ij] << " ";
	    }
	    cout << " (tot : "<<  this->boxes_A[k][l][n].v_triangles.size() << ")"<< endl;
      }
    }
  }
}

/* ****************************************************************************
   change the cell status due to spontaneous mutations TOMMASO PROJECT
   *****************************************************************************  */
void CoupledModel::cell_mutation(Cell& cell)
{
  double mutation = aleatorio();
  double lambda = params.mutation_probability;
  double nu = params.mutation_amount;
  double pR = 0.5;
  //double pR = cell.O2/60.0; // need to put 60 in as max cell.O2 somehow
  //lambda/2.0;//+1.0/2.0*(1-cell.O2/60.0); // taken from Stace et al. paper
  double leftside = 0.0;
  double rightside = 1.0;

      // mutate cells left or right or not at all
    if(mutation<lambda) {
        double pR_mutation = aleatorio();
        if(pR_mutation<pR) cell.cont_pheno = cell.cont_pheno+nu;
        else cell.cont_pheno = cell.cont_pheno-nu;
    }
    else cell.cont_pheno = cell.cont_pheno; //cell.type = O2_type;

    if(cell.cont_pheno<=0.0) cell.cont_pheno = leftside;
    if(cell.cont_pheno>=1.0) cell.cont_pheno = rightside;

  // =================================================================
  // write cell phenotypes to file
  // =================================================================
  //ofstream cellphenotype;
  //string cell_pheno =  "cell_phenotypes.txt";
  //cellphenotype.open(cell_pheno.c_str(),ios::app);
  //cellphenotype << cell.cont_pheno << " " ;
  //cellphenotype.close();
  //********************************   
}

/* ***************************************************************************
   change the status of cell according to O2 concentration
   *************************************************************************** */
/// @todo add state change in the Cell variable
void CoupledModel::oxy_in_cell(Cell& cell)
{
  // if cell is alive
  if(cell.type != 3) {

    // normoxic cell if oxygen is above cell.phenotype
    // !!!! CICELY QU. - Alfonso possible change for master code? Otherwise we could get a reverse phenotype !!!!
    if(cell.O2 >= cell.phenotype && cell.type == 1) {
      cell.type = 1;
      cell.hypoxic_count = 0;
    }
    
    // normoxic -> hypoxic if oxygen is below cell.phenotype
    if(cell.O2 >= params.threshold_death && cell.O2 < cell.phenotype) { 
      cell.type = 2;
      cell.hypoxic_count++;
    }

    // death if oxygen is below params.threshold_death
    if(cell.O2 < params.threshold_death) {
      cell.type = 3;
      //cout << "!!! cell " << cell.name << " has died due to low oxygen !!!" << endl;
    }	

    // death if cell stays too long in hypoxic state
    if(cell.hypoxic_count > params.time_death) {
      cell.type = 3;
      //cout << "!!! cell " << cell.name << " has died due to being hypoxic for too long!!!" << endl;
    }
  }		
}

/* ***************************************************************************
   change the status of cell according to phenotype
   *************************************************************************** */
/*
@brief function for the TOMMAS0 Project rather than oxy_in_cell
       cell considered normoxic if phenotype >=0.5
       cell considered hypoxic if phenotype <0.5
*/
void CoupledModel::phenotype_of_cell(Cell& cell)
{
    // providing cell alive
    if(cell.type != 3) {
        // normoxic if phenotype >=0.5
        if(cell.cont_pheno >= 0.5) cell.type = 1;
        // hypoxic if phenotype <0.5
        if(cell.cont_pheno < 0.5) cell.type = 2;
    }
}

/* ******************************************************************************
   revert cell phenotype hypoxic->normoxic if cell has enough oxygen
   ****************************************************************************** */
void CoupledModel::reverse_phenotype(Cell& cell)
{
  /// @todo add a parameter for regulating reversion // CICELYQU. - does this need refining???
  double time_steps_per_day = 60.*24./(params.time_step+0.0);
  float valor= this->aleatorio();  
  if( valor < 1./(2.*time_steps_per_day)) {
    cell.type = 1;
    cout << " ** reverse_phenotype: reverting cell " << cell.name << " to normoxic ** " << endl;
  } 
}

/* ********************************************************************************
   Calculates the new cells in the system for different phenotypes  
   ****************************************************************************** */
void CoupledModel::cell_birth(Cell& cell)
{
  //this->b_energy = params.birth_energy[cell.interaction_phenotype];

  double E, nu; // !!!! ALFONSO !!!!! CHANGE TO BE UPDATED TO MASTER
  if(params.n_phenotypes==1){
    E = params.alpha_YoungM[0];
    nu = params.alpha_PoissonNo[0];
  }
  else{
    E = fmax(params.alpha_YoungM[0],params.alpha_YoungM[1]);
    nu = fmax(params.alpha_PoissonNo[0],params.alpha_PoissonNo[1]);
  }
  
  double eff_modulus = E/(2.0*(1-nu*nu));
  double eff_radius = this->Ro*this->Ro/(2*this->Ro);
  double be_d = params.be_displacement;
  double force_rep = 4./3. * eff_modulus * sqrt(eff_radius) * pow(be_d,be_d);
  this->birth_energy = params.be_multiplier * force_rep / (2 * PIG * be_d);
  
  // coordinate of the box of the cell //TOMMASO
  //int bx=(int)(floor(cell.position[0]/this->box_sizex));
  //int by=(int)(floor(cell.position[1]/this->box_sizey));
  //int bz=(int)(floor(cell.position[2]/this->box_sizez));

  // CICELYQU. are we happy with all of these values?
  double max_contacts = params.contact_inhibition;
  bool birthconds = false;
  if (cell.contacts<=max_contacts &&
        cell.type!=3 &&
        //cell.type==1 && // type = 1: normoxic cells
        //this->boxes_A[bx][by][bz].cells.size()<=max_cell_in_box && //competition for space //TOMMASO
        cell.radius>0.99*max_radius_cell &&
        cell.energy<=this->birth_energy) {
      //cout << " the birth conditions are met " << endl;
      birthconds = true;
  }

  if(birthconds){
      
        double birth_probability = aleatorio();
    
        // TOMMASO FUNCTION
        double x = 1.0-cell.cont_pheno;
        double f = params.hypoxic_birth*(1.0-(1.0-x)*(1.0-x));
        double g = params.normoxic_birth*(cell.O2/(cell.O2+params.alpha_s))*(1-x*x);
        double birth_threshold = f+g;
        double death_threshold = params.death;//_2.0*params.hypoxic_birth;//*this->boxes_A[bx][by][bz].cells.size();
        // ===================
        //cout << "DEATH THRESHOLD " << death_threshold << endl;

        //double birth_threshold = params.alpha_birthrate[cell.interaction_phenotype];

        // TOMMASO Proj.
        if (birth_probability < death_threshold) {
            cell.type=3;
            if (params.verbose>2) {
                cout << "!!!!!!!!! WARNING: at time " << this->reloj << " cell " << cell.name <<
                     " has died and will be removed " << endl;
            }
        } else if(birth_probability < birth_threshold+death_threshold) {

            if (params.verbose > 2) {
                cout << " ** new birth at time " << this->reloj
                     << " mother: " << cell.name
                     << ", new number of cells " << this->total_no_of_cells + cells_counter + 1
                     << " ** " << endl;
            }

            // 3D: 2.rnew^3 = rold^3; 2D: 2.rnew^2 = rold^2
            double new_radius = cell.radius / pow(2., 1. / params.dimension);
            //We reset the radius and intracellular concentrations of the mother cell
            //this->boxes_A[cell.box[0]][cell.box[1]][cell.box[2]].cells[position_in_box_array].radius=new_radius;
            cell.radius = new_radius;

            // compute the position of the daugther cell
            double newpositionx, newpositiony, newpositionz;

            if (cell.contacts != 0) {
                //We take the preferred position calculated by the neighbours
                //Noise is necessary in order not to repeat births.
                // todo...
                //newpositionx=newpositionx + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
                //newpositiony=newpositiony + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
                //newpositionz=newpositionz + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
            }

            //If no cells around then we choose a position at random
            double phi = box_muller(0, 6.28);
            double theta = box_muller(0, 3.14);
            if (params.dimension == 2) theta = acos(-1.0) / 2.;
            newpositionx = cell.position[0] + this->birth_step * sin(theta) * cos(phi);
            newpositiony = cell.position[1] + this->birth_step * sin(theta) * sin(phi);
            newpositionz = cell.position[2] + this->birth_step * cos(theta);

            //check errors:
            if (newpositionx == cell.position[0] &&
                newpositiony == cell.position[1] &&
                newpositionz == cell.position[2]) {
                cout << " *** error in CoupledModel::birthday_phenotype : "
                     << " same position for daughter cell, file "
                     << __FILE__ << " line " << __LINE__ << endl;
                exit(1);
            }

            /* use periodic BC
                if (newpositionx<0)
                newpositionx=params.lattice_length_x + newpositionx;
                if (newpositionx>params.lattice_length_x)
                newpositionx=newpositionx - params.lattice_length_x;
                if (newpositiony<0)
                newpositiony=params.lattice_length_y + newpositiony;
                if (newpositiony>params.lattice_length_y)
                newpositiony=newpositiony - params.lattice_length_y;
                if (newpositionz<0)
                newpositionz=params.lattice_length_z + newpositionz;
                if (newpositionz>params.lattice_length_z)
                newpositionz=newpositionz - params.lattice_length_z;
            */

            if ((newpositionx < 0) || (newpositionx > params.lattice_length_x) ||
                (newpositiony < 0) || (newpositiony > params.lattice_length_y) ||
                (newpositionz < 0) || (newpositionz > params.lattice_length_z)) {
                /*cout << " ** (Time " << reloj
                << " ) ** !!!!!!!!! WARNING: the new cell created from " << cell.name
                << " is moving out of the domain (not supported) "
                << " this cell will no longer be taken into account " << endl;*/
                //if (cell.interaction_phenotype==0) this->phenotype1_count--;
                //else this->phenotype2_count--;
            } else {

                // coordinate of the box of the daughter cell
                int c1 = (int) (floor(newpositionx / this->box_sizex));
                int c2 = (int) (floor(newpositiony / this->box_sizey));
                int c3 = (int) (floor(newpositionz / this->box_sizez));

                if (this->boxes_A[c1][c2][c3].cells.size() >= max_cell_in_box) {
                    cout << " *** WARNING (time step " << reloj
                         << "): !! too many cells !! *** " << endl;
                    cout << " *** in box " << c1 << "," << c2 << "," << c3 << endl;
                    cout << " *** I found " << this->boxes_A[c1][c2][c3].cells.size() << " cells " << endl;
                    cout << " *** (max_cell_in_box = " << max_cell_in_box << ")" << endl;
                    //this->end();
                    //exit(1);
                }

                // ============================
                // create a new cell
                Cell newcell;
                // set the properties of the new cell
                newcell.birthday = this->reloj;
                newcell.name = this->total_no_of_cells + cells_counter + this->total_no_of_removed_cells;
                // counts the number of newborn cells
                cells_counter++;
                newcell.mother_name = cell.name;

                newcell.position[0] = newpositionx;
                newcell.position[1] = newpositiony;
                newcell.position[2] = newpositionz;

                // old position: position of mother cell
                newcell.position_old[0] = cell.position[0];
                newcell.position_old[1] = cell.position[1];
                newcell.position_old[2] = cell.position[2];

                newcell.box[0] = c1;
                newcell.box[1] = c2;
                newcell.box[2] = c3;

                // cell velocity (polarity) CHECK THIS
                newcell.polarity[0] = cell.polarity[0];
                newcell.polarity[1] = cell.polarity[1];
                newcell.polarity[2] = cell.polarity[2];

                // oxygen concentration
                newcell.type = cell.type; ///@attention we take type of mother
                newcell.O2 = cell.O2;
                newcell.dxO2 = cell.dxO2;
                newcell.dyO2 = cell.dyO2;
                newcell.dzO2 = cell.dzO2;

                // adhesion constant
                newcell.adhesion = cell.adhesion;

                newcell.radius = new_radius;
                newcell.clear_contacts();
                newcell.hypoxic_count = 0;
                // initialize counter of hypoxic state
                // newcell.phenotype_counter = 0;

                // follower-leader behavior
                newcell.is_follower = cell.is_follower;

                // CHECK 25/6/19 CICELYQU. - should this not be linked to position?
                // determine oxygen phenotype of new cell
                double alea = aleatorio();
                float random_threshold = 0.5;
                if (alea > random_threshold) {
                    double valor = params.threshold_hypo +
                                   params.variance_phenotype * (1. - 2. * rand()) / (RAND_MAX + 0.0);
                    if (valor > params.threshold_death) {
                        newcell.phenotype = valor;
                    } else {
                        newcell.phenotype = cell.phenotype;
                    }
                } else {
                    newcell.phenotype = cell.phenotype;
                }
                newcell.cont_pheno = cell.cont_pheno;
                //cout << "new cell has phenotype " << newcell.cont_pheno << endl;

                // determine interaction phenotype of the new cell
                newcell.interaction_phenotype = cell.interaction_phenotype;

                // determine polarisation of new cell
                newcell.polarised = 0;

                if (newcell.interaction_phenotype == 0) this->phenotype1_count++;
                else this->phenotype2_count++;

                // add the cell to the corresponding box
                boxes_A[c1][c2][c3].cells.push_back(newcell);

                // updating the box domain
                if (this->maxx < c1) {
                    this->new_maxx = c1;

                    if (c1 == (int) params.boxesx) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_maxx=0;
                        //c1=0;
                    }
                }
                if (this->maxy < c2) {
                    this->new_maxy = c2;

                    if (c2 == (int) params.boxesy) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_maxy=0;
                        //c2=0;
                    }
                }
                if (this->maxz < c3) {
                    this->new_maxz = c3;

                    if (c3 == (int) params.boxesz) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_maxz=0;
                        //c3=0;
                    }
                }

                if (this->minx > c1) {
                    this->new_minx = c1;

                    if (c1 < 0) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_minx=params.boxesx;
                        //c1=params.boxesx-1;
                    }
                }
                if (this->miny > c2) {
                    this->new_miny = c2;

                    if (c2 < 0) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_miny=params.boxesy;
                        //c2=params.boxesy-1;
                    }
                }
                if (this->minz > c3) {
                    this->new_minz = c3;

                    if (c3 < 0) {
                        cout << " ** error: Cell " << newcell.name
                             << " is moving out of the domain (not supported) " << endl;
                        exit(1);
                        //this->new_minz=params.boxesz;
                        //c3=params.boxesz-1;
                    }
                }
            }
        }
  }
        
  //grow if still possible
  double addgrowth = params.growth_rate[cell.interaction_phenotype]*params.time_step;
  if(cell.type!=3){
      if(cell.radius<(max_radius_cell)){
          cell.radius += addgrowth;
      }
  }
  /*else {
     if(cell.radius>=0.01){
          cell.radius -= addgrowth;
     }
  }*/
}

/* ********************************************************************************
   Removes any dead cells - currently not in use
   ****************************************************************************** */
void CoupledModel::cell_death(Cell& cell)
{
    cout << "cell " << cell.name << " is of type " << cell.type << " and will now be removed " << endl;
    this->total_no_of_removed_cells += 1;

    // remove the dead cell from its box
    int bx=(int)(floor(cell.position[0]/this->box_sizex));
    int by=(int)(floor(cell.position[1]/this->box_sizey));
    int bz=(int)(floor(cell.position[2]/this->box_sizez));

    unsigned int index = 0;
    for(unsigned int i=0; i<this->boxes_A[bx][by][bz].cells.size(); i++) {
        if (boxes_A[bx][by][bz].cells[i].name == cell.name) index = i;
    }
    //const unsigned int cell_death_status = 3;
    /*boxes_A[bx][by][bz].cells.erase(remove_if(boxes_A[bx][by][bz].cells.begin(), boxes_A[bx][by][bz].cells.end(),
                                              [&cell_death_status](const Cell &cell) -> bool {
                                                  return cell.type == cell_death_status;
                                              }), boxes_A[bx][by][bz].cells.end());*/

    // move last item in list to the position of index and then reduce list size by 1
    boxes_A[bx][by][bz].cells[index] = boxes_A[bx][by][bz].cells[ boxes_A[bx][by][bz].cells.size()-1 ];
    boxes_A[bx][by][bz].cells.pop_back();
    //if (boxes_A[bx][by][bz].cells.size()==0){break;}

    //clear memory
    //vector<Cell> swap(this->boxes_A[bx][by][bz].cells);

}


/* ****************************************************************************
   Cell velocity correction due to fibres                                        
   **************************************************************************** */
void CoupledModel::cell_fibres_interaction(Cell& cell)
{
  vector<double> fibre_adhesion, fibre_repulsion;
  fibre_adhesion.resize(params.dimension);
  fibre_adhesion.clear();
  fibre_repulsion.resize(params.dimension);
  fibre_repulsion.clear();
  unsigned int n_fibres = cell.contact_fibres.size();

  /*
  ofstream fibint;
  string fibint_list =  this->params.casedirectory + this->params.casename + "_fibre_interaction_info.txt";
  fibint.open(fibint_list.c_str(),ios::app);
  fibint << "    cell " << cell.name << " is in interacting with fibre " << endl;
  */

  // additional model parameters
  double cell_velocity_max = 10.;// [micron/s]
  double p_exponent = 1.; 
  double q_exponent = 1.;
  
  // compute current cell velocity
  double cell_velocity = 0.;
  for (unsigned int j=0; j<params.dimension; j++) {
    cell_velocity += cell.polarity[j]*cell.polarity[j];
  }
  cell_velocity = sqrt(cell_velocity);

  // for all the fibres in contact with the cell
  for(unsigned int k=0; k<n_fibres; k++) {

    //fibint << "    * " << cell.contact_fibres[k]->fibre_name << endl;

    // velocity_dot_direction = (cell.vel,fibre.direction)
    double velocity_dot_direction = 0.;  
    for (unsigned int j=0; j<params.dimension; j++) {
      velocity_dot_direction += cell.contact_fibres[k]->direction[j]*cell.polarity[j];
    }

    for (unsigned int j=0; j<params.dimension; j++) {
      double xi = fabs(velocity_dot_direction)/(cell_velocity+1e-8); // = |(v,f)|/|v|, |f| = 1
      double xip = pow(xi,p_exponent);
      double xiq = pow((1-xi*xi),q_exponent);

      // fibre_adhesion = alpha*|(v,f)| f   (f = fibre direction)
      fibre_adhesion[j] += params.vel_adhesion * xip
	  * (1 - cell_velocity/cell_velocity_max) * (cell.contact_fibres[k]->direction[j]);						   

      // fibre_repulsion = -beta*(1-|(v,f)|^2/|v|^2)*v
      fibre_repulsion[j] += -params.vel_contact * xiq * (cell.vel[j]);      
    }  
  }

  //fibint.close();
  
  // compute new velocity
  // used in the OLD version of contact/forces/velocity update
  //for (unsigned int j=0; j<params.dimension; j++) {
  //cell.vel[j] += fibre_adhesion[j] + fibre_repulsion[j]; //+ cell.polarity[j]; 
  //}

  // add force
  // used in the NEW version of contact/forces/velocity update
  for (unsigned int j=0; j<params.dimension; j++) {
    cell.force[j] += params.Gcm*(fibre_adhesion[j] + fibre_repulsion[j]); //+ cell.polarity[j];
    // note: we multiply by Gcm in order to maintain the old behavior
  }  
}

/* ***************************************************************************
   Adhesion and Repulsion forces between cells                               
   *************************************************************************** */
void CoupledModel::cell_cell_interaction(Cell& cell)
{
  cell.energy = 0;

  // forces for cell migration (leader-follower)
  double migration_force_x=0.;
  double migration_force_y=0.;
  double migration_force_z=0.;

  double fx = 0, fy = 0, fz = 0;
  for(unsigned int k=0; k<cell.contacts; k++) {

    // K = (1-nu^2)/E (for considered cell and neighbor)
    double nu = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- nu*nu)/params.alpha_YoungM[cell.interaction_phenotype];
    nu = params.alpha_PoissonNo[cell.neighbors[k]->interaction_phenotype];
    double neighbour_K = (1.-nu*nu)/params.alpha_YoungM[cell.neighbors[k]->interaction_phenotype];
    // sum of Ks
    double eff_K = (cell_K + neighbour_K);

    // effective radius
    double eff_radius = cell.radius * cell.neighbors[k]->radius /
      (cell.radius + cell.neighbors[k]->radius);

    double adhesion_coeff = fmin(cell.adhesion,cell.neighbors[k]->adhesion);
    
    // distance between centers
    double p_x = cell.neighbors[k]->position[0] - cell.position[0];
    double p_y = cell.neighbors[k]->position[1] - cell.position[1];
    double p_z = cell.neighbors[k]->position[2] - cell.position[2];
    double dist_cells = sqrt(p_x*p_x + p_y*p_y + p_z*p_z);
    
    double e_x = p_x/dist_cells;
    double e_y = p_y/dist_cells;
    double e_z = p_z/dist_cells;
   
    // distance between surfaces
    double d_ij = cell.radius + cell.neighbors[k]->radius - dist_cells;

    /*
       Theoretically d_ij>0 when cells are in contact. However, in some cases,
       (e.g. follow-leader) we also track the neighbors at time n-1 that
       are no longer neighbors at time n. In this case, we have to set manually
       d_ij = 0 so that the following terms vanished
    */
    if (d_ij<1e-8) {
      d_ij = 0;
    }
    
    // compute force
    double f_ij =
      // repulsion
      4./3. * (1./eff_K) * sqrt(eff_radius) * pow(d_ij,1.5) // NEW
      // adhesion
      //-cell.adhesion * (cell.radius*d_ij - d_ij*d_ij/4.);
      -adhesion_coeff * (cell.radius*d_ij - d_ij*d_ij/4.);
      //-adhesion_coeff * surface_area_new; // IGNACIO WONDERED ABOUT CHANGING THE ADHESION FORCE
    double surface = 2 * PIG * d_ij;
    cell.energy += fabs(f_ij)/surface;
     
    //Projection of the force over the three dimensions
    fx = fx + f_ij * (-e_x);
    fy = fy + f_ij * (-e_y);
    fz = fz + f_ij * (-e_z);

    /* ********************************************
       migration: based on variation in contact area
       ********************************************
    */
    ///@todo need to handle the case when the cells where not in contact before?
    double p_x_old = cell.neighbors[k]->position_old[0] - cell.position_old[0];
    double p_y_old = cell.neighbors[k]->position_old[1] - cell.position_old[1];
    double p_z_old = cell.neighbors[k]->position_old[2] - cell.position_old[2];
    double dist_cells_old = sqrt(p_x_old*p_x_old + p_y_old*p_y_old + p_z_old*p_z_old);
    double d_ij_old = cell.radius + cell.neighbors[k]->radius - dist_cells_old;
    
    //double contact_area_new = 2 * PIG * cell.radius * (d_ij)/2.0; // CICELYQU. is this more accurate or just PIG * cell.radius * (d_ij)
    double contact_area_new = 2 * PIG * this->Ro * (d_ij)/2.0;
    double contact_area_old = 2 * PIG * this->Ro * (d_ij_old)/2.0;
    
    double variation_area = contact_area_new-contact_area_old;

    if ((variation_area<0)&&(dist_cells>7)) {
        // intensity = f* a^2/(1+a^2)
        double add_denominator = params.follower_denominator;
        double intensity_migration_force = params.follower_force*
	    variation_area*variation_area/(add_denominator+variation_area*variation_area);
      
        migration_force_x = migration_force_x + intensity_migration_force*e_x;
        migration_force_y = migration_force_y + intensity_migration_force*e_y;
        migration_force_z = migration_force_z + intensity_migration_force*e_z;
    }  
  }

  // update total force on cell
  cell.force[0] += fx;
  cell.force[1] += fy;
  cell.force[2] += fz;

  if(cell.is_follower==1){
    // only follower feels migration forces
    cell.force[0] += migration_force_x;
    cell.force[1] += migration_force_y;
    cell.force[2] += migration_force_z;
  }  
}

/* ***************************************************************************
   Solve equation of motion for velocity                                     
   *************************************************************************** */
void CoupledModel::update_cell_velocity(Cell& cell)
{
  /* Possible change
    double force_value = cell.force[0]+cell.force[1]+cell.force[2];
    if(force_value==0){
        for (unsigned int j=0; j<params.dimension; j++) {
        cell.force[j]=params.Gcm*cell.polarity[j];
        }
    }
  */

  // compute cell velocity
  //double normgrad, diff_multiplier;
  double diff_multiplier;
  //double friction = params.alpha_gcm[cell.interaction_phenotype]*(params.Gcm);
  // TOMMASO - incorporating different phenotypes having different friction
  double friction_modifier = 1.0;
  if (cell.cont_pheno < 0.5){
      friction_modifier = params.hypoxic_friction;
  }
  double friction = params.Gcm * friction_modifier;

  switch (cell.type)
  {
    
    case 3:
    // dead cells: move only according to mechanical forces
    // TOMMASO dead cells no longer exist so we don't need to worry about this
    cell.vel[0] = 0.0;//cell.force[0]/(params.Gcm);
    cell.vel[1] = 0.0;//cell.force[1]/(params.Gcm);
    cell.vel[2] = 0.0;//cell.force[2]/(params.Gcm);
    break;
    
    case 2:
    // hypoxic cells
    // - might respond to oxygen gradient (params.oxygen_response>0)
    // - higher diffusion coefficient (diff_multiplier)
    diff_multiplier = 1.0; // increase variance for hypoxic cells
    cell.vel[0] = (cell.force[0] + box_muller(0,vf*diff_multiplier))/friction;
    cell.vel[1] = (cell.force[1] + box_muller(0,vf*diff_multiplier))/friction;
    cell.vel[2] = (cell.force[2] + box_muller(0,vf*diff_multiplier))/friction;
    //normgrad = sqrt(cell.dxO2*cell.dxO2 + cell.dyO2*cell.dyO2 + cell.dzO2*cell.dzO2)+1;
    // renormalized chemotaxis term: psi* gradc/(1+psi*|gradc|)  
    //cell.vel[0] = (cell.force[0]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dxO2/(1 + normgrad*params.oxygen_response)/friction;
    //cell.vel[1] = (cell.force[1]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dyO2/(1 + normgrad*params.oxygen_response)/friction;
    //cell.vel[2] = (cell.force[2]+ box_muller(0,vf*diff_multiplier))/friction +
    //params.oxygen_response * cell.dzO2/(1 + normgrad*params.oxygen_response)/friction;
    break;
    
    case 1:
    // normoxic cells
    cell.vel[0] = (cell.force[0] + box_muller(0,vf))/friction;
    cell.vel[1] = (cell.force[1] + box_muller(0,vf))/friction;
    cell.vel[2] = (cell.force[2] + box_muller(0,vf))/friction;
    break;
    
    default:
    cout << " *** ERROR in file " << __FILE__ << ", line " << __LINE__
    << " cell " << cell.name << " has type: " << cell.type  << ", unknown. " << endl;
    exit(1);
    break;  
  }

  // handle the 2-dimensional case
  if (params.dimension == 2)  cell.vel[2] = 0.;				
}

/* ***************************************************************************
   Adhesion and Repulsion between cells and vessels                          
   *************************************************************************** */
void CoupledModel::cell_vessel_interaction(Cell& cell)
{

  // attractive-repulsive forces
  double  c_v_fx=0.;
  double  c_v_fy=0.;
  double  c_v_fz=0.;

  for(unsigned int kk=0; kk<cell.contact_vessels.size(); kk++) {

    double cell_pois = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- cell_pois*cell_pois)/params.alpha_YoungM[cell.interaction_phenotype];
    double ves_pois = params.vessel_PoissonNo;
    double ves_K =  (1.- ves_pois*ves_pois)/params.vessel_YoungM;
    double eff_mod = 1.0/(cell_K + ves_K);

    double c_v_disp = 0.0;
    for (unsigned int j=0; j < vessels.size(); j++) {

      if (vessels[j].vessel_name == cell.contact_vessels[kk]->vessel_name) {
        double cell_vessel_dist = DISTANCE(cell,this->vessels[j]);
	    // displacement
	    c_v_disp = cell.radius + cell.contact_vessels[kk]->ves_radius - cell_vessel_dist;

	    this->VESSELPOINT(cell,this->vessels[j]);	
      }
    }

    double c_v_adh_coeff = fmin(cell.adhesion,params.vessel_adhesion);
    double c_v_f_ij =
      //c_v_repulsion
      4./3. * eff_mod * sqrt(cell.radius) * pow(c_v_disp,1.5)
      //c_v_adhesion
      -c_v_adh_coeff * (cell.radius*c_v_disp - c_v_disp*c_v_disp/4.);
    // NOT SURE ABOUT THIS BIT 
    double surf1 = 2 * PIG * c_v_disp;
    cell.energy += fabs(c_v_f_ij)/surf1;

    // co-ordinate distances between cell centre and vessel
    double c_v_x = this->vessel_point[0] - cell.position[0];
    double c_v_y = this->vessel_point[1] - cell.position[1];
    double c_v_z = this->vessel_point[2] - cell.position[2];
    double dist_cell_vessel = sqrt(c_v_x*c_v_x + c_v_y*c_v_y + c_v_z*c_v_z);
    
    c_v_x = c_v_x/dist_cell_vessel;
    c_v_y = c_v_y/dist_cell_vessel;
    c_v_z = c_v_z/dist_cell_vessel;

    c_v_fx = c_v_fx + c_v_f_ij * (-c_v_x);
    c_v_fy = c_v_fy + c_v_f_ij * (-c_v_y);
    c_v_fz = c_v_fz + c_v_f_ij * (-c_v_z);
  }

  cell.force[0] += c_v_fx;
  cell.force[1] += c_v_fy;
  cell.force[2] += c_v_fz;
}

/* ***************************************************************************
   hertz: gives the force from the hertz model  OLD VERSION COMMENTED OUT                          *
   *************************************************************************** */
/*void CoupledModel::hertz(Cell& cell)
{
  if (params.verbose > 3) {
    if (cell.neighbors.size()) {
      cout << " compute contact forces for cell " << cell.name << endl;
      cout << " contact with: ";
      for (unsigned int j=0; j < cell.neighbors.size(); j++) {
	    cout << cell.neighbors[j]->name << " ";
      }
      cout << endl;
    }
  }

  cell.energy = 0;

  //cell-vessel forces
  double  c_v_fx=0.;
  double  c_v_fy=0.;
  double  c_v_fz=0.;

  for(unsigned int kk=0; kk<cell.contact_vessels.size(); kk++) {

    double cell_pois = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- cell_pois*cell_pois)/params.alpha_YoungM[cell.interaction_phenotype];
    double ves_pois = params.vessel_PoissonNo;
    double ves_K =  (1.- ves_pois*ves_pois)/params.vessel_YoungM;
    double eff_mod = 1.0/(cell_K = ves_K);

    double c_v_disp = 0.0;
    for (unsigned int j=0; j < vessels.size(); j++) {

      if (vessels[j].vessel_name == cell.contact_vessels[kk]->vessel_name) {
        double cell_vessel_dist = DISTANCE(cell,this->vessels[j]);
	    c_v_disp = cell.radius + cell.contact_vessels[kk]->ves_radius - cell_vessel_dist;

	    this->VESSELPOINT(cell,this->vessels[j]);	
      }
    }

    double c_v_adh_coeff = fmin(cell.adhesion,params.vessel_adhesion);
    double c_v_f_ij =
      //c_v_repulsion
      4./3. * eff_mod * sqrt(cell.radius) * pow(c_v_disp,1.5)
      //c_v_adhesion
      -c_v_adh_coeff * (cell.radius*c_v_disp - c_v_disp*c_v_disp/4.);
    // NOT SURE ABOUT THIS BIT 
    double surf1 = 2 * PIG * c_v_disp;
    cell.energy += fabs(c_v_f_ij)/surf1;

    // co-ordinate distances between cell centre and vessel
    double c_v_x = this->vessel_point[0] - cell.position[0];
    double c_v_y = this->vessel_point[1] - cell.position[1];
    double c_v_z = this->vessel_point[2] - cell.position[2];
    double dist_cell_vessel = sqrt(c_v_x*c_v_x + c_v_y*c_v_y + c_v_z*c_v_z);
    
    c_v_x = c_v_x/dist_cell_vessel;
    c_v_y = c_v_y/dist_cell_vessel;
    c_v_z = c_v_z/dist_cell_vessel;

    c_v_fx = c_v_fx + c_v_f_ij * (-c_v_x);
    c_v_fy = c_v_fy + c_v_f_ij * (-c_v_y);
    c_v_fz = c_v_fz + c_v_f_ij * (-c_v_z);
  }

  // cell-cell forces
  double  fx=c_v_fx;
  double  fy=c_v_fy;
  double  fz=c_v_fz;

  // forces for cell migration (leader-follower)
  double migration_force_x=0.;
  double migration_force_y=0.;
  double migration_force_z=0.;
  
  for(unsigned int k=0; k<cell.contacts; k++) {

    // original version (Caiazzo & Ramis, JTB 2015)
    //double K_tap= 3./4.* 2 * (1.-params.PoissonNo*params.PoissonNo) / params.YoungM;
    // = 3./4. (1.-params.PoissonNo*params.PoissonNo) / (params.YoungM)
    //+ (1.-params.PoissonNo*params.PoissonNo) / (params.YoungM)
 
    // K = (1-nu^2)/E (for considered cell and neighbor)
    double nu = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- nu*nu)/params.alpha_YoungM[cell.interaction_phenotype];
    nu = params.alpha_PoissonNo[cell.neighbors[k]->interaction_phenotype];
    double neighbour_K = (1.-nu*nu)/params.alpha_YoungM[cell.neighbors[k]->interaction_phenotype];
    // sum of Ks
    double eff_K = (cell_K + neighbour_K);

    // effective radius
    double eff_radius = cell.radius * cell.neighbors[k]->radius /
      (cell.radius + cell.neighbors[k]->radius);

    double adhesion_coeff = fmin(cell.adhesion,cell.neighbors[k]->adhesion);
    
    // distance between centers
    double p_x = cell.neighbors[k]->position[0] - cell.position[0];
    double p_y = cell.neighbors[k]->position[1] - cell.position[1];
    double p_z = cell.neighbors[k]->position[2] - cell.position[2];
    double dist_cells = sqrt(p_x*p_x + p_y*p_y + p_z*p_z);
    
    double e_x = p_x/dist_cells;
    double e_y = p_y/dist_cells;
    double e_z = p_z/dist_cells;
   
    // distance between surfaces
    double d_ij = cell.radius + cell.neighbors[k]->radius - dist_cells;

    // Theoretically d_ij>0 when cells are in contact.
    // However, in some cases (e.g. follow-leader), we also track the neighbors
    // at time n-1 that are no longer neighbors at time n.
    // In this case, we have to set manually d_ij = 0
    if (d_ij<1e-8) {
      d_ij = 0;
    }
    
    // compute force
    double f_ij =
      // repulsion
      4./3. * (1./eff_K) * sqrt(eff_radius) * pow(d_ij,1.5) // NEW
      // adhesion
      //-cell.adhesion * (cell.radius*d_ij - d_ij*d_ij/4.);
      -adhesion_coeff * (cell.radius*d_ij - d_ij*d_ij/4.);
      //-adhesion_coeff * surface_area_new; // IGNACIO WONDERED ABOUT CHANGING THE ADHESION FORCE
    double surface = 2 * PIG * d_ij;
    cell.energy += fabs(f_ij)/surface;
    
    // // =================================================================
    // // write cell energies to data file
    // // =================================================================
    // ofstream cellenergy;
    // string ce_list =  "cell_energies.txt";
    // cellenergy.open(ce_list.c_str(),ios::app);
    // cellenergy << cell.energy << " " ;
    // cellenergy.close();
    // // ********************************

    //Projection of the force over the three dimensions
    fx = fx + f_ij * (-e_x);
    fy = fy + f_ij * (-e_y);
    fz = fz + f_ij * (-e_z);
    
    // migration: based on variation in contact area
    ///@todo need to handle the case when the cells where not in contact before?
    double p_x_old = cell.neighbors[k]->position_old[0] - cell.position_old[0];
    double p_y_old = cell.neighbors[k]->position_old[1] - cell.position_old[1];
    double p_z_old = cell.neighbors[k]->position_old[2] - cell.position_old[2];
    double dist_cells_old = sqrt(p_x_old*p_x_old + p_y_old*p_y_old + p_z_old*p_z_old);
    double d_ij_old = cell.radius + cell.neighbors[k]->radius - dist_cells_old;
    
    //double contact_area_new = 2 * PIG * cell.radius * (d_ij)/2.0;
    double contact_area_new = 2 * PIG * this->Ro * (d_ij)/2.0;
    double contact_area_old = 2 * PIG * this->Ro * (d_ij_old)/2.0;
    
    double variation_area = contact_area_new-contact_area_old;

    // The additional force is only applied only if:
    //  the cells are greater than 7 units apart (otherwise cells envelop other cells)
    //  variation_area < 0
    
    if ((variation_area<0)&&(dist_cells>7)) {
      // intensity = f* a^2/(1+a^2)
      double add_denominator = params.follower_denominator;
      double intensity_migration_force = params.follower_force*
	  variation_area*variation_area/(add_denominator+variation_area*variation_area);
      
      migration_force_x = migration_force_x + intensity_migration_force*e_x;
      migration_force_y = migration_force_y + intensity_migration_force*e_y;
      migration_force_z = migration_force_z + intensity_migration_force*e_z;
    } 
  }


  // normalize migration forces
  //double normalise_force = sqrt(migration_force_x*migration_force_x +
				//migration_force_y*migration_force_y +
				//migration_force_z*migration_force_z);
  //if (normalise_force!=0){
    //migration_force_x = migration_force_x/normalise_force;
    //migration_force_y = migration_force_y/normalise_force;
    //migration_force_z = migration_force_z/normalise_force;
  //}

  // double d_x = fabs(cell.position[0] - this->params.lattice_length_x/2.);
  // if (d_x<1.5*cell.radius) {
  //   double nu = params.alpha_PoissonNo[cell.interaction_phenotype];
  //   double cell_K =  (1.- nu*nu)/params.alpha_YoungM[cell.interaction_phenotype];
  //   fx += 4./3. * (1./cell_K) * cell.radius * pow(d_x,1.5) *
  //     (cell.position[0]-this->params.lattice_length_x/2.)/d_x;
  // }
  // if (d_x<cell.radius/2.) fx = 0.;
  
  // compute cell velocity
  //double normgrad, diff_multiplier;
  double diff_multiplier;
  double friction = params.alpha_gcm[cell.interaction_phenotype]*(params.Gcm);
  switch (cell.type) {
    
    case 3:
    // dead cells
    cell.vel[0] = fx/(params.Gcm);
    cell.vel[1] = fy/(params.Gcm);
    cell.vel[2] = fz/(params.Gcm);
    break;
    
    case 2:
    // hypoxic cells
    // - might respond to oxygem gradient (params.oxygen_response>0)
    // - higher diffusion coefficient (diff_multiplier)
    //normgrad = sqrt(cell.dxO2*cell.dxO2 + cell.dyO2*cell.dyO2 + cell.dzO2*cell.dzO2)+1;
    //diff_multiplier = 10.0; // increase variance for hypoxic cells
    diff_multiplier = 1.0;
    // renormalized chemotaxis term: psi* gradc/(1+psi*|gradc|)  
    cell.vel[0] = (fx+ box_muller(0,vf*diff_multiplier))/friction; //+
        //params.oxygen_response * cell.dxO2/(1 + normgrad*params.oxygen_response)/params.Gcm;
    cell.vel[1] = (fy+ box_muller(0,vf*diff_multiplier))/friction; //+
    //params.oxygen_response * cell.dyO2/(1 + normgrad*params.oxygen_response)/params.Gcm;
    cell.vel[2] = (fz+ box_muller(0,vf*diff_multiplier))/friction;
    //+ params.oxygen_response * cell.dzO2/(1 + normgrad*params.oxygen_response)/params.Gcm;	      
    break;
    
    case 1: 
    // normoxic      
    //double friction = params.alpha_gcm[cell.interaction_phenotype]*
    //(params.Gcm + 0.05*(1.+tanh(cell.position[0]-300)));
    //diff_multiplier = cell.type;
    diff_multiplier = 1.0;
    
    // FOLLOWER CELLS
    if(cell.is_follower==1){
        cell.vel[0] = (fx + box_muller(0,vf*diff_multiplier) + migration_force_x)/friction;
        cell.vel[1] = (fy + box_muller(0,vf*diff_multiplier) + migration_force_y)/friction;
        cell.vel[2] = (fz + box_muller(0,vf*diff_multiplier) + migration_force_z)/friction;
    }
    // LEADER CELLS
    else {
        cell.vel[0] = (fx + box_muller(0,vf*diff_multiplier))/friction;
        cell.vel[1] = (fy + box_muller(0,vf*diff_multiplier))/friction;
        cell.vel[2] = (fz + box_muller(0,vf*diff_multiplier))/friction;
    }
    break;
    
    default:
    cout << " *** ERROR in file " << __FILE__ << ", line " << __LINE__
	    << " cell " << cell.name << " has type: " << cell.type  << ", unknown. " << endl;
    exit(1);
    break;  
  }

  // handle the 2-dimensional case
  if (params.dimension == 2)  cell.vel[2] = 0.;  
}*/

/* ************************************************************************
   moves the cells according to the velocity previously computed
   ************************************************************************ */
void CoupledModel::movement(const Cell& cell, 
	      const int u, const int v, const int w, 
	      const unsigned int cont_cell)
{
  // copy the cell
  Cell celula_nueva = cell;
  // save the old position
  for (unsigned int j=0; j<3; j++) {
    celula_nueva.position_old[j] = cell.position[j];
  }

  // compute new position

  //LEADER CELL
  //the force-based movement of cells now includes an additional term for a basic linear trajectory of leader cells.
  double A = 5.0; // to get actual video I changed this to 1.0 
  if (cell.is_follower==0){
    celula_nueva.position[0]= this->boxes_A[u][v][w].cells[cont_cell].position[0]  + A*params.time_step +
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[0];
    celula_nueva.position[1]= this->boxes_A[u][v][w].cells[cont_cell].position[1] + A*params.time_step + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[1];
    celula_nueva.position[2]= this->boxes_A[u][v][w].cells[cont_cell].position[2] + 
      this->movez * params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[2];
  }
  // FOLLOWER CELLS
  else {
    celula_nueva.position[0]= this->boxes_A[u][v][w].cells[cont_cell].position[0] + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[0];
    celula_nueva.position[1]= this->boxes_A[u][v][w].cells[cont_cell].position[1] + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[1];
    celula_nueva.position[2]= this->boxes_A[u][v][w].cells[cont_cell].position[2] + 
      this->movez * params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[2];
  }

  /* periodic boundary
      if (celula_nueva.position[0]<0) {
        celula_nueva.position[0]=params.lattice_length_x+celula_nueva.position[0];
      } 
      if (celula_nueva.position[0]>params.lattice_length_x) {
        celula_nueva.position[0]=celula_nueva.position[0]-params.lattice_length_x;
      }
      if (celula_nueva.position[1]<0) {
        celula_nueva.position[1]=params.lattice_length_y+celula_nueva.position[1];
      } 
      if (celula_nueva.position[1]>params.lattice_length_y) {
        celula_nueva.position[1]=celula_nueva.position[1]-params.lattice_length_y;
      }
      if (celula_nueva.position[2]<0) {
        celula_nueva.position[2]=params.lattice_length_z+celula_nueva.position[2];
      } 
      if (celula_nueva.position[2]>params.lattice_length_z) {
        celula_nueva.position[2]=celula_nueva.position[2]-params.lattice_length_z;
      }
  */

  if ((celula_nueva.position[0]<0)||(celula_nueva.position[0]>params.lattice_length_x)||
      (celula_nueva.position[1]<0)||(celula_nueva.position[1]>params.lattice_length_y)||
      (celula_nueva.position[2]<0)||(celula_nueva.position[2]>params.lattice_length_z)) {
      /*cout << " ** (Time " << reloj
	  << " ) ** WARNING: the new cell created from " << cell.name
	  << " is moving out of the domain (not supported) "
	  << " this cell will no longer be taken into account " << endl;*/
      if (cell.interaction_phenotype==0) this->phenotype1_count--;
      else this->phenotype2_count--;
  } else {

    //nueva caja
    int c1=(int)(floor(celula_nueva.position[0]/this->box_sizex)); 
    int c2=(int)(floor(celula_nueva.position[1]/this->box_sizey)); 
    int c3=(int)(floor(celula_nueva.position[2]/this->box_sizez)); 
  
    celula_nueva.box[0]=c1;
    celula_nueva.box[1]=c2;
    celula_nueva.box[2]=c3;

    //cout << " replace the cell in the box " << endl;
    /// @todo this should be done without copying the cell
    this->boxes_new_A[c1][c2][c3].cells.push_back(celula_nueva);
    
    // update maximum and minumum boxes if necessary
    if(this->maxx<c1) this->new_maxx=c1;
    if(this->maxy<c2) this->new_maxy=c2;	
    if(this->maxz<c3) this->new_maxz=c3;	
    if(this->minx>c1) this->new_minx=c1;
    if(this->miny>c2) this->new_miny=c2;
    if(this->minz>c3) this->new_minz=c3;
  } 
}

/* ************************************************************************
   moves the cells due to presence and position of fibres COMMENTED OUT               
   ************************************************************************ */
/*void CoupledModel::fibre_induced_movement(Cell& cell)
{
    //double pol_mag = sqrt(pow(cell.pol_axis[0],2)+pow(cell.pol_axis[1],2)+pow(cell.pol_axis[2],2));   // TODO include pol_axis[3] as a cell variable in Cell.h
    double pol_axis[3];
    double prob[fibres.size()];
    double rand_num = aleatorio(); /// TODO check this generates the right type of random number
    double sum = 0;
    double size = 0;
    double mtmx, mtmy, mtmz;

    for (unsigned int j=0; j<fibres.size(); j++){
        for (unsigned int i=0; i<=2; i++){
            pol_axis[i]=1.0;
            cell.pol_axis[i]/pol_mag;   /// TODO include pol_axis[3] as a cell variable in Cell.h
            prob[j] = prob[j] + pol_axis[i]*fibres[j].direction[i];
            size = size + prob[j]; 
        }
    
        for (unsigned int i=0; i<=2; i++){
            prob[j] = prob[j]/size;
        }

        if (rand_num > sum && rand_num <= sum+prob[j]){
            mtmx = fibres[j].direction[0];
            mtmy = fibres[j].direction[1];
            mtmz = fibres[j].direction[2];
        }
        else{
            sum = sum + prob[j];
        }
    }

    // double num_cells = this->total_no_of_cells;
    // double fibspeed[numcells];
    // for (unsigned int k=0; k<num_cells; k++){
    //   fibspeed[k] = sqrt(pow(mtmx,2)+pow(mtmy,2)+pow(mtmz,2));
    //   mtmx = (mtmx/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
    //   mtmy = (mtmy/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
    //   mtmz = (mtmz/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
    //} 
}*/

/* ************************************************************************
   computes cells and fibres in contact               
   ************************************************************************ */
void CoupledModel::compute_cell_fibres_contact(const int u,
					       const int v,const int w,
					       const int j)
{
  unsigned int total_cells_in_box=this->boxes_A[u][v][w].cells.size();

  for(unsigned int i=0; i<total_cells_in_box; i++){
    double cell_fibre_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->fibres[j]);
	  
    if(cell_fibre_min_dist < 2.0*(fibres[j].fradius + this->boxes_A[u][v][w].cells[i].radius)){
      if(params.fib_deg==1){
	    fibre_degradation(this->boxes_A[u][v][w].cells[i],fibres[j],cell_fibre_min_dist);
      }
	  this->boxes_A[u][v][w].cells[i].contact_fibres.push_back(&this->fibres[j]);

      //cout << " cell: " << this->boxes_A[u][v][w].cells[i].name
	  //	 << " in contact with fibre: " << j << endl;
    }
  }
}

/* *******************************************************************************
   Calculate the local distance between cells and fibres OLD VERSION COMMENTED OUT
   ******************************************************************************* */
/*void CoupledModel::compute_cell_fibres_contact(const int u,
					       const int v,const int w)
{
  unsigned int total_cells_in_box=this->boxes_A[u][v][w].cells.size();

  // loop over cells, find all the fibres in contact
  for(unsigned int i=0; i<total_cells_in_box; i++){
	    
    for(unsigned int j=0; j<fibres.size(); j++){

      // if a fibre exists 
      if (fibres[j].fibre_exists>0.5){
	    double dist_ref = boxes_A[u][v][w].cells[i].radius + fibres[j].length;
	    bool go_ahead = false;

        if ((boxes_A[u][v][w].cells[i].position[0]-
	    (fibres[j].start[0]+0.5*fibres[j].length*fibres[j].direction[0])<dist_ref)) {
	        if ((boxes_A[u][v][w].cells[i].position[1]-
	        (fibres[j].start[1]+0.5*fibres[j].length*fibres[j].direction[1])<dist_ref)) {
	            if ((boxes_A[u][v][w].cells[i].position[2]-
			    (fibres[j].start[2]+0.5*fibres[j].length*fibres[j].direction[2])<dist_ref)) {
			        go_ahead=true;
			    }
	        }
	    }
	
	    //double fibre_centre[3];
	    //double c_f_centre=0;
	    //for (unsigned int jj=0; jj<=2; jj++){
	    //fibre_centre[jj] = fibres[j].start[jj]+0.5*fibres[j].length*fibres[j].direction[jj];
	    //c_f_centre = c_f_centre+(boxes_A[u][v][w].cells[i].position[jj]-fibre_centre[jj])*
	    (boxes_A[u][v][w].cells[i].position[jj]-fibre_centre[jj]);
	  }
	
	  double cell_fibre_dist_check = sqrt(c_f_centre);

      // and the cell and fibre are minimally separated
	  //if(cell_fibre_dist_check < dist_ref){
		    if (go_ahead) {
	            // compute distance between cells and fibres
	            double cell_fibre_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->fibres[j]);
	            double cell_rad = this->boxes_A[u][v][w].cells[i].radius;
	            double fibre_rad = fibres[j].fradius;

	            // if the cell-fibre distance is less than the double cell radius + fibre radius
	            if(cell_fibre_min_dist < 2.0*(fibre_rad+cell_rad)){
	
	                ofstream fibint;
	                string fibint_list =  this->params.casedirectory + this->params.casename 
	                //+ "_fibre_interaction_info.txt";
	                fibint.open(fibint_list.c_str(),ios::app);
	                // fibint << " at time " << reloj << " cell " << this->boxes_A[u][v][w].cells[i].name 
	                // << " is close to fibre " << this->fibres[j].fibre_name 
	                //<< " the distance between the cell and fibre is " 
	                //<< cell_fibre_min_dist << endl;

                    // potentially degrade/remove fibres
	                if(params.fib_deg==1){
	                 fibre_degradation(this->boxes_A[u][v][w].cells[i],fibres[j],cell_fibre_min_dist);
	                }

	                if(fibres[j].fibre_exists>0.5){
		                this->boxes_A[u][v][w].cells[i].contact_fibres.push_back(&this->fibres[j]);
	                } else{
	                    //fibint << "    fibre " << fibres[j].fibre_name << " does not affect polarity " << endl;
	                }
	                //fibint << "    ************************** " << endl;
	                //fibint.close(); 	    
	            }
	        }
      }
    }
  }
}*/

/**********************************************************************************************/
/* calculate the local distance between cells and vessels */
/**********************************************************************************************/
void CoupledModel::compute_cell_vessels_contact(const int u,
		 const int v,const int w) 
{
  unsigned int total_cells_in_box = this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<total_cells_in_box; i++){

      this->boxes_A[u][v][w].cells[i].vessel_interaction = 0.0;
	    
      for(unsigned int j=0; j<vessels.size(); j++){

	    double vessel_centre[3];
	    double c_to_v_centre=0;
	    for (unsigned int jj=0; jj<=2; jj++){
	      vessel_centre[jj] = vessels[j].ves_start[jj]+0.5*vessels[j].ves_length*vessels[j].ves_direction[jj];
	      c_to_v_centre = c_to_v_centre+pow(boxes_A[u][v][w].cells[i].position[jj]-vessel_centre[jj],2);
	    }
      
	    double cell_vessel_dist_check = sqrt(c_to_v_centre);
	    double dist_ref = boxes_A[u][v][w].cells[i].radius + vessels[j].ves_length;
	    if(cell_vessel_dist_check < dist_ref){
	      // compute distance between cells and vessels
	      double cell_vessel_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->vessels[j]);
	      double cell_rad = this->boxes_A[u][v][w].cells[i].radius;
	      double vessel_rad = vessels[j].ves_radius;
	
	      if(cell_vessel_min_dist < (cell_rad+vessel_rad)){
		    if (params.verbose>1) {
		        cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in contact with vessel " << this->vessels[j].vessel_name << endl;
		        cout << " distance between the cell and vessel is " << cell_vessel_min_dist-(cell_rad+vessel_rad) << endl;
		    }
		    this->boxes_A[u][v][w].cells[i].contact_vessels.push_back(&this->vessels[j]);

		    this->boxes_A[u][v][w].cells[i].vessel_interaction = 1.0;
	      }
	    }
      }
  }
}

/* *****************************************************************************
   for each in the box u,v,w, compute all contacts with other cells
   ****************************************************************************** */
void CoupledModel::compute_cell_cell_contact(const int u, const int v,const int w)
{
  for(unsigned int i=0; i<this->boxes_A[u][v][w].cells.size(); i++) {
    // loop in neighbor boxes to find cells in contact with cell i
    // we consider a layer of [-1,1]x[-1,1]x[-1,1] 
    for(int h=-1;h<2;h++) {
      int borderx=u+h;
      for(int l=-1;l<2;l++) {
	    int bordery=v+l;
	    for(int m=-1;m<2;m++) {
	        int borderz=w+m;		   		   
	  
	        // check that we are looking inside the computational domain
	        if( borderx<=(int) params.boxesx-1 && borderx>=0 && 
	             bordery<=(int) params.boxesy-1 && bordery>=0 && 
	             borderz<=(int) params.boxesz-1 && borderz>=0 ) {

	             for(unsigned int j=0; j<this->boxes_A[u+h][v+l][w+m].cells.size(); j++) {
	             // check that we are not looking at the same cell
	                if(this->boxes_A[u][v][w].cells[i].name !=this->boxes_A[u+h][v+l][w+m].cells[j].name) {	
		                // compute distance between cells
		                double cell_cell_min = DISTANCE(this->boxes_A[u][v][w].cells[i],
					    this->boxes_A[u+h][v+l][w+m].cells[j]); 
		                double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
		                this->boxes_A[u][v][w].cells[i].radius;
		
		                // cells in contact
		                if(cell_cell_min < cell_cell_center){
                            // store the pointer to the neighbor cell in the cell. vector
		                    this->boxes_A[u][v][w].cells[i].neighbors.push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
		                    // increase number of contacts
		                    this->boxes_A[u][v][w].cells[i].contacts++;

		                    if (params.verbose>3) {
		                        cout << " cell " << this->boxes_A[u][v][w].cells[i].name
			                    << " in  contact with "
			                    << this->boxes_A[u+h][v+l][w+m].cells[j].name << endl;
		                    }
		                } else {
		                    if (params.follower_force>0) {
		                        if ( (this->boxes_A[u][v][w].cells[i].birthday < reloj)&&
			                        (this->boxes_A[u+h][v+l][w+m].cells[j].birthday < reloj)) {
		                            // check if cells were in contact before
		                            // compute distance between cells
		                            double cell_cell_dist_old = 
			                           pow(this->boxes_A[u][v][w].cells[i].position_old[0] -
			                            this->boxes_A[u+h][v+l][w+m].cells[j].position_old[0],2)+
			                           pow(this->boxes_A[u][v][w].cells[i].position_old[1] -
			                            this->boxes_A[u+h][v+l][w+m].cells[j].position_old[1],2)+
			                           pow(this->boxes_A[u][v][w].cells[i].position_old[2] -
			                            this->boxes_A[u+h][v+l][w+m].cells[j].position_old[2],2);
		                            cell_cell_dist_old = sqrt(cell_cell_dist_old);
		      
		                            double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
			                        this->boxes_A[u][v][w].cells[i].radius;
		      
		                            if( cell_cell_dist_old < cell_cell_center){
			                            /*cout << " time step: " << this->reloj
			                            << ", cells " << this->boxes_A[u][v][w].cells[i].name
			                            << " and " << this->boxes_A[u+h][v+l][w+m].cells[j].name
			                            << " were in contact at the previous time step " << endl;
			                            */
			                            this->boxes_A[u][v][w].cells[i].neighbors.
			                            push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
			                            // increase number of contacts
			                            this->boxes_A[u][v][w].cells[i].contacts++;
		                            }
		                        }
		                    }
		                }//end if
	                }
	             } // for(unsigned int j=0;j<neighbors box cells;j++)    
	        } //end if (inside domain) 
	    }  
      }
    } //end for loops 
  } // loop of cells in box
}

/* *******************************************************************************
   calculate all forces acting on cells in the box u,v,w
   ****************************************************************************** */
void CoupledModel::compute_all_forces(const int u, const int v,const int w) 
{
  unsigned int n_cells = this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<n_cells; i++) {
    // set forces to 0
    for (unsigned int j=0; j<3; j++){
      this->boxes_A[u][v][w].cells[i].force[j]=0.;
    }
    // ******* CELL-CELLS ********
    this->cell_cell_interaction(this->boxes_A[u][v][w].cells[i]);

    // ******* CELL-VESSELS ********
    this->cell_vessel_interaction(this->boxes_A[u][v][w].cells[i]);

    // ******* CELL-FIBRES ********
    // set polarity (needed in the interaction with fibres)
    // note: this might not be needed now since we do not update velocity yet
    for (unsigned int k=0; k<3; k++) {
      this->boxes_A[u][v][w].cells[i].polarity[k] = this->boxes_A[u][v][w].cells[i].vel[k];
    } 
    //velocity correction (fibres and polarization)
    if (this->boxes_A[u][v][w].cells[i].contact_fibres.size()){      
      /*ofstream fibint;
      string fibint_list =  this->params.casedirectory + this->params.casename + "_fibre_interaction_info.txt";
      fibint.open(fibint_list.c_str(),ios::app);
      fibint << " number of fibres in contact is: " << this->boxes_A[u][v][w].cells[i].contact_fibres.size()
	     << endl;
      //fibint << "    CELL FIBRE INTERACTION ACTIVATED " << endl;
       fibint << "    the old cell polarity is " << this->boxes_A[u][v][w].cells[i].polarity[0]
	     << " " << this->boxes_A[u][v][w].cells[i].polarity[1] << " "
	     << this->boxes_A[u][v][w].cells[i].polarity[2] << endl;*/
          
      this->cell_fibres_interaction(this->boxes_A[u][v][w].cells[i]);
      // is this needed here?
      /*for (int k=0; k<3; k++) {
    	this->boxes_A[u][v][w].cells[i].polarity[k] = this->boxes_A[u][v][w].cells[i].vel[k];
	  }*/
      
      /*fibint << "    the new cell polarity is " << this->boxes_A[u][v][w].cells[i].polarity[0]
	     << " " << this->boxes_A[u][v][w].cells[i].polarity[1] << " "
	     << this->boxes_A[u][v][w].cells[i].polarity[2] << endl;
      fibint << " --------------------" << endl;
      fibint.close();*/     
    }
    if (params.alpha_birthrate[0]==0){
      ofstream cellforcesFile;
      string forces = this->params.casedirectory + this->params.casename + "_forces.txt";
      cellforcesFile.open(forces.c_str(),ios::app);
	      cellforcesFile << this->boxes_A[u][v][w].cells[i].force[0] << " "
			       << this->boxes_A[u][v][w].cells[i].force[1] << " "
			       << this->boxes_A[u][v][w].cells[i].force[2] << " "
			     << this->boxes_A[u][v][w].cells[i].contact_fibres.size() << " " << endl;
      cellforcesFile.close();
    }
  }//End cells loop
}

/* *******************************************************************************
   contact forces CURRENTLY UNUSED SO COMMENTED OUT
   ***************************************************************************** */
/*void CoupledModel::contact_forces(const int u,
		 const int v,const int w,
		 const unsigned int cells_number) 
{
  //cout << "contact_forces function is ACTIVATED" << endl;
  int borderx,bordery,borderz;
  
  for(unsigned int i=0; i<cells_number; i++) {

    // loop in neighbor boxes to find cells in contact with cell i
    // we consider a layer of [-1,1]x[-1,1]x[-1,1] 
    for(int h=-1;h<2;h++) {
      borderx=u+h;
      for(int l=-1;l<2;l++) {
	    bordery=v+l;
	    for(int m=-1;m<2;m++) {
	        borderz=w+m;		   		   
	  
	        // check that we are looking inside the computational domain
	        if( borderx<=(int) params.boxesx-1 && borderx>=0 && 
	            bordery<=(int) params.boxesy-1 && bordery>=0 && 
	            borderz<=(int) params.boxesz-1 && borderz>=0 ) {

	            for(unsigned int j=0; j<this->boxes_A[u+h][v+l][w+m].cells.size(); j++) {
	                // check that we are not looking at the same cell
	                if(this->boxes_A[u][v][w].cells[i].name != this->boxes_A[u+h][v+l][w+m].cells[j].name) {
		                // compute distance between cells
		                double cell_cell_min = DISTANCE(this->boxes_A[u][v][w].cells[i],
					    this->boxes_A[u+h][v+l][w+m].cells[j]); 
		                double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
		                this->boxes_A[u][v][w].cells[i].radius;
		
		                // cells in contact
		                if(cell_cell_min < cell_cell_center){
		                    // store the pointer to the neighbor cell in the cell. vector
		                    this->boxes_A[u][v][w].cells[i].neighbors.push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
		                    // increase number of contacts
		                    this->boxes_A[u][v][w].cells[i].contacts++;

		                    if (params.verbose>3) {
		                        cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in  contact with "
			                    << this->boxes_A[u+h][v+l][w+m].cells[j].name << endl;
		                    }
		                } else {

		                    if (params.follower_force>0) {
		                        if ( (this->boxes_A[u][v][w].cells[i].birthday < reloj)&&
			                        (this->boxes_A[u+h][v+l][w+m].cells[j].birthday < reloj)) {
		                            // check if cells were in contact before
		                            // compute distance between cells
		                            double cell_cell_dist_old = 
			                            pow(this->boxes_A[u][v][w].cells[i].position_old[0] -
			                                this->boxes_A[u+h][v+l][w+m].cells[j].position_old[0],2)+
			                            pow(this->boxes_A[u][v][w].cells[i].position_old[1] -
			                                this->boxes_A[u+h][v+l][w+m].cells[j].position_old[1],2)+
			                            pow(this->boxes_A[u][v][w].cells[i].position_old[2] -
			                                this->boxes_A[u+h][v+l][w+m].cells[j].position_old[2],2);
		                            cell_cell_dist_old = sqrt(cell_cell_dist_old);
		      
		                            double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
			                        this->boxes_A[u][v][w].cells[i].radius;
		      
		                            if( cell_cell_dist_old < cell_cell_center){
			                            //cout << " time step: " << this->reloj
			                            //<< ", cells " << this->boxes_A[u][v][w].cells[i].name
			                            //<< " and " << this->boxes_A[u+h][v+l][w+m].cells[j].name
			                            //<< " were in contact at the previous time step " << endl;
			                 
			                            this->boxes_A[u][v][w].cells[i].neighbors.
			                            push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
			                            // increase number of contacts
			                            this->boxes_A[u][v][w].cells[i].contacts++;
		                            }
		                        }
		                    }
		                }//end if
	                }
	            } // for(unsigned int j=0;j<neighbors box cells;j++)
	        }//end if (inside domain)
	    }  
      }
    }//end for loops
    
    // for debugging: print all contacts
    if (params.verbose > 3) {
      if (this->boxes_A[u][v][w].cells[i].contacts) {
 	cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " has " 
	     << this->boxes_A[u][v][w].cells[i].contacts << " contacts with ";
      }
      for (unsigned int j=0; j<this->boxes_A[u][v][w].cells[i].neighbors.size(); j++) {
	cout << this->boxes_A[u][v][w].cells[i].neighbors[j]->name << " ";
      }
      cout << endl;
    }

    //polarise cells - CICELY NEW
    //polarise(this->boxes_A[u][v][w].cells[i]);
    //compute cell-cell forces and new cell velocity
    // Equation (for the velocity)
    //x_dot = 1/tissue_friction * sum_{cells in contact} f_contact
      //       + 1/n_fibres_in_contact * sum_{fibres in contact} direction
      //     + polarity_velocity

    // contact forces
    hertz(this->boxes_A[u][v][w].cells[i]);
    for (int k=0; k<3; k++) {
      this->boxes_A[u][v][w].cells[i].polarity[k] = this->boxes_A[u][v][w].cells[i].vel[k];
    }
       
    //velocity correction (fibres and polarization)
    if (0) { //this->boxes_A[u][v][w].cells[i].contact_fibres.size()){
      ofstream fibint;
      string fibint_list =  this->params.casedirectory + this->params.casename + "_fibre_interaction_info.txt";
      fibint.open(fibint_list.c_str(),ios::app);
      fibint << " number of fibres in contact is: " << this->boxes_A[u][v][w].cells[i].contact_fibres.size() << endl;
      fibint << "    CELL FIBRE INTERACTION ACTIVATED " << endl;
      fibint << "    the old cell polarity is " << this->boxes_A[u][v][w].cells[i].polarity[0] << " " << this->boxes_A[u][v][w].cells[i].polarity[1] << " " << this->boxes_A[u][v][w].cells[i].polarity[2] << endl;
      cell_fibres_interaction(this->boxes_A[u][v][w].cells[i]);
      for (int k=0; k<3; k++) {
    	this->boxes_A[u][v][w].cells[i].polarity[k] = this->boxes_A[u][v][w].cells[i].vel[k];
      }
      fibint << "    the new cell polarity is " << this->boxes_A[u][v][w].cells[i].polarity[0] << " " << this->boxes_A[u][v][w].cells[i].polarity[1] << " " << this->boxes_A[u][v][w].cells[i].polarity[2] << endl;
      fibint << " --------------------" << endl;
      fibint.close();
    }

    
    // check if hypoxic cells can become normoxic again
    if (this->boxes_A[u][v][w].cells[i].type==2) {
        reverse_phenotype(this->boxes_A[u][v][w].cells[i]);
    }
  }//End cells loop
}*/

/* ****************************************************************************
   assign each cell in box[u][v][w] to the element with closest baricentre
   **************************************************************************** */  
void CoupledModel::compare_elements(int u, int v, int w,
				    Mesh& _mesh) 
{
  int borderx,bordery,borderz;
  
  //Now we run over the number of cells in box
  for(unsigned int i=0;i<this->boxes_A[u][v][w].cells.size();i++) {
    float minimaDistanciaElem=200;
    int minimoElem=-1;
    
    for (unsigned int ii = 0; ii < this->boxes_A[u][v][w].v_triangles.size(); ii++) {
      double x,y,z,dist_tri;
      unsigned int it = this->boxes_A[u][v][w].v_triangles[ii];

      x = (_mesh.xp[_mesh.tetra[4*it]-1] +  _mesh.xp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.xp[_mesh.tetra[4*it+2]-1] + _mesh.xp[_mesh.tetra[4*it+3]-1])/4.;
      
      y = (_mesh.yp[_mesh.tetra[4*it]-1] +  _mesh.yp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.yp[_mesh.tetra[4*it+2]-1] + _mesh.yp[_mesh.tetra[4*it+3]-1])/4.;
      
      z = (_mesh.zp[_mesh.tetra[4*it]-1] +  _mesh.zp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.zp[_mesh.tetra[4*it+2]-1] + _mesh.zp[_mesh.tetra[4*it+3]-1])/4.;

      dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
		    pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
		    +pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
      
      if(dist_tri<minimaDistanciaElem){
	        minimaDistanciaElem=dist_tri;
	        minimoElem = it;
      }  
    }

    if (minimoElem>=0) {
      _mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1;
      // ------------------------
      // fill array per type
      int cellType =  this->boxes_A[u][v][w].cells[i].type;
      if ( cellType==1) {
	    _mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1;
      } else if (cellType==2) {
	    _mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
      } else {
	    _mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1;
      }
    }

    //If there is no triangle we look in the nearby boxes
    if(minimaDistanciaElem==200){
        for(int h=-1;h<2;h++) {
	        borderx=u+h;
	        for(int l=-1;l<2;l++) {
	            bordery=v+l;							   
	            for(int m=-1;m<2;m++) {
	                borderz=w+m;		   		   
	                if(borderx<=(int) params.boxesx-1 && borderx>=0 && 
	                    bordery<=(int) params.boxesy-1 && bordery>=0 && 
	                    borderz<=(int) params.boxesz-1 && borderz>=0){
                        //We are inside of the boxes domain
 
		                for (unsigned int ii = 0; ii < this->boxes_A[u+h][v+l][w+m].v_triangles.size(); ii++) {
		                    double x,y,z,dist_tri;
		                    unsigned int it2 = this->boxes_A[u+h][v+l][w+m].v_triangles[ii];
		                    x = (_mesh.xp[_mesh.tetra[4*it2]-1] +  _mesh.xp[_mesh.tetra[4*it2+1]-1] 
		                     +  _mesh.xp[_mesh.tetra[4*it2+2]-1] + _mesh.xp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		                    y = (_mesh.yp[_mesh.tetra[4*it2]-1] +  _mesh.yp[_mesh.tetra[4*it2+1]-1] 
		                      +  _mesh.yp[_mesh.tetra[4*it2+2]-1] + _mesh.yp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		                    z = (_mesh.zp[_mesh.tetra[4*it2]-1] +  _mesh.zp[_mesh.tetra[4*it2+1]-1] 
		                      +  _mesh.zp[_mesh.tetra[4*it2+2]-1] + _mesh.zp[_mesh.tetra[4*it2+3]-1])/4.;

                            /// @attention there was a bug here: boxes_A[u+h][v+l][w+m] instead of boxes_A[u][v][w]
		                    dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
				                pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
				                +pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
	  
	                        if(dist_tri<minimaDistanciaElem){
		                      minimaDistanciaElem=dist_tri;
		                      minimoElem = it2;
		                    }
		                }
	                }
	            }
	        }
        }

        if (minimoElem>=0) {
	        _mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1; 
    
	        // ------------------------
	        // fill array per type
            int cellType =  this->boxes_A[u][v][w].cells[i].type;
            if ( cellType==1) {
	            _mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1; 
	        } else if (cellType==2) {
	            _mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
	        } else {
	            _mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1; 
	        }
	        // ------------------------
        }
    }
  } 
}

/* ****************************************************************************
   search in box
   **************************************************************************** */
/// @todo write a version search_in_box(Box b, Cell c)
int CoupledModel::search_in_box(const vector<int>& box_number,
				const unsigned int cell_name)
{

  int u=box_number[0];
  int v=box_number[1];
  int w=box_number[2];
  
  for(unsigned int i=0;i<this->boxes_A[u][v][w].cells.size();i++) {
    if(cell_name == this->boxes_A[u][v][w].cells[i].name) return i;
  }
  cout << "WARNING search_in_box(): cell " << cell_name << " not found in box "
       << u << "," << v << "," << w << endl;
  return -1;
}

/* ****************************************************************************
   MAIN LOOP
   **************************************************************************** */
void CoupledModel::loop()
{
  /*
     @todo move this to the init function
     for this, we need to replace fe_mesh with oxy_diff.mesh everywhere
  */
  Mesh fe_mesh;    
  if (this->params.femSolverType>0) {
    /// @todo move this part in init() function
    // ======================
    // read FE-mesh from file
    fe_mesh.read(this->params.meshdir + this->params.meshname);
    fe_mesh.info();
    // set which elements are contained in which box
    setElementsInBox(fe_mesh);
    // initialize arrays containing cell densities
    fe_mesh.initCellsArrays();
    // check boxes-tetra initialization
    if (this->params.verbose>2) {
      checkElementsInBoxes();
    }
    // ======================
    
    // set the dimension of the problem
    // according to the mesh
    if (fe_mesh.dim != (int) this->params.dimension) {
      cout << " !WARNING! input dimension = " << this->params.dimension
	   << " and mesh dimension = " << fe_mesh.dim << endl;
      cout << " I am changing the input dimension to " <<  fe_mesh.dim 
	   << " (file " << __FILE__ << ", line " 
	   << __LINE__ << ")" << endl;
    }
    
    // initialize PDE solver
    this->oxy_diff.init(this->params, fe_mesh);
    // ======================
  } else {
    // initialize empty PDE solver
    this->oxy_diff.init();
  }

  // initialize cell type counters
  this->totNorm.resize(0);
  this->totHypo.resize(0);
  this->totDead.resize(0);

  cout << " ========= " << endl;
  system(" date\n");
  cout << " starting main loop, end time: " << params.n_steps << endl;
  cout << " ========= " << endl;
  this->reloj=0; // initial time step
  
  // ====================
  // MAIN LOOP
  // ====================
  while(this->total_no_of_cells < this->max_cell && 
	    this->reloj < params.n_steps) {

        ///@todo remove when we have a better measure of density change
        //int n_cells_old = this->total_no_of_cells;

        if ((this->reloj%100)==0) {
	        cout << " *** time step: " << this->reloj << " (of " << params.n_steps << ") " 
	        << " *** n. of cells: " << this->total_no_of_cells
            //<< " *** n. of dead cells: " << this->total_no_of_removed_cells
            << endl;
            cout << "     count per type " ;
            this->count_cells_per_type();
        }
    
        if (params.verbose>0) {
            if ((this->reloj%500)==0) {
	            cout << " *** time step: " << this->reloj << " (of " << params.n_steps << ") " 
	            << " *** n. of cells: " << this->total_no_of_cells
                << " *** n. of dead cells: " << this->total_no_of_removed_cells << endl;
            }
            // count how many cells per type (and save on file)
            if ((this->reloj%500)==0) {
	        cout << " count per type " << endl;
	        this->count_cells_per_type();
       	    } 
        }
    
        // increase coupled model iteration number
        this->reloj++;
        // increase PDE iteration number
        this->oxy_diff.time = this->oxy_diff.time + 1.;

        // ********************************
        // Solve for nutrient concentration
        // ********************************
        if(this->params.femSolverType>0) {
   
            // compute number of cells/element
            fe_mesh.initCellsArrays(); // set to zero the counter of cells in each tria
            for(int k=this->minx; k<=this->maxx; k++) {
	            for(int l=this->miny; l<=this->maxy; l++) {
	                for(int n=this->minz; n<=this->maxz; n++) {
	                    this->compare_elements(k,l,n,fe_mesh);
	                }
	            }
            }
      
            double density_change = fe_mesh.cellDensityChange();

            // check if the solver has to be called
            /// @todo generalize with density change
            double relative_density_change = density_change/this->total_no_of_cells;

            if (total_no_of_cells<200) {
	            this->oxy_diff.launch = (this->reloj==1) || (this->reloj%100==0);
	            this->oxy_diff.launch = (this->oxy_diff.launch || (relative_density_change>0.1) );
            }
            else if (total_no_of_cells<2000) {
	            this->oxy_diff.launch = ( (relative_density_change>0.1) || (this->reloj%400==0) );
            }
            else {
	            this->oxy_diff.launch = ( (relative_density_change>0.1) || (this->reloj%1000==0) );
            }
        
            /*
            } else if (total_no_of_cells<5000) {
	            this->oxy_diff.launch = ((density_change>80) || (this->reloj%500==0) );
	
            } else if (total_no_of_cells<10000) {
                this->oxy_diff.launch = ((density_change>100) ||  (this->reloj%3000==0)) ;
	            this->oxy_diff.launch = this->oxy_diff.launch &&  (relative_density_change>0.08);
	
            } else {
	            this->oxy_diff.launch = (density_change>200) && (relative_density_change>0.1);	
            }
            */

            // We now run the fem every iteration
            if (this->oxy_diff.launch) {
                // REMOVED COUT
                //cout << " n cells: " << total_no_of_cells << " --> density change: " << density_change << ", relative: " << relative_density_change << endl;
                // REMOVED COUT
                //cout <<  "write piecewise constant cell density (-> FreeFem)" << endl;
                ofstream o_file;
                o_file.open(this->params.fileCellsDensity2FEM.c_str(),ios::out);
                for (unsigned int i=0; i<fe_mesh.nElem; i++) {	  
	                o_file << fe_mesh.cellsInTria[i] << " "
		            << fe_mesh.cellsInTriaNorm[i] << " " // normoxic
		            << fe_mesh.cellsInTriaHypo[i] << " " // hypoxic
		            << fe_mesh.cellsInTriaDead[i] << " " // dead
		            << endl;
                }
                o_file.close();
     
                // store current cell density
                fe_mesh.storeCellDensity();
	
                // write file containing cell information
                ///@todo change interface: provide only cell list and positions
                ofstream outCellFile;
                outCellFile.open(this->params.fileCells2FEM.c_str(),ios::out);
                if (!outCellFile) {
	                cerr << " *** ERROR, file " << __FILE__ << ", line " << __LINE__ 
	                << " *** could not open file " << this->params.fileCells2FEM << endl;
	                exit(1);
                }
	
                // header: length, iteration (pde), iteration (total)
                outCellFile << this->total_no_of_cells << " "
                        << this->oxy_diff.iter << " "
		                << this->reloj << " "
		                << endl;
	            // body: x,y,z,type,r,phenotype,adhesion coeff,id,energy
	            for(int k=this->minx; k<=this->maxx; k++) {
	                for(int l=this->miny; l<=this->maxy; l++){
	                    for(int n=this->minz; n<=this->maxz; n++) {
	                        for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
		                        outCellFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
			                    << this->boxes_A[k][l][n].cells[i].position[1] << " "
			                    << this->boxes_A[k][l][n].cells[i].position[2] << " "
			                    << this->boxes_A[k][l][n].cells[i].type << " "
			                    << this->boxes_A[k][l][n].cells[i].radius << " "
		                        // << this->boxes_A[k][l][n].cells[i].phenotype << " " //REMOVED 25/6/19 TOMMASO
			                    << this->boxes_A[k][l][n].cells[i].cont_pheno << " " //ADDED 25/6/19 TOMMASO
			                    << this->boxes_A[k][l][n].cells[i].adhesion << " "
			                    << this->boxes_A[k][l][n].cells[i].name << " "
			                    << this->boxes_A[k][l][n].cells[i].energy << endl;
	                        }
                        }
	                }
                }
                outCellFile.close();

                // solve PDE (and write new concentration)
                // REMOVED COUT
                cout << " solve PDE at time " << reloj << endl;
                this->oxy_diff.solve(reloj);

                // read the new concentration file
                ifstream o2_conc_file;
                string ifname = this->params.fileFEM2Cells;
                o2_conc_file.open(ifname.c_str(),ios::in);

                unsigned int ic = 0;
                for(int k=this->minx; k<=this->maxx; k++) {
	              for(int l=this->miny; l<=this->maxy; l++){
	                for(int n=this->minz; n<=this->maxz; n++) {
	                  for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].O2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dxO2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dyO2;
		                o2_conc_file >> this->boxes_A[k][l][n].cells[i].dzO2;
		                ic++;

		                // change cell status depending on new concentration
		                //oxy_in_cell(this->boxes_A[k][l][n].cells[i]);
                        // change cell status depending on the O2 phenotype
                        // phenotype_of_cell(this->boxes_A[k][l][n].cells[i]); TOMMASO 21/4/22 removed from here since we want this calculated everytime!
		                //cell_mutation(this->boxes_A[k][l][n].cells[i]); // TOMMASO 6/08/19 - removed function here.
	                  }
	                }
	              }
                }
                o2_conc_file.close();
            }
        }
    
        // ********************************
        // Loop of mutation, birth and death
        // ********************************
        cells_counter=0;  // total number of new cells      
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
	                this->compute_cell_cell_contact(k,l,n); // TOMMASO 7/08/19 - need this as otherwise number of contacts is zero. !!! ALFONSO !!! 
	                for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                        this->cell_mutation(this->boxes_A[k][l][n].cells[i]); // TOMMASO 6/08/19 - moved function here.
                        this->phenotype_of_cell(this->boxes_A[k][l][n].cells[i]); // TOMMASO 21/4/22 - moved function here.
                        this->cell_birth(this->boxes_A[k][l][n].cells[i]);
                        const unsigned int cell_death_status = 3;
                        if (this->boxes_A[k][l][n].cells[i].type == cell_death_status) {
                            //this->cell_death(this->boxes_A[k][l][n].cells[i]);
                            /*cout << "DEATH at time " << reloj << endl;
                            cout << "cell " << this->boxes_A[k][l][n].cells[i].name << " is of type "
                                 << this->boxes_A[k][l][n].cells[i].type << " and will now be removed " << endl;*/
                            this->total_no_of_removed_cells += 1;
                            this->boxes_A[k][l][n].cells[i] = this->boxes_A[k][l][n].cells[
                                    this->boxes_A[k][l][n].cells.size() - 1];
                            this->boxes_A[k][l][n].cells.pop_back();
                            if (boxes_A[k][l][n].cells.size() == 0) break;
                        }
                    }
                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                        this->boxes_A[k][l][n].cells[i].clear_contacts(); // TOMMASO 23/08/19 - needed not to double number of contacts used below. !!! ALFONSO !!
                        //vector<Cell> swap(this->boxes_A[k][l][n].cells);
	                }				
	            }       
            }
        }

        // update maximum number of occupied boxes after births
        this->update_maximum();

        // update number of cells
        this->total_no_of_cells=this->total_no_of_cells+cells_counter;
        if (cells_counter) {
            if (params.verbose>2) {
                cout << " *** time step " << this->reloj << ", newborn: " << cells_counter
	            << " total number of cells: " << this->total_no_of_cells << endl;
	            cout << " phenotype 1 - " << this->phenotype1_count <<  " cells "
	            << " phenotype 2 - " << this->phenotype2_count << " cells " << endl;
            }
        }

        // double help = this->total_no_of_cells-this->phenotype1_count;
        //   if (help>0){
        //     cout << " !!!!!!! HELP !!!!!!!!!!!!!" << endl;
        //     	cout << " total number of cells: " << this->total_no_of_cells << endl;
        // 	cout << " phenotype 1 - " << this->phenotype1_count <<  " cells "
        // 	     << " phenotype 2 - " << this->phenotype2_count << " cells " << endl;
        // }

        double dist_centre_x = 0.0;
        double dist_centre_y = 0.0;
        double dist_centre_z = 0.0;
        double dist_centre = 0.0;
        double mean_dist = 0.0;
        double variance_dist = 0.0;

        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++){
	            for(int n=this->minz; n<=this->maxz; n++) {
                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	                    dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	                    dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	                    dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	                    dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);
	                    mean_dist += dist_centre;
	                }
                }
            }
        }
        mean_dist = mean_dist/this->total_no_of_cells;

        for(int k=this->minx; k<=this->maxx; k++) {
          for(int l=this->miny; l<=this->maxy; l++){
	        for(int n=this->minz; n<=this->maxz; n++) {
	          for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	              dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	              dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	              dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	              dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);
	              variance_dist += (dist_centre-mean_dist)*(dist_centre-mean_dist);
	          }
            }
          }
        }
        variance_dist = variance_dist/this->total_no_of_cells;

        if (this->params.verbose>1){
          // =================================================================
          // write mean position (relative to centre of domain) to data file
          // =================================================================
          ofstream meandist;
          string meandist_list =  this->params.casedirectory + this->params.casename + "_mean_distance.txt";
          meandist.open(meandist_list.c_str(),ios::app);
          meandist << mean_dist << " " ;
          meandist.close();
          // ********************************
       
          // ============================================
          // write variance of above to data file
          // ============================================
          ofstream vardist;
          string vardist_list =  this->params.casedirectory + this->params.casename + "_variance_distance.txt";
          vardist.open(vardist_list.c_str(),ios::app);
          vardist << variance_dist << " " ;
          vardist.close();
          // ********************************
        }
    
        // ============================================
        // write intermediary cell numbers to data file
        // ============================================
        ofstream pheno1File;
        string full_list =  this->params.casedirectory + this->params.casename + "_full.txt";
        pheno1File.open(full_list.c_str(),ios::app);
        pheno1File << this->phenotype1_count << " " ;
        pheno1File.close();
        // if (params.ic_phenotype.size()==2){
        //   ofstream pheno2File;
        //   pheno2File.open("phenotype2_cellnumbers.txt",ios::app);
        //   pheno2File << this->phenotype2_count << " " ;
        //   pheno2File.close();
        // }
        // ********************************

        // ********************************
        // Loop of forces and movement
        // ********************************
        // count cells for each box (before movement)
        if (this->params.verbose>0) {
          this->count_cells_per_box();
        }
    
        // ============================================
        // NEW VERSION
        // ============================================
        // compute contacts and forces
        //int count_f = 0;

        for(unsigned int j=0; j<fibres.size(); j++){

            //calculate fibre crosslinks
            if (fibres[j].fibre_exists>0.5) {
                //this->compute_fibre_crosslinks(k,l,n,j);
                this->compute_fibre_crosslinks(j);
                crosslink_count = 0;
            }


            // check if the fibre intersects the cells
            int min_fibre_x = (int) fibres[j].start[0]/this->box_sizex;
            int max_fibre_x = (int) fibres[j].end[0]/this->box_sizex + 1;
            if (min_fibre_x>max_fibre_x) {
	            min_fibre_x = max_fibre_x;
	            max_fibre_x = (int) fibres[j].start[0]/this->box_sizex;
            }
            int min_fibre_y = (int) fibres[j].start[1]/this->box_sizey;
            int max_fibre_y = (int) fibres[j].end[1]/this->box_sizey + 1;
            if (min_fibre_y>max_fibre_y) {
                min_fibre_y = max_fibre_y;
	                max_fibre_y = (int) fibres[j].start[1]/this->box_sizey;
            }
            int min_fibre_z = (int) fibres[j].start[2]/this->box_sizey;
            int max_fibre_z = (int) fibres[j].end[2]/this->box_sizey + 1;
            if (min_fibre_z>max_fibre_z) {
	            min_fibre_z = max_fibre_z;
	            max_fibre_z = (int) fibres[j].start[2]/this->box_sizez;
            }
      
            // take the largest 'min'
            if (min_fibre_x<this->minx) min_fibre_x = this->minx;
            if (min_fibre_y<this->miny) min_fibre_y = this->miny;
            if (min_fibre_z<this->minz) min_fibre_z = this->minz;
      
            // take the smallest 'max'
            if (max_fibre_x>this->maxx) max_fibre_x = this->maxx;
            if (max_fibre_y>this->maxy) max_fibre_y = this->maxy;
            if (max_fibre_z>this->maxz) max_fibre_z = this->maxz;

            //bool tmp_bool = true;
            for(int k=min_fibre_x; k<=max_fibre_x; k++) {
	            for(int l=min_fibre_y; l<=max_fibre_y; l++) {
	                for(int n=min_fibre_z; n<=max_fibre_z; n++) {
	                    /*if (tmp_bool) {
	                        count_f++;
	                        tmp_bool = false;
	                    }*/
	    
	                    if (fibres[j].fibre_exists>0.5) {
	                        this->compute_cell_fibres_contact(k,l,n,j);
	                    }
	                }
	            }
            }
        }
        //cout << " I have checked " << count_f << " fibres." << endl;
    
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
    	            // contacts
    	            //this->compute_cell_fibres_contact(k,l,n);
    	            this->compute_cell_vessels_contact(k,l,n);
    	            this->compute_cell_cell_contact(k,l,n);									  
    	            //forces
    	            this->compute_all_forces(k,l,n);	  
    	        }	
            }
        }

        // update velocity (and other operations, e.g., oxygen/dependent)
        ///@todo merge this loop with the movement loop
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
                for(int n=this->minz; n<=this->maxz; n++) {
    	            // update velocity
    	            for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
    	                this->update_cell_velocity(this->boxes_A[k][l][n].cells[i]);

                        // POSSIBLE CHANGE ???
	                    /*for (int j=0; j<3; j++) {
	                        this->boxes_A[k][l][n].cells[i].polarity[j] = this->boxes_A[k][l][n].cells[i].vel[j];
	                    }*/

	                    /*fibint << "    the new cell polarity is " << this->boxes_A[k][l][n].cells[i].polarity[0]
	                    << " " << this->boxes_A[k][l][n].cells[i].polarity[1] << " "
	                    << this->boxes_A[k][l][n].cells[i].polarity[2] << endl;
                        fibint << " --------------------" << endl;
                        fibint.close();*/
	    
    	                /* check if hypoxic cells can become normoxic again
                        switched this off for now [TOMMASO]
    	                if (this->boxes_A[k][l][n].cells[i].type==2) {
	                        //this->reverse_phenotype(this->boxes_A[k][l][n].cells[i]);
    	                }
                        */
    	            }
    	        }
            }
        }

        // ============================================
        // OLD VERSION
        // ============================================
        // Compute contacts and  forces
        /*for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
                for(int n=this->minz; n<=this->maxz; n++) {
     	            this->compute_cell_fibres_contact(k,l,n);
     	            this->compute_cell_vessels_contact(k,l,n);    	  
     	            // the number of cells in the box must be stored since it
     	            // might change during the contact force computation ?
     	            initial_cells_in_box=this->boxes_A[k][l][n].cells.size();
     	            this->contact_forces(k,l,n,initial_cells_in_box);
     	        }
            }
        }*/

        // move cells
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++) {
	            for(int n=this->minz; n<=this->maxz; n++) {
	                for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	                    this->movement(this->boxes_A[k][l][n].cells[i],k,l,n,i);
	                }           
	            }
            }
        } 
   
        // print all cells infos (if required)
        if (this->params.verbose>4) {
            for(int k=this->minx; k<=this->maxx; k++) {
	            for(int l=this->miny; l<=this->maxy; l++){
	                for(int n=this->minz; n<=this->maxz; n++) {
	                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	                        this->boxes_A[k][l][n].cells[i].printInfo();
	                    }
	                }
	            }
            }
        }
        this->update_maximum();
        this->update_box();
    

        // output routines
        // write output list of cells
        std::stringstream outputFileName;
        std::string s0;
    
        if (params.writeVtkCells) {
            s0 = this->params.outputDirectory + this->params.testcase + "_cells.";
            outputFileName << s0  << reloj << ".vtk";
            this->writeVtk(outputFileName.str());      
        }
    
        if (params.writeVtkFibres) {
            s0 = this->params.outputDirectory + this->params.testcase + "_fibres.";
            std::stringstream fibres_outputFileName;
            fibres_outputFileName << s0  << reloj << ".vtk";
            this->writeFibresVtk(fibres_outputFileName.str());
        }

        if (params.writeVtkVessels) {
            s0 = this->params.outputDirectory + this->params.testcase + "_vessels.";
            std::stringstream vessels_outputFileName;
            vessels_outputFileName << s0  << reloj << ".vtk";
            this->writeVesselsVtk(vessels_outputFileName.str());
        }

        //if (params.alpha_birthrate[0]==0){
        if (reloj == params.n_steps){
            ofstream cellpositionFile;
            //string final_pos = this->params.casedirectory + this->params.casename + "_cell_positions.txt";
            string final_pos = this->params.casedirectory + "cell_positions.txt";
            cellpositionFile.open(final_pos.c_str(),ios::app);
            for(int k=this->minx; k<=this->maxx; k++) {
	            for(int l=this->miny; l<=this->maxy; l++){
	                for(int n=this->minz; n<=this->maxz; n++) {
	                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	                        cellpositionFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
			                << this->boxes_A[k][l][n].cells[i].position[1] << " "
			                << this->boxes_A[k][l][n].cells[i].position[2] << " " << endl;
	                    }
                    }
                }
            }
            cellpositionFile.close();

            ofstream cellphenotypeFile;
            string phenotypes = this->params.casedirectory + "cell_phenotypes.txt";
            cellphenotypeFile.open(phenotypes.c_str(),ios::app);
            for(int k=this->minx; k<=this->maxx; k++) {
                for(int l=this->miny; l<=this->maxy; l++){
                    for(int n=this->minz; n<=this->maxz; n++) {
                        for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                            cellphenotypeFile << this->boxes_A[k][l][n].cells[i].cont_pheno << endl;
                        }
                    }
                }
            }
            cellphenotypeFile.close();
        }
  } // end of main loop while

  cout << " count per type " << endl;
  this->count_cells_per_type();
  //cout << " number of removed cells: " << total_no_of_removed_cells << endl;

  writeParameterList();        

  if (this->params.verbose>1){
    // ============================================
    // write mean position to data file
    // ============================================
    ofstream meandist;
    string meandist_list =  this->params.casedirectory + this->params.casename + "_mean_distance.txt";
    meandist.open(meandist_list.c_str(),ios::app);
    meandist << endl;
    meandist.close();
    // ********************************

    // ============================================
    // write variance of above to data file
    // ============================================
    ofstream vardist;
    string vardist_list =  this->params.casedirectory + this->params.casename + "_variance_distance.txt";
    vardist.open(vardist_list.c_str(),ios::app);
    vardist << endl;
    vardist.close();
    // ********************************
  }

  // write statistics
  if (params.writeStatistics) {
    cout << " write Stats " << endl;
    ofstream cellStatsFile;
    string statfilename = this->params.testcase + "_stats.txt";
    // max - min of occupied boxes
    if (simulation_id>-1) {
      cellStatsFile.open(statfilename.c_str(),ios::app);
    } else {
      cellStatsFile.open(statfilename.c_str());
    }
    cellStatsFile << this->total_no_of_cells << " "
		  << this->reloj << " "
		  << this->maxx-this->minx << " " << this->maxy-this->miny << " "
		  << this->maxz-this->minz << endl;
    cellStatsFile.close(); 
  }
 
  // ============================================
  // intermediary cell numbers end line
  // ============================================
  /*
    ofstream pheno1File;
    string full_list = this->params.casedirectory + this->params.casename + "_full.txt";
    pheno1File.open(full_list.c_str(),ios::app);
    pheno1File << endl;
    pheno1File.close();
    if (params.ic_phenotype.size()==2){
      ofstream pheno2File;
      pheno2File.open("phenotype2_cellnumbers.txt",ios::app);
      pheno2File << endl;
      pheno2File.close();
    }
  */

  // ======================================
  // write final cell numbers to data file
  // ======================================
  /*
     ofstream cellnumbersFile;
     string final_list = this->params.casedirectory + this->params.casename + "_final_timestep.txt";
     cellnumbersFile.open(final_list.c_str(),ios::app);
     cellnumbersFile << this->total_no_of_cells << " "
     << this->reloj << " "
     << this->maxx-this->minx << " " << this->maxy-this->miny << " " << this->maxz-this->minz << " "
     << this->phenotype1_count << " " << this->phenotype2_count << " "
     << (this->maxx-this->minx)*(this->maxy - this->miny)*(this->maxz - this->minz) << endl;
     cellnumbersFile.close();
  */

  // ======================================
  // write single cell positions to file
  // NOTE: This is not by the input flag writeCellList
  // ======================================
  /*
    if (params.alpha_birthrate[0]==0){
        ofstream cellpositionFile;
        string final_pos = this->params.casedirectory + this->params.casename + "_cell_positions.txt";
        cellpositionFile.open(final_pos.c_str(),ios::app);
        for(int k=this->minx; k<=this->maxx; k++) {
            for(int l=this->miny; l<=this->maxy; l++){
                for(int n=this->minz; n<=this->maxz; n++) {
                    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
                        cellpositionFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
                        << this->boxes_A[k][l][n].cells[i].position[1] << " "
                        << this->boxes_A[k][l][n].cells[i].position[2] << " " << endl;
                    }
                }
            }
        }
        cellpositionFile.close();
    }
    -params.ic_cell_x[i]-params.ic_cell_y[i]-params.ic_cell_z[i];
  */ 
}

/* ****************************************************************************
   count cells per box
   **************************************************************************** */
void CoupledModel::count_cells_per_box()
{
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    if (this->boxes_A[k][l][n].cells.size()) {
	        if (params.verbose>1) {
	            cout << " box [" << k << "," << l << "," << n << "] has " 
		        << this->boxes_A[k][l][n].cells.size() << " cells " << endl;
	        }
	    }
      }
    }
  } 
}

/* ****************************************************************************
   count cells per type
   **************************************************************************** */
void CoupledModel::count_cells_per_type()
{
  unsigned int countNorm = 0;
  unsigned int countHypo = 0;
  unsigned int countDead = 0;
  unsigned int countPhenotypeOxygen_0_10 = 0;
  unsigned int countPhenotypeOxygen_0_12 = 0;
  unsigned int countPhenotypeOxygen_14_24 = 0;
  unsigned int countPhenotypeOxygen_12_24 = 0;
  
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        unsigned int cellType = this->boxes_A[k][l][n].cells[i].type;
	        countNorm = countNorm + (cellType==1);
	        countHypo = countHypo + (cellType==2);
	        //countDead = countDead + (cellType==3);
	        countDead = total_no_of_removed_cells; // TOMMASO
	        countDead = total_no_of_removed_cells; // TOMMASO
	        if (this->boxes_A[k][l][n].cells[i].phenotype < 10) {
	            countPhenotypeOxygen_0_10 ++;
	        } else if ( (this->boxes_A[k][l][n].cells[i].phenotype > 14) &&
		        (this->boxes_A[k][l][n].cells[i].phenotype < 24) ){
	            countPhenotypeOxygen_14_24 ++;
	        }
	  
	        if (this->boxes_A[k][l][n].cells[i].phenotype < 12) {
	            countPhenotypeOxygen_0_12 ++;
	        } else if ( (this->boxes_A[k][l][n].cells[i].phenotype > 12) &&
		        (this->boxes_A[k][l][n].cells[i].phenotype < 24) ){
	            countPhenotypeOxygen_12_24 ++;
            }
	    }
      }
    }
  }
  this->totNorm.push_back(countNorm);
  this->totHypo.push_back(countHypo);
  this->totDead.push_back(countDead);
  
  ofstream nCellFile;
  nCellFile.open("cell_counter.txt",ios::app);
  if (!nCellFile) {
    cerr << " *** ERROR *** could not open file *** " << endl;
    exit(1);
  }
  
  unsigned int nt = this->totNorm.size()-1;
  nCellFile << this->totNorm[nt] << " "
  << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  cout << " cells (norm,hypo,dead): " 
  << this->totNorm[nt] << " " << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  if ( (this->totNorm[nt] + this->totHypo[nt] + this->totDead[nt]) < this->total_no_of_cells ) {
     cout << " ERROR in CoupledModel::count_cells_per_type(): not all cells found ! " << endl;
     this->end();
     exit(1);
  }
}
  
/* ****************************************************************************
   end
   **************************************************************************** */
void CoupledModel::end()
{

  if (params.getGenealogy) {
    this->gather_all_births();
  }
  // write boxes, fibres and vessels on file (if not written before)
  std::string s0;
  if (params.writeVtkFibres==0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_fibres.";
    std::stringstream fibres_outputFileName;
    fibres_outputFileName << s0  << reloj << ".vtk";
    this->writeFibresVtk(fibres_outputFileName.str());
  }

  if (params.writeVtkVessels==0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_vessels.";
    std::stringstream vessels_outputFileName;
    vessels_outputFileName << s0  << reloj << ".vtk";
    this->writeVesselsVtk(vessels_outputFileName.str());
  }

  if (params.writeVtkCells==0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_cells.";
    std::stringstream cells_outputFileName;
    cells_outputFileName << s0  << reloj << ".vtk";
    this->writeVtk(cells_outputFileName.str());
  }
  
  s0 = this->params.outputDirectory + this->params.testcase + "_boxes.vtk";
  this->writeBoxesVtk(s0);

  if (params.writeCellList) {
    ofstream outCellFile;
    outCellFile.open(this->params.fileCells.c_str(),ios::out);
    if (!outCellFile) {
      cerr << " *** ERROR, file " << __FILE__ << ", line " << __LINE__ 
	   << " *** could not open file " << this->params.fileCells2FEM << endl;
      exit(1);
    }
	
    // header: length, iteration (pde), iteration (total)
    cout << " writing on " << this->params.fileCells.c_str() << endl;
    outCellFile << this->total_no_of_cells << " "
		<< this->oxy_diff.iter << " "
		<< this->reloj << " "
		<< endl;
    // body: x,y,z,type,r,phenotype,adhesion coeff,id,energy
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	            outCellFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
			    << this->boxes_A[k][l][n].cells[i].position[1] << " "
			    << this->boxes_A[k][l][n].cells[i].position[2] << " "
			    << this->boxes_A[k][l][n].cells[i].type << " "
			    << this->boxes_A[k][l][n].cells[i].radius;
	            if (params.writeCellList>1) {
		            outCellFile << " "
		            //<< this->boxes_A[k][l][n].cells[i].phenotype << " " //REMOVED 25/6/19 TOMMASO
			        << this->boxes_A[k][l][n].cells[i].cont_pheno << " " //ADDED 25/6/19 TOMMASO
			        << this->boxes_A[k][l][n].cells[i].adhesion << " "
			        << this->boxes_A[k][l][n].cells[i].name << " "
			        << this->boxes_A[k][l][n].cells[i].energy;
	            }
	            outCellFile << endl;
	        }
	    }
      }
    }
    outCellFile.close();
  }   
}

/* ****************************************************************************
   WRITE CELLS TO VTK
   **************************************************************************** */
void CoupledModel::writeVtk(string filename,unsigned int onlyCoord)
{
  ofstream outfile(filename.c_str());

  // rename (locally) the total number of cells
  unsigned int nCells = this->total_no_of_cells;
  
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Unstructured grid legacy vtk file with point scalar data" << endl;
  outfile << "ASCII\n\n";
  outfile << "DATASET UNSTRUCTURED_GRID\n";
  outfile << "POINTS " << nCells << " double\n";
  
  // write positions: we loop on all boxes
  //@warning the 2-dimensional output is not supported
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        outfile << this->boxes_A[k][l][n].cells[i].position[0] << " "
		    << this->boxes_A[k][l][n].cells[i].position[1] << " "
		    << this->boxes_A[k][l][n].cells[i].position[2] << endl;
	    }
      }
    }
  }
  outfile << endl;

  outfile << "POINT_DATA " << nCells << endl;
  // write radii
  outfile << "SCALARS radius double" << endl;
  outfile << "LOOKUP_TABLE default" << endl;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        outfile << this->boxes_A[k][l][n].cells[i].radius << endl;
	    }
      }
    }
  }
  outfile << endl;

  if (onlyCoord==0) {
    // write cell type
    outfile << "SCALARS status double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].type << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    // write oxygen concentration
    outfile << "SCALARS concentration double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].O2 << endl;
	        }
	    }
      }
    }
    outfile << endl;

    //write fibre contacts
    outfile << "SCALARS fibres double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].contact_fibres.size() << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    //write phenoype
    outfile << "SCALARS phenot double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].phenotype << endl;
	        }
	    }
      }
    }
    outfile << endl;

    //write phenoype //ADDED 25/6/19 TOMMASO
    outfile << "SCALARS cont_pheno double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].cont_pheno << endl;
	        }
	    }
      }
    }
    outfile << endl;
    
    //write adhesion
    outfile << "SCALARS name double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].name << endl;
	        }
	    }
      }
    }
    outfile << endl;

     //write interaction phenotype
    outfile << "SCALARS interaction_phenotype double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].interaction_phenotype << endl;
	        }
	    }
      }
    }
    outfile << endl;

    //write vessel interaction
    outfile << "SCALARS vessel_interaction double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].vessel_interaction << endl;
	        }
	    }
      }
    }
    outfile << endl;

    //write interaction phenotype
    outfile << "SCALARS is_follower double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].is_follower << endl;
	        }
	    }
      }
    }
    outfile << endl;

     //write polarised - CICELY NEW
    outfile << "SCALARS polarised double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	        for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	            outfile << this->boxes_A[k][l][n].cells[i].polarised << endl;
	        }
	    }
      }
    }
    outfile << endl;  
  }
}

/* ****************************************************************************
   WRITE BOXES TO VTK
   **************************************************************************** */
void CoupledModel::writeBoxesVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  // header
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Boxes output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  // nodes
  int nNodesX = params.boxesx+1;
  int nNodesY = params.boxesy+1;
  int nNodesZ = params.boxesz+1;
  int nNodes = nNodesX*nNodesY*nNodesZ;
  ofile << "POINTS " << nNodes << " float" << endl;
  for (int k=0; k<nNodesZ; k++){
    for (int j=0; j<nNodesY; j++){
      for (int i=0; i<nNodesX; i++){
	    ofile << i*this->box_sizex << " " << j*this->box_sizey << " "
	    <<  k*this->box_sizez << endl;
      }
    }  
  }
  // cells
  int nCells = (params.boxesx)*(params.boxesy)*(params.boxesz);
  ofile << "CELLS " << nCells << " " << 9*nCells << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    // coordinates of box_(i,j,k)
	    // the points must be ordered according to vtk data structures
	    // see e.g. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
	    ofile << 8 << " " <<  i + j*nNodesX + k*nNodesX*nNodesY << " "
	    << (i+1) + j*nNodesX + k*nNodesX*nNodesY << " "
	    << i + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	    << (i+1) + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	    <<  i + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << (i+1) + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << i + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY << " "
	    << (i+1) + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY
	    << endl;
      }
    }  
  }
  // nodes
  ofile << "CELL_TYPES " << nCells  << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    ofile << 11 << endl;
      }
    }  
  }
  ofile << "CELL_DATA  " << nCells << endl;
  ofile << "SCALARS nCells float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	    ofile << this->boxes_A[i][j][k].cells.size()<< endl;
      }
    }
  }
  ofile.close();
}

/* ****************************************************************************
   WRITE FIBRES TO VTK
   **************************************************************************** */
void CoupledModel::writeFibresVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  // file version and identifier in the form "# vtk DataFile Version x.x"
  ofile << "# vtk DataFile Version 3.0" << endl;
  // header used to describe the data
  ofile << "Fibres output - vtk " << endl;
  // file format either "ASCII" or "BINARY"
  ofile << "ASCII" << endl;
  // dataset structure in the form "DATASET_type"
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  // dataset attributes using "POINTS " or "CELLS " followed by an integer number specifying the number of points or cells.
  ofile << "POINTS " << 2*fibres.size() << " double" << endl; // Each fibre has two points (start and end)
  for(unsigned int ff=0; ff < fibres.size(); ff++){
    ofile << this->fibres[ff].start[0] << " " << this->fibres[ff].start[1] << " " << this->fibres[ff].start[2] << endl;
    ofile << this->fibres[ff].start[0]+fibres[ff].length*fibres[ff].direction[0] << " " << this->fibres[ff].start[1]+fibres[ff].length*fibres[ff].direction[1] << " " << this->fibres[ff].start[2]+fibres[ff].length*fibres[ff].direction[2] << endl;
  }
  ofile << "CELLS " << fibres.size() << " " << 3*fibres.size() << endl; // Each Fibre element has to detail its type and its points
  // 	// coordinates of each fibre
  // 	// the points must be ordered according to vtk data structures
  // 	// see e.g. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  for (unsigned int ff=0; ff<fibres.size(); ff++){
    ofile << 2 << " " << 2*ff << " " << 2*ff + 1 << endl;
  }
  ofile << "CELL_TYPES " << fibres.size() << endl;
  for(unsigned int ff=0; ff < fibres.size(); ff++){
    ofile << 3 << endl;
  }
  ofile << "CELL_DATA  " << fibres.size() << endl;
  ofile << "SCALARS length float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
    ofile << fibres[ff].length << endl;
  }
  ofile << "SCALARS fibre_name double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
   ofile << fibres[ff].fibre_name << endl;
  }
  ofile << "SCALARS fibre_crosslinks double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
      ofile << fibres[ff].fibre_crosslinks << endl;
  }
  ofile << "SCALARS fibre_exists double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
   ofile << fibres[ff].fibre_exists << endl;
  }
  ofile << "VECTORS orientation double" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
   ofile << fibres[ff].direction[0] << " " << fibres[ff].direction[1] << " " << fibres[ff].direction[2] << endl;
  }
  ofile << "SCALARS fibre_radius double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++){
    ofile << fibres[ff].fradius << endl;
  }
  ofile.close();
}

/* ****************************************************************************
   WRITE VESSELS TO VTK
   **************************************************************************** */
void CoupledModel::writeVesselsVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Vessels output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  ofile << "POINTS " << 2*vessels.size() << " double" << endl; // Each vessel has two points (start and end)
  for(unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << this->vessels[ff].ves_start[0] << " " << this->vessels[ff].ves_start[1] << " " << this->vessels[ff].ves_start[2] << endl;
    ofile << this->vessels[ff].ves_start[0]+vessels[ff].ves_length*vessels[ff].ves_direction[0] << " " << this->vessels[ff].ves_start[1]+vessels[ff].ves_length*vessels[ff].ves_direction[1] << " " << this->vessels[ff].ves_start[2]+vessels[ff].ves_length*vessels[ff].ves_direction[2] << endl;
  }
  ofile << "CELLS " << vessels.size() << " " << 3*vessels.size() << endl; // Each Vessel element has to detail its type and its points
  for (unsigned int ff=0; ff<vessels.size(); ff++){
    ofile << 2 << " " << 2*ff << " " << 2*ff + 1 << endl;
  }
  ofile << "CELL_TYPES " << vessels.size() << endl;
  for(unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << 3 << endl;
  }
  ofile << "CELL_DATA  " << vessels.size() << endl;
  ofile << "SCALARS length float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << vessels[ff].ves_length << endl;
  }
  ofile << "SCALARS vessel_name double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
   ofile << vessels[ff].vessel_name << endl;
  }
  ofile << "SCALARS vessel_radius double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++){
    ofile << vessels[ff].ves_radius << endl;
  }
  ofile.close();
}

/* ****************************************************************************
   FIBRE DEGRADATION
   **************************************************************************** */
/* @brief
   For fibres located within a cell diameter away from cell
   -- degrade due to contact based on a probability of 0.5
            (cell must be moving, not located at either end point of the fibre
            and the direction of movement must be directed towards the fibre)
   -- degrade due to diffusion based on a probability of 0.1
*/  
void CoupledModel::fibre_degradation(Cell& cell,
                        Fibre& fibre, double cell_fibre_min_dist)
{
  /*
  ofstream fibint;
  string fibint_list =  this->params.casedirectory + this->params.casename + "_fibre_interaction_info.txt";
  fibint.open(fibint_list.c_str(),ios::app);
  */
  
  double cf1[3], cf2[3];
  double n[3], ncf1[3], ncf2[3];
  double vncoeff = 0.0;
  double vfcoeff = 0.0;
  double v_n[3], v_p[3], v_f[3];
  double test1 = 0.0;
  double test2 = 0.0;

  double prob_degradation_cell = params.prob_cell_deg;
  double prob_degradation_diffuse = params.prob_diff_deg;
  
  for (unsigned int l=0; l<3; l++) {
    cf1[l] = fibre.start[l]-cell.position[l];  // start->cell
    cf2[l] = fibre.start[l]+fibre.length*fibre.direction[l]-cell.position[l]; //end->cell
  }

  // cross product of cf1 and cf2 giving the normal to the plane containing these two vectors
  n[0] = cf1[1]*cf2[2]-cf1[2]*cf2[1];
  n[1] = cf1[2]*cf2[0]-cf1[0]*cf2[2];
  n[2] = cf1[0]*cf2[1]-cf1[1]*cf2[0];
  // magnitude of normal such that n/nnorm is the unit normal
  double nnorm = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

  // dot product of cell polarity and unit normal
  for (unsigned int l=0; l<3; l++) {
    vncoeff += cell.polarity[l]*n[l]/nnorm;
  }

  // projections of the cell velocity onto the normal and in the plane
  for (unsigned int l=0; l<3; l++) {
    v_n[l] = vncoeff*n[l]/nnorm;
    v_p[l] = cell.polarity[l] - v_n[l];
    vfcoeff += v_p[l]*fibre.direction[l];
  }

  // magnitude of normal cell velocity component
  double vnnorm = sqrt(v_n[0]*v_n[0]+v_n[1]*v_n[1]+v_n[2]*v_n[2]);

  // projection of the cell velocity towards the fibre
  for (unsigned int l=0; l<3; l++) {
    v_f[l] = v_p[l] - vfcoeff*fibre.direction[l];
  }

  // magnitude of cell velocity component directed towards the fibre
  double vfnorm = sqrt(v_f[0]*v_f[0]+v_f[1]*v_f[1]+v_f[2]*v_f[2]);

  // cross product of cf1 and cell velocity in plane
  ncf1[0] = cf1[1]*v_p[2]-cf1[2]*v_p[1];
  ncf1[1] = cf1[2]*v_p[0]-cf1[0]*v_p[2];
  ncf1[2] = cf1[0]*v_p[1]-cf1[1]*v_p[0];
  // cross product of cf2 and cell velocity in plane
  ncf2[0] = cf2[1]*v_p[2]-cf2[2]*v_p[1];
  ncf2[1] = cf2[2]*v_p[0]-cf2[0]*v_p[2];
  ncf2[2] = cf2[0]*v_p[1]-cf2[1]*v_p[0];

  for (unsigned int l=0; l<3; l++) {
    test1 += ncf1[l]*n[l];
    test2 += ncf2[l]*-n[l];
  }

  // check if velocity of cell is such that cell hits the fibre
  if (test1>0.0 || test2>0.0) {

    double alpha = (cell.radius+fibre.fradius)/cell_fibre_min_dist;
    
    if (abs(vnnorm)<alpha*abs(vfnorm)) {

      // fibint << "   ** cell " << cell.name << " is moving towards fibre " << fibre.fibre_name << " ** " << endl;

      double rand_contact_degrade = aleatorio();
      double prob_degradation = prob_degradation_cell*(1. - abs(vnnorm)/(alpha*abs(vfnorm)+1e-8));
      if (rand_contact_degrade<prob_degradation){

	        //fibint << "    rand_contact_degrade is " << rand_contact_degrade
	        //	       << " -> decision to degrade by contact" << endl;   
	        fibre.fibre_exists=0.0;
	
	        if (params.verbose>0) {
	            cout << "at time " << reloj << " fibre " << fibre.fibre_name
	            << " eliminated due to contact with cell " << cell.name << endl;
	            cout << "cell: " << cell.position[0] << " " << cell.position[1] << " " << cell.position[2] << endl;
	            cout << "polarty: " << cell.polarity[0] << " " << cell.polarity[1] << " " << cell.polarity[2] << endl;
	            cout << "fibre: " << fibre.start[0] << " " << fibre.start[1] << " " << fibre.start[2] << endl;
	            cout << "       " << fibre.start[0]+fibre.length*fibre.direction[0] << " "
	            << fibre.start[1]+fibre.length*fibre.direction[1] << " "
	            << fibre.start[2]+fibre.length*fibre.direction[2] << endl;
	        }
      }
    }
  }
    
  if (fibre.fibre_exists!=0.0){
    double rand_diffuse_degrade = aleatorio();
    
    if (rand_diffuse_degrade<prob_degradation_diffuse){

      //fibint << "    rand_diffuse_degrade is " << rand_diffuse_degrade << " -> decision to degrade by diffusion " << endl;
      //fibint.close();
      
      fibre.fibre_exists=0.5;
      if (params.verbose>0) {
        cout << "at time " << reloj << " fibre " << fibre.fibre_name << " eliminated by diffusion from cell " << cell.name << endl;
      }
    }
  } 
}		    

/* *******************************************************************************
   get sub domain
   ***************************************************************************** */
int CoupledModel::get_sub_domain_id(double x,double y,double z)
{
  int sub_domain_model_type = 0;
  
  switch (sub_domain_model_type) {
    case 0:
    {  
      if (x <= params.lattice_length_x/2.0) 
	    return 0;
      else
	    return 1;
      break;
    }
    default:
    {
      cout << " ** ERROR in CoupledModel::get_sub_domain_id: sub_domain_model_type = "
	   << sub_domain_model_type << " not implemented yet" << endl;
      exit(1);
      return -1;
      break;
    }
    return -1;
  }
}

/* ************************************************************************
   collect children            
   ************************************************************************ */
void CoupledModel::collect_all_children(unsigned int mother_id,
					std::vector<int>& children,
					std::vector<int>& children_times)
{
  children.resize(0);
  children_times.resize(0);

  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	        if (boxes_A[k][l][n].cells[i].name != mother_id) {
	            if (boxes_A[k][l][n].cells[i].mother_name == mother_id) {
	                children.push_back(boxes_A[k][l][n].cells[i].name);
	                children_times.push_back(boxes_A[k][l][n].cells[i].birthday);
	                cout << mother_id << " --> " << mother_id << ","
		            << boxes_A[k][l][n].cells[i].name << " at time "
		            << boxes_A[k][l][n].cells[i].birthday << endl;
	            }
	        }
	    }
      }
    }
  }
}

/* ************************************************************************
   gather births             
   ************************************************************************ */
void CoupledModel::gather_all_births()
{
  cout << " find births" << endl;
  std::vector<int> children,children_2,children_times;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {	  
	        //if (this->boxes_A[k][l][n].cells[i].birthday == 0) {
	        if (this->boxes_A[k][l][n].cells[i].name == 0) {
	            cout << " ... for cell " << this->boxes_A[k][l][n].cells[i].name
		        << " (birthday: " << this->boxes_A[k][l][n].cells[i].birthday << ")" << endl;
	            collect_all_children(this->boxes_A[k][l][n].cells[i].name,children,children_times);
	            while (children.size()>0) {
	                vector<int> save_children(children);
	                for (unsigned int l=0; l<children.size(); l++) {
		                save_children[l] = children[l];
	                }
	                for (unsigned int l=0; l<children.size(); l++) {
		                collect_all_children(save_children[l],children,children_times);
	                }
	            }
	        }
	    }
      }
    }
  }
} 

/* ************************************************************************
   parameter list              
   ************************************************************************ */
void CoupledModel::writeParameterList()
{
  int sum_phenotype = 0.0;
  unsigned int nt = this->totNorm.size()-1;
  for (int i=0; i<params.n_initial_cells; i++) {
    sum_phenotype += params.ic_phenotype[i];
  }
  
  ofstream parameters;
  string paramlist = this->params.casedirectory + this->params.casename + "_parameters.txt";
  parameters.open(paramlist.c_str(),ios::out);
  parameters << "The parameters for the case **" << this->params.casename.c_str() << "** are:" << endl;
  parameters << endl;
  parameters << "There are " << params.n_initial_cells << " cell(s) at the start of the simulation: " << endl;
  // << params.n_initial_cells-sum_phenotype << " of phenotype 0 and " << sum_phenotype << " of phenotype 1" << endl; TOMMASO
  parameters << "After " << reloj << " timesteps there are " << this->total_no_of_cells << " cell(s) at the end of the simulation " << endl;
  parameters << "(norm,hypo,dead): " << this->totNorm[nt] << " " << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  //<< this->phenotype1_count << " of phenotype 0 and " << this->phenotype2_count << " of phenotype 1" << endl; TOMMASO
  parameters << endl;
  parameters << "[cells]" << endl;
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "growth rate = " << params.growth_rate[i] << " ";
  }
  parameters << endl;
  /* TOMMASO for now just commented out
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "birth rate = " << params.alpha_birthrate[i] << " ";
  }
  parameters << endl;*/ 
  parameters << "normoxic birth rate = " << params.normoxic_birth << endl;
  parameters << "hypoxic birth rate = " << params.hypoxic_birth << endl;
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "Youngs Modulus = " << params.alpha_YoungM[i] << " ";
  }
  parameters << endl;
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "Poisson number = " << params.alpha_PoissonNo[i] << " ";
  }
  parameters << endl;
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "GCM = " << params.alpha_gcm[i]*params.Gcm << " ";
  }
  parameters << endl;
  for (int i=0; i<params.n_phenotypes; i++) {
    parameters << "adhesion value = " << params.adhesion_value[i] << " ";
  }
  parameters << endl;
  parameters << "random motion = " << params.variance_motion << endl;
  parameters << "birth energy = " << this->birth_energy << endl;
  parameters << endl;
  parameters << "[vessels]" << endl;
  if (params.n_initial_vessels==0){
    parameters << "There are no vessels in this simulation" << endl;
  }
  if (params.n_initial_vessels!=0){
    parameters << "Youngs Modulus = " << params.vessel_YoungM << endl <<
      "Poisson number = " << params.vessel_PoissonNo << endl <<
      "adhesion value = " << params.vessel_adhesion << endl;
  }
  parameters << endl;
  parameters << "[mutations]" << endl;
  parameters << "mutation_amount = " << params.mutation_amount << endl;
  parameters << "mutation_probability = " << params.mutation_probability << endl;  
  parameters << endl;
  parameters << "[oxygen]" << endl;
  parameters << "threshold_hypo = " << params.threshold_hypo << endl;
  parameters << "threshold_death = " << params.threshold_death << endl;  
  parameters << endl;
  parameters << "[fibres]" << endl;
  if (params.n_initial_fibres[0]==0 && params.n_initial_fibres[1]==0){
    parameters << "There are no fibres in this simulation" << endl;
  }
  if (params.n_initial_fibres[0]!=0 || params.n_initial_fibres[0]!=0){
      parameters << "n_initial_fibres = " << params.n_initial_fibres[0] << " " << params.n_initial_fibres[1] << endl <<
	  "fibre_orientation_distribution = " << params.fibre_orientation_distribution[0] << " " << params.fibre_orientation_distribution[1] << endl <<
	  "fibre_orientation_mean_phi = " << params.fibre_orientation_mean_phi[0] << " " << params.fibre_orientation_mean_phi[1] << endl <<
	  "fibre_orientation_variance_phi = " << params.fibre_orientation_variance_phi[0] << " " << params.fibre_orientation_variance_phi[1] << endl <<
	  "fibre_orientation_mean_theta = " << params.fibre_orientation_mean_theta[0] << " " << params.fibre_orientation_mean_theta[1] << endl <<
	  "fibre_orientation_variance_theta = " << params.fibre_orientation_variance_theta[0] << " " << params.fibre_orientation_variance_theta[1] << endl <<
	  "fibre_length_mean = " << params.fibre_length_mean[0] << " " << params.fibre_length_mean[1] << endl <<
	  "fibre_length_variance = " << params.fibre_length_variance[0] << " " << params.fibre_length_variance[1] << endl <<
	  "fibre_radius = " << params.fibre_radius << endl <<
	  "vel_adhesion = " << params.vel_adhesion << endl <<
	  "vel_contact = " << params.vel_contact << endl;
      if (params.fib_deg==1){
	    parameters << "Fibre degradation switched on" << endl;
      }
      if (params.fib_deg==0){
	    parameters << "Fibre degradation switched off" << endl;
      }
  }
  parameters.close();
}