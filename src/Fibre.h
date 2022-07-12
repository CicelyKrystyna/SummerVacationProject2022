#ifndef _FIBRE_H_
#define _FIBRE_H_

#include <vector>

/** ************************************************************************ 
*
* @class      Fibre
* @brief      Store and evolve single fibre 
* 
* 
* @author     Ignacio Ramis & Alfonso Caiazzo & Cicely Macnamara
* @date       01.08.2016
*
****************************************************************************/
class Fibre {

  public:
  Fibre();
  ~Fibre(){};
  
  ///@brief ID in the whole list
  unsigned int fibre_name;

  ///@brief current box
  //int fibre_box[3]; 
  ///@brief box after position update
  //int fibre_new_box[3];

  ///@brief start position
  float start[3];
  ///@brief end position
  float end[3];
  ///@brief elevation angle
  double theta;
  ///@brief azimuth angle
  double phi;
  ///@brief length
  double length;
  ///@brief direction
  float direction[3];
  ///@brief fibre radius
  double fradius;

  ///@brief fibre existance
  double fibre_exists;

  ///@brief fibre crosslinks
  double fibre_crosslinks;

  ///@brief adhesion forces constant
  //double fibre_adhesion;
  
  ///@brief n. of contacts
  // unsigned int fibre_contacts;
  
  /// @brief store the neighbors (in contact)
  //std::vector<Fibre*> fibre_neighbors;

  // methods
  void fibre_clear_contacts();

  /// @brief print fibre information
  void fibre_printInfo() const;


};


#endif
