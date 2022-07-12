#ifndef _VESSEL_H_
#define _VESSEL_H_

#include <vector>

class Vessel {

  public:
  Vessel();
  ~Vessel(){};
  
  ///@brief ID in the whole list
  unsigned int vessel_name;

  ///@brief start position
  float ves_start[3];
  ///@brief length
  double ves_length;
  ///@brief direction
  float ves_direction[3];
  ///@brief fibre radius
  double ves_radius;

  /// @brief store the neighbors (in contact)
  //std::vector<Vessel*> vessel_neighbors;

  // methods
  void vessel_clear_contacts();

  /// @brief print fibre information
  void vessel_printInfo() const;

};


#endif
