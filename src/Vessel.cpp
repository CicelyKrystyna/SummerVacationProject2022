#include "Vessel.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

Vessel::Vessel()
{
  
}

void Vessel::vessel_printInfo() const 
{
  
  cout << " *** print vessel " << vessel_name << " info *** " << endl;
  cout << " - start position: " << ves_start[0] << "," 
        << ves_start[1] << "," << ves_start[2] << endl;
  cout << " - length: " << ves_length << endl;
  cout << " - direction: " << ves_direction[0] << "," 
        << ves_direction[1] << "," << ves_direction[2] << endl;
  cout << " - radius: " << ves_radius << endl;
  cout << " ***********************\n " << endl;
  
}

void Vessel::vessel_clear_contacts()
{
  
}

