#ifndef _PDE_H_
#define _PDE_H_

#include <string.h>
#include "Param.h"
#include "Mesh.h"

/** ************************************************************************ 
*
* @class      PDE
* @brief      class for PDE solver (continuous) - FreeFem based
* 
* 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       12.2014
*
****************************************************************************/
 

class PDE
{
 public:
  PDE(){};
  ~PDE(){};

  /// @brief the real time of simulation
  double time; 
  
  /// @brief how many times we call the PDE solver
  unsigned int iter; 

  /// @brief  how to call the solver (with system(FreeFemCall))
  std::string FreeFemCall; 

  /// @brief whether we run the solver or not
  bool launch; 

  Mesh mesh;

  // methods

  /// @brief initialize class
  void init(const Param& _params, const Mesh& _m); 
  void init();

  ///@brief solve the PDE
  void solve(double t); 



};

#endif
