#ifndef _MESH_H
#define _MESH_H 1

#include <stdlib.h>
#include <vector>
#include <string>

/** ************************************************************************ 
*
* @class      Mesh
* @brief      Store a mesh geometry (tetra in 3D)
* 
* 
* @author     Ignacio Ramis & Alfonso Caiazzo 
* @date       06.06.2014
*
****************************************************************************/
class Mesh
{

 public:

  int nNodes;
  int dim;
  std::string meshFile;
  std::vector<double> xp,yp,zp;
  std::vector<int> node_ind;

  unsigned int nEdges,nTria,nTetra,nElem;
  std::vector<int> edge;
  std::vector<int> tria;
  std::vector<int> quad;
  std::vector<int> tetra;
  std::vector<int> hexa;
  std::vector<int> edgeRef,triaRef,tetraRef,pRef;
  std::vector<int> cellsInTria;
  std::vector<int> cellsInTriaOld;
  std::vector<int> cellsInTriaNorm;
  std::vector<int> cellsInTriaHypo;
  std::vector<int> cellsInTriaDead;
  std::vector<int> elementVolume;

  ///@brief empty constructor
  Mesh(){dim=0;nElem=0;};

  /// @brief read Mesh from a file (.mesh)
  void read(std::string _f); 

  void computeVolumes();
  
  /// @brief write mesh (.mesh)
  void write(std::string _f); 

  /// @brief allocate and clear fem-cell coupling arrays
  /// @todo move into another function
  void initCellsArrays();

  /// @brief evaluate the density change
  void storeCellDensity();
  double cellDensityChange();
  
  /// @brief print info about the mesh
  void info();

};


#endif


