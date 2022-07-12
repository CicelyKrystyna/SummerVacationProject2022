/* 
read and write Utilities for ensight cases   
*/

#include <math.h>
#include <iostream>
#include <fstream>
#include "Mesh.h"

using namespace std;

void stripSpace(string &str) 
{
  for (size_t i=0;i<str.length();i++)
    {
      if (str[i]==' ') 
	{
	  str.erase(i,1);
	  i--;
	}
    }
}

void Mesh::read(string _file)
{
  
  cout << " -- Mesh:: reading mesh from : " << _file << endl;

  ifstream input_file;
  string _line;
  int current_line_number;

  meshFile = _file;
  input_file.open(_file.c_str(),ios::in);
  if (!input_file) {
    cerr << " *** ERROR *** could not open file " << _file << endl;
    exit(1);
  }
  // read header
  do {
    getline(input_file,_line,'\n');
    stripSpace(_line); 
  } while (_line.find("Dimension")==string::npos);
  _line.erase(0,9);

  if (_line.length()>0) {
    sscanf(_line.c_str(),"%d",&dim);
  }
  else {
    input_file >> dim; 
  }

  nNodes = 0;
  do {
    getline(input_file,_line,'\n');
    stripSpace(_line); // borra espacios en blanco
  } while (_line.find("Vertices")==string::npos);
  _line.erase(0,8);
  
  if (_line.length()>0) {
    sscanf(_line.c_str(),"%d",&nNodes);
  }
  input_file >> nNodes;
  
  if (nNodes) {
    // read vertices
    node_ind.resize(nNodes);
    xp.resize(nNodes);
    yp.resize(nNodes);
    pRef.resize(nNodes);
    if (dim==3) {
      zp.resize(nNodes);
    }
  } else {
    cerr << " *** ERROR: mesh has no nodes! " << endl;
  }
  current_line_number = input_file.tellg(); 
  
  for (int i=0;i<nNodes;i++) {
    input_file >> xp[i] >> yp[i];
    if (dim>2) {
      input_file >> zp[i];
    }
    input_file>> pRef[i];
    node_ind[i] = i+1;
  }


  // Edges
  nEdges = 0;
  do {
    getline(input_file,_line,'\n');
    stripSpace(_line); 
  } 
  while (_line.find("Edges")==string::npos && !input_file.eof());
  
  if (_line.find("Edges")!=string::npos) {
    input_file >> nEdges;
  }
  if (nEdges) {
    edge.resize(2*nEdges);
    edgeRef.resize(nEdges);
    for (int i=0;i<(int)nEdges;i++) {
      input_file >> edge[2*i] >> edge[2*i+1] >> edgeRef[i];
    }
  }

  input_file.clear();
  input_file.seekg(current_line_number);
  
  nTria = 0;
  // Triangles
  do {
    getline(input_file,_line,'\n');
    stripSpace(_line); // borra espacios en blanco
    // debido a que el programa NB2MESH usa en comentario "Triangles" ==> *****!!
    if (_line.find("#")!=std::string::npos)
      _line="####COMENTARIO###";
  } while (_line.find("Triangles")==std::string::npos && !input_file.eof());
  
  if (_line.find("Triangles")!=std::string::npos) {
    input_file >> nTria;
    if (nTria) {
      tria.resize(3*nTria);
      triaRef.resize(nTria);
      for (int i=0;i<(int)nTria;i++) {
	input_file >> tria[3*i] >> tria[3*i+1] >> tria[3*i+2] >> triaRef[i];
      }
    }
    
  }
  input_file.clear();
  input_file.seekg(current_line_number);

  nTetra = 0;
  // Tetra
  do {
    getline(input_file,_line,'\n');
    stripSpace(_line); // borra espacios en blanco
    // debido a que el programa NB2MESH usa en comentario "Triangles" ==> *****!!
    if (_line.find("#")!=std::string::npos)
      _line="####COMENTARIO###";
  } while (_line.find("Tetrahedra")==std::string::npos && !input_file.eof());
  
  if (_line.find("Tetrahedra")!=std::string::npos) {
    input_file >> nTetra;
    if (nTetra) {
      tetra.resize(4*nTetra);
      tetraRef.resize(nTetra);
      for (int i=0;i<(int)nTetra;i++) {
	input_file >> tetra[4*i] >> tetra[4*i+1] 
		   >> tetra[4*i+2] >> tetra[4*i+3] >> tetraRef[i];
      }
    }
    
  }


  if (dim==2){
    nElem = nTria;
  } else {
    nElem = nTetra;
  }
  
  cellsInTriaOld.resize(nElem);

  computeVolumes();

}

void Mesh::computeVolumes()
{
  cout << " computing volumes " << endl;
  //double total_volume = 0;
  if (dim==2){
    // todo...
    cout << " function not yet implemented " << endl;
    exit(1);
  } else {
    elementVolume.resize(nTetra);
    for (unsigned int i=0; i<nTetra; i++) {
      double v11 = xp[tetra[4*i+1]-1] - xp[tetra[4*i]-1];
      double v12 = yp[tetra[4*i+1]-1] - yp[tetra[4*i]-1];
      double v13 = zp[tetra[4*i+1]-1] - zp[tetra[4*i]-1];
      
      double v21 = xp[tetra[4*i+2]-1] - xp[tetra[4*i]-1];
      double v22 = yp[tetra[4*i+2]-1] - yp[tetra[4*i]-1];
      double v23 = zp[tetra[4*i+2]-1] - zp[tetra[4*i]-1];
      
      double v31 = xp[tetra[4*i+3]-1] - xp[tetra[4*i]-1];
      double v32 = yp[tetra[4*i+3]-1] - yp[tetra[4*i]-1];
      double v33 = zp[tetra[4*i+3]-1] - zp[tetra[4*i]-1];

      elementVolume[i] = 1./6.*(v11*(v22*v33-v23*v32) -
				v12*(v21*v33-v23*v31) +
				v13*(v21*v32-v22*v31));
      //cout << "tetrahedra " << i << " has volume " << elementVolume[i] << endl;
      //total_volume += fabs(elementVolume[i]);
      
    }
    
  }
  //cout << " total volume = " << total_volume << endl;
}

void Mesh::info(){

  if (dim>0) {
    cout << " --- dimension: " << dim << endl;
    cout << " --- n. of nodes: " << nNodes << endl;
    cout << " --- n. of edges: " << nEdges << endl;
    cout << " --- n. of triangles: " << nTria << endl;
    cout << " --- n. of tetrahedra: " << nTetra << endl;
  }
}



void Mesh:: write(string _file)
{

  ofstream o_file;
  o_file.open(_file.c_str(),ios::out);
  if (!o_file) {
    cerr << " Mesh::write(): *** ERROR *** could not open file " 
	 << _file << endl;
    exit(1);
  }

  cout << " ** Mesh:: writing on " << _file << endl;
  info();


  o_file << "MeshVersionFormatted 1" << endl;
  o_file << "Dimension" << endl;
  o_file << dim << endl;
  o_file << "Vertices" << endl;
  o_file << nNodes << endl;

  for (int i=0; i<nNodes; i++){
    o_file << xp[i] << " " << yp[i] << " ";
    if (dim==3) {
      o_file << zp[i] << " ";
    }
    o_file << pRef[i] << endl;
  }

  if (dim==2){

    // edges
    o_file << endl;
    o_file << "Edges" << endl;
    o_file << nEdges << endl;
    for (unsigned int i=0; i<edgeRef.size(); i++) {
      o_file << edge[2*i] << " " << edge[2*i+1] 
	     << " " << edgeRef[i] << endl;
    }
    
  }
  
  // triangles
  o_file << endl;
  o_file << "Triangles" << endl;
  o_file << nTria << endl;
  for (unsigned int i=0; i<triaRef.size(); i++) {
    o_file << tria[3*i] << " " << tria[3*i+1] << " " << tria[3*i+2]
	   << " " << triaRef[i] << endl;
  }


  if (dim==3){
    // tetra
    o_file << endl;
    o_file << "Tetrahedra" << endl;
    o_file << nTetra << endl;
    for (unsigned int i=0; i<tetraRef.size(); i++) {
      o_file << tetra[4*i] << " " << tetra[4*i+1] << " "  
	     << tetra[4*i+2] << " " << tetra[4*i+3]
	     << " " << tetraRef[i] << endl;
    }
  }


  cout << " ... done " << endl;

}





 
// ---- fem-cells coupling
void Mesh::initCellsArrays()
{
  
  cellsInTria.resize(nElem);
  cellsInTria.clear();

  cellsInTriaNorm.resize(nElem);
  cellsInTriaNorm.clear();

  cellsInTriaHypo.resize(nElem);
  cellsInTriaHypo.clear();

  cellsInTriaDead.resize(nElem);
  cellsInTriaDead.clear();

}

void Mesh::storeCellDensity()
{
  cellsInTriaOld.resize(nElem);
  for (unsigned int i=0; i<nElem; i++){
    cellsInTriaOld[i] = cellsInTria[i];
  }
}

double Mesh::cellDensityChange()
{
  ///@todo should we divide by element size?
  double difference=0;
  for (unsigned int i=0; i<nElem; i++){
      difference += (cellsInTria[i]-cellsInTriaOld[i])*(cellsInTria[i]-cellsInTriaOld[i]);
  }
  //cellsInTriaOld.clear();
  return sqrt(difference);
}
