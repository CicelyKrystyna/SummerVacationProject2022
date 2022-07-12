/*
 *  PDE.cpp
 *  class for continous PDE
 *
 */
#include "PDE.h"
#include <iostream>
#include <sstream>



// initialize class
void PDE::init(const Param &_params, const Mesh& _m)
{

  string FreeFemMACOS = "FreeFem++"; // for MacOS
  string FreeFemLinux = "FreeFem++3.20-x11"; // for linux
  string FreeFemOS = _params.FreeFemCall; // passed as a parameter

  FreeFemCall =  FreeFemOS + " " + _params.FreeFemFile;
  cout << " Freefem will be called with: " << endl;
  cout << ">> " << FreeFemCall.c_str()  << endl;


  mesh = _m;
  time = 0.;
  iter = 0;
  launch = false;


}

void PDE::init()
{
  cout << " PDE.cpp:: no FEM solver specified" << endl;
  time = 0.;
  iter = 0;
  launch = false;
}

// call the solver
void PDE::solve(double _time)
{

  // solve diffusion
  //cout << " Start Freefem " << endl;
  iter++;
  ostringstream inp;
  inp <<  FreeFemCall << " -v 0 -iter " << iter << " -time " << _time;
  //string FreeFemCall =  FreeFemLinux + " " + FreeFemFile + " -iter "+ i 
  system(inp.str().c_str());
  //cout << " End Freefem " << endl;

}
  
