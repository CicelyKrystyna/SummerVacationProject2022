/************************************************************************/
/*                              MAIN FUNCTION                           */
/*    authors: A. Caiazzo & I. Ramis-Conde                              */
/*    date: 12.2014                                                      */
/************************************************************************/

#include <iostream>
#include "CoupledModel.h"
#include <sstream>

using namespace std;

int main(int argc, char** const argv)
{
  clock_t total_time = clock();
  CoupledModel m;
  
  if (argc==1) {
    cout << " ERROR in main.cpp: input file missing " << endl;
    exit(1);
  }

  int tot_n_of_simulations=1;
  if (argc==3) {
    istringstream ss(argv[2]);
    ss >> tot_n_of_simulations;
  }
  cout << " === I will run  " << tot_n_of_simulations << " simulations" << endl;

  for (int i_sim = 0; i_sim < tot_n_of_simulations; i_sim++) {
    cout << " ====== START SIMULATION " << i_sim << endl;
    system("date\n");
    cout << " =============================== " << endl;
    
    // --- READ PARAMETERS FROM INPUT FILE ---
    // --- INITIALIZE SIMULATION ---
    m.init(argv[1]);
    m.simulation_id = i_sim;
    // --- MAIN LOOP
    m.loop();
  
    // -- END SIMULATION
    m.end();
  }
  total_time = clock() - total_time;
  cout << " # of cells total: " << m.total_no_of_cells + m.total_no_of_removed_cells << endl;
  cout << " total Wall Clock time:" << ((float)total_time)/CLOCKS_PER_SEC << endl;
  cout << " ========= " << endl;
  cout << " program finished (Normal end)\n";
  system("date\n");
  cout << " ========= " << endl;
  
  return 0;
  
}
