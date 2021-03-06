# SummerVacationProject2022
In this project we will use a previously developed c++ 3D off-lattice agent-based multi-scale model which simulates the behaviour of, and spatio-temporal interactions between various tumour agents (namely cells and blood vessels) to investigate how tumour growth patterns and inherent spatial heterogeneity are driven by oxygen phenotypes of cells. Within the existing model mechanical interactions between agents occur through repulsion and adhesive forces and these are coupled to a finite element solver which solves a reaction-diffusion equation for chemical substances (e.g. oxygen) which diffuse throughout the 3D tissue domain from sources (e.g. blood vessels) and are consumed by agents (e.g. cells). Agent behaviour is governed both by the mechanical interactions and changes to their phenotype. In this project we will allow cells to sample across a discrete range of phenotypes representing a continuum between normoxic and hypoxic behaviour in order to investigate how cell phenotype is structured by the local environment and how in turn cell phenotype affects patterns of tumour growth. Thus allowing us to better understand the specific mechanisms that occur in the tumour microenvironment and underpin tumour development. This project seeks to apply an agent-based approach to compare to an earlier continuous PDE formalism of phenotypic heterogeneity in vascularised tumours.

# Downloading the source code:

Step 1: become familiar with the windows command line (CMD) if you are not already https://www.makeuseof.com/tag/a-beginners-guide-to-the-windows-command-line/

Step 2: within the CMD navigate to your preferred directory (folder e.g. I have my code sitting within my Desktop) and type git clone https://github.com/CicelyKrystyna/SummerVacationProject2022.git to download the code (you may need to install git https://git-scm.com/download/win). You should then have a new directory which contains the source code (src), this README file, an examples directory plus some additional documentation.

Step 3: In order to run the code you will need to have FreeFem installed (download from https://freefem.org/) and a c++ compiler (this video may be helpful https://www.youtube.com/watch?v=ml3xutoE5JM). To view the code (optional) you will need something like visual studio (https://visualstudio.microsoft.com/vs/features/cplusplus/)

Step 4: Navigate to src and type make (again for windows you may need to install make - see https://www.technewstoday.com/install-and-use-make-in-windows/) in CMD which will compile the code (fingers crossed for no errors). 

Step 5: To run the code navigate to example/first_simulation and type ../../src/cell++ input.dat in CMD hopefully the code will run :)

Step 6: In order to store results from the simulation create directories called vtk; output and fem_results within the first_simulation directory.
