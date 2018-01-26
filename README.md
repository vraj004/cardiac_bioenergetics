cardiac_bioenergetics
===========

A set of codes to simulate mitochondrial OXPHOS and energy metabolism in cardiac cells using FE methods. 
The code uses the OpenCMISS libraries (www.opencmiss.org) to set up the simulation framework. The meshes used in these simulations were generated from
2D cross sectional images of cardiomyocytes. 

**required Applications/Packages**
----------------------------------
openCMISS and associated libraries - www.opencmiss.org

**BRANCH INFORMAITON**
----------------------
master:
 - contains the base codes for simulating mitochondrial OXPHOS and subsequent metabolite exchange between mitochondria and myofibrils. 
The simulation uses reaction diffusion equations to model the diffusion of metabolites over a realistic FE mesh dervied from electron 
microsocpy images. 

**FOLDER INFORMATION**
----------------------
 src/:
 - directory containing the fortran 90 program routine that uses opencmiss libraries to simulate mitochondrial OXHPOS and diffusion of metabolites in a cross section of a cell.
 build/:
  - directory containing all the necessary inputs that are required for the simulations
    - inputs.txt containts the diffusvity values and intial micromolar concentrations of the metabolites simulated in the model 
    MESH/:
        - filenames containing .1.node/ele/face : trinagle generated mesh files of a 2D cross section from SBM derived rat ventricular myocyte
 


RUNNING THE PROGRAM
-------------------
1. Interested users need to install compiled versions of the opencmiss libraries. Instructions and support are available at www.opencmiss.org. 
2. Once installed, generate a Makefile in the cardiac_bioenergetics directory using using CMake: "cmake -DOpenCMISSLibs_DIR=YOUR_OPENCMISS_INSTALL_DIR_HERE ."
3. Build the Fortran source file with the generated Makefile: "make"
4. File 'inputs.txt' contains default settings to run a simulation. Comments (prefixed by #) explain the different input variables.
5. Run the executable: "./src/cardiac_bioenergetics"
