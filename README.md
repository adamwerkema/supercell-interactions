# supercell-interactions
Code and data for Interactions Between Supercells in Multistorm Simulations

This repository contains the custom files needed to reproduce the CM1 simulations used in "Interactions Between Supercells in Multistorm Simulations". The files can be found in the /cm1_input directory:

1) init3d.F files
- There is a separate init3d.F file for each of the 16 multistorm simulations in the study (four separation distances x four separation angles).
- The same init3d_isolated file can be used for each of the 8 isolated supercell simulations.
- Remember to rename each file as init3d.F when placing it in the /src directory in CM1.

2) input soundings
- The same input sounding is used for all of the 16 multistorm simulations. It is simply titled "input_sounding".
- For the 8 isolated supercell ensemble, there are 7 perturbed input soundings. Use one for each supercell, along with the original input_sounding file for the 8th supercell.
- Remember to rename each file as input_sounding when using it in CM1.

3) namelists
- There are five namelist files: one for the isolated ensemble and one for each of the four separation distances. Use each namelist according which simulation you would like to reproduce.
- Reminder to rename the namelist file as namelist.input

The track_mesocyclones.py script tracks the mesocyclone(s) from a CM1 simulation. It outputs the location of the mesocyclone(s) at each CM1 output time, allowing statistics about a supercell/mesocyclone to be calculated from a mesocyclone-relative perspective.

Direct any questions to adam.werkema@noaa.gov
