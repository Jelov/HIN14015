* Codes for Lxy efficiency for a given y, pT, cent region
* Need 3DEff root files as an input
  - ex) .../Rap0.0-1.6_Pt6.5-30.0/NPMC_eff.root
* Arrays for bin boundaries MUST follow "binArrays.h" from 3DEff, not binArrays3D.h
  - binArrays.h contains bin boundaries for finer ctau efficiency 
  - binArrays3D.h contains bin boundaries for 3DEff
* Compile and run: ./compEffUnfoldingProfile
