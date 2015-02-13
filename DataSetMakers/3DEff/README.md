* Codes for efficiency for each y, pT, cent bins
* ./drawing.sh compiles and runs jobs with options
* binArrays.h, binArrays3D.h: Define y, pT, cent arrays, have to be put under a given y-pT directory
  - ex) Rap0.0-1.6_Pt6.5-30.0, notAbs_Rap-2.4--1.6_Pt3.0-30.0
  - binArrays.h: Finer ctau efficiency for a given y,pT,cent region
  - binArrays3D.h: Finer 1D efficiency of y,pT,cent (array[]), 2D efficiency of y,pT,cent (array2[])
* Below files need to be symbolic-linked to a given y-pT directory
  - ex) lJpsiEff.h, 3DEff.cpp, 3DEff_draw.cpp, lJpsiEff.cpp, lJpsiEff_draw.cpp, LxyEff_draw.cpp, LxyTrueReco.cpp
* Then ./drawing.sh visits all sub-directories and fill up efficiency histograms
* ...3DAnaBins_eff.root: Histograms at differential regions (array2[] is used) will be drawn as a function of pT, plotted at \<pT\>. "array[]" used histograms will give 1D efficiency as a function of y,pT,cent, respectively.
* ..._eff.root: All histograms will be drawn as a function of Lxy, plotted at \<Lxy\>
