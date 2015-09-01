* Input directories are used in the following way
  * _3DEff: eff (y, pT, cent) from PRMC. Always measured at Lxy=0.
  * _4DEff: eff (y, pT, cent) at Lxy=0 - eff (y, pT, cent, Lxy). Obtained from NPMC
  * Final eff = 3D Eff - 4D Eff

* checkDimuons.cpp
  * Require 2 RooDataSets and produce various comparison plots

* CompDimuon.cpp
  * Dimuon mass distributions will be compared from 3 different RooDataSets
