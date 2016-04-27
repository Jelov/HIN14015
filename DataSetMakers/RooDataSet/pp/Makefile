ROOTCFLAGS = $(shell root-config --cflags)
ROOTGLIBS = $(shell root-config --glibs)

CPP = g++
CPPFLAGS = -g -fPIC -Wno-deprecated -O2 -ansi
LD = g++
LDFLAGS = -g
SOFLAGS = -shared

CPPFLAGS += $(ROOTCFLAGS)
NGLIBS = $(ROOTGLIBS)
NGLIBS += -lMathMore -lMinuit -lRooFit -lRooFitCore -lFoam
GLIBS = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR = ./
CPP += -I$(INCLUDEDIR)
OUTLIB = ./

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./

RooHistPdfConv.o:	$(INCLUDEDIR)/RooHistPdfConv.cpp
	$(CPP) $(CPPFLAGS) -c -o $(OUTLIB)/libRooHistPdfConv.o $(NGLIBS) $<

Tree2Datasets:	$(INCLUDEDIR)tree2Datasets.cpp
	$(CPP) $(CPPFLAGS) -o Tree2Datasets $(GLIBS) $ $<

Tree2DatasetsMC:	$(INCLUDEDIR)tree2Datasets_mc.cpp
	$(CPP) $(CPPFLAGS) -o Tree2DatasetsMC $(GLIBS) $ $<

compEffUnfoldingProfile:	$(INCLUDEDIR)compEffUnfoldingProfile.cpp
	$(CPP) $(CPPFLAGS) -o compEffUnfoldingProfile $(GLIBS) $ $<

Fit2DDataPbPb:	$(INCLUDEDIR)fit2DData_pbpb.cpp
	$(CPP) $(CPPFLAGS) -o Fit2DDataPbPb $(OUTLIB)/*.o $(GLIBS) $ $<

RandomSets:	$(INCLUDEDIR)makeRandomSet.cpp
	$(CPP) $(CPPFLAGS) -o RandomSets $(GLIBS) $ $<

