include ../../pythPath
CC=g++

#for ROOT
#ROOTCFLAGS = $(shell root-config --cflags) 
#ROOTLIBS   = $(shell root-config --libs)

mcExclusiveInc = ../../inc/
mcExclusiveLib = ../../lib/

#Possible additional Libraries  -lpythia8lhapdf6  -lfastjet
dijets: dijets.cpp
	$(CXX) -I./usrInc -I$(mcExclusiveInc) -I$(pythInc) $(ROOTCFLAGS)   $^ -o $@   -L$(mcExclusiveLib) -lExclusive  -Wl,-R$(mcExclusiveLib)  -L$(pythLib) -lpythia8  -Wl,-R$(pythLib)   $(ROOTLIBS) -ldl
