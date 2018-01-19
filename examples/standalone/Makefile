include ../../pythPath
CC=g++

#ROOTCFLAGS = $(shell root-config --cflags)
#ROOTLIBS   = $(shell root-config --libs)

mcExclusiveInc = ../../inc/
mcExclusiveLib = ../../lib/


dijets: dijets.cpp
	$(CXX) -I./usrInc -I$(mcExclusiveInc) -I$(pythInc) $(ROOTCFLAGS)   $^ -o $@   -L$(mcExclusiveLib) -lExclusive  -Wl,-R$(mcExclusiveLib)  -L$(pythLib) -lpythia8   -lpythia8lhapdf6 -Wl,-R$(pythLib)  -lfastjet  $(ROOTLIBS) -ldl
