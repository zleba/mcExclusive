# mcExclusive

This Add-On to Pythia 8 allows to generate the Central Exclusive Production.
It is based on UserHooks machinery.
The installation instructions are described in the INSTALL file, at the current state, one needs to add 2 lines to the Pythia source code and recompile it before compiling this package.

For more information about the physics behind, see [arXiv:1608.03765], where the method is explained and some results are shown.
In case of usage, please cite this paper.

Very briefly, the program is based on the KMR model extended for possible ISR.
The S2 can be calculated using Pythia's MPI or take from other program/paper.
Beware of high sensitivity of the cross sections to the PDF behaviour in the low-scale region.

The list of currently implemented processes is given in examples/standalone/dijets.cpp
Note, that the produced events are weighted.
