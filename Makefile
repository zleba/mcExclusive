include pythPath


CC=g++
CFLAGS=-g -MMD -MP -I./inc  -I$(pythInc)  \
                      -ansi -pedantic -W  -fPIC 


DEPS = inc/ColorTensor.h inc/ColorData.h inc/SpinTensor.h inc/ExclusiveHooks.h  \
       inc/Sigma2Process.h inc/SigmaQCD.h inc/SigmaElectroweak.h \
       inc/KMRlum.h inc/KMRlumi.h inc/ColorReconection.h inc/Basics.h inc/TopoSelector.h \

SRCS = src/ColorData.cpp src/ColorDataGG.cpp  src/ColorDataQQ.cpp  src/ColorTensor.cpp  \
       src/ColorReconection.cpp src/SpinTensor.cpp  \
       src/Sigma2Process.cpp src/SigmaQCD.cpp src/SigmaElectroweak.cpp  \
       src/ExclusiveHooks.cpp  src/KMRlum.cpp src/KMRlumi.cpp 


OBJS = obj/ColorData.o obj/ColorDataGG.o  obj/ColorDataQQ.o  obj/ColorTensor.o  \
       obj/ColorReconection.o obj/SpinTensor.o  \
       obj/Sigma2Process.o obj/SigmaQCD.o obj/SigmaElectroweak.o \
       obj/ExclusiveHooks.o  obj/KMRlum.o obj/KMRlumi.o 

DEP=$(OBJS:.o=.d) $(USROBJ:.o=.d)


obj/%.o: src/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

lib/libExclusive.so: $(OBJS)
	$(CXX) -shared $^ -o $@ -Wl,--whole-archive  -Wl,--no-whole-archive



-include $(DEP)
