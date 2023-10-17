device=cpu

#configure machine 
include ../../config/default.mk
include ../../config/*.mk 
include ../../config/devices/devices.mk
include ../../config/version.mk

#INCLUDE+= -I../         # other src libraries
INCLUDE+= -I../../inc   # other project libraries

all: hw mima

toefl: toefl.cpp toefl.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -DWITH_GLFW -g

hw: hw.cu hw.cuh 
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(JSONLIB) -DDG_BENCHMARK 

mima: mima.cu mima.cuh 
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(JSONLIB) -DDG_BENCHMARK

.PHONY: clean

clean:
	rm -f hw mima toefl
