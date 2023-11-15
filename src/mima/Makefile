device=omp

#configure machine
include ../../config/default.mk
include ../../config/version.mk
include ../../config/*.mk
include ../../config/devices/devices.mk

INCLUDE+= -I../../inc   # other project libraries

all: mima mima_hpc mima_mpi

mima: mima.cpp mima.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -DWITH_GLFW -g

mima_hpc: mima.cpp mima.h parameters.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -g

mima_mpi: mima.cpp mima.h parameters.h diag.h
	$(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS)

doc:
	pdflatex -shell-escape ./mima.tex;
	bibtex mima.aux;
	pdflatex -shell-escape ./mima.tex;
	pdflatex -shell-escape ./mima.tex;


.PHONY: clean doc

clean:
	rm -rf mima mima_hpc mima_mpi mima.aux mima.log mima.out mima.pyg mima.pdf mima.bbl mima.blg
