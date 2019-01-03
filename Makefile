LIBNAME:=libaction.a
lib_reduction:=modules/reduction
lib_core:=modules/core
libraries:=$(lib_reduction)	$(lib_core)
library_files:=modules/core/libactioncore.a modules/reduction/libactionred.a

MATLAB_SRC:=$(wildcard	modules/reduction/*.cc	modules/core/*.cc)
MATLAB_MEX=$(MATLAB_SRC:.cc=.mexa64)
MATLAB_FLAGS=-DUSE_BLAS_LIB -DAXPBY -DINT_64BITS -DNDEBUG -largeArrayDims

MATLAB = $(shell matlab -e | sed -n 's/MATLAB=//p')
MEX = $(MATLAB)/bin/mex

.PHONY: all matlab R clean $(libraries)
all: $(LIBNAME)

$(libraries):
	$(MAKE) --directory=$@

$(LIBNAME): $(libraries)
	ar -rcs $@ $(library_files)

matlab:
	$(MAKE) matlab --directory=$(lib_reduction)
	$(MAKE) matlab --directory=$(lib_core)
    
R:
	$(MAKE) R --directory=$(lib_reduction)
	$(MAKE) R --directory=$(lib_core)	

clean:
	rm -f $(LIBNAME)
	
