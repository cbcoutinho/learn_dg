# makefile: makes the learn_dg program

##############################
#### Project Directories #####
##############################

SRC=./src
OBJ=./obj
BIN=./bin
DOC=./docs
TEST=./test
FLIB_SRC=./src/fortranlib/src

##############################
###### Compiler options ######
##############################

# Use debug flags unless otherwise stated
ifndef DEBUG
DEBUG="yes"
endif

FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
ifeq ($(DEBUG),"yes")
# Debug flags:
# $(info DEBUG is $(DEBUG))
$(info Building in Debug mode)
FFLAGS += -O0 -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
else
# Release flags:
# $(info DEBUG is $(DEBUG))
FFLAGS += -O3 -march=native -ffast-math -funroll-loops
$(info Building in Release mode)
endif

FLIBS = -lblas -llapack
# FLIBS += -fopenmp

##############################
######## FORD options ########
##############################

ifneq ("$(wildcard $(HOME)/Projects/ford/ford.py)","")
# FILE_EXISTS = 1
FORD = $(HOME)/Projects/ford/ford.py
else
# FILE_EXISTS = 0
FORD = ford
endif

FORD_FLAGS = -d $(SRC) \
	-p $(DOC)/user-guide \
	-o $(DOC)/html

##############################
#### Project Dependencies ####
##############################

# Dependencies of main program
objects=$(OBJ)/lib_array.o \
	$(OBJ)/integration.o \
	$(OBJ)/misc.o \
	$(OBJ)/legendre.o \
	$(OBJ)/io.o \
	$(OBJ)/assembly.o \
	$(OBJ)/linalg.o



##############################
####### Build Recepies #######
##############################

default: build

# Fortran library
$(OBJ)/lib_array.o: $(FLIB_SRC)/lib_array.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/misc.o: $(SRC)/misc.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/assembly.o: $(SRC)/assembly.f90 $(OBJ)/legendre.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/linalg.o: $(SRC)/linalg.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/legendre.o: $(SRC)/legendre.f90 $(OBJ)/misc.o $(OBJ)/lib_array.o $(OBJ)/integration.o $(OBJ)/linalg.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/integration.o: $(SRC)/integration.f90 $(OBJ)/lib_array.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/io.o: $(SRC)/io.f90 $(OBJ)/lib_array.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Main program
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $< $(FLIBS)
$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFLAGS) -o $@ $+ $(FLIBS)

submodules:
	git submodule init
	git submodule update

build: submodules $(BIN)/main

mesh: $(TEST)/test1D.geo $(TEST)/test2D.geo
	gmsh $(TEST)/test1D.geo -order 1 -1 -o $(TEST)/test1D.msh > /dev/null 2>&1
	gmsh $(TEST)/test2D.geo -order 1 -2 -o $(TEST)/test2D.msh > /dev/null 2>&1

run: build mesh
	$(BIN)/main $(TEST)/test1D.msh

plot: run
	python plotter.py

.PHONY: docs
docs: submodules $(DOC)/learn_dg.md README.md
	cp README.md $(DOC)/README.md
	$(FORD) $(FORD_FLAGS) $(DOC)/learn_dg.md
	rm $(DOC)/README.md

debug: clean build
	/usr/bin/valgrind --track-origins=yes --leak-check=full $(BIN)/main $(TEST)/test1D.msh

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main
	rm -f $(TEST)/test1D.msh $(TEST)/test2D.msh
