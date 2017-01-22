# makefile: makes the learn_dg program

SRC=$(shell pwd)/src
OBJ=$(shell pwd)/obj
BIN=$(shell pwd)/bin
DOCS=$(shell pwd)/docs
FORTRANLIB_SRC=$(shell pwd)/src/fortranlib/src

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -O0 -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS = -lblas -llapack
# FLIBS += -fopenmp

FORD_FLAGS = -d $(SRC) \
	-p $(DOCS)/user-guide \
	-o $(DOCS)/html

# Dependencies of main program
objects=$(OBJ)/lib_array.o \
	$(OBJ)/lib_algebra.o \
	$(OBJ)/integration.o \
	$(OBJ)/misc.o \
	$(OBJ)/legendre.o \
	$(OBJ)/io.o \
	$(OBJ)/assembly.o \
	$(OBJ)/linalg.o

default: all

# Fortran library
$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/lib_algebra.o: $(FORTRANLIB_SRC)/lib_algebra.f90
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

all: $(BIN)/main

mesh: test1D.geo test2D.geo
	gmsh test1D.geo -order 1 -1 -o test1D.msh > /dev/null 2>&1
	gmsh test2D.geo -order 1 -2 -o test2D.msh > /dev/null 2>&1

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main test1D.msh test2D.msh

run: all mesh
	$(BIN)/main test1D.msh

debug: clean all
	/usr/bin/valgrind --track-origins=yes --leak-check=full $(BIN)/main test1D.msh

plot: run
	python plotter.py

.PHONY: docs
docs: $(DOCS)/learn_dg.md
	cp README.md $(DOCS)/README.md
	ford $(FORD_FLAGS) $(DOCS)/learn_dg.md
	rm $(DOCS)/README.md
