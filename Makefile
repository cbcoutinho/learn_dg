# makefile: makes the learn_dg program

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin
FORTRANLIB_SRC=$(current_dir)/src/fortranlib/src

# Compiler
FF = gfortran
FFLAGS = -Wall -std=f2008 -Wextra -fPIC -fmax-errors=1 -Wimplicit-interface
# Debug flags:
FFLAGS += -O0 -g -fcheck=all -fbacktrace #-ffpe-trap=zero,overflow,underflow
# Release flags:
# FFLAGS += -O3 -march=native -ffast-math -funroll-loops

FLIBS = -lblas -llapack

.DEFAULT_GOAL := $(BIN)/main

# Dependencies of main program
objects=$(OBJ)/lib_array.o \
	$(OBJ)/integration.o \
	$(OBJ)/misc.o \
	$(OBJ)/legendre.o \
	$(OBJ)/io.o

# Fortran library
$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/misc.o: $(SRC)/misc.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/legendre.o: $(SRC)/legendre.f90 $(OBJ)/misc.o $(OBJ)/lib_array.o $(OBJ)/integration.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/integration.o: $(SRC)/integration.f90 $(OBJ)/lib_array.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/io.o: $(SRC)/io.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Main program
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $<
$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFLAGS) -o $@ $+ $(FLIBS)

mesh: $(BIN)/test1D.geo
	gmsh $(BIN)/test1D.geo -order 1 -1 $(BIN)/test1D.msh > /dev/null 2>&1

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main $(BIN)/test1D.msh

run: $(BIN)/main mesh
	$(BIN)/main $(BIN)/test1D.msh

debug: clean $(BIN)/main mesh
	/usr/bin/valgrind --track-origins=yes --leak-check=full $(BIN)/main $(BIN)/test1D.msh
