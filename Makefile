# makefile: makes the learn_dg program

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin
FORTRANLIB_SRC=$(current_dir)/src/fortranlib/src

# Compiler
FF = gfortran
RM = rm -f
FFlags = -Wall -fbounds-check
FLIBS = -lblas -llapack

.DEFAULT_GOAL := $(BIN)/main

# Dependencies of main program
objects=$(OBJ)/base_types.o \
	$(OBJ)/lib_array.o \
	$(OBJ)/integration.o \
	$(OBJ)/misc.o \
	$(OBJ)/legendre.o

# Fortran library
$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/base_types.o: $(FORTRANLIB_SRC)/base_types.f90
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/misc.o: $(SRC)/misc.f90 $(OBJ)/base_types.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/legendre.o: $(SRC)/legendre.f90 $(OBJ)/misc.o $(OBJ)/lib_array.o $(OBJ)/base_types.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/integration.o: $(SRC)/integration.f90 $(OBJ)/lib_array.o $(OBJ)/base_types.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<

# Main program
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFlags) -I$(OBJ) -c -o $@ $<
$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFlags) -o $@ $+ $(FLIBS)

clean:
	$(RM) $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main

# run: $(BIN)/main
#		mpiexec $(BIN)/main
run: $(BIN)/main
	$(BIN)/main
