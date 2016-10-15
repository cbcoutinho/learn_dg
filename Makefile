# makefile: makes the learn_dg program

current_dir = $(shell pwd)
SRC=$(current_dir)/src
OBJ=$(current_dir)/obj
BIN=$(current_dir)/bin
FORTRANLIB_SRC=$(current_dir)/fortranlib/src

# Compiler
FF = gfortran
#FFlags = -Wall -fbounds-check
#FLIBS = -lblas -llapack

# Dependencies of main program
objects=$(OBJ)/base_types.o $(OBJ)/integration.o

$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFlags) -o $@ $+ $(FLIBS)
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFlags) -I$(OBJ) -c -o $@ $<
$(OBJ)/integration.o: $(SRC)/integration.f90 $(OBJ)/lib_array.o $(OBJ)/base_types.o
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<

$(OBJ)/lib_array.o: $(FORTRANLIB_SRC)/lib_array.f90
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<
$(OBJ)/base_types.o: $(FORTRANLIB_SRC)/base_types.f90
	$(FF) $(FFlags) -J$(OBJ) -c -o $@ $<

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(BIN)/main
# run: $(BIN)/main
	# mpiexec $(BIN)/main
run: $(BIN)/main
	$(BIN)/main
