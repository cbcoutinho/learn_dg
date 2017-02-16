# makefile: makes the learn_dg program

##############################
###### Windows vs Linux ######
##############################



##############################
#### Project Directories #####
##############################

SRC=./src
OBJ=./obj
BIN=./bin
DOC=./docs
TEST=./test
BLD=./build
FLIB_SRC=./src/fortranlib/src

##############################
###### Compiler options ######
##############################

ifeq ($(OS),Windows_NT)
RM = del /Q
FixPath = $(subst /,\,$1)
CMFLAGS = -G "MinGW Makefiles"
MAIN = main.exe
DOUBLEINT = doubleint.exe
else
#  ifeq ($(shell uname), Linux)
RM = rm -rf
FixPath = $1
MAIN = main
DOUBLEINT = doubleint
#  endif
endif

# Use debug flags unless otherwise stated
ifndef DEBUG
DEBUG:=1
endif

ifndef FF
FF:=gfortran
endif

# FFVERSIONGTEQ6 := $(shell expr `$(FF) -dumpversion | cut -f1 -d.` \>= 6)
# ifeq "$(FFVERSIONGTEQ6)" "0"
# 	@echo "Fortran Compiler must be at least version 6"
# endif

FFLAGS := -std=f2008 -fPIC -fmax-errors=1
ifeq ($(DEBUG),1)
# if [[ "$(DEBUG)" =~ ^[Yy]$ ]]; then
# Debug flags:
$(info DEBUG is $(DEBUG))
$(info Building in Debug mode)
FFLAGS += -Og -g -fcheck=all -fbacktrace -Wimplicit-interface -Wall -Wextra #-ffpe-trap=zero,overflow,underflow
else
# Release flags:
$(info DEBUG is $(DEBUG))
FFLAGS += -Ofast -march=native -ffast-math -funroll-loops
$(info Building in Release mode)
endif

FLIBS := -lblas -llapack
# FLIBS += -fopenmp

##############################
######## FORD options ########
##############################

ifneq ("$(wildcard $(HOME)/Projects/ford/ford.py)","")
# FILE_EXISTS = 1
FORD := $(HOME)/Projects/ford/ford.py
else
# FILE_EXISTS = 0
FORD := ford
endif

FORD_FLAGS := -d $(SRC) \
	-p $(DOC)/user-guide \
	-o $(DOC)/html

##############################
#### Project Dependencies ####
##############################

# Dependencies of main program
objects:=$(OBJ)/lib_array.o \
	$(OBJ)/integration.o \
	$(OBJ)/misc.o \
	$(OBJ)/legendre.o \
	$(OBJ)/pascal_1D.o \
	$(OBJ)/pascal_2D.o \
	$(OBJ)/io.o \
	$(OBJ)/assembly.o \
	$(OBJ)/linalg.o

##############################
####### Build Recepies #######
##############################

default: all

# Fortran library
$(OBJ)/lib_array.o: $(FLIB_SRC)/lib_array.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Modules
$(OBJ)/misc.o: $(SRC)/misc_mod.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/assembly.o: $(SRC)/assembly_mod.f90 $(OBJ)/legendre.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<
$(OBJ)/linalg.o: $(SRC)/linalg_mod.f90
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/pascal_1D.o: $(SRC)/pascal_1D_smod.f90 $(OBJ)/legendre.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/pascal_2D.o: $(SRC)/pascal_2D_smod.f90 $(OBJ)/legendre.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/legendre.o: $(SRC)/legendre_mod.f90 $(OBJ)/misc.o $(OBJ)/lib_array.o $(OBJ)/integration.o $(OBJ)/linalg.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/integration.o: $(SRC)/integration_mod.f90 $(OBJ)/lib_array.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $< $(FLIBS)
$(OBJ)/io.o: $(SRC)/io_mod.f90 $(OBJ)/lib_array.o
	$(FF) $(FFLAGS) -J$(OBJ) -c -o $@ $<

# Main program
$(OBJ)/main.o: $(SRC)/main.f90 $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $< $(FLIBS)
$(BIN)/main: $(OBJ)/main.o $(objects)
	$(FF) $(FFLAGS) -o $@ $^ $(FLIBS)

$(OBJ)/doubleint.o: $(SRC)/doubleint.f90 $(objects)
	$(FF) $(FFLAGS) -I$(OBJ) -c -o $@ $< $(FLIBS)
$(BIN)/doubleint: $(OBJ)/doubleint.o $(objects)
	$(FF) $(FFLAGS) -o $@ $^ $(FLIBS)

submodules:
	git submodule init
	git submodule update

all: submodules $(BIN)/main $(BIN)/doubleint

mesh: $(call FixPath,$(TEST)/test1D.geo) $(call FixPath,$(TEST)/test2D.geo)
	gmsh $(call FixPath,$(TEST)/test1D.geo) -order 1 -1 -o $(call FixPath,$(TEST)/test1D.msh)
	gmsh $(call FixPath,$(TEST)/test2D.geo) -order 1 -2 -o $(call FixPath,$(TEST)/test2D.msh)

run1: all mesh
	$(call FixPath,$(BIN)/$(MAIN)) $(call FixPath,$(TEST)/test1D.msh)

run2: all
	$(call FixPath,$(BIN)/$(DOUBLEINT))

run: run1 run2

plot: cmake
	python plotter.py

.PHONY: docs
docs: submodules $(DOC)/learn_dg.md README.md
	cp README.md $(DOC)/README.md
	$(FORD) $(FORD_FLAGS) $(DOC)/learn_dg.md
	$(RM) $(DOC)/README.md

debug: clean all mesh
	valgrind --track-origins=yes --leak-check=full $(call FixPath,$(BIN)/$(MAIN)) $(call FixPath,$(TEST)/test1D.msh)
	valgrind --track-origins=yes --leak-check=full $(call FixPath,$(BIN)/$(DOUBLEINT))

# .PHONY: cmake
cmake: submodules mesh | $(BLD)
	cd $(call FixPath,$(BLD)) && cmake .. $(CMFLAGS) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) && cd ..
	$(MAKE) -C $(call FixPath,$(BLD))
	$(call FixPath,$(BLD)/$(MAIN)) $(call FixPath,$(TEST)/test1D.msh)

.ONESHELL:
$(BLD):
ifeq ($(OS),Windows_NT)
	-@ if not exist "$(BLD)" mkdir $(BLD)
else
	test -d $(BLD) || mkdir $(BLD)
endif

clean:
	$(RM) $(call FixPath,$(OBJ)/*.o) $(call FixPath,$(OBJ)/*.mod) $(call FixPath,$(OBJ)/*.smod) $(call FixPath,$(BIN)/$(MAIN)) $(call FixPath,$(OBJ)/doubleint)
	$(RM) $(call FixPath,$(TEST)/test1D.msh) $(call FixPath,$(TEST)/test2D.msh)
	$(RM) $(call FixPath,$(BLD))

# $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.smod $(BIN)/main
