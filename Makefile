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

# Use debug flags unless otherwise stated
ifndef DEBUG
DEBUG:=1
endif

FF:=gfortran
RM:=rm -rf

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
####### Build Recepies #######
##############################

default: all

submodules:
	git submodule init
	git submodule update

all: cmake # $(BIN)/main $(BIN)/doubleint

mesh: $(TEST)/test1D.geo $(TEST)/test2D.geo
	gmsh $(TEST)/test1D.geo -order 1 -1 -o $(TEST)/test1D.msh
	gmsh $(TEST)/test2D.geo -order 1 -2 -o $(TEST)/test2D.msh

run1: all mesh
	$(BIN)/$(MAIN) $(TEST)/test1D.msh

run2: all
	$(BIN)/$(DOUBLEINT)

run: run1 run2

plot: cmake
	python plotter.py

.PHONY: docs
docs: submodules $(DOC)/learn_dg.md README.md
	cp README.md $(DOC)/README.md
	$(FORD) $(FORD_FLAGS) $(DOC)/learn_dg.md
	$(RM) $(DOC)/README.md

debug: clean all mesh
	valgrind --track-origins=yes --leak-check=full $(BIN)/$(MAIN)) $(TEST)/test1D.msh
	valgrind --track-origins=yes --leak-check=full $(BIN)/$(DOUBLEINT)

# .PHONY: cmake
cmake: submodules mesh | $(BLD)
	cd $(BLD) && cmake .. $(CMFLAGS) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) && cd ..
	$(MAKE) -C $(BLD)

cmake_win: submodules mesh | $(BLD)
	cd $(BLD)) && cmake -DCMAKE_TOOLCHAIN_FILE:STRING=../cmake/Toolchain-x64-mingw32.cmake .. $(CMFLAGS) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) && cd ..
	$(MAKE) -C $(BLD)

test: cmake
	$(BLD)/bin/driverA
	$(BLD)/bin/driverA $(TEST)/test1D.msh

.ONESHELL:
$(BLD):
	test -d $(BLD) || mkdir $(BLD)

clean:
	$(RM) $(OBJ)/*.o $(OBJ)/*.*mod
	$(RM) $(TEST)/test1D.msh $(TEST)/test2D.msh
	$(RM) $(BLD)
