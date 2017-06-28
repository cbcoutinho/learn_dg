# makefile: makes the learn_dg program

##############################
#### Project Directories #####
##############################

SRC_DIR=./src
DOC_DIR=./docs
TEST_DIR=./test
BLD_DIR=./build
FLIB_SRC_DIR=./src/fortranlib/src

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

FORD_FLAGS := -d $(SRC_DIR) \
	-p $(DOC_DIR)/user-guide \
	-o $(DOC_DIR)/html

##############################
####### Build Recepies #######
##############################

default: all

submodules:
	git submodule init
	git submodule update

all: cmake

mesh: $(TEST_DIR)/test1D.geo $(TEST_DIR)/test2D.geo
	gmsh $(TEST_DIR)/test1D.geo -order 1 -1 -o $(TEST_DIR)/test1D.msh
	gmsh $(TEST_DIR)/test2D.geo -order 1 -2 -o $(TEST_DIR)/test2D.msh


# Build and test the project
cmake: submodules mesh | $(BLD_DIR)
	cd $(BLD_DIR) && cmake .. $(CMFLAGS) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) && cd ..
	$(MAKE) -C $(BLD_DIR)

cmake_win: submodules mesh | $(BLD_DIR)
	cd $(BLD_DIR) && cmake -DCMAKE_TOOLCHAIN_FILE:STRING=../cmake/Toolchain-x64-mingw32.cmake .. $(CMFLAGS) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) && cd ..
	$(MAKE) -C $(BLD_DIR)

tests: cmake
	$(BLD_DIR)/bin/driverA
	$(BLD_DIR)/bin/driverA $(TEST_DIR)/test1D.msh

# After running one of the tests, plot the output
plot: cmake tests
	python test/plotter.py
	$(RM) data.out



# Other

.PHONY: docs
docs: submodules $(DOC_DIR)/learn_dg.md README.md
	cp README.md $(DOC_DIR)/README.md
	$(FORD) $(FORD_FLAGS) $(DOC_DIR)/learn_dg.md
	$(RM) $(DOC_DIR)/README.md

.ONESHELL:
$(BLD_DIR):
	test -d $(BLD_DIR) || mkdir $(BLD_DIR)

clean:
	$(RM) $(TEST_DIR)/test1D.msh $(TEST_DIR)/test2D.msh
	$(RM) $(BLD_DIR)
