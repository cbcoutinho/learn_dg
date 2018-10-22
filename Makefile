# makefile: makes the learn_dg program

##############################
#### Project Directories #####
##############################

SRC_DIR=src
DOC_DIR=docs
TEST_DIR=tests
BLD_DIR=build
FLIB_SRC_DIR=src/fortranlib/src

##############################
###### Compiler options ######
##############################

FC         := gfortran
RM         := rm -rf

BUILD_TYPE ?= Debug
PROFILE    ?= OFF
USE_OPENMP ?= OFF

CMFLAGS= -B$(BLD_DIR) -H. \
	-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
	-DUSE_OPENMP=$(USE_OPENMP) \
	-DPROFILE=$(PROFILE) \
	-DCMAKE_Fortran_COMPILER=$(FC)

##############################
######## FORD options ########
##############################

FORD_FLAGS := -d ../$(SRC_DIR) \
	-p ../$(DOC_DIR)/user-guide \
	-o ../$(DOC_DIR)/html

##############################
####### Build Recepies #######
##############################

.PHONY: clean docs tags

all: cmake

mesh: $(TEST_DIR)/test1D.geo $(TEST_DIR)/test2D.geo
	gmsh -format msh2 -order 1 -1 -o $(TEST_DIR)/test1D.msh $(TEST_DIR)/test1D.geo
	gmsh -format msh2 -order 2 -2 -o $(TEST_DIR)/test2D.msh $(TEST_DIR)/test2D.geo

cleantest: clean
	${MAKE} test

cleanslowtest: clean
	${MAKE} slowtest

# Build and test the project
cmake: $(BLD_DIR)
	cmake $(CMFLAGS)
	cmake --build $(BLD_DIR)

cmake_win: $(BLD_DIR)
	cmake $(CMFLAGS) -DCMAKE_TOOLCHAIN_FILE:STRING=cmake/Toolchain-x86_64-w64-mingw32.cmake
	cmake --build $(BLD_DIR)

cmake_win32: $(BLD_DIR)
	cmake $(CMFLAGS) -DCMAKE_TOOLCHAIN_FILE:STRING=cmake/Toolchain-i686-w64-mingw32.cmake
	cmake --build $(BLD_DIR)

driver: cmake mesh $(BLD_DIR)
	$(BLD_DIR)/bin/driver1D
	$(BLD_DIR)/bin/driver1D $(TEST_DIR)/test1D.msh

sync: Pipfile.lock
	pipenv sync

test: cmake mesh driver sync
	pipenv run pytest -vs --cache-clear --duration=5 -m 'not slowtest'

slowtest: cmake mesh driver sync
	pipenv run pytest -vs --cache-clear --duration=5

# After running one of the tests, plot the output
plot: cmake driver
	python3 test/plotter.py

# Other
docs: $(DOC_DIR)/learn_dg.md README.md
	cp README.md $(DOC_DIR)/README.md
	@ford $(FORD_FLAGS) $(DOC_DIR)/learn_dg.md
	$(RM) $(DOC_DIR)/README.md

# Requires ctags
tags:
	ctags --exclude=build --exclude=templates --exclude=venv --exclude=.git -R .

.ONESHELL:
$(BLD_DIR):
	test -d $(BLD_DIR) || mkdir $(BLD_DIR)

clean:
	$(RM) data.out gmon.out
	$(RM) $(TEST_DIR)/test1D.msh $(TEST_DIR)/test2D.msh
	$(RM) .cache $(TEST_DIR)/__pycache__ $(TEST_DIR)/helpers.pyc
	$(RM) $(BLD_DIR) .pytest_cache
	$(RM) docs/html
