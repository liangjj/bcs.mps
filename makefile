###############################################################################
#
#  makefile_on_disk template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = gps-mps.x
CPPSRC	= main.cpp \
		  mps_gen.cpp \
		  SpinHamiltonian.cpp \
		  schmidt.cpp \
		  mps_op.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------
BRIGHT_ROOT=.

BOOSTDIR=/home/boxiao/usr
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization -lboost_filesystem -lboost_system -lboost_mpi

BTASINC=-I$(HOME)/mps/btas-master/include
BTASLIB=-L$(HOME)/mps/btas-master/lib
MPSXXINC=-I$(HOME)/mps/mpsxx-master

NEWMATLIB=-L./newmat10 -lnewmat

INCLUDE = $(BTASINC) $(MPSXXINC) $(BOOSTINC)

LIBS=-lpthread -lmkl_intel_lp64 -lmkl_sequential -lmkl_core $(BOOSTLIB) $(BTASLIB) $(NEWMATLIB) -lbtas

#CC	= gcc
#CXX	= g++
CXX = mpic++
CC = mpicc

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= $(INCLUDE) -O3 -std=c++11 -D_HAS_CBLAS -D_HAS_INTEL_MKL

LDFLAGS	= -O3 -std=c++11

# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) -f makefile $(BRIGHT_ROOT)/$(BINNAME)
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@

# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning executable ... '
	@echo $(BINNAME)
	@rm -f $(BINNAME)
	@echo -n '  +++ Cleaning all object files ... '
	@echo $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

# ====================== End of file 'makefile_on_disk.in' ========================== #
