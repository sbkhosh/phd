###########################################################################
### template for architecture dependent INCOMPACT3D makefiles           ###
### <insert machine name>     		     				###
### <insert a web link to the machine description if available>         ###
###########################################################################

###########################################################################
###                            SWITCHES                                 ###
###########################################################################

# FFTLIB - needed in /src/files.mk
# set to: mkl, fftw, essl or generic
FFTLIB = generic

# double precision should be default; single precision, however, might 
# be sufficient and much faster
PRECISION = double

# switch on DEBUG mode in case you observed errors while running INCOMPACT3D
DEBUG = no

# performance profiling library 
# * switch to none for normal simulations (!)
# * possible choice distributed with INCOMPACT3D is USE_PERFLIB=FR
USE_PERFLIB = none

###########################################################################
#   COMPULSORY LIBRARIES						  #
###########################################################################

#INCLUDE PATHS
INCPATHS = -I$(OBJDIR) -I. -I$(SRCDIR)

ifeq ($(FFTLIB),mkl)
#note: check the existing makefiles as a hint on how to arrange the order
#of the mkl libraries if you need to link them, e.g., for lapack or blas
#routines, as well
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),essl)
 INCPATHS +=
 LIBS +=
 PREPROC += -DWITHESSL
endif
ifeq ($(FFTLIB),fftw)
 INCPATHS +=
 LIBS +=
endif
ifeq ($(FFTLIB),generic)
 INCPATHS +=
 LIBS +=
endif

###########################################################################
###                     COMPILER & COMPILER FLAGS     		  	###
###########################################################################

#FORTRAN COMPILER
FC = mpif90

#ARCHIVE command
ARCHIVE = ar r

#FORTAN COMPILER FLAGS
FFLAGS = -cpp

#LINKER (usually same as FC)
LD =$(FC)

#LDFLAGS
LDFLAGS =

###########################################################################
# ADDITIONAL COMPILER FLAGS (set via SWITCH in header)			  #
###########################################################################
ifeq ($(DEBUG),yes)
#switch off any optimization 
 OPT=	-O0
#and set some flags for backtracing, e.g.:
OPT+= -g -cpp -traceback -warn all
else
 OPT= -O2
#check the compiler documentation for optimization flags
endif

FFLAGS += $(OPT)
LDFLAGS += $(OPT)

ifeq ($(PRECISION),double)
#Mandatory: specify compiler flag for changing the default
#from single to double precision (real = 8 Byte)
FFLAGS += #-r8 (ifort), -qrealsize=8 (xlf), -fdefault-real-8 (gfortran) etc.
PREPROC +=-DDOUBLE_PREC
endif

ifneq ($(USE_PERFLIB),none)
 PREPROC += -DWITHPERF
endif

###########################################################################
### Linking                                                             ###
###########################################################################

include $(FFILES)

#------- mpirun
run:	$(EXECDIR)/$(EXEC)
	ulimit -s unlimited;\
	cd $(RUNDIR);\
	OMP_NUM_THREADS=$(OMP_NUM_THREADS);\
	(mpirun command with $(N_PES) MPI processes) $(EXEC)

#e.g. mpiexec -n $N_PES $(EXEC)

