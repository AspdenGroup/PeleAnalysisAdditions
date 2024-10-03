PELE_ANALYSIS_HOME = ../../PeleAnalysis

# AMReX
DEBUG	      = FALSE
DIM	      = 3
PRECISION     = DOUBLE

# Compilation
COMP          = gnu
USE_MPI       = FALSE
USE_OMP       = FALSE
USE_CUDA      = FALSE

EBASE = partStream

CEXE_sources += $(EBASE).cpp

ifeq ($(EBASE),partStream)
  USE_PARTICLES = TRUE
  CEXE_headers += StreamPC.H
  CEXE_sources += StreamPC.cpp
endif

Blocs := .
include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.Analysis
