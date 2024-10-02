PELE_ANALYSIS_HOME = ../../../PeleAnalysis

# AMReX
DIM        = 3
PRECISION  = DOUBLE
PROFILE    = FALSE
VERBOSE    = TRUE
DEBUG      = FALSE

# Compiler
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE
TINY_PROFILE = FALSE

CEXE_sources += $(EBASE).cpp

Eos_Model       = Fuego
Chemistry_Model = LiDryer
Transport_Model = Simple

Blocs   := .

include $(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.ModelSpecific
