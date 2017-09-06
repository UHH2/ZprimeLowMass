LIBRARY := SUHH2KinFitTest
USERLDFLAGS := -lSUHH2core -lSUHH2common -lSUHH2JetMETObjects -lSUHH2ZprimeSemiLeptonic -lSUHH2KinFitter -lGenVector
#-lSUHH2KinFitter
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
include ../Makefile.local
