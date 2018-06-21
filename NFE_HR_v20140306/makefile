#FC= gfortran
#FFLAGS= -Wall -pedantic -std=f95 -fbounds-check -O \
#  -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
#  -g -free -m64
FC=ifort
FFLAGS= -O2 -traceback
#FFLAGS= -O2 -axAVX2,AVX,SSE4.2,SSE4.1,SSE3,SSSE3,SSE2
#FFLAGS= -O2 -axAVX2,AVX,SSE4.2,SSE4.1,SSE3,SSSE3,SSE2 -parallel -ipo
#FFLAGS= -fast -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -parallel
#FFLAGS= -check all -warn all -std -fpe0 -traceback -g
#LIBS = -L/usr/lib -lblas -llapack -llapack95
#LIBS=
#INCLUDE= 
#OBJS= makethhrklist.o makethhrinput.o calcthhr.o calcnfehr.o calcff_from_hr.o calcff.o
EXE= makethhrklist makethhrinput calcthhr calcnfehr calcff_from_hr ecalcff ecalcff_old

all: makethhrklist makethhrinput calcthhr calcnfehr calcff_from_hr ecalcff ecalcff_old
makethhrklist: makethhrklist.f90
	$(FC) $(FFLAGS) -o makethhrklist makethhrklist.f90
#makethhrklist.o : makethhrklist.f90
#	$(FC) $(FFLAGS) -c makethhrklist.f90

makethhrinput: makethhrinput.f90
	$(FC) $(FFLAGS) -o makethhrinput makethhrinput.f90
#makethhrinput.o : makethhrinput.f90
#	$(FC) $(FFLAGS) -c makethhrinput.f90

calcthhr: calcthhr.f90
	$(FC) $(FFLAGS) -o calcthhr calcthhr.f90
#calcthhr.o : calcthhr.f90
#	$(FC) $(FFLAGS) -c calcthhr.f90

calcnfehr: calcnfehr.f90
	$(FC) $(FFLAGS) -o calcnfehr calcnfehr.f90
#calcnfehr.o : calcnfehr.f90
#	$(FC) $(FFLAGS) -c calcnfehr.f90

calcff_from_hr: calcff_from_hr.f90
	$(FC) $(FFLAGS) -o calcff_from_hr calcff_from_hr.f90
#calcff.o : calcff.f90
#	$(FC) $(FFLAGS) -c calcff.f90

ecalcff: calcff.f90
	$(FC) $(FFLAGS) -o ecalcff calcff.f90
#calcff.o : calcff.f90
#	$(FC) $(FFLAGS) -c calcff.f90

ecalcff_old: calcff_old.f90
	$(FC) $(FFLAGS) -o ecalcff_old calcff_old.f90
#calcff.o : calcff_old.f90
#	$(FC) $(FFLAGS) -c calcff_old.f90

clean:
	 -rm -f $(EXE) *.o *.mod

